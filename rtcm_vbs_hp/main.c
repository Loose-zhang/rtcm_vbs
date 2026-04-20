/*------------------------------------------------------------------------------
 * rtcm_vbs_hp : High-precision RTCM3 Virtual Base Station
 *
 * Single-base pipeline:
 *   TCP src   -> input_rtcm3 -> ret dispatch
 *     ret==1 (MSM obs) -> vbs_correct_obs -> gen_rtcm3(MSM) -> broadcast
 *     ret==2 (eph)     -> mirror nav + optional gen_rtcm3(eph)
 *     ret==5 (sta)     -> rewrite sta.pos -> gen_rtcm3(1005/6)
 *     ret==10 (SSR)    -> raw pass-through (VBS-agnostic)
 *
 * Dual-base pipeline (--dual):
 *   Two src channels (A, B). Obs epochs are stashed into merger_t which
 *   emits the CNR-priority merged epoch; the merged set is then VBS-
 *   corrected and MSM-encoded as in single-base mode.
 *
 * Optional secondary streams:
 *   --eph-port : ephemeris-only stream (1019/1020/1042/1045/1046)
 *   --ssr-port : SSR-only stream (any SSR MT). Passed through untouched
 *                to the output listeners.
 *----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>

#include "rtklib.h"
#include "vbs_core.h"
#include "vbs_merge.h"
#include "tcp_io.h"

/* ---- CLI config ---------------------------------------------------------- */
typedef struct {
    /* primary / A source */
    const char *src_host;
    int         src_port;
    /* dual-base B source */
    int         dual;
    const char *b_host;
    int         b_port;
    int         epoch_window_ms;

    /* optional feeders */
    const char *eph_host;
    int         eph_port;
    const char *ssr_host;
    int         ssr_port;

    /* output */
    int         out_port;

    /* VBS geometry */
    int         use_llh_offset;
    double      dn, de, du;
    double      vbs_lat, vbs_lon, vbs_alt;

    int         apply_trop;
    int         msm_out_type;
    int         eph_out_enable;
    int         ssr_passthrough;

    int         trace_level;
    const char *trace_file;
} config_t;

static volatile int g_stop = 0;
static void on_sigint(int sig) { (void)sig; g_stop = 1; }

/* ---- SSR message-type classification ------------------------------------- */
static int is_ssr_msg(int mtype)
{
    /* Standard + draft SSR ranges per RTCM3 / IGS SSR */
    if (mtype >= 1057 && mtype <= 1068) return 1;
    if (mtype >= 1240 && mtype <= 1254) return 1;
    if (mtype >= 1258 && mtype <= 1270) return 1;
    if (mtype == 11 || mtype == 12 || mtype == 13 || mtype == 14) return 1;
    return 0;
}

/* Parse message type from a complete RTCM3 frame (buff[0]==0xD3). */
static int frame_mtype(const uint8_t *buff, int len)
{
    if (len < 6 || buff[0] != 0xD3) return -1;
    /* 12 bits starting at byte 3 (bit 24) */
    return (int)getbitu(buff, 24, 12);
}

/* ---- helpers: output encoding ------------------------------------------- */
static int msm_type_for_sys(int sys, int subtype)
{
    int base;
    switch (sys) {
        case SYS_GPS: base = 1070; break;
        case SYS_GLO: base = 1080; break;
        case SYS_GAL: base = 1090; break;
        case SYS_SBS: base = 1100; break;
        case SYS_QZS: base = 1110; break;
        case SYS_CMP: base = 1120; break;
        case SYS_IRN: base = 1130; break;
        default:      return 0;
    }
    return base + subtype;
}

static int eph_msg_for_sys(int sys)
{
    switch (sys) {
        case SYS_GPS: return 1019;
        case SYS_GLO: return 1020;
        case SYS_QZS: return 1044;
        case SYS_GAL: return 1046;
        case SYS_CMP: return 1042;
        case SYS_IRN: return 1041;
        default:      return 0;
    }
}

/* Emit one or more MSM frames per constellation from rtcm_out->obs.
 * RTCM MSM sync: sync=1 means more MSM for this epoch follow; sync=0 ends
 * the epoch. RTKLIB (RTKNAVI) only delivers obs to positioning on sync=0, and
 * merges all prior sync=1 MSM into one epoch — so only the final frame of
 * this emit may use sync=0 (not once per GNSS). */
static void emit_msm(rtcm_t *out, tcp_server_t *srv, int subtype)
{
    obsd_t *data = out->obs.data;
    int nobs = out->obs.n;
    if (nobs <= 0) return;

    obsd_t master[MAXOBS];
    int master_n = nobs < MAXOBS ? nobs : MAXOBS;
    memcpy(master, data, master_n * sizeof(obsd_t));

    static const int systems[] = {
        SYS_GPS, SYS_GLO, SYS_GAL, SYS_CMP, SYS_QZS, SYS_IRN, SYS_SBS
    };
    obsd_t sys_buf[MAXOBS];

    int total = 0;
    for (size_t k = 0; k < sizeof(systems)/sizeof(systems[0]); k++) {
        int sys = systems[k];
        int msg = msm_type_for_sys(sys, subtype);
        if (!msg) continue;

        int nsat = 0, nsig = 0, mask[MAXCODE] = {0};
        for (int i = 0; i < master_n; i++) {
            if (satsys(master[i].sat, NULL) != sys) continue;
            for (int j = 0; j < NFREQ+NEXOBS; j++) {
                int code = master[i].code[j];
                if (!code || mask[code-1]) continue;
                mask[code-1] = 1; nsig++;
            }
            nsat++;
        }
        if (nsat == 0 || nsig == 0) continue;
        if (nsig > 64) continue;

        int ns   = 64 / nsig;
        int nmsg = (nsat - 1) / ns + 1;
        total += nmsg;
    }
    if (total <= 0) return;

    int sent = 0;
    for (size_t k = 0; k < sizeof(systems)/sizeof(systems[0]); k++) {
        int sys = systems[k];
        int msg = msm_type_for_sys(sys, subtype);
        if (!msg) continue;

        int nsat = 0, nsig = 0, mask[MAXCODE] = {0};
        for (int i = 0; i < master_n; i++) {
            if (satsys(master[i].sat, NULL) != sys) continue;
            sys_buf[nsat++] = master[i];
            for (int j = 0; j < NFREQ+NEXOBS; j++) {
                int code = master[i].code[j];
                if (!code || mask[code-1]) continue;
                mask[code-1] = 1; nsig++;
            }
        }
        if (nsat == 0 || nsig == 0) continue;
        if (nsig > 64) continue;

        int ns   = 64 / nsig;
        int nmsg = (nsat - 1) / ns + 1;

        out->obs.data = sys_buf;
        int offset = 0;
        for (int m = 0; m < nmsg; m++) {
            int chunk = (nsat - offset) < ns ? (nsat - offset) : ns;
            if (offset > 0) {
                memmove(sys_buf, sys_buf + offset, chunk * sizeof(obsd_t));
            }
            out->obs.n = chunk;
            int sync = (sent < total - 1) ? 1 : 0;
            if (gen_rtcm3(out, msg, 0, sync)) {
                tcps_broadcast(srv, out->buff, out->nbyte);
            }
            offset += chunk;
            sent++;
        }
        out->obs.data = data;
        out->obs.n    = nobs;
    }
}

static void emit_eph(rtcm_t *out, tcp_server_t *srv, int sat)
{
    int sys = satsys(sat, NULL);
    int msg = eph_msg_for_sys(sys);
    if (!msg) return;
    if (sys == SYS_GAL) {
        if (gen_rtcm3(out, 1046, 0, 0)) tcps_broadcast(srv, out->buff, out->nbyte);
        if (gen_rtcm3(out, 1045, 0, 0)) tcps_broadcast(srv, out->buff, out->nbyte);
    } else {
        if (gen_rtcm3(out, msg, 0, 0)) tcps_broadcast(srv, out->buff, out->nbyte);
    }
}

static void emit_station(rtcm_t *out, tcp_server_t *srv)
{
    if (gen_rtcm3(out, 1006, 0, 0)) {
        tcps_broadcast(srv, out->buff, out->nbyte);
    }
}

/* ---- monotonic milliseconds --------------------------------------------- */
static long now_ms_mono(void)
{
#ifdef _WIN32
    return (long)GetTickCount64();
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (long)ts.tv_sec * 1000 + ts.tv_nsec / 1000000;
#endif
}

/* ---- per-channel RTCM state + input processing -------------------------- */
typedef struct {
    tcp_client_t tcp;
    rtcm_t       rtcm;
    const char  *host;
    int          port;
    const char  *tag;
} channel_t;

static void channel_init(channel_t *c, const char *host, int port, const char *tag)
{
    c->tcp.fd = -1;
    c->host = host;
    c->port = port;
    c->tag  = tag;
    init_rtcm(&c->rtcm);
    c->rtcm.time.time = time(NULL);
    c->rtcm.time.sec  = 0.0;
}

static void channel_close(channel_t *c)
{
    tcpc_close(&c->tcp);
    free_rtcm(&c->rtcm);
}

/* Attempt to connect if not yet. Non-blocking one-shot; caller is the loop. */
static void channel_try_connect(channel_t *c)
{
    if (c->tcp.fd >= 0) return;
    if (tcpc_connect(&c->tcp, c->host, c->port) == 0) {
        fprintf(stderr, "[%s] connected %s:%d\n", c->tag, c->host, c->port);
    }
}

/* ---- nav mirror helpers ------------------------------------------------- */
static void mirror_eph(rtcm_t *out, const rtcm_t *in)
{
    int sat = in->ephsat, set = in->ephset;
    int sys = satsys(sat, NULL), prn = 0;
    satsys(sat, &prn);
    if (sys == SYS_GLO) {
        out->nav.geph[prn-1] = in->nav.geph[prn-1];
    } else if (sys != 0) {
        out->nav.eph[sat-1 + MAXSAT*set] = in->nav.eph[sat-1 + MAXSAT*set];
    }
    out->ephsat = sat;
    out->ephset = set;
}

/* Mirror GLONASS FCN (frequency channel number) table used by signal
 * wavelength lookups. */
static void mirror_glo_fcn(rtcm_t *out, const rtcm_t *in)
{
    memcpy(out->nav.glo_fcn, in->nav.glo_fcn, sizeof(in->nav.glo_fcn));
}

/* ---- main loop ----------------------------------------------------------- */
static void run(const config_t *cfg)
{
    if (cfg->trace_file) traceopen(cfg->trace_file);
    tracelevel(cfg->trace_level);

    channel_t chA, chB, chEph, chSsr;
    channel_init(&chA,   cfg->src_host, cfg->src_port, "A");
    if (cfg->dual)       channel_init(&chB,   cfg->b_host,   cfg->b_port,   "B");
    if (cfg->eph_port)   channel_init(&chEph, cfg->eph_host, cfg->eph_port, "eph");
    if (cfg->ssr_port)   channel_init(&chSsr, cfg->ssr_host, cfg->ssr_port, "ssr");

    /* Shared output RTCM state (encoder + mirrored nav/sta) */
    rtcm_t rout;
    init_rtcm(&rout);
    rout.time.time = time(NULL);

    vbs_ctx_t vbs;
    double vbs_llh[3] = { cfg->vbs_lat, cfg->vbs_lon, cfg->vbs_alt };
    vbs_init(&vbs, cfg->use_llh_offset, cfg->dn, cfg->de, cfg->du,
             cfg->use_llh_offset ? NULL : vbs_llh, cfg->apply_trop);

    merger_t mrg;
    merger_init(&mrg, cfg->epoch_window_ms);

    tcp_server_t *srv = tcps_create(cfg->out_port, 16);
    if (!srv) {
        fprintf(stderr, "[error] cannot bind output port %d\n", cfg->out_port);
        return;
    }

    int subtype = cfg->msm_out_type ? cfg->msm_out_type : 7;
    long last_stat_t = (long)time(NULL);
    long ssr_fwd = 0;

    uint8_t buf[4096];
    while (!g_stop) {
        channel_try_connect(&chA);
        if (cfg->dual)     channel_try_connect(&chB);
        if (cfg->eph_port) channel_try_connect(&chEph);
        if (cfg->ssr_port) channel_try_connect(&chSsr);

        int any_data = 0;

        /* ---- iterate over all input channels ---- */
        channel_t *channels[4] = { &chA, NULL, NULL, NULL };
        int nch = 1;
        if (cfg->dual)     channels[nch++] = &chB;
        if (cfg->eph_port) channels[nch++] = &chEph;
        if (cfg->ssr_port) channels[nch++] = &chSsr;

        for (int ci = 0; ci < nch; ci++) {
            channel_t *ch = channels[ci];
            if (ch->tcp.fd < 0) continue;

            int r = tcpc_read(&ch->tcp, buf, sizeof(buf));
            if (r == 0 || r < 0) {
                if (r != 0 || ch->tcp.fd >= 0) {
                    fprintf(stderr, "[%s] disconnected\n", ch->tag);
                }
                tcpc_close(&ch->tcp);
                continue;
            }
            any_data = 1;

            for (int i = 0; i < r; i++) {
                int ret = input_rtcm3(&ch->rtcm, buf[i]);
                if (ret == 0) continue;

                /* At this point ch->rtcm.buff[0..rtcm.len+3) holds the full
                   frame that was just consumed (preamble..CRC). */
                int mtype = frame_mtype(ch->rtcm.buff, ch->rtcm.len + 3);

                /* --- SSR: raw pass-through regardless of channel --- */
                if (is_ssr_msg(mtype)) {
                    if (cfg->ssr_passthrough) {
                        tcps_broadcast(srv, ch->rtcm.buff, ch->rtcm.len + 3);
                        ssr_fwd++;
                    }
                    continue;
                }

                /* From here on, dispatch by decoded return code. */
                rout.time  = ch->rtcm.time;
                rout.staid = ch->rtcm.staid;

                if (ret == 1) {
                    /* observation epoch */
                    mirror_glo_fcn(&rout, &ch->rtcm);

                    if (cfg->dual && (ch == &chA || ch == &chB)) {
                        int side = (ch == &chA) ? 0 : 1;
                        merger_stash(&mrg, side, ch->rtcm.obs.data,
                                     ch->rtcm.obs.n, ch->rtcm.time,
                                     now_ms_mono());
                    } else if (ch == &chA) {
                        /* single-base: correct + emit immediately; all obs
                           come from base A so side array is NULL/implicit. */
                        int n = ch->rtcm.obs.n;
                        if (n > MAXOBS) n = MAXOBS;
                        memcpy(rout.obs.data, ch->rtcm.obs.data,
                               n * sizeof(obsd_t));
                        rout.obs.n = n;
                        vbs_correct_obs(&vbs, &rout, NULL);
                        if (rout.obs.n > 0) emit_msm(&rout, srv, subtype);
                    }
                    /* eph/ssr channels typically don't carry obs; ignore */
                }
                else if (ret == 2) {
                    mirror_eph(&rout, &ch->rtcm);
                    if (cfg->eph_out_enable) {
                        emit_eph(&rout, srv, ch->rtcm.ephsat);
                    }
                }
                else if (ret == 5) {
                    /* Learn whichever real base this channel represents.
                       In dual-base mode BOTH bases need to be known so
                       per-sat correction can pick the right one.
                       The outgoing 1005/1006 frame always carries the VBS
                       coordinate (from channel A's frame to preserve the
                       original staid / antenna info). */
                    int side = -1;
                    if (ch == &chA) side = 0;
                    else if (cfg->dual && ch == &chB) side = 1;

                    if (side >= 0) {
                        vbs_update_base(&vbs, side, &ch->rtcm.sta);
                    }
                    if (ch == &chA) {
                        rout.sta = ch->rtcm.sta;
                        vbs_rewrite_station(&vbs, &rout);
                        emit_station(&rout, srv);
                    }
                }
                /* ret==10 handled by SSR passthrough above; skip here */
            }
        }

        /* ---- dual-base: poll merger & emit if ready ---- */
        if (cfg->dual) {
            obsd_t merged[MAXOBS];
            unsigned char side_arr[MAXOBS];
            int    nm = 0;
            gtime_t mt;
            if (merger_poll(&mrg, now_ms_mono(),
                            merged, side_arr, &nm, &mt) && nm > 0) {
                memcpy(rout.obs.data, merged, nm * sizeof(obsd_t));
                rout.obs.n = nm;
                rout.time  = mt;
                vbs_correct_obs(&vbs, &rout, side_arr);
                if (rout.obs.n > 0) emit_msm(&rout, srv, subtype);
            }
        }

        /* brief sleep if nothing happened, to avoid busy-spin */
        if (!any_data) {
#ifdef _WIN32
            Sleep(5);
#else
            struct timespec ts = { 0, 5*1000*1000 };  /* 5 ms */
            nanosleep(&ts, NULL);
#endif
        }

        /* ---- periodic stats ---- */
        long now_s = (long)time(NULL);
        if (now_s - last_stat_t >= 15) {
            if (cfg->dual) {
                fprintf(stderr,
                    "[stat] 1005/6=%ld MSM=%ld | obs in=%ld out=%ld | "
                    "sat corr=%ld no_eph=%ld | epochs merged=%ld A-only=%ld "
                    "B-only=%ld | CNR both_sat=%ld chose_A=%ld chose_B=%ld "
                    "| SSR fwd=%ld | clients=%d\n",
                    vbs.n_frames_1005_6, vbs.n_frames_msm,
                    vbs.n_obs_in, vbs.n_obs_out,
                    vbs.n_sat_corrected, vbs.n_sat_no_eph,
                    mrg.n_merged, mrg.n_a_only, mrg.n_b_only,
                    mrg.n_both, mrg.n_chose_a, mrg.n_chose_b,
                    ssr_fwd, tcps_num_clients(srv));
            } else {
                fprintf(stderr,
                    "[stat] 1005/6=%ld MSM=%ld | obs in=%ld out=%ld | "
                    "sat corr=%ld no_eph=%ld | SSR fwd=%ld | clients=%d\n",
                    vbs.n_frames_1005_6, vbs.n_frames_msm,
                    vbs.n_obs_in, vbs.n_obs_out,
                    vbs.n_sat_corrected, vbs.n_sat_no_eph,
                    ssr_fwd, tcps_num_clients(srv));
            }
            last_stat_t = now_s;
        }
    }

    channel_close(&chA);
    if (cfg->dual)     channel_close(&chB);
    if (cfg->eph_port) channel_close(&chEph);
    if (cfg->ssr_port) channel_close(&chSsr);
    tcps_destroy(srv);
    free_rtcm(&rout);
}

/* ---- CLI parser ---------------------------------------------------------- */
static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Source (primary / base A):\n"
        "  --src-host HOST       (default 127.0.0.1)\n"
        "  --src-port PORT       (default 50001)\n"
        "\n"
        "Dual-base mode:\n"
        "  --dual                enable dual-base CNR-merge mode\n"
        "  --b-host HOST         base B host (default 127.0.0.1)\n"
        "  --b-port PORT         base B TCP port\n"
        "  --window MS           epoch alignment window ms (default 40)\n"
        "\n"
        "Optional feeder streams:\n"
        "  --eph-host HOST / --eph-port PORT  ephemeris-only stream\n"
        "  --ssr-host HOST / --ssr-port PORT  SSR-only stream (passthrough)\n"
        "\n"
        "Output:\n"
        "  --out-port  PORT      VBS broadcast port (default 50002)\n"
        "\n"
        "VBS position (pick ONE mode):\n"
        "  --offset N E U        metres NED offset from learned base\n"
        "  --vbs LAT LON ALT     absolute VBS (deg, deg, m)\n"
        "\n"
        "Corrections:\n"
        "  --trop                enable differential Saastamoinen troposphere\n"
        "  --msm N               output MSM subtype 4 or 7 (default 7)\n"
        "  --no-eph              do not re-emit ephemeris messages\n"
        "  --no-ssr              do not forward SSR messages from any channel\n"
        "\n"
        "Diagnostics:\n"
        "  --trace FILE / --level N\n"
        "\n", prog);
}

int main(int argc, char **argv)
{
    config_t cfg = {
        .src_host = "127.0.0.1", .src_port = 50001,
        .dual = 0,
        .b_host = "127.0.0.1", .b_port = 50003,
        .epoch_window_ms = 40,
        .eph_host = "127.0.0.1", .eph_port = 0,
        .ssr_host = "127.0.0.1", .ssr_port = 0,
        .out_port = 50002,
        .use_llh_offset = 1,
        .dn = 0.0, .de = 0.0, .du = 0.0,
        .vbs_lat = 0.0, .vbs_lon = 0.0, .vbs_alt = 0.0,
        .apply_trop = 0,
        .msm_out_type = 7,
        .eph_out_enable = 1,
        .ssr_passthrough = 1,
        .trace_level = 0, .trace_file = NULL,
    };

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--src-host") && i+1<argc) cfg.src_host = argv[++i];
        else if (!strcmp(argv[i], "--src-port") && i+1<argc) cfg.src_port = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--dual"))                  cfg.dual = 1;
        else if (!strcmp(argv[i], "--b-host")   && i+1<argc) cfg.b_host = argv[++i];
        else if (!strcmp(argv[i], "--b-port")   && i+1<argc) cfg.b_port = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--window")   && i+1<argc) cfg.epoch_window_ms = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--eph-host") && i+1<argc) cfg.eph_host = argv[++i];
        else if (!strcmp(argv[i], "--eph-port") && i+1<argc) cfg.eph_port = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--ssr-host") && i+1<argc) cfg.ssr_host = argv[++i];
        else if (!strcmp(argv[i], "--ssr-port") && i+1<argc) cfg.ssr_port = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--out-port") && i+1<argc) cfg.out_port = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--offset")   && i+3<argc) {
            cfg.use_llh_offset = 1;
            cfg.dn = atof(argv[++i]);
            cfg.de = atof(argv[++i]);
            cfg.du = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--vbs") && i+3<argc) {
            cfg.use_llh_offset = 0;
            cfg.vbs_lat = atof(argv[++i]);
            cfg.vbs_lon = atof(argv[++i]);
            cfg.vbs_alt = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--trop"))    cfg.apply_trop = 1;
        else if (!strcmp(argv[i], "--msm") && i+1<argc) cfg.msm_out_type = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--no-eph"))  cfg.eph_out_enable = 0;
        else if (!strcmp(argv[i], "--no-ssr"))  cfg.ssr_passthrough = 0;
        else if (!strcmp(argv[i], "--trace") && i+1<argc) cfg.trace_file = argv[++i];
        else if (!strcmp(argv[i], "--level") && i+1<argc) cfg.trace_level = atoi(argv[++i]);
        else { usage(argv[0]); return 1; }
    }

    if (cfg.msm_out_type != 4 && cfg.msm_out_type != 7) {
        fprintf(stderr, "[cfg] --msm must be 4 or 7\n");
        return 1;
    }

    signal(SIGINT,  on_sigint);
    signal(SIGTERM, on_sigint);
#ifndef _WIN32
    signal(SIGPIPE, SIG_IGN);
#endif

    fprintf(stderr,
        "rtcm_vbs_hp: A=%s:%d -> out=:%d | mode=%s | trop=%d | MSM=%d "
        "| SSR fwd=%d\n",
        cfg.src_host, cfg.src_port, cfg.out_port,
        cfg.use_llh_offset ? "NED-offset" : "absolute-LLH",
        cfg.apply_trop, cfg.msm_out_type, cfg.ssr_passthrough);
    if (cfg.dual) fprintf(stderr, "        B=%s:%d  epoch window=%d ms\n",
                          cfg.b_host, cfg.b_port, cfg.epoch_window_ms);
    if (cfg.eph_port) fprintf(stderr, "        eph stream %s:%d\n",
                              cfg.eph_host, cfg.eph_port);
    if (cfg.ssr_port) fprintf(stderr, "        ssr stream %s:%d\n",
                              cfg.ssr_host, cfg.ssr_port);
    if (cfg.use_llh_offset) {
        fprintf(stderr, "        offset N=%.3f E=%.3f U=%.3f m\n",
                cfg.dn, cfg.de, cfg.du);
    } else {
        fprintf(stderr, "        VBS LAT=%.8f LON=%.8f ALT=%.3f\n",
                cfg.vbs_lat, cfg.vbs_lon, cfg.vbs_alt);
    }
    run(&cfg);
    return 0;
}
