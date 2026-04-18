/*------------------------------------------------------------------------------
 * rtcm_vbs_hp : High-precision RTCM3 Virtual Base Station
 *
 * Pipeline:
 *   TCP source --> input_rtcm3() --> rtcm_in (obs + eph + sta)
 *                                  |
 *                                  +--> rtcm_out (mirrored state)
 *                                        |
 *                                        +-- ret==1 : VBS-correct obs,
 *                                        |           then gen_rtcm3(MSM4/MSM7)
 *                                        +-- ret==2 : gen_rtcm3(eph msg)
 *                                        +-- ret==5 : rewrite sta.pos,
 *                                                    gen_rtcm3(1005/1006)
 *
 * Built against RTKLIB (../src + ../include).
 *----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>

#include "rtklib.h"
#include "vbs_core.h"
#include "tcp_io.h"

/* ---- CLI config ---------------------------------------------------------- */
typedef struct {
    const char *src_host;
    int         src_port;
    const char *eph_host;    /* optional secondary stream just for eph */
    int         eph_port;
    int         out_port;

    int         use_llh_offset;
    double      dn, de, du;
    double      vbs_lat, vbs_lon, vbs_alt;

    int         apply_trop;

    int         msm_out_type;   /* 7 -> MSM7 (default), 4 -> MSM4, 0 -> same as input */
    int         eph_out_enable; /* 1: emit 1019/1020/1042/1045/1046 */

    int         trace_level;
    const char *trace_file;
} config_t;

static volatile int g_stop = 0;
static void on_sigint(int sig) { (void)sig; g_stop = 1; }

/* ---- helpers ------------------------------------------------------------- */
static int msm_type_for_sys(int sys, int subtype)
{
    /* subtype: 4 or 7 */
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
        case SYS_GAL: return 1046;  /* I/NAV (1045 is F/NAV) */
        case SYS_CMP: return 1042;
        case SYS_IRN: return 1041;
        default:      return 0;
    }
}

/* Emit one MSM frame per constellation from rtcm_out->obs.
   Mirrors the essential logic of streamsvr.c::write_rtcm3_msm() but without
   stream/time-interval concerns (we emit whenever upstream gave us obs). */
static void emit_msm(rtcm_t *out, tcp_server_t *srv, int subtype)
{
    obsd_t *data = out->obs.data;
    int nobs = out->obs.n;
    if (nobs <= 0) return;

    /* Keep a master copy; we'll temporarily restamp out->obs per-system */
    obsd_t master[MAXOBS];
    int master_n = nobs < MAXOBS ? nobs : MAXOBS;
    memcpy(master, data, master_n * sizeof(obsd_t));

    static const int systems[] = {
        SYS_GPS, SYS_GLO, SYS_GAL, SYS_CMP, SYS_QZS, SYS_IRN, SYS_SBS
    };

    obsd_t sys_buf[MAXOBS];

    for (size_t k = 0; k < sizeof(systems)/sizeof(systems[0]); k++) {
        int sys = systems[k];
        int msg = msm_type_for_sys(sys, subtype);
        if (!msg) continue;

        /* Gather sats of this sys, count sigs */
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

        /* Pack: nsat * nsig <= 64 */
        int ns   = 64 / nsig;
        int nmsg = (nsat - 1) / ns + 1;

        out->obs.data = sys_buf;   /* RTKLIB encoder reads from here */
        int offset = 0;
        for (int m = 0; m < nmsg; m++) {
            int chunk = (nsat - offset) < ns ? (nsat - offset) : ns;
            /* move chunk to front of sys_buf */
            if (offset > 0) {
                memmove(sys_buf, sys_buf + offset, chunk * sizeof(obsd_t));
            }
            out->obs.n = chunk;
            int sync = (m < nmsg - 1) ? 1 : 0;
            if (gen_rtcm3(out, msg, 0, sync)) {
                tcps_broadcast(srv, out->buff, out->nbyte);
            }
            offset += chunk;
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
        /* emit both 1045 (F/NAV, set=1) and 1046 (I/NAV, set=0) as available */
        if (gen_rtcm3(out, 1046, 0, 0)) tcps_broadcast(srv, out->buff, out->nbyte);
        if (gen_rtcm3(out, 1045, 0, 0)) tcps_broadcast(srv, out->buff, out->nbyte);
    } else {
        if (gen_rtcm3(out, msg, 0, 0)) tcps_broadcast(srv, out->buff, out->nbyte);
    }
}

static void emit_station(rtcm_t *out, tcp_server_t *srv)
{
    /* Prefer 1006 (has antenna height); if deltas are zero encode 1005. */
    int msg = 1006;
    if (gen_rtcm3(out, msg, 0, 0)) {
        tcps_broadcast(srv, out->buff, out->nbyte);
    }
}

/* ---- main loop ----------------------------------------------------------- */
static void run(const config_t *cfg)
{
    if (cfg->trace_file) traceopen(cfg->trace_file);
    tracelevel(cfg->trace_level);

    rtcm_t rin, rout;
    init_rtcm(&rin);
    init_rtcm(&rout);
    /* Approximate current GPST so RTCM time disambiguation works at startup */
    time_t now = time(NULL);
    gtime_t t0 = { now, 0.0 };
    rin.time  = t0;
    rout.time = t0;

    vbs_ctx_t vbs;
    double vbs_llh[3] = { cfg->vbs_lat, cfg->vbs_lon, cfg->vbs_alt };
    vbs_init(&vbs,
             cfg->use_llh_offset,
             cfg->dn, cfg->de, cfg->du,
             cfg->use_llh_offset ? NULL : vbs_llh,
             cfg->apply_trop);

    tcp_server_t *srv = tcps_create(cfg->out_port, 16);
    if (!srv) {
        fprintf(stderr, "[error] cannot bind output port %d\n", cfg->out_port);
        return;
    }

    tcp_client_t src;  src.fd  = -1;
    tcp_client_t ephc; ephc.fd = -1;

    int subtype = cfg->msm_out_type ? cfg->msm_out_type : 7;
    long last_stat_t = (long)time(NULL);

    uint8_t buf[4096];
    while (!g_stop) {
        /* ---- (re)connect source ---- */
        if (src.fd < 0) {
            if (tcpc_connect(&src, cfg->src_host, cfg->src_port) == 0) {
                fprintf(stderr, "[src] connected %s:%d\n",
                        cfg->src_host, cfg->src_port);
            } else {
                fprintf(stderr, "[src] connect %s:%d failed, retry in 3 s\n",
                        cfg->src_host, cfg->src_port);
                for (int i = 0; i < 30 && !g_stop; i++) {
                    struct timespec ts = {0, 100*1000*1000};
                    nanosleep(&ts, NULL);
                }
                continue;
            }
        }

        if (cfg->eph_port > 0 && ephc.fd < 0) {
            if (tcpc_connect(&ephc, cfg->eph_host, cfg->eph_port) == 0) {
                fprintf(stderr, "[eph] connected %s:%d\n",
                        cfg->eph_host, cfg->eph_port);
            }
        }

        /* ---- drain source ---- */
        int r = tcpc_read(&src, buf, sizeof(buf));
        if (r <= 0) {
            fprintf(stderr, "[src] disconnected, will reconnect\n");
            tcpc_close(&src);
            continue;
        }

        for (int i = 0; i < r; i++) {
            int ret = input_rtcm3(&rin, buf[i]);
            if (ret == 0 || ret < 0) continue;

            /* mirror time/staid */
            rout.time  = rin.time;
            rout.staid = rin.staid;

            if (ret == 1) {
                /* --- observation frame --- */
                /* copy obs into rout, then correct */
                int n = rin.obs.n;
                if (n > MAXOBS) n = MAXOBS;
                memcpy(rout.obs.data, rin.obs.data, n * sizeof(obsd_t));
                rout.obs.n = n;
                /* relay GLONASS FCNs */
                memcpy(rout.nav.glo_fcn, rin.nav.glo_fcn,
                       sizeof(rin.nav.glo_fcn));

                vbs_correct_obs(&vbs, &rout);
                if (rout.obs.n > 0) {
                    emit_msm(&rout, srv, subtype);
                }
            }
            else if (ret == 2) {
                /* --- ephemeris frame --- */
                int sat = rin.ephsat, set = rin.ephset;
                int sys = satsys(sat, NULL);
                if (sys == SYS_GLO) {
                    int prn; satsys(sat, &prn);
                    rout.nav.geph[prn-1] = rin.nav.geph[prn-1];
                } else if (sys != 0) {
                    rout.nav.eph[sat-1 + MAXSAT*set] =
                        rin.nav.eph[sat-1 + MAXSAT*set];
                }
                rout.ephsat = sat;
                rout.ephset = set;

                /* also propagate into rin's copy used by satposs of rout */
                /* (already copied above) */

                /* keep eph for own computations even if not emitting */
                if (cfg->eph_out_enable) emit_eph(&rout, srv, sat);

                /* We must also feed the *input* nav into rout.nav fully,
                   because satposs reads from rtcm->nav. */
            }
            else if (ret == 5) {
                /* --- station params --- */
                rout.sta = rin.sta;
                vbs_update_base_from_sta(&vbs, &rin.sta);
                vbs_rewrite_station(&vbs, &rout);
                emit_station(&rout, srv);
            }
            else {
                /* 10: SSR; we don't VBS-correct them, but pass-through is
                   non-trivial with gen_rtcm3. Skip for now. */
            }
        }

        /* ---- drain optional ephemeris-only stream ---- */
        if (ephc.fd >= 0) {
            int re = tcpc_read(&ephc, buf, sizeof(buf));
            if (re > 0) {
                for (int i = 0; i < re; i++) {
                    int ret = input_rtcm3(&rin, buf[i]);
                    if (ret == 2) {
                        int sat = rin.ephsat, set = rin.ephset;
                        int sys = satsys(sat, NULL);
                        if (sys == SYS_GLO) {
                            int prn; satsys(sat, &prn);
                            rout.nav.geph[prn-1] = rin.nav.geph[prn-1];
                        } else if (sys != 0) {
                            rout.nav.eph[sat-1 + MAXSAT*set] =
                                rin.nav.eph[sat-1 + MAXSAT*set];
                        }
                        if (cfg->eph_out_enable) emit_eph(&rout, srv, sat);
                    }
                }
            } else if (re == 0 || re < 0) {
                fprintf(stderr, "[eph] disconnected\n");
                tcpc_close(&ephc);
            }
        }

        /* ---- periodic stats ---- */
        long now_s = (long)time(NULL);
        if (now_s - last_stat_t >= 15) {
            fprintf(stderr,
                "[stat] frames: 1005/6=%ld MSM=%ld | obs in=%ld out=%ld "
                "| sat corr=%ld no_eph=%ld | clients=%d\n",
                vbs.n_frames_1005_6, vbs.n_frames_msm,
                vbs.n_obs_in, vbs.n_obs_out,
                vbs.n_sat_corrected, vbs.n_sat_no_eph,
                tcps_num_clients(srv));
            last_stat_t = now_s;
        }
    }

    tcpc_close(&src);
    tcpc_close(&ephc);
    tcps_destroy(srv);
    free_rtcm(&rin);
    free_rtcm(&rout);
}

/* ---- CLI parser ---------------------------------------------------------- */
static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Source:\n"
        "  --src-host HOST       source RTCM TCP host (default 127.0.0.1)\n"
        "  --src-port PORT       source RTCM TCP port (default 50001)\n"
        "  --eph-host HOST       optional ephemeris-only stream host\n"
        "  --eph-port PORT       optional ephemeris-only stream port\n"
        "  --out-port  PORT      VBS broadcast output port (default 50002)\n"
        "\n"
        "VBS position (pick ONE mode):\n"
        "  --offset N E U        metres, relative to the learned base (NED)\n"
        "  --vbs LAT LON ALT     absolute VBS, lat/lon in degrees, alt in m\n"
        "\n"
        "Corrections:\n"
        "  --trop                enable differential Saastamoinen troposphere\n"
        "  --msm N               output MSM subtype (4 or 7, default 7)\n"
        "  --no-eph              do not re-emit ephemeris messages\n"
        "\n"
        "Diagnostics:\n"
        "  --trace FILE          open RTKLIB trace log\n"
        "  --level  N            RTKLIB trace level (1..5)\n"
        "\n", prog);
}

int main(int argc, char **argv)
{
    config_t cfg = {
        .src_host = "127.0.0.1",
        .src_port = 50001,
        .eph_host = "127.0.0.1",
        .eph_port = 0,
        .out_port = 50002,
        .use_llh_offset = 1,
        .dn = 0.0, .de = 0.0, .du = 0.0,
        .vbs_lat = 0.0, .vbs_lon = 0.0, .vbs_alt = 0.0,
        .apply_trop = 0,
        .msm_out_type = 7,
        .eph_out_enable = 1,
        .trace_level = 0,
        .trace_file = NULL,
    };

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--src-host") && i+1<argc) cfg.src_host = argv[++i];
        else if (!strcmp(argv[i], "--src-port") && i+1<argc) cfg.src_port = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--eph-host") && i+1<argc) cfg.eph_host = argv[++i];
        else if (!strcmp(argv[i], "--eph-port") && i+1<argc) cfg.eph_port = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--out-port") && i+1<argc) cfg.out_port = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--offset")   && i+3<argc) {
            cfg.use_llh_offset = 1;
            cfg.dn = atof(argv[++i]);
            cfg.de = atof(argv[++i]);
            cfg.du = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--vbs")      && i+3<argc) {
            cfg.use_llh_offset = 0;
            cfg.vbs_lat = atof(argv[++i]);
            cfg.vbs_lon = atof(argv[++i]);
            cfg.vbs_alt = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--trop"))   cfg.apply_trop = 1;
        else if (!strcmp(argv[i], "--msm") && i+1<argc) cfg.msm_out_type = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--no-eph")) cfg.eph_out_enable = 0;
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
        "rtcm_vbs_hp: src=%s:%d -> out=:%d | mode=%s | trop=%d | MSM=%d\n",
        cfg.src_host, cfg.src_port, cfg.out_port,
        cfg.use_llh_offset ? "NED-offset" : "absolute-LLH",
        cfg.apply_trop, cfg.msm_out_type);
    if (cfg.use_llh_offset) {
        fprintf(stderr, "        offset N=%.3f E=%.3f U=%.3f m\n",
                cfg.dn, cfg.de, cfg.du);
    } else {
        fprintf(stderr, "        VBS  LAT=%.8f LON=%.8f ALT=%.3f\n",
                cfg.vbs_lat, cfg.vbs_lon, cfg.vbs_alt);
    }

    run(&cfg);
    return 0;
}
