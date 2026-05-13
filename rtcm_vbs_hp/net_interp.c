/*------------------------------------------------------------------------------
 * net_interp.c : lightweight dual-base network/ionosphere interpolation
 *----------------------------------------------------------------------------*/
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "net_interp.h"

#define NET_MIN_ELEV (10.0 * D2R)

static int usable_satpos(const double *rs, int svh)
{
    return svh >= 0 && (rs[0] != 0.0 || rs[1] != 0.0 || rs[2] != 0.0);
}

static void satposs_with_fallback(gtime_t time, const obsd_t *obs, int n,
                                  const nav_t *nav, int ephopt,
                                  double *rs, double *dts, double *var,
                                  int *svh)
{
    static double rs_brdc[6*MAXOBS];
    static double dts_brdc[2*MAXOBS];
    static double var_brdc[MAXOBS];
    static int svh_brdc[MAXOBS];

    satposs(time, obs, n, nav, ephopt, rs, dts, var, svh);
    if (ephopt == EPHOPT_BRDC) return;

    satposs(time, obs, n, nav, EPHOPT_BRDC, rs_brdc, dts_brdc, var_brdc, svh_brdc);
    for (int i = 0; i < n && i < MAXOBS; i++) {
        if (usable_satpos(rs + 6*i, svh[i])) continue;
        if (!usable_satpos(rs_brdc + 6*i, svh_brdc[i])) continue;
        memcpy(rs + 6*i, rs_brdc + 6*i, sizeof(double) * 6);
        memcpy(dts + 2*i, dts_brdc + 2*i, sizeof(double) * 2);
        var[i] = var_brdc[i];
        svh[i] = svh_brdc[i];
    }
}

static const obsd_t *find_sat(const obsd_t *obs, int n, int sat)
{
    for (int i = 0; i < n; i++) {
        if (obs[i].sat == sat) return &obs[i];
    }
    return NULL;
}

static double first_common_p(const obsd_t *a, const obsd_t *b,
                             double *pa, double *pb)
{
    for (int ia = 0; ia < NFREQ + NEXOBS; ia++) {
        if (a->P[ia] == 0.0 || !a->code[ia]) continue;
        for (int ib = 0; ib < NFREQ + NEXOBS; ib++) {
            if (b->P[ib] == 0.0 || b->code[ib] != a->code[ia]) continue;
            *pa = a->P[ia];
            *pb = b->P[ib];
            return 1.0;
        }
    }
    return 0.0;
}

static int cmp_double(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

static double median_of(double *v, int n)
{
    if (n <= 0) return 0.0;
    qsort(v, (size_t)n, sizeof(double), cmp_double);
    return (n & 1) ? v[n/2] : 0.5 * (v[n/2 - 1] + v[n/2]);
}

/* Fix 2: IDW weights based on 3-D distance from VBS to each base.
 * wa + wb = 1.0; no extrapolation needed since IDW is naturally bounded. */
static void idw_weights(const double *vbs_ecef,
                        const double *a_ecef, const double *b_ecef,
                        double *wa, double *wb)
{
    double da = 0.0, db = 0.0;
    for (int k = 0; k < 3; k++) {
        double ta = vbs_ecef[k] - a_ecef[k];
        double tb = vbs_ecef[k] - b_ecef[k];
        da += ta * ta;
        db += tb * tb;
    }
    da = sqrt(da);
    db = sqrt(db);

    if (da < 1.0) { *wa = 1.0; *wb = 0.0; return; }
    if (db < 1.0) { *wa = 0.0; *wb = 1.0; return; }

    double ra = 1.0 / da;
    double rb = 1.0 / db;
    double rs = ra + rb;
    *wa = ra / rs;
    *wb = rb / rs;
}

void net_interp_init(net_interp_t *ni, int enable)
{
    memset(ni, 0, sizeof(*ni));
    ni->enable = enable;
    ni->max_abs_corr_m = 30.0;
}

int net_interp_build_dual(net_interp_t *ni,
                          const vbs_ctx_t *vbs,
                          const rtcm_t *rtcm_nav_src,
                          const obsd_t *obs_a, int na,
                          const obsd_t *obs_b, int nb,
                          double *sat_extra_m)
{
    if (!ni || !sat_extra_m) return 0;
    memset(sat_extra_m, 0, sizeof(double) * MAXSAT);

    if (!ni->enable || !vbs || !rtcm_nav_src) return 0;
    if (!vbs->base_valid[0] || !vbs->base_valid[1]) return 0;
    if (!obs_a || !obs_b || na <= 0 || nb <= 0) return 0;

    /* Fix 2: IDW weights replace the old 1-D baseline projection. */
    double wa, wb;
    idw_weights(vbs->vbs_ecef, vbs->base_ecef[0], vbs->base_ecef[1], &wa, &wb);

    /* Pre-compute LLH and troposphere inputs for Fix 1 & Fix 3. */
    double pos_a[3], pos_b[3];
    ecef2pos(vbs->base_ecef[0], pos_a);
    ecef2pos(vbs->base_ecef[1], pos_b);

    int n = na < nb ? na : nb;
    if (n > MAXOBS) n = MAXOBS;

    static obsd_t pair_a[MAXOBS];
    static obsd_t pair_b[MAXOBS];
    static double rs_a[6*MAXOBS], rs_b[6*MAXOBS];
    static double dts_a[2*MAXOBS], dts_b[2*MAXOBS];
    static double var_a[MAXOBS], var_b[MAXOBS];
    static int svh_a[MAXOBS], svh_b[MAXOBS];
    static double raw_dres[MAXOBS];
    static int raw_sys[MAXOBS];

    int np = 0;
    for (int i = 0; i < na && np < MAXOBS; i++) {
        const obsd_t *ob = find_sat(obs_b, nb, obs_a[i].sat);
        double pa, pb;
        if (!ob) continue;
        if (!first_common_p(&obs_a[i], ob, &pa, &pb)) continue;
        pair_a[np] = obs_a[i];
        pair_b[np] = *ob;
        np++;
    }
    if (np <= 0) return 0;

    satposs_with_fallback(pair_a[0].time, pair_a, np, &rtcm_nav_src->nav,
            vbs->sat_ephopt, rs_a, dts_a, var_a, svh_a);
    satposs_with_fallback(pair_b[0].time, pair_b, np, &rtcm_nav_src->nav,
            vbs->sat_ephopt, rs_b, dts_b, var_b, svh_b);

    for (int i = 0; i < np; i++) {
        raw_dres[i] = NAN;
        raw_sys[i] = 0;
    }

    for (int i = 0; i < np; i++) {
        const obsd_t *oa = &pair_a[i];
        const obsd_t *ob = &pair_b[i];
        double pa, pb;
        if (oa->sat < 1 || oa->sat > MAXSAT) continue;
        if (!usable_satpos(rs_a + 6*i, svh_a[i]) ||
            !usable_satpos(rs_b + 6*i, svh_b[i])) continue;
        if (!first_common_p(oa, ob, &pa, &pb)) continue;

        double ea[3], eb[3];
        double rho_a = geodist(rs_a + 6*i, vbs->base_ecef[0], ea);
        double rho_b = geodist(rs_b + 6*i, vbs->base_ecef[1], eb);
        if (rho_a <= 0.0 || rho_b <= 0.0) continue;

        double dres = (pb - pa) - (rho_b - rho_a);
        double azel_a[2] = {0.0, 0.0};
        double azel_b[2] = {0.0, 0.0};
        satazel(pos_a, ea, azel_a);
        satazel(pos_b, eb, azel_b);
        if (azel_a[1] < NET_MIN_ELEV || azel_b[1] < NET_MIN_ELEV) {
            ni->n_sat_reject++;
            continue;
        }
        if (vbs->apply_trop) {
            double ta = tropmodel(oa->time, pos_a, azel_a, vbs->humi);
            double tb = tropmodel(ob->time, pos_b, azel_b, vbs->humi);
            dres -= (tb - ta);
        }
        if (!isfinite(dres)) continue;

        raw_dres[i] = dres;
        raw_sys[i] = satsys(oa->sat, NULL);
        ni->n_sat_pairs++;
    }

    static const int syses[] = {SYS_GPS, SYS_GLO, SYS_GAL, SYS_CMP,
                                SYS_QZS, SYS_IRN, SYS_SBS};
    double sys_bias[sizeof(syses)/sizeof(syses[0])];
    int sys_nvalid[sizeof(syses)/sizeof(syses[0])];
    for (size_t s = 0; s < sizeof(syses)/sizeof(syses[0]); s++) {
        double vals[MAXOBS];
        int nv = 0;
        for (int i = 0; i < np; i++) {
            if (raw_sys[i] == syses[s] && isfinite(raw_dres[i])) {
                vals[nv++] = raw_dres[i];
            }
        }
        sys_nvalid[s] = nv;
        sys_bias[s] = nv >= 3 ? median_of(vals, nv) : 0.0;
    }

    int out = 0;
    for (int i = 0; i < np; i++) {
        const obsd_t *oa = &pair_a[i];
        if (!isfinite(raw_dres[i])) continue;

        double bias = 0.0;
        int nvalid = 0;
        for (size_t s = 0; s < sizeof(syses)/sizeof(syses[0]); s++) {
            if (raw_sys[i] == syses[s]) {
                bias = sys_bias[s];
                nvalid = sys_nvalid[s];
                break;
            }
        }
        if (nvalid < 1) continue;

        double ea[3], eb[3];
        double rho_a = geodist(rs_a + 6*i, vbs->base_ecef[0], ea);
        double rho_b = geodist(rs_b + 6*i, vbs->base_ecef[1], eb);
        if (rho_a <= 0.0 || rho_b <= 0.0) continue;

        double dres = raw_dres[i] - bias;

        /* dres has already had inter-base troposphere and the per-system
         * receiver clock/code common mode removed. The remaining term is the
         * spatial residual used for the lightweight ionosphere interpolation. */
        double azel_a[2] = {0.0, 0.0};
        double azel_b[2] = {0.0, 0.0};
        satazel(pos_a, ea, azel_a);
        satazel(pos_b, eb, azel_b);
        if (azel_a[1] < NET_MIN_ELEV || azel_b[1] < NET_MIN_ELEV) {
            ni->n_sat_reject++;
            continue;
        }

        if (!isfinite(dres) || fabs(dres) > ni->max_abs_corr_m) {
            ni->n_sat_reject++;
            continue;
        }

        /* Fix 3: normalise dres to the zenith (VTEC) domain using the mean
         * mapping function of the two bases, then re-apply the VBS-side mapping
         * function. This removes the bias introduced when the two bases observe
         * the same satellite at noticeably different elevations (worst case for
         * low-elevation satellites on long baselines). VBS mapping function is
         * approximated by the IDW-weighted average of the two bases. */
        double mf_a = ionmapf(pos_a, azel_a);
        double mf_b = ionmapf(pos_b, azel_b);
        double mf_mean = 0.5 * (mf_a + mf_b);

        double corr;
        if (mf_mean > 1e-3 && mf_a > 1e-3 && mf_b > 1e-3 &&
            mf_a < 5.0 && mf_b < 5.0) {
            double dres_vtec = dres / mf_mean;
            double mf_v = wa * mf_a + wb * mf_b;
            corr = wb * dres_vtec * mf_v;
        } else {
            /* Fix 2 only (fallback for extreme geometry): IDW without
             * ionosphere normalisation. */
            corr = wb * dres;
        }

        if (!isfinite(corr) || fabs(corr) < 1e-4) continue;
        sat_extra_m[oa->sat - 1] = corr;
        ni->n_sat_interp++;
        out++;
    }
    if (out > 0) ni->n_epochs_used++;
    return out;
}
