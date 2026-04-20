/*------------------------------------------------------------------------------
 * net_interp.c : lightweight dual-base network/ionosphere interpolation
 *----------------------------------------------------------------------------*/
#include <string.h>
#include <math.h>
#include "net_interp.h"

static const obsd_t *find_sat(const obsd_t *obs, int n, int sat)
{
    for (int i = 0; i < n; i++) {
        if (obs[i].sat == sat) return &obs[i];
    }
    return NULL;
}

static double first_valid_p(const obsd_t *o)
{
    for (int j = 0; j < NFREQ + NEXOBS; j++) {
        if (o->P[j] != 0.0) return o->P[j];
    }
    return 0.0;
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

    int np = 0;
    for (int i = 0; i < na && np < MAXOBS; i++) {
        const obsd_t *ob = find_sat(obs_b, nb, obs_a[i].sat);
        if (!ob) continue;
        if (first_valid_p(&obs_a[i]) == 0.0 || first_valid_p(ob) == 0.0) continue;
        pair_a[np] = obs_a[i];
        pair_b[np] = *ob;
        np++;
    }
    if (np <= 0) return 0;

    satposs(pair_a[0].time, pair_a, np, (nav_t *)&rtcm_nav_src->nav,
            EPHOPT_BRDC, rs_a, dts_a, var_a, svh_a);
    satposs(pair_b[0].time, pair_b, np, (nav_t *)&rtcm_nav_src->nav,
            EPHOPT_BRDC, rs_b, dts_b, var_b, svh_b);

    int out = 0;
    for (int i = 0; i < np; i++) {
        const obsd_t *oa = &pair_a[i];
        const obsd_t *ob = &pair_b[i];
        if (oa->sat < 1 || oa->sat > MAXSAT) continue;
        if (svh_a[i] < 0 || svh_b[i] < 0) continue;

        double pa = first_valid_p(oa);
        double pb = first_valid_p(ob);
        if (pa == 0.0 || pb == 0.0) continue;

        double ea[3], eb[3];
        double rho_a = geodist(rs_a + 6*i, vbs->base_ecef[0], ea);
        double rho_b = geodist(rs_b + 6*i, vbs->base_ecef[1], eb);
        if (rho_a <= 0.0 || rho_b <= 0.0) continue;

        double dgeom = rho_b - rho_a;
        double dres  = (pb - pa) - dgeom;

        /* Fix 1: when trop correction is also active, vbs_core already applies
         * dtrop = trop(vbs) - trop(base) for each satellite. The dres computed
         * above still contains the inter-base troposphere difference
         * (trop_B - trop_A). Subtract it here so the residual passed to the
         * interpolator is dominated by ionosphere only, preventing the trop
         * component from being corrected twice. */
        double azel_a[2] = {0.0, 0.0};
        double azel_b[2] = {0.0, 0.0};
        satazel(pos_a, ea, azel_a);
        satazel(pos_b, eb, azel_b);
        if (vbs->apply_trop) {
            double ta = tropmodel(oa->time, pos_a, azel_a, vbs->humi);
            double tb = tropmodel(ob->time, pos_b, azel_b, vbs->humi);
            dres -= (tb - ta);
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

        sat_extra_m[oa->sat - 1] = corr;
        ni->n_sat_pairs++;
        ni->n_sat_interp++;
        out++;
    }
    if (out > 0) ni->n_epochs_used++;
    return out;
}
