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

static double baseline_fraction(const double *a, const double *b, const double *v)
{
    double ab[3], av[3];
    for (int i = 0; i < 3; i++) {
        ab[i] = b[i] - a[i];
        av[i] = v[i] - a[i];
    }
    double den = dot(ab, ab, 3);
    if (den <= 1e-6) return 0.0;
    return dot(av, ab, 3) / den;
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

    double alpha = baseline_fraction(vbs->base_ecef[0], vbs->base_ecef[1],
                                     vbs->vbs_ecef);
    if (alpha < -1.0) alpha = -1.0;
    if (alpha >  2.0) alpha =  2.0;

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

        if (!isfinite(dres) || fabs(dres) > ni->max_abs_corr_m) {
            ni->n_sat_reject++;
            continue;
        }

        double corr = alpha * dres;
        sat_extra_m[oa->sat - 1] = corr;
        ni->n_sat_pairs++;
        ni->n_sat_interp++;
        out++;
    }
    if (out > 0) ni->n_epochs_used++;
    return out;
}
