/*------------------------------------------------------------------------------
 * net_interp.h : lightweight dual-base network/ionosphere interpolation
 *
 * First-stage model for long-baseline VBS enhancement:
 *   - Works with the existing dual-base A/B architecture.
 *   - Builds a per-satellite residual difference between base A and base B.
 *   - Interpolates that residual to the VBS position along the A-B baseline.
 *
 * The residual is formed in measurement space after subtracting pure geometry:
 *
 *   res = (P_B - P_A) - (rho_B - rho_A)
 *
 * so the remaining term is dominated by spatially varying atmosphere,
 * especially ionosphere over longer baselines. The interpolated residual is
 * then added on top of the existing geometric/trop VBS correction.
 *----------------------------------------------------------------------------*/
#ifndef NET_INTERP_H
#define NET_INTERP_H

#include "rtklib.h"
#include "vbs_core.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int    enable;
    double max_abs_corr_m;
    long   n_epochs_used;
    long   n_sat_pairs;
    long   n_sat_interp;
    long   n_sat_reject;
} net_interp_t;

void net_interp_init(net_interp_t *ni, int enable);

/* Build additive per-satellite corrections in metres, indexed by sat-1.
 * Returns number of satellites for which a correction was produced. */
int net_interp_build_dual(net_interp_t *ni,
                          const vbs_ctx_t *vbs,
                          const rtcm_t *rtcm_nav_src,
                          const obsd_t *obs_a, int na,
                          const obsd_t *obs_b, int nb,
                          double *sat_extra_m /* size MAXSAT */);

#ifdef __cplusplus
}
#endif
#endif /* NET_INTERP_H */
