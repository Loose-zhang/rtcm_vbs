/*------------------------------------------------------------------------------
 * net_interp.h : lightweight dual-base network/ionosphere interpolation
 *
 * First-stage model for long-baseline VBS enhancement:
 *   - Works with the existing dual-base A/B architecture.
 *   - Builds a per-satellite residual difference between base A and base B.
 *   - Interpolates that residual to the VBS position with inverse-distance
 *     weights from VBS to A/B. When geometry is usable, the residual is first
 *     normalised by an ionosphere mapping function and then mapped back to the
 *     VBS direction.
 *
 * The residual is formed in measurement space after subtracting pure geometry:
 *
 *   res = (P_B - P_A) - (rho_B - rho_A)
 *
 * so the remaining term is dominated by spatially varying atmosphere,
 * especially ionosphere over longer baselines. If troposphere correction is
 * also active, the inter-base troposphere component is removed before
 * interpolation so it is not applied twice. The residual field is then
 * interpolated to the VBS position and converted into side-specific additive
 * corrections:
 *
 *   extra_A = residual(VBS) - residual(A)
 *   extra_B = residual(VBS) - residual(B)
 *
 * This lets A and B observations be corrected into one virtual-reference
 * measurement space instead of treating B only as a diagnostic source.
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
 * sat_extra_a is applied to A-side observations, sat_extra_b to B-side
 * observations. Either output pointer may be NULL if that side is not needed.
 * Returns number of satellites for which a correction was produced. */
int net_interp_build_dual(net_interp_t *ni,
                          const vbs_ctx_t *vbs,
                          const rtcm_t *rtcm_nav_src,
                          const obsd_t *obs_a, int na,
                          const obsd_t *obs_b, int nb,
                          double *sat_extra_a /* size MAXSAT */,
                          double *sat_extra_b /* size MAXSAT */);

#ifdef __cplusplus
}
#endif
#endif /* NET_INTERP_H */
