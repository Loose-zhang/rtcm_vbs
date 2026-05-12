/*------------------------------------------------------------------------------
 * vbs_core.h : High-precision Virtual Base Station core corrections
 *
 * Depends on RTKLIB (include/rtklib.h).
 *
 * The VBS correction replaces a reference station's observations with those
 * that *would have been* measured at a different (virtual) location, by using
 * the geometric range difference between the real base and the virtual base
 * for every satellite in view.
 *
 *   P_vbs = P_base  +  [ |r_sat - r_vbs|  -  |r_sat - r_base| ]
 *   L_vbs = L_base  +  [ |r_sat - r_vbs|  -  |r_sat - r_base| ] / lambda
 *
 * Satellite positions are computed at signal-transmission time via
 * RTKLIB's satposs() which already accounts for:
 *   - signal propagation time iteration (t - rho/c)
 *   - satellite clock bias (af0 + af1*dt + af2*dt^2)
 *   - relativistic clock correction (-2 sqrt(mu a) e sinE / c^2)
 *   - BDS GEO orientation special case
 *
 * Ranges use RTKLIB's geodist() which includes the Sagnac (Earth-rotation)
 * correction, removing the ~30 m systematic bias present in the old Python
 * implementation.
 *
 * Dual-base correctness:
 *   The context can hold TWO real-base positions (A and B). In single-base
 *   mode only index 0 (A) is used. In dual-base mode each satellite's
 *   observation is corrected from whichever of A or B actually measured it
 *   (the merger tags this per satellite). This eliminates the error that
 *   arises when |r_base_A - r_base_B| is large (km scale) and a sat's obs
 *   is from B but gets corrected as if it came from A.
 *----------------------------------------------------------------------------*/
#ifndef VBS_CORE_H
#define VBS_CORE_H

#include "rtklib.h"

#ifdef __cplusplus
extern "C" {
#endif

#define VBS_MAX_BASES 2  /* 0 = A, 1 = B */

typedef struct {
    double base_ecef[VBS_MAX_BASES][3];  /* real base ECEF(s) in metres */
    int    base_valid[VBS_MAX_BASES];    /* per-base learned flag */

    double vbs_ecef[3];                  /* virtual base target ECEF (m) */

    /* Offset mode: if use_llh_offset != 0, vbs_ecef is derived from the
       *A* base's first learned position. Otherwise vbs_ecef is an absolute
       coordinate set at init time. */
    int    use_llh_offset;
    double dn, de, du;                   /* North/East/Up offset (m) */

    int    apply_trop;
    double humi;                         /* relative humidity 0..1 */
    int    sat_ephopt;                   /* EPHOPT_BRDC or EPHOPT_SSR??? */

    /* Statistics */
    long   n_obs_in, n_obs_out;
    long   n_sat_corrected, n_sat_no_eph;
    long   n_frames_1005_6, n_frames_msm;
    long   n_sat_netcorr, n_sat_netcorr_seen;
    long   n_sat_ssr_used, n_sat_ssr_fallback;
} vbs_ctx_t;

/* ---- lifecycle ---- */
void vbs_init(vbs_ctx_t *ctx,
              int use_llh_offset,
              double dn, double de, double du,
              const double *vbs_llh /* NULL if offset mode */,
              int apply_trop);

/* Update base[side] from a just-decoded sta_t (1005/1006 payload).
 * side = 0 (A) or 1 (B). Returns 1 if the base ecef changed. Side 0 also
 * drives VBS recomputation in offset mode. */
int  vbs_update_base(vbs_ctx_t *ctx, int side, const sta_t *sta);

/* Rewrite rtcm->sta.pos to the VBS coordinate so a subsequent
 * gen_rtcm3(1005/1006) encodes the virtual ARP. */
void vbs_rewrite_station(vbs_ctx_t *ctx, rtcm_t *rtcm);

/* Select satellite position source for subsequent VBS corrections.
 * If an SSR option is selected, corrections fall back per-satellite to BRDC
 * whenever orbit/clock SSR is unavailable or stale. */
void vbs_set_ephopt(vbs_ctx_t *ctx, int sat_ephopt);

/* Apply VBS correction to rtcm->obs in-place.
 *
 *   side_per_obs:
 *     NULL        -> all obs use base 0 (A).  Single-base path.
 *     non-NULL    -> array of length rtcm->obs.n; each element is 0 or 1
 *                    indicating which real base produced that obsd_t.
 *
 *   sat_extra_m:
 *     NULL        -> no extra correction.
 *     non-NULL    -> additive per-satellite correction in metres, indexed by
 *                    sat-1. Intended for network/ionosphere interpolation.
 *
 * Satellites without usable ephemeris, or whose tagged base has not yet
 * been learned, are dropped from rtcm->obs so downstream encoding only
 * contains geometrically consistent measurements.
 *
 * Returns the number of satellites successfully corrected. */
int  vbs_correct_obs_ex(vbs_ctx_t *ctx, rtcm_t *rtcm,
                        const unsigned char *side_per_obs,
                        const double *sat_extra_m);

/* Backward-compatible wrapper without extra network correction. */
int  vbs_correct_obs(vbs_ctx_t *ctx, rtcm_t *rtcm,
                     const unsigned char *side_per_obs);

#ifdef __cplusplus
}
#endif
#endif /* VBS_CORE_H */
