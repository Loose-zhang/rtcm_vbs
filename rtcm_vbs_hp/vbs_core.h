/*------------------------------------------------------------------------------
 * vbs_core.h : High-precision Virtual Base Station core corrections
 *
 * Depends on RTKLIB (include/rtklib.h).
 *
 * The VBS correction replaces the reference station's observations with those
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
 * Optionally a differential troposphere correction (Saastamoinen) can be
 * applied when the VBS vertical offset is significant.
 *----------------------------------------------------------------------------*/
#ifndef VBS_CORE_H
#define VBS_CORE_H

#include "rtklib.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double base_ecef[3];       /* real base station ECEF (m), updated from 1005/1006 */
    int    base_valid;         /* true when base_ecef has been set */
    double vbs_ecef[3];        /* virtual base station target ECEF (m) */

    /* Offset mode: if use_llh_offset != 0, vbs_ecef is derived each time the
       real-base ECEF is (re)learned by applying (dn,de,du) in the local
       geodetic frame of the real base. Otherwise vbs_ecef is fixed. */
    int    use_llh_offset;
    double dn, de, du;         /* North/East/Up offset (m) */

    int    apply_trop;         /* 1: differential Saastamoinen troposphere */
    double humi;               /* relative humidity 0..1 for Saastamoinen */

    /* Statistics */
    long   n_obs_in, n_obs_out;
    long   n_sat_corrected, n_sat_no_eph;
    long   n_frames_1005_6, n_frames_msm;
} vbs_ctx_t;

/* Initialize context. If use_llh_offset=1, (dn,de,du) are used to derive the
 * VBS from the *first* learned base position. If use_llh_offset=0, vbs_llh
 * is used as the absolute VBS coordinate (lat deg, lon deg, h m). */
void vbs_init(vbs_ctx_t *ctx,
              int use_llh_offset,
              double dn, double de, double du,
              const double *vbs_llh /* may be NULL if use_llh_offset=1 */,
              int apply_trop);

/* Called after input_rtcm3 returns 5 (station params). Updates base_ecef and,
 * if offset mode is active and VBS not yet resolved, computes vbs_ecef. */
void vbs_update_base_from_sta(vbs_ctx_t *ctx, const sta_t *sta);

/* Rewrite sta->pos in the rtcm struct to the VBS coordinate so subsequent
 * gen_rtcm3(1005/1006) encodes the virtual position. */
void vbs_rewrite_station(vbs_ctx_t *ctx, rtcm_t *rtcm);

/* Apply VBS correction to rtcm->obs in-place.
 * Returns number of satellites successfully corrected. Satellites without
 * usable ephemeris are removed from rtcm->obs so downstream encoding only
 * contains geometrically-consistent measurements. */
int  vbs_correct_obs(vbs_ctx_t *ctx, rtcm_t *rtcm);

#ifdef __cplusplus
}
#endif
#endif /* VBS_CORE_H */
