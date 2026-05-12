/*------------------------------------------------------------------------------
 * vbs_core.c : high-precision VBS correction engine
 *----------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "vbs_core.h"

#define VBS_MIN_ELEV (10.0 * D2R)

static int usable_satpos(const double *rs, int svh)
{
    return svh >= 0 && (rs[0] != 0.0 || rs[1] != 0.0 || rs[2] != 0.0);
}

static void satposs_with_fallback(gtime_t time, const obsd_t *obs, int n,
                                  const nav_t *nav, int ephopt,
                                  double *rs, double *dts, double *var,
                                  int *svh, int *used_ssr)
{
    static double rs_brdc[6*MAXOBS];
    static double dts_brdc[2*MAXOBS];
    static double var_brdc[MAXOBS];
    static int svh_brdc[MAXOBS];

    if (used_ssr) memset(used_ssr, 0, sizeof(int) * (n > MAXOBS ? MAXOBS : n));
    satposs(time, obs, n, nav, ephopt, rs, dts, var, svh);
    if (ephopt == EPHOPT_BRDC) return;

    satposs(time, obs, n, nav, EPHOPT_BRDC, rs_brdc, dts_brdc, var_brdc, svh_brdc);
    for (int i = 0; i < n && i < MAXOBS; i++) {
        if (usable_satpos(rs + 6*i, svh[i])) {
            if (used_ssr) used_ssr[i] = 1;
            continue;
        }
        if (!usable_satpos(rs_brdc + 6*i, svh_brdc[i])) continue;
        memcpy(rs + 6*i, rs_brdc + 6*i, sizeof(double) * 6);
        memcpy(dts + 2*i, dts_brdc + 2*i, sizeof(double) * 2);
        var[i] = var_brdc[i];
        svh[i] = svh_brdc[i];
    }
}

/* local NED -> ECEF offset -------------------------------------------------- */
static void apply_ned(const double *base_ecef, double dn, double de, double du,
                      double *out_ecef)
{
    double pos[3], denu[3] = {de, dn, du}, dxyz[3];
    ecef2pos(base_ecef, pos);
    enu2ecef(pos, denu, dxyz);
    out_ecef[0] = base_ecef[0] + dxyz[0];
    out_ecef[1] = base_ecef[1] + dxyz[1];
    out_ecef[2] = base_ecef[2] + dxyz[2];
}

void vbs_init(vbs_ctx_t *ctx,
              int use_llh_offset,
              double dn, double de, double du,
              const double *vbs_llh,
              int apply_trop)
{
    memset(ctx, 0, sizeof(*ctx));
    ctx->use_llh_offset = use_llh_offset;
    ctx->dn = dn; ctx->de = de; ctx->du = du;
    ctx->apply_trop = apply_trop;
    ctx->humi = 0.7;
    ctx->sat_ephopt = EPHOPT_BRDC;

    if (!use_llh_offset && vbs_llh) {
        double pos[3] = { vbs_llh[0]*D2R, vbs_llh[1]*D2R, vbs_llh[2] };
        pos2ecef(pos, ctx->vbs_ecef);
    }
}

void vbs_set_ephopt(vbs_ctx_t *ctx, int sat_ephopt)
{
    if (sat_ephopt == EPHOPT_SSRAPC || sat_ephopt == EPHOPT_SSRCOM ||
        sat_ephopt == EPHOPT_BRDC) {
        ctx->sat_ephopt = sat_ephopt;
    }
}

int vbs_update_base(vbs_ctx_t *ctx, int side, const sta_t *sta)
{
    if (side < 0 || side >= VBS_MAX_BASES) return 0;
    if (norm(sta->pos, 3) < 1e3) return 0;  /* not yet populated */

    double *be = ctx->base_ecef[side];
    int changed = (fabs(sta->pos[0]-be[0]) > 1e-3 ||
                   fabs(sta->pos[1]-be[1]) > 1e-3 ||
                   fabs(sta->pos[2]-be[2]) > 1e-3);

    be[0] = sta->pos[0];
    be[1] = sta->pos[1];
    be[2] = sta->pos[2];
    ctx->base_valid[side] = 1;

    double pos_llh[3];
    ecef2pos(be, pos_llh);
    if (changed) {
        fprintf(stderr,
            "[vbs] base %c ECEF = (%.3f, %.3f, %.3f)  LLH = (%.8f, %.8f, %.3f)\n",
            'A' + side, be[0], be[1], be[2],
            pos_llh[0]*R2D, pos_llh[1]*R2D, pos_llh[2]);
    }

    /* In offset mode, VBS is anchored to base A: recompute once A is known. */
    if (side == 0 && ctx->use_llh_offset &&
        (changed || norm(ctx->vbs_ecef, 3) < 1e3)) {
        apply_ned(be, ctx->dn, ctx->de, ctx->du, ctx->vbs_ecef);
        double vpos[3];
        ecef2pos(ctx->vbs_ecef, vpos);
        fprintf(stderr,
            "[vbs] VBS   ECEF = (%.3f, %.3f, %.3f)  LLH = (%.8f, %.8f, %.3f)\n",
            ctx->vbs_ecef[0], ctx->vbs_ecef[1], ctx->vbs_ecef[2],
            vpos[0]*R2D, vpos[1]*R2D, vpos[2]);
    }
    return changed;
}

void vbs_rewrite_station(vbs_ctx_t *ctx, rtcm_t *rtcm)
{
    if (norm(ctx->vbs_ecef, 3) < 1e3) return;
    rtcm->sta.pos[0] = ctx->vbs_ecef[0];
    rtcm->sta.pos[1] = ctx->vbs_ecef[1];
    rtcm->sta.pos[2] = ctx->vbs_ecef[2];
    ctx->n_frames_1005_6++;
}

/* differential troposphere (Saastamoinen): returns trop(vbs) - trop(base) */
static double delta_trop(const double *base_pos_llh, const double *azel_base,
                         const double *vbs_pos_llh, const double *azel_vbs,
                         double humi, gtime_t time)
{
    double tb = tropmodel(time, base_pos_llh, azel_base, humi);
    double tv = tropmodel(time, vbs_pos_llh,  azel_vbs,  humi);
    return tv - tb;
}

int vbs_correct_obs_ex(vbs_ctx_t *ctx, rtcm_t *rtcm,
                       const unsigned char *side_per_obs,
                       const double *sat_extra_m)
{
    if (norm(ctx->vbs_ecef, 3) < 1e3) {
        /* VBS not resolved yet: drop obs to avoid sending anything
           geometrically inconsistent. */
        rtcm->obs.n = 0;
        return 0;
    }

    int n = rtcm->obs.n;
    if (n <= 0) return 0;
    if (n > MAXOBS) n = MAXOBS;

    static double rs[6*MAXOBS];
    static double dts[2*MAXOBS];
    static double var[MAXOBS];
    static int    svh[MAXOBS];
    static int    used_ssr[MAXOBS];

    satposs_with_fallback(rtcm->obs.data[0].time, rtcm->obs.data, n,
            &rtcm->nav, ctx->sat_ephopt, rs, dts, var, svh, used_ssr);

    /* Pre-compute LLH for each real base we'll possibly use, for the trop
       differential. */
    double base_pos_llh[VBS_MAX_BASES][3];
    for (int s = 0; s < VBS_MAX_BASES; s++) {
        if (ctx->base_valid[s]) ecef2pos(ctx->base_ecef[s], base_pos_llh[s]);
    }
    double vbs_pos[3];
    ecef2pos(ctx->vbs_ecef, vbs_pos);

    int n_out = 0, n_noeph = 0;
    for (int i = 0; i < n; i++) {
        obsd_t *o = &rtcm->obs.data[i];
        double *rs_i = rs + 6*i;

        if (sat_extra_m && o->sat >= 1 && o->sat <= MAXSAT &&
            fabs(sat_extra_m[o->sat - 1]) > 0.0) {
            ctx->n_sat_netcorr_seen++;
        }

        if (!usable_satpos(rs_i, svh[i])) {
            n_noeph++;
            continue;
        }
        if (ctx->sat_ephopt != EPHOPT_BRDC) {
            if (used_ssr[i]) ctx->n_sat_ssr_used++;
            else ctx->n_sat_ssr_fallback++;
        }

        int side = side_per_obs ? side_per_obs[i] : 0;
        if (side < 0 || side >= VBS_MAX_BASES || !ctx->base_valid[side]) {
            n_noeph++;
            continue;
        }
        const double *base = ctx->base_ecef[side];

        double e_base[3], e_vbs[3];
        double r_base = geodist(rs_i, base,           e_base);
        double r_vbs  = geodist(rs_i, ctx->vbs_ecef,  e_vbs);
        if (r_base <= 0.0 || r_vbs <= 0.0) { n_noeph++; continue; }

        double azel_b[2], azel_v[2];
        satazel(base_pos_llh[side], e_base, azel_b);
        satazel(vbs_pos,             e_vbs,  azel_v);
        if (azel_b[1] < VBS_MIN_ELEV || azel_v[1] < VBS_MIN_ELEV) {
            n_noeph++;
            continue;
        }

        double drho = r_vbs - r_base;

        double dtrop = 0.0;
        if (ctx->apply_trop) {
            dtrop = delta_trop(base_pos_llh[side], azel_b, vbs_pos, azel_v,
                               ctx->humi, o->time);
        }

        double dnet = 0.0;
        if (sat_extra_m && o->sat >= 1 && o->sat <= MAXSAT) {
            dnet = sat_extra_m[o->sat - 1];
            if (fabs(dnet) > 0.0) ctx->n_sat_netcorr++;
        }

        for (int j = 0; j < NFREQ+NEXOBS; j++) {
            if (o->P[j] != 0.0) {
                o->P[j] += drho + dtrop + dnet;
            }
            if (o->L[j] != 0.0) {
                double freq = sat2freq(o->sat, o->code[j], &rtcm->nav);
                if (freq > 0.0) {
                    double lam = CLIGHT / freq;
                    o->L[j] += (drho + dtrop + dnet) / lam;
                }
            }
            /* Doppler unchanged: VBS is static wrt the real bases. */
        }

        if (n_out != i) rtcm->obs.data[n_out] = *o;
        n_out++;
    }

    rtcm->obs.n = n_out;
    ctx->n_obs_in        += n;
    ctx->n_obs_out       += n_out;
    ctx->n_sat_corrected += n_out;
    ctx->n_sat_no_eph    += n_noeph;
    ctx->n_frames_msm++;
    return n_out;
}

int vbs_correct_obs(vbs_ctx_t *ctx, rtcm_t *rtcm,
                    const unsigned char *side_per_obs)
{
    return vbs_correct_obs_ex(ctx, rtcm, side_per_obs, NULL);
}
