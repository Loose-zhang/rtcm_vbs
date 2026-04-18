/*------------------------------------------------------------------------------
 * vbs_core.c : high-precision VBS correction engine
 *----------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "vbs_core.h"

#define CLIGHT_MS   (CLIGHT*1e-3)   /* m per ms (unused, keep for clarity) */

/* local NED -> ECEF offset -------------------------------------------------- */
static void apply_ned(const double *base_ecef, double dn, double de, double du,
                      double *out_ecef)
{
    double pos[3], E[9], dxyz[3];
    double denu[3] = {de, dn, du};  /* RTKLIB enu2ecef expects E,N,U */

    ecef2pos(base_ecef, pos);       /* pos = {lat, lon, h} */

    /* RTKLIB provides enu2ecef which rotates ENU->ECEF at a given lat/lon */
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

    if (!use_llh_offset && vbs_llh) {
        double pos[3] = { vbs_llh[0]*D2R, vbs_llh[1]*D2R, vbs_llh[2] };
        pos2ecef(pos, ctx->vbs_ecef);
    }
}

void vbs_update_base_from_sta(vbs_ctx_t *ctx, const sta_t *sta)
{
    if (norm(sta->pos, 3) < 1e3) return;  /* not yet learned */

    int changed = (fabs(sta->pos[0]-ctx->base_ecef[0]) > 1e-3 ||
                   fabs(sta->pos[1]-ctx->base_ecef[1]) > 1e-3 ||
                   fabs(sta->pos[2]-ctx->base_ecef[2]) > 1e-3);

    ctx->base_ecef[0] = sta->pos[0];
    ctx->base_ecef[1] = sta->pos[1];
    ctx->base_ecef[2] = sta->pos[2];
    ctx->base_valid = 1;

    if (ctx->use_llh_offset && (changed || norm(ctx->vbs_ecef, 3) < 1e3)) {
        apply_ned(ctx->base_ecef, ctx->dn, ctx->de, ctx->du, ctx->vbs_ecef);

        double bpos[3], vpos[3];
        ecef2pos(ctx->base_ecef, bpos);
        ecef2pos(ctx->vbs_ecef,  vpos);
        fprintf(stderr,
            "[vbs] base ECEF = (%.3f, %.3f, %.3f)  LLH = (%.8f, %.8f, %.3f)\n",
            ctx->base_ecef[0], ctx->base_ecef[1], ctx->base_ecef[2],
            bpos[0]*R2D, bpos[1]*R2D, bpos[2]);
        fprintf(stderr,
            "[vbs] VBS  ECEF = (%.3f, %.3f, %.3f)  LLH = (%.8f, %.8f, %.3f)\n",
            ctx->vbs_ecef[0], ctx->vbs_ecef[1], ctx->vbs_ecef[2],
            vpos[0]*R2D, vpos[1]*R2D, vpos[2]);
    }
}

void vbs_rewrite_station(vbs_ctx_t *ctx, rtcm_t *rtcm)
{
    if (!ctx->base_valid || norm(ctx->vbs_ecef, 3) < 1e3) return;
    rtcm->sta.pos[0] = ctx->vbs_ecef[0];
    rtcm->sta.pos[1] = ctx->vbs_ecef[1];
    rtcm->sta.pos[2] = ctx->vbs_ecef[2];
    /* deltas/hgt for 1006 stay as-is (ARP offset at the station); they are
       interpreted relative to the reported position. */
    ctx->n_frames_1005_6++;
}

/* differential troposphere (Saastamoinen): returns trop(vbs) - trop(base) ---- */
static double delta_trop(const double *base_pos_llh,
                         const double *vbs_pos_llh,
                         const double *azel,
                         double humi, gtime_t time)
{
    double tb = tropmodel(time, base_pos_llh, azel, humi);
    double tv = tropmodel(time, vbs_pos_llh,  azel, humi);
    return tv - tb;
}

int vbs_correct_obs(vbs_ctx_t *ctx, rtcm_t *rtcm)
{
    if (!ctx->base_valid || norm(ctx->vbs_ecef, 3) < 1e3) {
        /* Can't correct yet: drop obs to avoid geometric inconsistency */
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

    /* satellite positions at signal transmission time, with all the subtle
       corrections (Sagnac-unaware; Sagnac is added in geodist()) */
    satposs(rtcm->obs.data[0].time, rtcm->obs.data, n, &rtcm->nav,
            EPHOPT_BRDC, rs, dts, var, svh);

    double base_pos[3], vbs_pos[3];
    ecef2pos(ctx->base_ecef, base_pos);
    ecef2pos(ctx->vbs_ecef,  vbs_pos);

    int n_out = 0, n_noeph = 0;
    for (int i = 0; i < n; i++) {
        obsd_t *o = &rtcm->obs.data[i];
        double *rs_i = rs + 6*i;

        if (svh[i] < 0 || (rs_i[0]==0.0 && rs_i[1]==0.0 && rs_i[2]==0.0)) {
            n_noeph++;
            continue;   /* drop this sat */
        }

        /* ranges with Sagnac correction */
        double e_base[3], e_vbs[3];
        double r_base = geodist(rs_i, ctx->base_ecef, e_base);
        double r_vbs  = geodist(rs_i, ctx->vbs_ecef,  e_vbs);
        if (r_base <= 0.0 || r_vbs <= 0.0) { n_noeph++; continue; }

        /* geometric range correction (m) */
        double drho = r_vbs - r_base;

        /* optional differential troposphere (uses VBS elevation) */
        double dtrop = 0.0;
        if (ctx->apply_trop) {
            double azel_b[2], azel_v[2];
            satazel(base_pos, e_base, azel_b);
            satazel(vbs_pos,  e_vbs,  azel_v);
            /* use average elevation for a single mapping */
            double azel_m[2] = {
                0.5*(azel_b[0]+azel_v[0]),
                0.5*(azel_b[1]+azel_v[1])
            };
            dtrop = delta_trop(base_pos, vbs_pos, azel_m, ctx->humi,
                               o->time);
        }

        /* apply to every frequency channel of this satellite */
        for (int j = 0; j < NFREQ+NEXOBS; j++) {
            if (o->P[j] != 0.0) {
                o->P[j] += drho + dtrop;
            }
            if (o->L[j] != 0.0) {
                double freq = sat2freq(o->sat, o->code[j], &rtcm->nav);
                if (freq > 0.0) {
                    double lam = CLIGHT / freq;
                    o->L[j] += (drho + dtrop) / lam;
                }
            }
            /* Doppler unchanged: relative velocity of VBS vs base is 0. */
        }

        /* compact array */
        if (n_out != i) rtcm->obs.data[n_out] = *o;
        n_out++;
    }

    rtcm->obs.n = n_out;
    ctx->n_obs_in  += n;
    ctx->n_obs_out += n_out;
    ctx->n_sat_corrected += n_out;
    ctx->n_sat_no_eph    += n_noeph;
    ctx->n_frames_msm++;
    return n_out;
}
