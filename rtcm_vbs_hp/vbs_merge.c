/*------------------------------------------------------------------------------
 * vbs_merge.c : dual-base merger with per-(sat,code) carrier-phase alignment
 *----------------------------------------------------------------------------*/
#include "vbs_merge.h"
#include <string.h>
#include <math.h>

/* Tunables ----------------------------------------------------------------- */
#define INIT_TOL_CYC      0.20   /* per-sample stability during init     */
#define UPDATE_TOL_CYC    0.30   /* IIR-update gating tolerance          */
#define IIR_ALPHA         0.10   /* bias = (1-a)*bias + a*diff  (faster) */
#define INIT_NEEDED       5      /* consecutive stable samples to lock   */
#define RESET_RESID_CYC   3.00   /* abs(diff - bias) larger than this    */
#define RESET_GAP_S       10.0   /* pair gap larger than this resets     */
#define PROP_DT_THRESH_S  0.001  /* skip Doppler propagation below 1 ms  */
#define PROP_DT_MAX_S     0.300  /* refuse propagation beyond 300 ms     */

/* Propagate one B obs epoch from its own timestamp to target time t_ref.
 *
 * For each signal slot with valid Doppler:
 *   L(t_ref) = L(t_B) + D * dt          [cycles, D in Hz = cycles/s]
 *   P(t_ref) = P(t_B) - D*(c/f) * dt    [metres]
 *
 * RTKLIB Doppler sign: D>0 means satellite is approaching (decreasing range),
 * so carrier phase increases and pseudorange decreases with time.
 *
 * If Doppler is zero (MSM4 omits Doppler; some receivers leave D=0) the
 * rigorous L/P step cannot be applied.  We keep the measured L and P
 * (same as leaving P unchanged in older versions) so small |dt| from
 * A/B epoch skew does not wipe all B-side carrier, which would destroy
 * dual-base AR.  n_b_prop_nodopp counts these no-Doppler propagation skips.
 *
 * The obs time field is set to t_ref after propagation. */
static void propagate_b_to_ref(merger_t *m, obsd_t *ob, gtime_t t_ref)
{
    double dt = timediff(t_ref, ob->time);  /* t_ref - t_B, seconds */
    ob->time = t_ref;

    if (fabs(dt) < PROP_DT_THRESH_S) return;  /* negligible: skip */

    if (fabs(dt) > PROP_DT_MAX_S) {
        /* Gap too large: suppress all carrier signals to avoid
         * false cycle-slip storms after long timing discontinuities. */
        for (int j = 0; j < NFREQ+NEXOBS; j++) {
            if (ob->L[j] != 0.0) {
                ob->L[j]  = 0.0;
                ob->LLI[j] = 0;
            }
        }
        return;
    }

    int sys = satsys(ob->sat, NULL);
    /* GLO FCN is stored as freq+7 (0-13, offset 7 → FCN -7..+6) */
    int fcn = (sys == SYS_GLO) ? (int)ob->freq - 7 : 0;

    for (int j = 0; j < NFREQ+NEXOBS; j++) {
        unsigned char code = ob->code[j];
        if (!code) continue;

        double freq = code2freq(sys, code, fcn);
        if (freq <= 0.0) continue;

        double D = (double)ob->D[j];  /* Hz (cycles/s) */

        if (D == 0.0) {
            /* No Doppler: cannot apply the L/P differential correction.
             * Keeping L and P avoids wiping B-side phase when streams are
             * MSM4 or D is missing; inter-base time skew is bounded by
             * pair_tol (default 50 ms) so residual phase error is usually
             * far smaller than losing carrier entirely. */
            m->n_b_prop_nodopp++;
            continue;
        }

        if (ob->L[j] != 0.0)
            ob->L[j] += D * dt;
        if (ob->P[j] != 0.0)
            ob->P[j] -= D * (CLIGHT / freq) * dt;
    }
}

void merger_init(merger_t *m, int window_ms, int pair_tol_ms)
{
    memset(m, 0, sizeof(*m));
    m->window_ms = window_ms > 0 ? window_ms : 40;
    if (pair_tol_ms < 1) pair_tol_ms = 50;
    if (pair_tol_ms > 250) pair_tol_ms = 250;
    m->pair_tol_s = pair_tol_ms / 1000.0;
}

void merger_stash(merger_t *m, int side,
                  const obsd_t *obs, int n, gtime_t t, long now_ms)
{
    epoch_buf_t *buf = side == 0 ? &m->a : &m->b;
    if (n > MAXOBS) n = MAXOBS;
    memcpy(buf->data, obs, n * sizeof(obsd_t));
    buf->n       = n;
    buf->time    = t;
    buf->valid   = 1;
    buf->ms_recv = now_ms;
}

static int copy_buf(const epoch_buf_t *src, obsd_t *out)
{
    int n = src->n > MAXOBS ? MAXOBS : src->n;
    memcpy(out, src->data, n * sizeof(obsd_t));
    return n;
}

int merger_poll_pair(merger_t *m, long now_ms,
                     obsd_t *out_a, int *na, int *have_a,
                     obsd_t *out_b, int *nb, int *have_b,
                     gtime_t *out_time)
{
    *na = *nb = 0;
    *have_a = *have_b = 0;

    if (m->a.valid && m->b.valid) {
        double dt = fabs(timediff(m->a.time, m->b.time));
        if (dt <= m->pair_tol_s) {
            *na = copy_buf(&m->a, out_a); *have_a = 1;
            *nb = copy_buf(&m->b, out_b); *have_b = 1;
            *out_time = m->a.time;
            m->n_merged++;
            m->a.valid = m->b.valid = 0;
            return 1;
        }
        /* Different epoch times: flush the older side as single. */
        if (timediff(m->a.time, m->b.time) < 0) {
            *na = copy_buf(&m->a, out_a); *have_a = 1;
            *out_time = m->a.time;
            m->n_a_only++;
            m->a.valid = 0;
            return 1;
        }
        *nb = copy_buf(&m->b, out_b); *have_b = 1;
        *out_time = m->b.time;
        m->n_b_only++;
        m->b.valid = 0;
        return 1;
    }

    if (m->a.valid && !m->b.valid && now_ms - m->a.ms_recv >= m->window_ms) {
        *na = copy_buf(&m->a, out_a); *have_a = 1;
        *out_time = m->a.time;
        m->n_a_only++;
        m->a.valid = 0;
        return 1;
    }
    if (m->b.valid && !m->a.valid && now_ms - m->b.ms_recv >= m->window_ms) {
        *nb = copy_buf(&m->b, out_b); *have_b = 1;
        *out_time = m->b.time;
        m->n_b_only++;
        m->b.valid = 0;
        return 1;
    }
    return 0;
}

/* Find the slot index in obs[] for the same (sat,code). Returns -1 if not. */
static int find_slot(const obsd_t *o, unsigned char code)
{
    if (!code) return -1;
    for (int j = 0; j < NFREQ+NEXOBS; j++) {
        if (o->code[j] == code && o->L[j] != 0.0) return j;
    }
    return -1;
}

static void reset_align(align_state_t *st, unsigned char code)
{
    st->valid      = 0;
    st->code       = code;
    st->init_count = 0;
    st->init_sum   = 0.0;
    st->bias_cyc   = 0.0;
    st->last_diff  = 0.0;
}

/* Update one (sat,code) state from the current pair. Returns 1 if state is
 * valid for normalisation after this update; else 0. */
static int update_one(merger_t *m, int sat_idx, int slot_a,
                      const obsd_t *oa, const obsd_t *ob, gtime_t t)
{
    int code = oa->code[slot_a];
    int slot_b = find_slot(ob, (unsigned char)code);
    if (slot_b < 0) {
        /* B doesn't have this signal; we can still normalise B's other
         * signals independently, but this slot has no joint info. */
        return 0;
    }

    align_state_t *st = &m->align[sat_idx][slot_a];

    /* If we are tracking a different code on this slot, restart. */
    if (st->code != (unsigned char)code) {
        reset_align(st, (unsigned char)code);
        m->n_align_reset_code++;
    }

    /* Cycle-slip on either side -> reset. */
    if ((oa->LLI[slot_a] & LLI_SLIP) || (ob->LLI[slot_b] & LLI_SLIP)) {
        reset_align(st, (unsigned char)code);
        m->n_align_reset_slip++;
        return 0;
    }

    /* Long gap -> reset. */
    if (st->last_pair_t.time != 0 &&
        fabs(timediff(t, st->last_pair_t)) > RESET_GAP_S) {
        reset_align(st, (unsigned char)code);
        m->n_align_reset_gap++;
    }

    double diff = oa->L[slot_a] - ob->L[slot_b];
    st->last_diff = diff;

    if (!st->valid) {
        if (st->init_count == 0) {
            st->init_sum   = diff;
            st->init_count = 1;
        } else {
            double mean = st->init_sum / st->init_count;
            if (fabs(diff - mean) <= INIT_TOL_CYC) {
                st->init_sum   += diff;
                st->init_count += 1;
                if (st->init_count >= INIT_NEEDED) {
                    st->bias_cyc = st->init_sum / st->init_count;
                    st->valid    = 1;
                    m->n_align_init_ok++;
                }
            } else {
                /* Restart init around the new sample. */
                st->init_sum   = diff;
                st->init_count = 1;
            }
        }
    } else {
        double resid = diff - st->bias_cyc;
        if (fabs(resid) > RESET_RESID_CYC) {
            reset_align(st, (unsigned char)code);
            m->n_align_reset_resid++;
        } else if (fabs(resid) <= UPDATE_TOL_CYC) {
            st->bias_cyc += IIR_ALPHA * resid;
        }
    }

    st->last_pair_t = t;
    return st->valid ? 1 : 0;
}

/* Apply per-(sat,slot) bias to B's L[] when available.
 *
 * Policy for B-only signals:
 *   1) If a valid B->A bias exists, normalise into A frame.
 *   2) If this (sat,code) has been seen in A/B overlap before but bias is
 *      currently invalid, suppress L/D to avoid injecting unstable phase.
 *   3) If this (sat,code) has never had A/B overlap history, keep original
 *      B-side L/D so independent-constellation fixing (e.g. BDS-only on B)
 *      remains possible.
 */
static int normalise_b(merger_t *m, obsd_t *ob)
{
    int sat_idx = ob->sat - 1;
    if (sat_idx < 0 || sat_idx >= MAXSAT) return 0;
    int n_norm = 0;

    for (int j = 0; j < NFREQ+NEXOBS; j++) {
        if (ob->L[j] == 0.0) continue;
        unsigned char code = ob->code[j];
        if (!code) continue;

        align_state_t *match = NULL;
        int has_history = 0;
        for (int k = 0; k < NFREQ+NEXOBS; k++) {
            align_state_t *st = &m->align[sat_idx][k];
            if (st->code != code) continue;
            has_history = 1;
            if (!match) match = st;
            if (st->valid) { match = st; break; }
        }

        if (match && match->valid) {
            ob->L[j] += match->bias_cyc;
            n_norm++;
            m->n_b_norm_applied++;
        } else if (has_history) {
            /* Alignment known but not yet valid: suppress L silently.
             * RTK tolerates missing phase far better than repeated
             * LLI_SLIP which forces ambiguity resets every epoch. */
            ob->L[j]  = 0.0;
            ob->LLI[j] = 0;
            m->n_b_unaligned_lli++;
        }
    }
    return n_norm;
}

int merger_align_and_merge(merger_t *m,
                           const obsd_t *obs_a, int na, int have_a,
                           const obsd_t *obs_b, int nb, int have_b,
                           gtime_t pair_time,
                           obsd_t *out_obs)
{
    int out_n = 0;

    /* Phase 1: propagate B obs to A epoch, then update per-(sat,code) bias
     * state for every satellite present on both sides.
     *
     * We work on a propagated copy of obs_b so that the bias estimate is
     * computed between two obs sets that refer to the same epoch (pair_time),
     * matching what will actually be emitted in the output. */
    obsd_t obs_b_prop[MAXOBS];
    int    nb_prop = 0;
    if (have_b) {
        nb_prop = nb > MAXOBS ? MAXOBS : nb;
        memcpy(obs_b_prop, obs_b, nb_prop * sizeof(obsd_t));
        for (int i = 0; i < nb_prop; i++)
            propagate_b_to_ref(m, &obs_b_prop[i], pair_time);
    }

    if (have_a && have_b) {
        for (int i = 0; i < na; i++) {
            int sat = obs_a[i].sat;
            int sat_idx = sat - 1;
            if (sat_idx < 0 || sat_idx >= MAXSAT) continue;
            const obsd_t *ob = NULL;
            for (int j = 0; j < nb_prop; j++) {
                if (obs_b_prop[j].sat == sat) { ob = &obs_b_prop[j]; break; }
            }
            if (!ob) continue;
            for (int slot_a = 0; slot_a < NFREQ+NEXOBS; slot_a++) {
                if (obs_a[i].L[slot_a] == 0.0) continue;
                if (!obs_a[i].code[slot_a]) continue;
                update_one(m, sat_idx, slot_a, &obs_a[i], ob, pair_time);
            }
        }
    }

    /* Phase 2: assemble output. */
    if (have_a && !have_b) {
        for (int i = 0; i < na && out_n < MAXOBS; i++) {
            out_obs[out_n++] = obs_a[i];
        }
        return out_n;
    }
    if (have_b && !have_a) {
        /* B-only epoch: obs already propagated to pair_time in obs_b_prop. */
        for (int i = 0; i < nb_prop && out_n < MAXOBS; i++) {
            out_obs[out_n] = obs_b_prop[i];
            normalise_b(m, &out_obs[out_n]);
            out_n++;
        }
        return out_n;
    }

    /* Both present: A sats first (A is primary), then B-only supplemental sats.
     * B supplemental sats use the propagated copy (already at pair_time). */
    unsigned char used_b[MAXOBS] = {0};
    for (int i = 0; i < na && out_n < MAXOBS; i++) {
        int sat = obs_a[i].sat;
        int bi  = -1;
        for (int j = 0; j < nb_prop; j++) {
            if (obs_b_prop[j].sat == sat) { bi = j; break; }
        }
        if (bi < 0) {
            out_obs[out_n++] = obs_a[i];
            continue;
        }
        used_b[bi] = 1;
        m->n_both++;

        /* Shared satellite: always keep A. */
        out_obs[out_n] = obs_a[i];
        m->n_chose_a++;
        out_n++;
    }
    /* B-only supplemental sats: propagated to pair_time, then normalise. */
    for (int j = 0; j < nb_prop && out_n < MAXOBS; j++) {
        if (used_b[j]) continue;
        out_obs[out_n] = obs_b_prop[j];
        normalise_b(m, &out_obs[out_n]);
        m->n_chose_b++;
        out_n++;
    }
    return out_n;
}
