/*------------------------------------------------------------------------------
 * vbs_merge.c : dual-base merger with per-(sat,code) carrier-phase alignment
 *----------------------------------------------------------------------------*/
#include "vbs_merge.h"
#include <string.h>
#include <math.h>

/* Tunables ----------------------------------------------------------------- */
#define INIT_TOL_CYC      0.15   /* per-sample stability during init     */
#define UPDATE_TOL_CYC    0.15   /* IIR-update gating tolerance          */
#define IIR_ALPHA         0.05   /* bias = (1-a)*bias + a*diff           */
#define INIT_NEEDED       3      /* consecutive stable samples to lock   */
#define RESET_RESID_CYC   0.75   /* abs(diff - bias) larger than this    */
#define RESET_GAP_S       5.0    /* pair gap larger than this resets     */
#define SNR_HYSTERESIS    2000   /* obsd_t SNR is in 0.001 dBHz          */

void merger_init(merger_t *m, int window_ms)
{
    memset(m, 0, sizeof(*m));
    m->window_ms = window_ms > 0 ? window_ms : 40;
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

static void set_obs_time(obsd_t *obs, int n, gtime_t t)
{
    for (int i = 0; i < n; i++) obs[i].time = t;
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
        double pair_window_s = (double)m->window_ms / 1000.0;
        if (dt <= pair_window_s) {
            *na = copy_buf(&m->a, out_a); *have_a = 1;
            *nb = copy_buf(&m->b, out_b); *have_b = 1;
            *out_time = m->a.time;
            set_obs_time(out_a, *na, *out_time);
            set_obs_time(out_b, *nb, *out_time);
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

/* avg SNR across non-empty signals */
static double avg_snr(const obsd_t *o)
{
    double sum = 0.0;
    int    cnt = 0;
    for (int j = 0; j < NFREQ+NEXOBS; j++) {
        if (o->SNR[j] > 0) { sum += (double)o->SNR[j]; cnt++; }
    }
    return cnt > 0 ? sum / cnt : -1.0;
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
 * B observations are not appended to the outbound MSM stream. This helper is
 * kept for diagnostics/experiments, but the production A-master policy below
 * uses B only internally for alignment state and network residual modelling.
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
        } else {
            ob->L[j] = 0.0;
            ob->D[j] = 0.0;
            ob->LLI[j] = 0;
            if (has_history) m->n_b_unaligned_lli++;
        }
    }
    return n_norm;
}

int merger_align_and_merge(merger_t *m,
                           const obsd_t *obs_a, int na, int have_a,
                           const obsd_t *obs_b, int nb, int have_b,
                           obsd_t *out_obs)
{
    int out_n = 0;

    /* Phase 1: update bias state for every sat present on both sides. */
    if (have_a && have_b) {
        for (int i = 0; i < na; i++) {
            int sat = obs_a[i].sat;
            int sat_idx = sat - 1;
            if (sat_idx < 0 || sat_idx >= MAXSAT) continue;
            const obsd_t *ob = NULL;
            for (int j = 0; j < nb; j++) {
                if (obs_b[j].sat == sat) { ob = &obs_b[j]; break; }
            }
            if (!ob) continue;
            for (int slot_a = 0; slot_a < NFREQ+NEXOBS; slot_a++) {
                if (obs_a[i].L[slot_a] == 0.0) continue;
                if (!obs_a[i].code[slot_a]) continue;
                update_one(m, sat_idx, slot_a, &obs_a[i], ob, obs_a[i].time);
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
        return 0;
    }

    /* Both present: keep A as the phase reference for shared satellites. */
    for (int i = 0; i < na && out_n < MAXOBS; i++) {
        int sat = obs_a[i].sat;
        int bi  = -1;
        for (int j = 0; j < nb; j++) {
            if (obs_b[j].sat == sat) { bi = j; break; }
        }
        if (bi >= 0) m->n_both++;
        out_obs[out_n] = obs_a[i];
        m->n_chose_a++;
        out_n++;
    }
    return out_n;
}
