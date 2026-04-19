/*------------------------------------------------------------------------------
 * vbs_merge.c : dual-base CNR-priority observation merger
 *----------------------------------------------------------------------------*/
#include "vbs_merge.h"
#include <string.h>
#include <math.h>
#include <stdio.h>

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

/* average SNR (0.001 dBHz units) across all non-empty signals of obsd_t. */
static double avg_snr(const obsd_t *o)
{
    double sum = 0.0;
    int    cnt = 0;
    for (int j = 0; j < NFREQ+NEXOBS; j++) {
        if (o->SNR[j] > 0) { sum += o->SNR[j]; cnt++; }
    }
    return cnt > 0 ? sum / cnt : -1.0;
}

/* merge two epochs. caller guarantees a/b valid. */
static int do_merge(const epoch_buf_t *a, const epoch_buf_t *b,
                    obsd_t *out, int *n_both)
{
    int out_n = 0;
    int used_b[MAXOBS] = {0};

    /* iterate A: either keep A, or replace with B if B has same sat & higher SNR */
    for (int i = 0; i < a->n && out_n < MAXOBS; i++) {
        int sat = a->data[i].sat;
        int bi  = -1;
        for (int j = 0; j < b->n; j++) {
            if (b->data[j].sat == sat) { bi = j; break; }
        }
        if (bi < 0) {
            out[out_n++] = a->data[i];
            continue;
        }
        used_b[bi] = 1;
        (*n_both)++;
        double sa = avg_snr(&a->data[i]);
        double sb = avg_snr(&b->data[bi]);
        out[out_n++] = (sb > sa) ? b->data[bi] : a->data[i];
    }
    /* append satellites that appear only in B */
    for (int j = 0; j < b->n && out_n < MAXOBS; j++) {
        if (!used_b[j]) out[out_n++] = b->data[j];
    }
    return out_n;
}

int merger_poll(merger_t *m, long now_ms,
                obsd_t *out_obs, int *out_n, gtime_t *out_time)
{
    /* Case 1: both sides have the same epoch (≤ 1 ms diff) → merge immediately */
    if (m->a.valid && m->b.valid) {
        double dt = fabs(timediff(m->a.time, m->b.time));
        if (dt <= 0.010) {
            int nb = 0;
            *out_n    = do_merge(&m->a, &m->b, out_obs, &nb);
            *out_time = m->a.time;
            m->n_both   += nb;
            m->n_merged += 1;
            m->a.valid = m->b.valid = 0;
            return 1;
        }
        /* different epochs: flush the older one right away */
        if (timediff(m->a.time, m->b.time) < 0) {
            /* a is older → emit A */
            memcpy(out_obs, m->a.data, m->a.n * sizeof(obsd_t));
            *out_n    = m->a.n;
            *out_time = m->a.time;
            m->n_a_only++;
            m->a.valid = 0;
            return 1;
        } else {
            memcpy(out_obs, m->b.data, m->b.n * sizeof(obsd_t));
            *out_n    = m->b.n;
            *out_time = m->b.time;
            m->n_b_only++;
            m->b.valid = 0;
            return 1;
        }
    }

    /* Case 2: only one side valid → wait window_ms then emit that side alone */
    if (m->a.valid && !m->b.valid && now_ms - m->a.ms_recv >= m->window_ms) {
        memcpy(out_obs, m->a.data, m->a.n * sizeof(obsd_t));
        *out_n    = m->a.n;
        *out_time = m->a.time;
        m->n_a_only++;
        m->a.valid = 0;
        return 1;
    }
    if (m->b.valid && !m->a.valid && now_ms - m->b.ms_recv >= m->window_ms) {
        memcpy(out_obs, m->b.data, m->b.n * sizeof(obsd_t));
        *out_n    = m->b.n;
        *out_time = m->b.time;
        m->n_b_only++;
        m->b.valid = 0;
        return 1;
    }

    return 0;
}
