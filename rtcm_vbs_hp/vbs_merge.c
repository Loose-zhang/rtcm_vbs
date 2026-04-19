/*------------------------------------------------------------------------------
 * vbs_merge.c : dual-base CNR-priority observation merger
 *----------------------------------------------------------------------------*/
#include "vbs_merge.h"
#include <string.h>
#include <math.h>

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

/* average SNR (0.001 dBHz units) across all non-empty signals of an obsd_t */
static double avg_snr(const obsd_t *o)
{
    double sum = 0.0;
    int    cnt = 0;
    for (int j = 0; j < NFREQ+NEXOBS; j++) {
        if (o->SNR[j] > 0) { sum += o->SNR[j]; cnt++; }
    }
    return cnt > 0 ? sum / cnt : -1.0;
}

static int fill_single(const epoch_buf_t *src, int side,
                       obsd_t *out, unsigned char *out_side)
{
    int n = src->n > MAXOBS ? MAXOBS : src->n;
    memcpy(out, src->data, n * sizeof(obsd_t));
    memset(out_side, (unsigned char)side, n);
    return n;
}

/* Merge two epochs; populate out[] + out_side[] and return length. */
static int do_merge(merger_t *m,
                    const epoch_buf_t *a, const epoch_buf_t *b,
                    obsd_t *out, unsigned char *out_side)
{
    int out_n = 0;
    unsigned char used_b[MAXOBS] = {0};

    for (int i = 0; i < a->n && out_n < MAXOBS; i++) {
        int sat = a->data[i].sat;
        int bi  = -1;
        for (int j = 0; j < b->n; j++) {
            if (b->data[j].sat == sat) { bi = j; break; }
        }
        if (bi < 0) {
            out[out_n]      = a->data[i];
            out_side[out_n] = 0;
            out_n++;
            continue;
        }
        used_b[bi] = 1;
        m->n_both++;
        double sa = avg_snr(&a->data[i]);
        double sb = avg_snr(&b->data[bi]);
        if (sb > sa) {
            out[out_n]      = b->data[bi];
            out_side[out_n] = 1;
            m->n_chose_b++;
        } else {
            out[out_n]      = a->data[i];
            out_side[out_n] = 0;
            m->n_chose_a++;
        }
        out_n++;
    }
    for (int j = 0; j < b->n && out_n < MAXOBS; j++) {
        if (!used_b[j]) {
            out[out_n]      = b->data[j];
            out_side[out_n] = 1;
            out_n++;
        }
    }
    return out_n;
}

int merger_poll(merger_t *m, long now_ms,
                obsd_t *out_obs, unsigned char *out_side,
                int *out_n, gtime_t *out_time)
{
    /* Case 1: both sides have observations at (approximately) the same
       epoch → merge immediately. */
    if (m->a.valid && m->b.valid) {
        double dt = fabs(timediff(m->a.time, m->b.time));
        if (dt <= 0.010) {
            *out_n    = do_merge(m, &m->a, &m->b, out_obs, out_side);
            *out_time = m->a.time;
            m->n_merged++;
            m->a.valid = m->b.valid = 0;
            return 1;
        }
        /* Different epochs: flush the older one unmerged so we don't lose
           it and so the newer stays pending for its twin. */
        if (timediff(m->a.time, m->b.time) < 0) {
            *out_n    = fill_single(&m->a, 0, out_obs, out_side);
            *out_time = m->a.time;
            m->n_a_only++;
            m->a.valid = 0;
            return 1;
        } else {
            *out_n    = fill_single(&m->b, 1, out_obs, out_side);
            *out_time = m->b.time;
            m->n_b_only++;
            m->b.valid = 0;
            return 1;
        }
    }

    /* Case 2: only one side valid and its wait window has elapsed. */
    if (m->a.valid && !m->b.valid && now_ms - m->a.ms_recv >= m->window_ms) {
        *out_n    = fill_single(&m->a, 0, out_obs, out_side);
        *out_time = m->a.time;
        m->n_a_only++;
        m->a.valid = 0;
        return 1;
    }
    if (m->b.valid && !m->a.valid && now_ms - m->b.ms_recv >= m->window_ms) {
        *out_n    = fill_single(&m->b, 1, out_obs, out_side);
        *out_time = m->b.time;
        m->n_b_only++;
        m->b.valid = 0;
        return 1;
    }

    return 0;
}
