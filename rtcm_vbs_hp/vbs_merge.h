/*------------------------------------------------------------------------------
 * vbs_merge.h : dual-base observation merger (CNR-priority)
 *
 * Holds the latest observation epoch from base A and base B, waits up to
 * window_ms ms for the twin, then produces a merged observation set with
 * per-observation provenance flags:
 *
 *   - sat only in A : take A, side = 0
 *   - sat only in B : take B, side = 1
 *   - sat in both   : take the one with higher average SNR across signals,
 *                     side = 0 or 1 accordingly
 *
 * The side array is consumed by vbs_correct_obs() so each satellite is
 * corrected against the real base that actually measured it — necessary
 * for correctness when the A/B baseline is non-trivial (km scale).
 *----------------------------------------------------------------------------*/
#ifndef VBS_MERGE_H
#define VBS_MERGE_H

#include "rtklib.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    obsd_t  data[MAXOBS];
    int     n;
    gtime_t time;
    int     valid;
    long    ms_recv;          /* monotonic ms when this epoch was stored */
} epoch_buf_t;

typedef struct {
    epoch_buf_t a, b;
    int         window_ms;

    long        n_merged;     /* merged epochs (both A and B present)     */
    long        n_a_only;     /* epochs emitted with A only               */
    long        n_b_only;     /* epochs emitted with B only               */
    long        n_both;       /* running count of satellites where both   */
                              /* A and B had observations and we picked   */
                              /* one by SNR                               */
    long        n_chose_a;
    long        n_chose_b;
} merger_t;

void merger_init(merger_t *m, int window_ms);

/* Stash the obs arriving on channel A (side==0) or B (side==1).
 * Overwrites any stale epoch without consuming. */
void merger_stash(merger_t *m, int side, const obsd_t *obs, int n,
                  gtime_t t, long now_ms);

/* Poll the merger. If an epoch is ready (twin arrived OR window elapsed,
 * OR A/B had mismatched epoch times so we flush the older), write the
 * merged set into out_obs[MAXOBS] and out_side[MAXOBS] (0=A, 1=B) and
 * return 1 with *out_n and *out_time set. Return 0 if nothing yet. */
int  merger_poll(merger_t *m, long now_ms,
                 obsd_t *out_obs, unsigned char *out_side,
                 int *out_n, gtime_t *out_time);

#ifdef __cplusplus
}
#endif
#endif
