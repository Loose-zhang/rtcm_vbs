/*------------------------------------------------------------------------------
 * vbs_merge.h : dual-base observation merger (CNR-priority)
 *
 * Holds the latest observation epoch from base A and base B, waits up to
 * EPOCH_WINDOW_MS ms for the twin, then produces a merged observation set:
 *
 *   - sat only in A : take A
 *   - sat only in B : take B
 *   - sat in both   : take the one with higher average SNR across signals
 *
 * The merged set is then passed back to VBS correction + MSM encoding.
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
    int         window_ms;    /* epoch alignment window */
    long        n_merged;
    long        n_a_only;
    long        n_b_only;
    long        n_both;
} merger_t;

void merger_init(merger_t *m, int window_ms);

/* Stash the obs coming from channel A (side==0) or B (side==1).
 * Overwrites any stale epoch without consuming. */
void merger_stash(merger_t *m, int side, const obsd_t *obs, int n,
                  gtime_t t, long now_ms);

/* If the pending epoch is ready (twin has arrived OR window elapsed),
 * write the merged set into out_obs[MAXOBS] and return its length.
 * Return 0 if nothing to emit yet. */
int  merger_poll(merger_t *m, long now_ms,
                 obsd_t *out_obs, int *out_n, gtime_t *out_time);

#ifdef __cplusplus
}
#endif
#endif
