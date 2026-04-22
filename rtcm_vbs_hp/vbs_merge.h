/*------------------------------------------------------------------------------
 * vbs_merge.h : dual-base observation merger with carrier-phase alignment
 *
 * Pipeline (intended caller order):
 *
 *   1. merger_stash()  : buffer raw obs from each side as they arrive.
 *   2. merger_poll_pair() : when both sides for the same epoch are ready (or
 *                           a wait window has elapsed and only one is left),
 *                           hand back the raw A and/or B obs sets so the
 *                           caller can run vbs_correct_obs_ex() on each
 *                           side independently.
 *   3. merger_align_and_merge() : given the two VBS-corrected obs sets,
 *                                 update the per-(sat,code) B->A carrier
 *                                 bias state, normalise B's carrier phases
 *                                 into A's reference frame, and merge them
 *                                 with continuity-aware policy. Sets LLI
 *                                 cycle-slip on signals whose alignment is
 *                                 not (yet) trusted.
 *
 * The bias state lives inside merger_t and is keyed by (sat-1, signal code)
 * so that different signal layouts on A and B do not get aligned to each
 * other accidentally.
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

/* Per-(satellite, signal-code) carrier-phase alignment state.
 *
 *   diff_cyc = L_A_vbs - L_B_vbs
 * is computed each epoch where both A and B see the same satellite + code.
 * Once `init_count` consecutive samples agree to within INIT_TOL_CYC the
 * mean of those samples is locked in as `bias_cyc` and `valid` becomes 1.
 * Subsequent epochs slowly track the bias with an IIR step. Any of: input
 * cycle-slip LLI, missing/changed code on either side, an inter-epoch gap
 * larger than RESET_GAP_S, or a residual larger than RESET_RESID_CYC,
 * resets the state and forces LLI_SLIP downstream until re-initialised. */
typedef struct {
    unsigned char valid;        /* bias_cyc usable for normalisation */
    unsigned char code;         /* CODE_??? this state is keyed to */
    unsigned char init_count;   /* consecutive stable samples */
    double        init_sum;     /* running sum during init phase */
    double        bias_cyc;     /* L_A - L_B once valid */
    double        last_diff;    /* most recent raw diff (debug) */
    gtime_t       last_pair_t;  /* time of last successful pair update */
} align_state_t;

typedef struct {
    epoch_buf_t a, b;
    int         window_ms;
    /* Max |time_A - time_B| (s) to treat A/B as the same epoch for dual merge.
     * Larger values reduce A-only/B-only ping-pong when streams have timing skew. */
    double      pair_tol_s;

    /* Carrier-alignment state per (sat, freq slot). The slot index is the
     * index inside obsd_t::L; we additionally validate `code` to make sure
     * A and B's slot j really refer to the same signal. */
    align_state_t align[MAXSAT][NFREQ+NEXOBS];

    /* Stats */
    long        n_merged;       /* epochs where both A and B were paired */
    long        n_a_only;       /* epochs emitted with A only            */
    long        n_b_only;       /* epochs emitted with B only            */
    long        n_both;         /* sat counts where both A and B had it  */
    long        n_chose_a;
    long        n_chose_b;

    long        n_align_init_ok;
    long        n_align_reset_slip;
    long        n_align_reset_code;
    long        n_align_reset_gap;
    long        n_align_reset_resid;

    long        n_b_norm_applied;
    long        n_b_unaligned_lli;
    long        n_b_prop_nodopp;   /* B signals skipped in propagation: D==0 & dt significant */
} merger_t;

/* pair_tol_ms: max A/B obs time difference for paired merge (1–250 ms, default 50). */
void merger_init(merger_t *m, int window_ms, int pair_tol_ms);

/* Stash the obs arriving on channel A (side==0) or B (side==1).
 * Overwrites any stale epoch without consuming. */
void merger_stash(merger_t *m, int side, const obsd_t *obs, int n,
                  gtime_t t, long now_ms);

/* Poll for a ready pair.
 *
 * On success returns 1 and sets:
 *   *have_a / *have_b    : flags
 *   out_a[*na] / out_b[*nb] : raw obs copies (callers must run their own
 *                              VBS correction on each side before merging).
 *   *out_time            : epoch timestamp (A's if both, otherwise the
 *                          present side's).
 *
 * Returns 0 if nothing is ready yet. */
int  merger_poll_pair(merger_t *m, long now_ms,
                      obsd_t *out_a, int *na, int *have_a,
                      obsd_t *out_b, int *nb, int *have_b,
                      gtime_t *out_time);

/* Take two already VBS-corrected obs sets (either may be empty) and
 * produce the final merged obs set. B's carrier phases are normalised
 * into A's reference frame using the per-(sat,code) bias state, which
 * is updated in-place. Sets LLI cycle-slip on signals whose alignment
 * cannot be trusted this epoch.
 *
 * Returns the number of obs written to out_obs. */
int  merger_align_and_merge(merger_t *m,
                            const obsd_t *obs_a, int na, int have_a,
                            const obsd_t *obs_b, int nb, int have_b,
                            gtime_t pair_time,
                            obsd_t *out_obs);

#ifdef __cplusplus
}
#endif
#endif
