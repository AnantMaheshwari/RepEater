/*
 * wfa_align.c — Gap-linear Wavefront Alignment (WFA); 
 * following Myers O(Nd) algorithm, simplified for short LTR-RT sequences
 *
 *
 * Algorithm
 * ---------
 * We use the formulation:
 *
 *   WF[d][k] = farthest row i in A reachable on diagonal k using exactly
 *              d edits, where diagonal k = i - j (j = column in B).
 *
 * Transitions from score d to d+1:
 *   Mismatch  (substitute): row advances by 1, diagonal unchanged
 *     candidate_mm = WF[d][k] + 1
 *   Deletion  (skip char in A, advance A → diagonal k increases by 1):
 *     candidate_del = WF[d][k-1] + 1   (from previous diagonal k-1, row +1)
 *   Insertion (skip char in B, advance B → diagonal k decreases by 1):
 *     candidate_ins = WF[d][k+1]        (from previous diagonal k+1, row unchanged)
 *
 *   WF[d+1][k] = max(candidate_mm, candidate_del, candidate_ins)
 *   Then extend: while A[WF[d+1][k]] == B[WF[d+1][k] - k]: advance row.
 *
 *
 * Identity: approximated as 1 - d / max(la, lb).
 */

#include "wfa_align.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Internal entry points (no longer exported; see wfa_align.h). */
static int wfa_align(const char *a, int la, const char *b, int lb,
                     int max_edit, float min_coverage,
                     float *out_identity, int *out_aln_len);
static int wfa_align_anchored(const char *a, int la, const char *b, int lb,
                              int anchor_a, int anchor_b,
                              int max_edit, float min_coverage,
                              float *out_identity, int *out_aln_len);

#define WFA_NEG  (-0x3FFFFFFF)   /* unreachable diagonal */

/* ── Inline extend ────────────────────────────────────────────────────────── */

static inline int extend_diag(
    const char *a, int la,
    const char *b, int lb,
    int k, int row_in
) {
    int i = row_in;
    int j = i - k;
    while (i < la && j >= 0 && j < lb && a[i] == b[j]) {
        i++; j++;
    }
    return i;   /* updated row */
}

/* ── wfa_align ───────────────────────────────────────────────────────────── */

static int wfa_align(
    const char *a, int la,
    const char *b, int lb,
    int         max_edit,
    float       min_coverage,
    float      *out_identity,
    int        *out_aln_len
) {
    *out_identity = -1.0f;
    *out_aln_len  =  0;

    if (la <= 0 || lb <= 0) return -1;

    /* Coverage pre-filter (Wicker 80/80/80 rule): to skip on alignment of lengths are too different. */
    int lmin = la < lb ? la : lb;
    int lmax = la > lb ? la : lb;
    if (min_coverage > 0.0f && (float)lmin / (float)lmax < min_coverage)
        return -2;

    /* Default edit budget: 25% of the longer sequence. */
    if (max_edit <= 0)
        max_edit = lmax / 4 + 1;

    int k_final = la - lb;

    int buf_size = 2 * max_edit + 3;
    int offset   = max_edit + 1;

    int *wf_cur  = (int *)malloc((size_t)buf_size * sizeof(int));
    int *wf_prev = (int *)malloc((size_t)buf_size * sizeof(int));
    if (!wf_cur || !wf_prev) { free(wf_cur); free(wf_prev); return -1; }

    for (int i = 0; i < buf_size; i++) { wf_cur[i] = WFA_NEG; wf_prev[i] = WFA_NEG; }

    /* ── Initialise score 0: extend diagonal 0 from (0,0) ── */
    {
        int row = extend_diag(a, la, b, lb, 0, 0);
        wf_cur[0 + offset] = row;
    }

    /* Check done at score 0 (sequences identical, or share common prefix
     * up to the full alignment). */
    {
        int kfi = k_final + offset;
        if (kfi >= 0 && kfi < buf_size && wf_cur[kfi] >= la) {
            free(wf_cur); free(wf_prev);
            *out_identity = 1.0f;
            *out_aln_len  = lmax;
            return 0;
        }
    }

    int result_d = -1;

    for (int d = 1; d <= max_edit; d++) {
        int *tmp = wf_prev; wf_prev = wf_cur; wf_cur = tmp;
        for (int i = 0; i < buf_size; i++) wf_cur[i] = WFA_NEG;

        /* Active diagonal range for this score */
        int k_lo = -d < -lb ? -lb : -d;
        int k_hi =  d >  la ?  la :  d;

        for (int k = k_lo; k <= k_hi; k++) {
            int ki = k + offset;
            if (ki < 0 || ki >= buf_size) continue;

            int best = WFA_NEG;

            /* Mismatch: from (d-1, k), advance both A and B */
            if (wf_prev[ki] != WFA_NEG) {
                int v = wf_prev[ki] + 1;
                if (v > best) best = v;
            }

            /* Deletion: from (d-1, k-1), advance A only  → diagonal k-1 → k */
            int ki_del = k - 1 + offset;
            if (ki_del >= 0 && ki_del < buf_size && wf_prev[ki_del] != WFA_NEG) {
                int v = wf_prev[ki_del] + 1;
                if (v > best) best = v;
            }

            /* Insertion: from (d-1, k+1), advance B only → diagonal k+1 → k */
            int ki_ins = k + 1 + offset;
            if (ki_ins >= 0 && ki_ins < buf_size && wf_prev[ki_ins] != WFA_NEG) {
                int v = wf_prev[ki_ins];
                if (v > best) best = v;
            }

            if (best == WFA_NEG) continue;

            int row_min = k > 0 ? k : 0;
            int row_max = la < lb + k ? la : lb + k;
            if (row_max < 0) continue;    /* diagonal unreachable */
            if (best < row_min) best = row_min;
            if (best > row_max) { wf_cur[ki] = row_max; continue; }

            /* Extend along matching characters */
            best = extend_diag(a, la, b, lb, k, best);
            wf_cur[ki] = best;
        }

        /* Check termination: reached (la, lb)? */
        int kfi = k_final + offset;
        if (kfi >= 0 && kfi < buf_size && wf_cur[kfi] >= la) {
            result_d = d;
            break;
        }
    }

    free(wf_cur);
    free(wf_prev);

    if (result_d < 0) return -1;   /* did not converge within max_edit */

    /* Approximate identity: 1 - d / lmax. */
    float identity = (float)(1.0 - (double)result_d / (double)lmax);
    if (identity < 0.0f) identity = 0.0f;

    *out_identity = identity;
    *out_aln_len  = lmax;
    return result_d;
}

/* ── wfa_align_reuse ─────────────────────────────────────────────────────
 *
 * Same algorithm as wfa_align but accepts pre-allocated wavefront buffers
 */
static int wfa_align_reuse(
    const char *a, int la,
    const char *b, int lb,
    int         max_edit,
    float       min_coverage,
    float      *out_identity,
    int        *out_aln_len,
    int        *wf_cur,
    int        *wf_prev,
    int         buf_cap
) {
    *out_identity = -1.0f;
    *out_aln_len  =  0;

    if (la <= 0 || lb <= 0) return -1;

    int lmin = la < lb ? la : lb;
    int lmax = la > lb ? la : lb;
    if (min_coverage > 0.0f && (float)lmin / (float)lmax < min_coverage)
        return -2;

    if (max_edit <= 0)
        max_edit = lmax / 4 + 1;

    int k_final  = la - lb;
    int buf_size = 2 * max_edit + 3;
    int offset   = max_edit + 1;

    if (buf_size > buf_cap)
        return wfa_align(a, la, b, lb, max_edit, min_coverage,
                         out_identity, out_aln_len);

    for (int i = 0; i < buf_size; i++) { wf_cur[i] = WFA_NEG; wf_prev[i] = WFA_NEG; }

    {
        int row = extend_diag(a, la, b, lb, 0, 0);
        wf_cur[0 + offset] = row;
    }
    {
        int kfi = k_final + offset;
        if (kfi >= 0 && kfi < buf_size && wf_cur[kfi] >= la) {
            *out_identity = 1.0f;
            *out_aln_len  = lmax;
            return 0;
        }
    }

    int result_d = -1;
    for (int d = 1; d <= max_edit; d++) {
        int *tmp = wf_prev; wf_prev = wf_cur; wf_cur = tmp;
        for (int i = 0; i < buf_size; i++) wf_cur[i] = WFA_NEG;

        int k_lo = -d < -lb ? -lb : -d;
        int k_hi =  d >  la ?  la :  d;

        for (int k = k_lo; k <= k_hi; k++) {
            int ki = k + offset;
            if (ki < 0 || ki >= buf_size) continue;

            int best = WFA_NEG;
            if (wf_prev[ki] != WFA_NEG) {
                int v = wf_prev[ki] + 1;
                if (v > best) best = v;
            }
            int ki_del = k - 1 + offset;
            if (ki_del >= 0 && ki_del < buf_size && wf_prev[ki_del] != WFA_NEG) {
                int v = wf_prev[ki_del] + 1;
                if (v > best) best = v;
            }
            int ki_ins = k + 1 + offset;
            if (ki_ins >= 0 && ki_ins < buf_size && wf_prev[ki_ins] != WFA_NEG) {
                int v = wf_prev[ki_ins];
                if (v > best) best = v;
            }

            if (best == WFA_NEG) continue;
            int row_min = k > 0 ? k : 0;
            int row_max = la < lb + k ? la : lb + k;
            if (row_max < 0) continue;
            if (best < row_min) best = row_min;
            if (best > row_max) { wf_cur[ki] = row_max; continue; }
            best = extend_diag(a, la, b, lb, k, best);
            wf_cur[ki] = best;
        }

        int kfi = k_final + offset;
        if (kfi >= 0 && kfi < buf_size && wf_cur[kfi] >= la) {
            result_d = d;
            break;
        }
    }

    if (result_d < 0) return -1;
    float identity = (float)(1.0 - (double)result_d / (double)lmax);
    if (identity < 0.0f) identity = 0.0f;
    *out_identity = identity;
    *out_aln_len  = lmax;
    return result_d;
}

/* ── wfa_align_anchored_reuse ──────────────────────────────────────────────
 */
static int wfa_align_anchored_reuse(
    const char *a, int la,
    const char *b, int lb,
    int         anchor_a,
    int         anchor_b,
    int         max_edit,
    float       min_coverage,
    float      *out_identity,
    int        *out_aln_len,
    int        *wf_cur,
    int        *wf_prev,
    int         wf_cap,
    char       *rev_buf_a,
    char       *rev_buf_b,
    int         rev_cap
) {
    *out_identity = -1.0f;
    *out_aln_len  =  0;

    if (la <= 0 || lb <= 0) return -1;

    if (anchor_a <= 0 || anchor_b <= 0)
        return wfa_align_reuse(a, la, b, lb, max_edit, min_coverage,
                               out_identity, out_aln_len,
                               wf_cur, wf_prev, wf_cap);

    if (anchor_a > la) anchor_a = la;
    if (anchor_b > lb) anchor_b = lb;

    int lmax = la > lb ? la : lb;
    {
        int lmin = la < lb ? la : lb;
        if (min_coverage > 0.0f && (float)lmin / (float)lmax < min_coverage)
            return -2;
    }

    if (max_edit <= 0)
        max_edit = lmax / 4 + 1;

    int lmax_right = (la - anchor_a) > (lb - anchor_b)
                     ? (la - anchor_a) : (lb - anchor_b);
    int lmax_left  = anchor_a > anchor_b ? anchor_a : anchor_b;
    int total_half = lmax_right + lmax_left;
    int budget_right = (total_half > 0)
        ? (int)((float)max_edit * (float)lmax_right / (float)total_half) + 1
        : max_edit / 2 + 1;
    int budget_left = max_edit - budget_right;
    if (budget_left  < 0) budget_left  = 0;
    if (budget_right < 0) budget_right = 0;

    /* Right half */
    float id_right = -1.0f; int al_right = 0;
    int d_right = wfa_align_reuse(a + anchor_a, la - anchor_a,
                                   b + anchor_b, lb - anchor_b,
                                   budget_right, 0.0f, &id_right, &al_right,
                                   wf_cur, wf_prev, wf_cap);

    /* Left half: reverse and align */
    if (anchor_a > rev_cap || anchor_b > rev_cap) {
        return wfa_align_anchored(a, la, b, lb, anchor_a, anchor_b,
                                   max_edit, min_coverage,
                                   out_identity, out_aln_len);
    }
    for (int i = 0; i < anchor_a; i++) rev_buf_a[i] = a[anchor_a - 1 - i];
    for (int i = 0; i < anchor_b; i++) rev_buf_b[i] = b[anchor_b - 1 - i];

    float id_left = -1.0f; int al_left = 0;
    int d_left = wfa_align_reuse(rev_buf_a, anchor_a, rev_buf_b, anchor_b,
                                  budget_left, 0.0f, &id_left, &al_left,
                                  wf_cur, wf_prev, wf_cap);

    if (d_right < 0 || d_left < 0)
        return -1;

    int total_d = d_right + d_left;
    float identity = (float)(1.0 - (double)total_d / (double)lmax);
    if (identity < 0.0f) identity = 0.0f;
    *out_identity = identity;
    *out_aln_len  = lmax;
    return total_d;
}

/* ── wfa_align_cigar and helpers ─────────────────────────────────────────── */

/*
 * CigarRun: a single run-length-encoded alignment operation.
 * Used by _wfa_build_runs, wfa_align_cigar, and wfa_align_cigar_anchored.
 */
typedef struct { int n; char op; } CigarRun;

#define EXT(d, k)  wf_ext[(size_t)(d) * (size_t)buf_size + (size_t)((k) + offset)]
#define OP(d, k)   wf_op [(size_t)(d) * (size_t)buf_size + (size_t)((k) + offset)]

/*
 * _wfa_build_runs() — internal WFA + traceback, returns CigarRun.
 */
static int _wfa_build_runs(
    const char *a, int la,
    const char *b, int lb,
    int         max_edit,
    float       min_coverage,
    float      *out_id,
    int        *out_al,
    CigarRun  **out_runs,
    int        *out_n_runs
) {
    *out_id    = -1.0f;
    *out_al    =  0;
    *out_runs  = NULL;
    *out_n_runs = 0;

    if (la <= 0 || lb <= 0) return -1;

    int lmin = la < lb ? la : lb;
    int lmax = la > lb ? la : lb;
    if (min_coverage > 0.0f && (float)lmin / (float)lmax < min_coverage)
        return -2;

    if (max_edit <= 0)
        max_edit = lmax / 4 + 1;

    int buf_size = 2 * max_edit + 3;
    int offset   = max_edit + 1;
    int k_final  = la - lb;

    size_t cells = (size_t)(max_edit + 1) * (size_t)buf_size;
    int  *wf_ext = (int  *)malloc(cells * sizeof(int));
    char *wf_op  = (char *)malloc(cells * sizeof(char));
    if (!wf_ext || !wf_op) { free(wf_ext); free(wf_op); return -1; }

    for (size_t i = 0; i < cells; i++) wf_ext[i] = WFA_NEG;
    for (size_t i = 0; i < cells; i++) wf_op[i]  = 'S';

    {
        int row = extend_diag(a, la, b, lb, 0, 0);
        EXT(0, 0) = row;
        OP(0, 0)  = 'S';
    }

    int result_d = -1;

    {
        int kfi = k_final + offset;
        if (kfi >= 0 && kfi < buf_size && EXT(0, k_final) >= la) {
            result_d = 0;
            goto tb;
        }
    }

    for (int d = 1; d <= max_edit; d++) {
        int k_lo = (-d > -lb) ? -d : -lb;
        int k_hi = ( d <  la) ?  d :  la;

        for (int k = k_lo; k <= k_hi; k++) {
            int ki = k + offset;
            if (ki < 0 || ki >= buf_size) continue;

            int best    = WFA_NEG;
            char best_op = 0;

            if (EXT(d-1, k) != WFA_NEG) {
                int v = EXT(d-1, k) + 1;
                if (v > best) { best = v; best_op = 'M'; }
            }
            if (k - 1 + offset >= 0 && k - 1 + offset < buf_size &&
                EXT(d-1, k-1) != WFA_NEG) {
                int v = EXT(d-1, k-1) + 1;
                if (v > best) { best = v; best_op = 'D'; }
            }
            if (k + 1 + offset >= 0 && k + 1 + offset < buf_size &&
                EXT(d-1, k+1) != WFA_NEG) {
                int v = EXT(d-1, k+1);
                if (v > best) { best = v; best_op = 'I'; }
            }

            if (best == WFA_NEG) continue;

            int row_min = k > 0 ? k : 0;
            int row_max = la < lb + k ? la : lb + k;
            if (row_max < 0) continue;
            if (best < row_min) best = row_min;

            OP(d, k) = best_op;

            if (best > row_max) { EXT(d, k) = row_max; continue; }

            best = extend_diag(a, la, b, lb, k, best);
            EXT(d, k) = best;
        }

        {
            int kfi = k_final + offset;
            if (kfi >= 0 && kfi < buf_size && EXT(d, k_final) >= la) {
                result_d = d;
                goto tb;
            }
        }
    }

    free(wf_ext); free(wf_op);
    return -1;

tb:
    {
        float identity = (float)(1.0 - (double)result_d / (double)lmax);
        if (identity < 0.0f) identity = 0.0f;
        *out_id = identity;
        *out_al = lmax;
    }

    int max_runs = 2 * result_d + 2;
    CigarRun *runs = (CigarRun *)malloc((size_t)(max_runs < 2 ? 2 : max_runs)
                                        * sizeof(CigarRun));
    if (!runs) { free(wf_ext); free(wf_op); return result_d; }

    int n_runs = 0;
    int d      = result_d;
    int k      = k_final;

    while (d > 0) {
        char op     = OP(d, k);
        int  k_prev = (op == 'D') ? (k - 1) : (op == 'I') ? (k + 1) : k;
        int  delta  = (op != 'I') ? 1 : 0;

        int prev_ext = EXT(d - 1, k_prev);
        int pre_raw  = (prev_ext != WFA_NEG) ? (prev_ext + delta) : delta;
        int ext_val  = EXT(d, k);
        int n_match  = ext_val - pre_raw;
        if (n_match < 0) n_match = 0;

        if (n_match > 0) { runs[n_runs].n = n_match; runs[n_runs].op = 'M'; n_runs++; }
        runs[n_runs].n = 1; runs[n_runs].op = op; n_runs++;

        k = k_prev;
        d--;
    }

    {
        int n_match = EXT(0, 0);
        if (n_match == WFA_NEG) n_match = 0;
        if (n_match > 0) { runs[n_runs].n = n_match; runs[n_runs].op = 'M'; n_runs++; }
    }

    /* Reverse to forward order. */
    for (int i = 0, j = n_runs - 1; i < j; i++, j--) {
        CigarRun tmp = runs[i]; runs[i] = runs[j]; runs[j] = tmp;
    }

    free(wf_ext); free(wf_op);
    *out_runs   = runs;
    *out_n_runs = n_runs;
    return result_d;
}

#undef EXT
#undef OP


static int _encode_runs(
    const CigarRun *runs, int n_runs,
    char *out, int cap,
    char prev_op, int prev_len
) {
    int pos = 0;
    for (int i = 0; i < n_runs; i++) {
        int  len = runs[i].n;
        char op  = runs[i].op;
        if (op == prev_op) {
            prev_len += len;
            continue;
        }
        if (prev_op && prev_len > 0) {
            int w = snprintf(out + pos, (size_t)(cap - pos), "%d%c", prev_len, prev_op);
            if (w <= 0 || w >= cap - pos) goto done;
            pos += w;
        }
        prev_op  = op;
        prev_len = len;
    }
    if (prev_op && prev_len > 0) {
        int w = snprintf(out + pos, (size_t)(cap - pos), "%d%c", prev_len, prev_op);
        if (w > 0 && w < cap - pos) pos += w;
    }
done:
    if (pos < cap) out[pos] = '\0'; else out[cap - 1] = '\0';
    return pos;
}

int wfa_align_cigar(
    const char *a, int la,
    const char *b, int lb,
    int         max_edit,
    float       min_coverage,
    float      *out_identity,
    int        *out_aln_len,
    char       *out_cigar,
    int         cigar_cap
) {
    *out_identity = -1.0f;
    *out_aln_len  =  0;
    if (out_cigar && cigar_cap > 0) out_cigar[0] = '\0';

    CigarRun *runs = NULL; int n_runs = 0;
    int rc = _wfa_build_runs(a, la, b, lb, max_edit, min_coverage,
                             out_identity, out_aln_len, &runs, &n_runs);
    if (rc < 0 || !out_cigar || cigar_cap < 1) { free(runs); return rc; }

    _encode_runs(runs, n_runs, out_cigar, cigar_cap, 0, 0);
    free(runs);
    return rc;
}

/* ── wfa_align_batch ─────────────────────────────────────────────────────── */

void wfa_align_batch(
    const char *const *seqs,
    const int         *lens,
    const int         *pair_i,
    const int         *pair_j,
    int                n_pairs,
    float              max_edit_frac,
    float              min_coverage,
    float             *out_identity,
    int               *out_aln_len
) {
    if (n_pairs <= 0) return;

    /* Pre-compute max wavefront buffer size across all pairs. */
    int global_max_edit = 0;
    for (int k = 0; k < n_pairs; k++) {
        int la = lens[pair_i[k]];
        int lb = lens[pair_j[k]];
        int lmax = la > lb ? la : lb;
        int me = (int)(lmax * max_edit_frac) + 1;
        if (me > global_max_edit) global_max_edit = me;
    }
    int wf_cap = 2 * global_max_edit + 3;

#ifdef _OPENMP
    #pragma omp parallel
    {
        /* Per-thread pre-allocated wavefront buffers. */
        int *wf_cur  = (int *)malloc((size_t)wf_cap * sizeof(int));
        int *wf_prev = (int *)malloc((size_t)wf_cap * sizeof(int));

        #pragma omp for schedule(dynamic, 64)
        for (int k = 0; k < n_pairs; k++) {
            int i  = pair_i[k];
            int j  = pair_j[k];
            int la = lens[i];
            int lb = lens[j];
            int lmax = la > lb ? la : lb;
            int max_edit = (int)(lmax * max_edit_frac) + 1;

            float id = -1.0f;
            int   al =  0;
            int rc = wfa_align_reuse(seqs[i], la, seqs[j], lb,
                                     max_edit, min_coverage, &id, &al,
                                     wf_cur, wf_prev, wf_cap);
            out_identity[k] = (rc >= 0) ? id : -1.0f;
            out_aln_len [k] = (rc >= 0) ? al :  0;
        }
        free(wf_cur); free(wf_prev);
    }
#else
    {
        int *wf_cur  = (int *)malloc((size_t)wf_cap * sizeof(int));
        int *wf_prev = (int *)malloc((size_t)wf_cap * sizeof(int));
        for (int k = 0; k < n_pairs; k++) {
            int i  = pair_i[k];
            int j  = pair_j[k];
            int la = lens[i];
            int lb = lens[j];
            int lmax = la > lb ? la : lb;
            int max_edit = (int)(lmax * max_edit_frac) + 1;

            float id = -1.0f;
            int   al =  0;
            int rc = wfa_align_reuse(seqs[i], la, seqs[j], lb,
                                     max_edit, min_coverage, &id, &al,
                                     wf_cur, wf_prev, wf_cap);
            out_identity[k] = (rc >= 0) ? id : -1.0f;
            out_aln_len [k] = (rc >= 0) ? al :  0;
        }
        free(wf_cur); free(wf_prev);
    }
#endif
}

/* ── wfa_align_anchored ──────────────────────────────────────────────────── */

static int wfa_align_anchored(
    const char *a, int la,
    const char *b, int lb,
    int         anchor_a,
    int         anchor_b,
    int         max_edit,
    float       min_coverage,
    float      *out_identity,
    int        *out_aln_len
) {
    *out_identity = -1.0f;
    *out_aln_len  =  0;

    if (la <= 0 || lb <= 0) return -1;

    /* Fall back to global if anchor is at or before the start. */
    if (anchor_a <= 0 || anchor_b <= 0)
        return wfa_align(a, la, b, lb, max_edit, min_coverage,
                         out_identity, out_aln_len);

    /* Clamp anchor to valid range. */
    if (anchor_a > la) anchor_a = la;
    if (anchor_b > lb) anchor_b = lb;

    int lmax = la > lb ? la : lb;

    /* Coverage pre-filter on full sequences. */
    {
        int lmin = la < lb ? la : lb;
        if (min_coverage > 0.0f && (float)lmin / (float)lmax < min_coverage)
            return -2;
    }

    if (max_edit <= 0)
        max_edit = lmax / 4 + 1;

    /* Split budget proportional to half-lengths. */
    int lmax_right = (la - anchor_a) > (lb - anchor_b)
                     ? (la - anchor_a) : (lb - anchor_b);
    int lmax_left  = anchor_a > anchor_b ? anchor_a : anchor_b;
    int total_half = lmax_right + lmax_left;
    int budget_right = (total_half > 0)
        ? (int)((float)max_edit * (float)lmax_right / (float)total_half) + 1
        : max_edit / 2 + 1;
    int budget_left = max_edit - budget_right;
    if (budget_left  < 0) budget_left  = 0;
    if (budget_right < 0) budget_right = 0;

    /* Right half: a[anchor_a..la) vs b[anchor_b..lb) */
    float id_right = -1.0f; int al_right = 0;
    int d_right = wfa_align(a + anchor_a, la - anchor_a,
                             b + anchor_b, lb - anchor_b,
                             budget_right, 0.0f, &id_right, &al_right);

    /* Left half: reverse a[0..anchor_a) and b[0..anchor_b), align. */
    char *rev_a = (char *)malloc((size_t)anchor_a);
    char *rev_b = (char *)malloc((size_t)anchor_b);
    if (!rev_a || !rev_b) { free(rev_a); free(rev_b); return -1; }
    for (int i = 0; i < anchor_a; i++) rev_a[i] = a[anchor_a - 1 - i];
    for (int i = 0; i < anchor_b; i++) rev_b[i] = b[anchor_b - 1 - i];

    float id_left = -1.0f; int al_left = 0;
    int d_left = wfa_align(rev_a, anchor_a, rev_b, anchor_b,
                            budget_left, 0.0f, &id_left, &al_left);
    free(rev_a); free(rev_b);

    /* If either half exceeded its budget */
    if (d_right < 0 || d_left < 0)
        return -1;

    int total_d = d_right + d_left;

    float identity = (float)(1.0 - (double)total_d / (double)lmax);
    if (identity < 0.0f) identity = 0.0f;
    *out_identity = identity;
    *out_aln_len  = lmax;
    return total_d;
}

/* ── wfa_align_batch_anchored ────────────────────────────────────────────── */

void wfa_align_batch_anchored(
    const char *const *seqs,
    const int         *lens,
    const int         *pair_i,
    const int         *pair_j,
    const int         *anchor_pi,
    const int         *anchor_pj,
    int                n_pairs,
    float              max_edit_frac,
    float              min_coverage,
    float             *out_identity,
    int               *out_aln_len
) {
    if (n_pairs <= 0) return;

    int global_max_edit = 0;
    int global_max_len  = 0;
    for (int k = 0; k < n_pairs; k++) {
        int la = lens[pair_i[k]];
        int lb = lens[pair_j[k]];
        int lmax = la > lb ? la : lb;
        int me = (int)(lmax * max_edit_frac) + 1;
        if (me > global_max_edit) global_max_edit = me;
        if (lmax > global_max_len) global_max_len = lmax;
    }
    int wf_cap = 2 * global_max_edit + 3;

#ifdef _OPENMP
    #pragma omp parallel
    {
        /* Per-thread pre-allocated buffers. */
        int  *wf_cur   = (int  *)malloc((size_t)wf_cap * sizeof(int));
        int  *wf_prev  = (int  *)malloc((size_t)wf_cap * sizeof(int));
        char *rev_a    = (char *)malloc((size_t)global_max_len);
        char *rev_b    = (char *)malloc((size_t)global_max_len);

        #pragma omp for schedule(dynamic, 64)
        for (int k = 0; k < n_pairs; k++) {
            int i  = pair_i[k];
            int j  = pair_j[k];
            int la = lens[i];
            int lb = lens[j];
            int lmax = la > lb ? la : lb;
            int max_edit = (int)(lmax * max_edit_frac) + 1;

            int anch_a = anchor_pi[k];
            int anch_b = anchor_pj[k];

            float id = -1.0f;
            int   al =  0;
            int rc = wfa_align_anchored_reuse(
                seqs[i], la, seqs[j], lb, anch_a, anch_b,
                max_edit, min_coverage, &id, &al,
                wf_cur, wf_prev, wf_cap,
                rev_a, rev_b, global_max_len);
            out_identity[k] = (rc >= 0) ? id : -1.0f;
            out_aln_len [k] = (rc >= 0) ? al :  0;
        }
        free(wf_cur); free(wf_prev);
        free(rev_a);  free(rev_b);
    }
#else
    {
        int  *wf_cur   = (int  *)malloc((size_t)wf_cap * sizeof(int));
        int  *wf_prev  = (int  *)malloc((size_t)wf_cap * sizeof(int));
        char *rev_a    = (char *)malloc((size_t)global_max_len);
        char *rev_b    = (char *)malloc((size_t)global_max_len);
        for (int k = 0; k < n_pairs; k++) {
            int i  = pair_i[k];
            int j  = pair_j[k];
            int la = lens[i];
            int lb = lens[j];
            int lmax = la > lb ? la : lb;
            int max_edit = (int)(lmax * max_edit_frac) + 1;

            int anch_a = anchor_pi[k];
            int anch_b = anchor_pj[k];

            float id = -1.0f;
            int   al =  0;
            int rc = wfa_align_anchored_reuse(
                seqs[i], la, seqs[j], lb, anch_a, anch_b,
                max_edit, min_coverage, &id, &al,
                wf_cur, wf_prev, wf_cap,
                rev_a, rev_b, global_max_len);
            out_identity[k] = (rc >= 0) ? id : -1.0f;
            out_aln_len [k] = (rc >= 0) ? al :  0;
        }
        free(wf_cur); free(wf_prev);
        free(rev_a);  free(rev_b);
    }
#endif
}

/* ── wfa_align_cigar_anchored ────────────────────────────────────────────── */

int wfa_align_cigar_anchored(
    const char *a, int la,
    const char *b, int lb,
    int         anchor_a,
    int         anchor_b,
    int         max_edit,
    float       min_coverage,
    float      *out_identity,
    int        *out_aln_len,
    char       *out_cigar,
    int         cigar_cap
) {
    *out_identity = -1.0f;
    *out_aln_len  =  0;
    if (out_cigar && cigar_cap > 0) out_cigar[0] = '\0';

    if (la <= 0 || lb <= 0) return -1;

    /* Fall back to global alignment if anchor is invalid. */
    if (anchor_a <= 0 || anchor_b <= 0)
        return wfa_align_cigar(a, la, b, lb, max_edit, min_coverage,
                               out_identity, out_aln_len, out_cigar, cigar_cap);

    if (anchor_a > la) anchor_a = la;
    if (anchor_b > lb) anchor_b = lb;

    int lmax = la > lb ? la : lb;

    {
        int lmin = la < lb ? la : lb;
        if (min_coverage > 0.0f && (float)lmin / (float)lmax < min_coverage)
            return -2;
    }

    if (max_edit <= 0)
        max_edit = lmax / 4 + 1;

    int lmax_right = (la - anchor_a) > (lb - anchor_b)
                     ? (la - anchor_a) : (lb - anchor_b);
    int lmax_left  = anchor_a > anchor_b ? anchor_a : anchor_b;
    int total_half = lmax_right + lmax_left;
    int budget_right = (total_half > 0)
        ? (int)((float)max_edit * (float)lmax_right / (float)total_half) + 1
        : max_edit / 2 + 1;
    int budget_left = max_edit - budget_right;
    if (budget_left  < 0) budget_left  = 0;
    if (budget_right < 0) budget_right = 0;

    /* Right half. */
    CigarRun *runs_right = NULL; int n_right = 0;
    float id_right = -1.0f; int al_right = 0;
    int d_right = _wfa_build_runs(a + anchor_a, la - anchor_a,
                                   b + anchor_b, lb - anchor_b,
                                   budget_right, 0.0f,
                                   &id_right, &al_right,
                                   &runs_right, &n_right);

    /* Left half on reversed prefixes. */
    char *rev_a = (char *)malloc((size_t)anchor_a);
    char *rev_b = (char *)malloc((size_t)anchor_b);
    if (!rev_a || !rev_b) {
        free(rev_a); free(rev_b); free(runs_right); return -1;
    }
    for (int i = 0; i < anchor_a; i++) rev_a[i] = a[anchor_a - 1 - i];
    for (int i = 0; i < anchor_b; i++) rev_b[i] = b[anchor_b - 1 - i];

    CigarRun *runs_left = NULL; int n_left = 0;
    float id_left = -1.0f; int al_left = 0;
    int d_left = _wfa_build_runs(rev_a, anchor_a, rev_b, anchor_b,
                                  budget_left, 0.0f,
                                  &id_left, &al_left,
                                  &runs_left, &n_left);
    free(rev_a); free(rev_b);

    /* If both halves failed, report failure. */
    if (d_right < 0 && d_left < 0) {
        free(runs_left); free(runs_right);
        return -1;
    }

    if (d_right < 0 || d_left < 0) {
        free(runs_left); free(runs_right);
        return -1;
    }

    /* Combined identity. */
    int total_d = d_right + d_left;
    float identity = (float)(1.0 - (double)total_d / (double)lmax);
    if (identity < 0.0f) identity = 0.0f;
    *out_identity = identity;
    *out_aln_len  = lmax;

    /* Encode CIGAR: reverse(runs_left) + runs_right. */
    if (out_cigar && cigar_cap > 0 && (n_left > 0 || n_right > 0)) {
        for (int i = 0, j = n_left - 1; i < j; i++, j--) {
            CigarRun tmp = runs_left[i]; runs_left[i] = runs_left[j]; runs_left[j] = tmp;
        }
        char prev_op = 0; int prev_len = 0;
        int pos = 0;
        for (int i = 0; i < n_left; i++) {
            char op  = runs_left[i].op;
            int  len = runs_left[i].n;
            if (op == prev_op) { prev_len += len; continue; }
            if (prev_op && prev_len > 0) {
                int w = snprintf(out_cigar + pos, (size_t)(cigar_cap - pos),
                                 "%d%c", prev_len, prev_op);
                if (w <= 0 || w >= cigar_cap - pos) goto enc_done;
                pos += w;
            }
            prev_op = op; prev_len = len;
        }
        for (int i = 0; i < n_right; i++) {
            char op  = runs_right[i].op;
            int  len = runs_right[i].n;
            if (op == prev_op) { prev_len += len; continue; }
            if (prev_op && prev_len > 0) {
                int w = snprintf(out_cigar + pos, (size_t)(cigar_cap - pos),
                                 "%d%c", prev_len, prev_op);
                if (w <= 0 || w >= cigar_cap - pos) goto enc_done;
                pos += w;
            }
            prev_op = op; prev_len = len;
        }
        /* Flush last run. */
        if (prev_op && prev_len > 0) {
            int w = snprintf(out_cigar + pos, (size_t)(cigar_cap - pos),
                             "%d%c", prev_len, prev_op);
            if (w > 0 && w < cigar_cap - pos) pos += w;
        }
enc_done:
        out_cigar[pos < cigar_cap ? pos : cigar_cap - 1] = '\0';
    }

    free(runs_left);
    free(runs_right);
    return total_d;
}


