/*
 * wfa_align.h — Gap-linear Wavefront Alignment (WFA) for 80/80/80 validation.
 *
 * Implements O(N·d) wavefront algorithm with linear gap
 * penalties. Uses ~80% identity thresholding to match the Wicker TE family definition. 
 *
 */
#pragma once

#ifdef __cplusplus
extern "C" {
#endif


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
);


int wfa_align_cigar(
    const char *a, int la,
    const char *b, int lb,
    int         max_edit,
    float       min_coverage,
    float      *out_identity,
    int        *out_aln_len,
    char       *out_cigar,
    int         cigar_cap
);


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
);


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
);

#ifdef __cplusplus
}
#endif
