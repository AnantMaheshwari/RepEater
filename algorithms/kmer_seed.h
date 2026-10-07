/*
 * kmer_seed.h — Canonical 12-mer seeding for LTR-RT candidate pair generation.
 */
#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * kmer_free()
*/
void kmer_free(void *ptr);

/*
 * kmer_emit_pairs_pos_alloc()
 */
int kmer_emit_pairs_pos_alloc(
    const char *const *seqs,
    const int         *lens,
    int                n_seqs,
    int                kmer_size,
    int                max_bucket,
    int              **out_i,
    int              **out_j,
    int              **out_pi,
    int              **out_pj
);

/*
 * kmer_seed_and_anchor_v2_alloc()
 * Memory-efficient k-mer seeding with anchor selection. 
 */
int kmer_seed_and_anchor_v2_alloc(
    const char *const *seqs,
    const int         *lens,
    int                n_seqs,
    int                kmer_size,
    int                max_bucket,
    int                min_shared,
    float              max_len_ratio,
    int                filter_n,
    int                filter_mode,
    int              **out_i,
    int              **out_j,
    int              **out_pi,
    int              **out_pj,
    int              **out_rc,
    long              *out_n_raw_events,
    long              *out_n_unique_pairs,
    int               *out_bucket_hist,
    int               *out_hist_overflow
);

#ifdef __cplusplus
}
#endif
