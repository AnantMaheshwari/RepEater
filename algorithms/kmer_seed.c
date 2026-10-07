/*
 * kmer_seed.c — Canonical 12-mer seeding for LTR-RT candidate pair generation.
 * This is to find candidate pairs with an orientation and anchor position.
 * Extra engineering here (inspired by FASTGA) is to avoid memory blow-up during kmer seeding
 * hence we take a 2-phase approach:
 * 1) use a compact pair to only compute key, count as a first pass
 * 2) pairs reaching min_shared (e.g. 3) we then compute the full PairStats to avoid computing on thrown away pairs
 * We use a standard hashing approach to store these open-addressed tables (from phase 1, 2) 
 * also handle reverse complements via heuristic: for shared kmer-hits compute variance of diagonal, anti-diagonal
 * and use that as a cheaper proxy for orientation.
 */

#include "kmer_seed.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ── Memory profiling helper ─────────────────────────────────────────────── */
static long _vmrss_mb(void)
{
    FILE *f = fopen("/proc/self/status", "r");
    if (!f) return -1;
    char line[256];
    long rss_kb = -1;
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            sscanf(line + 6, " %ld", &rss_kb);
            break;
        }
    }
    fclose(f);
    return (rss_kb > 0) ? rss_kb / 1024 : -1;
}

static int _kmer_seed_verbose(void)
{
    static int v = -1;
    if (v < 0) {
        const char *e = getenv("KMER_SEED_VERBOSE");
        v = (e && *e && *e != '0') ? 1 : 0;
    }
    return v;
}
#define KSLOG(...) do { if (_kmer_seed_verbose()) fprintf(stderr, __VA_ARGS__); } while (0)

/* ── Base encoding ──────────────────────────────────────────────────────────
 * A=0, C=1, G=2, T=3.  Any other character returns 4 (invalid).           */

static const uint8_t BASE_ENC[256] = {
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  /* 0x00 */
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  /* 0x10 */
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  /* 0x20 */
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  /* 0x30 */
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,  /* 0x40  A C G */
    4,4,4,4, 3,4,4,4, 4,4,4,4, 4,4,4,4,  /* 0x50  T */
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,  /* 0x60  a c g */
    4,4,4,4, 3,4,4,4, 4,4,4,4, 4,4,4,4,  /* 0x70  t */
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
};

/* Reverse-complement a k-mer encoded in the lower 2k bits of v. */
static inline uint32_t rc_kmer(uint32_t v, int k)
{
    uint32_t rc = 0;
    for (int i = 0; i < k; i++) {
        rc = (rc << 2) | (3u ^ (v & 3u));
        v >>= 2;
    }
    return rc;
}

/* ── kmer_encode_canonical ──────────────────────────────────────────────── */

static int kmer_encode_canonical(const char *seq, int len, int k,
                                uint32_t *out_kmers)
{
    if (len < k || k < 1 || k > 15) return 0;

    int out_n = 0;
    uint32_t fwd = 0;
    uint32_t mask = (k < 16) ? ((1u << (2 * k)) - 1u) : 0xFFFFFFFFu;
    int valid = 0;   /* number of consecutive valid bases in the window */

    for (int i = 0; i < len; i++) {
        uint8_t b = BASE_ENC[(unsigned char)seq[i]];
        if (b == 4) {
            /* Invalid base (N etc.): reset window */
            fwd   = 0;
            valid = 0;
            continue;
        }
        fwd = ((fwd << 2) | b) & mask;
        valid++;

        if (valid >= k) {
            uint32_t rev = rc_kmer(fwd, k);
            out_kmers[out_n++] = (fwd < rev) ? fwd : rev;
        }
    }
    return out_n;
}

/* ── Sorting helpers ────────────────────────────────────────────────────── */

/* Flat entry: (kmer, elem_idx). */
typedef uint64_t KmerEntry;   /* high 32 bits = kmer, low 32 bits = elem_idx */

#define ENTRY(kmer, idx)  (((uint64_t)(kmer) << 32) | (uint32_t)(idx))
#define ENTRY_KMER(e)     ((uint32_t)((e) >> 32))
#define ENTRY_IDX(e)      ((uint32_t)((e) & 0xFFFFFFFFu))

/* Position-aware flat entry for kmer_emit_pairs_pos_alloc().
 * Sorted on .entry (kmer, elem_idx); .pos carries the k-mer start position
 * in the source sequence and is not part of the sort key. */
typedef struct { uint64_t entry; uint32_t pos; } KmerEntryPos;

static int cmp_entry_pos(const void *a, const void *b)
{
    const KmerEntryPos *x = (const KmerEntryPos *)a;
    const KmerEntryPos *y = (const KmerEntryPos *)b;
    if (x->entry < y->entry) return -1;
    if (x->entry > y->entry) return  1;
    return 0;
}

#define PACK_KEP(kmer, idx, pos) \
    (((uint64_t)(kmer) << 40) | ((uint64_t)(idx) << 16) | (uint64_t)(pos))
#define KEP_KMER(e)  ((uint32_t)((e) >> 40))
#define KEP_IDX(e)   ((uint32_t)(((e) >> 16) & 0xFFFFFFu))
#define KEP_POS(e)   ((uint16_t)((e) & 0xFFFFu))

static int cmp_uint64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a;
    uint64_t y = *(const uint64_t *)b;
    if (x < y) return -1;
    if (x > y) return  1;
    return 0;
}

/* Comparator for output anchors -- sort by (i, j) for deterministic order. */
typedef struct { int i, j, pi, pj, rc; } OutAnchor;
static int cmp_out_anchor(const void *a, const void *b)
{
    const OutAnchor *x = (const OutAnchor *)a;
    const OutAnchor *y = (const OutAnchor *)b;
    if (x->i != y->i) return (x->i < y->i) ? -1 : 1;
    if (x->j != y->j) return (x->j < y->j) ? -1 : 1;
    return 0;
}

/* ── seeding ─────────────────────────────────
 *
 * wo passes over the sorted k-mer array:
 *
 *   Pass 1: bucket scan → for each pair, accumulate count and
 *           diagonal / anti-diagonal running sums in the hash table.
 *   Pass 2: re-scan buckets → for surviving pairs, pick the k-mer hit
 *           closest to the target diagonal as the anchor.

 */

/* ── PairStats: per-pair accumulator in the hash table ─────────────────── */
typedef struct {
    uint64_t key;       /* 8: (i+1) << 32 | (j+1); 0 = empty slot */
    uint32_t count;     /* 4: number of shared k-mers */
    uint32_t flags;     /* 4: bit 0 = alive, bit 1 = is_rc */
    union {
        struct { double mean_d, mean_a, M2_d, M2_a; } p1; 
        struct { double target, best_diff;
                 int32_t best_pi, best_pj;           } p2;  
    };
} PairStats;  /* total: 48 bytes */

#define HT_EMPTY_KEY 0ULL

static inline uint64_t _ht_pair_key(int i, int j)
{
    /* i < j guaranteed; add 1 so that (0,0) doesn't collide with empty */
    return ((uint64_t)(i + 1) << 32) | (uint64_t)(j + 1);
}

/* standard Fibonacci hashing for good distribution (open-addressed table) */
static inline size_t _ht_hash(uint64_t key, int shift)
{
    return (size_t)((key * 11400714819323198485ULL) >> shift);
}

/* Find or insert a key in the open-addressed table.
 * Returns pointer to the slot (existing or newly inserted). */
static inline PairStats *_ht_find_or_insert(
    PairStats *table, size_t capacity, int shift, uint64_t key, int *did_insert)
{
    size_t idx = _ht_hash(key, shift);
    size_t mask = capacity - 1; 
    *did_insert = 0;

    for (;;) {
        idx &= mask;
        if (table[idx].key == key)
            return &table[idx];
        if (table[idx].key == HT_EMPTY_KEY) {
            table[idx].key = key;
            *did_insert = 1;
            return &table[idx];
        }
        idx++; 
    }
}

/* Find a key; return NULL if not present. */
static inline PairStats *_ht_find(
    PairStats *table, size_t capacity, int shift, uint64_t key)
{
    size_t idx = _ht_hash(key, shift);
    size_t mask = capacity - 1;

    for (;;) {
        idx &= mask;
        if (table[idx].key == key)
            return &table[idx];
        if (table[idx].key == HT_EMPTY_KEY)
            return NULL;
        idx++;
    }
}


/* ── CompactPair: count-only HT for two-phase seeding ────────────────────  */
typedef struct {
    uint64_t key;    
    uint32_t count;
} CompactPair;

static inline CompactPair *_cht_find_or_insert(
    CompactPair *table, size_t capacity, int shift, uint64_t key, int *did_insert)
{
    size_t idx = _ht_hash(key, shift);
    size_t mask = capacity - 1;
    *did_insert = 0;

    for (;;) {
        idx &= mask;
        if (table[idx].key == key)
            return &table[idx];
        if (table[idx].key == HT_EMPTY_KEY) {
            table[idx].key = key;
            *did_insert = 1;
            return &table[idx];
        }
        idx++;
    }
}

static inline CompactPair *_cht_find(
    CompactPair *table, size_t capacity, int shift, uint64_t key)
{
    size_t idx = _ht_hash(key, shift);
    size_t mask = capacity - 1;

    for (;;) {
        idx &= mask;
        if (table[idx].key == key)
            return &table[idx];
        if (table[idx].key == HT_EMPTY_KEY)
            return NULL;
        idx++;
    }
}



/* ── kmer_free ──────────────────────────────────────────────────────────── */

void kmer_free(void *ptr) { free(ptr); }


/* ── kmer_emit_pairs_pos_alloc ─────────────────────────────────────────── */

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
) {
    *out_i = *out_j = *out_pi = *out_pj = NULL;
    if (n_seqs <= 0 || kmer_size < 4 || kmer_size > 15) return 0;

    /* 1. Estimate total k-mers and allocate flat array */
    size_t total_est = 0;
    for (int i = 0; i < n_seqs; i++)
        if (lens[i] >= kmer_size) total_est += (size_t)(lens[i] - kmer_size + 1);

    KmerEntryPos *flat = (KmerEntryPos *)malloc(total_est * sizeof(KmerEntryPos));
    if (!flat) return 0;

    int max_len = 0;
    for (int i = 0; i < n_seqs; i++)
        if (lens[i] > max_len) max_len = lens[i];

    uint32_t *tmp_kmers = (uint32_t *)malloc((size_t)max_len * sizeof(uint32_t));
    if (!tmp_kmers) { free(flat); return 0; }

    /* 2. Build flat array with positions */
    size_t flat_n = 0;
    for (int i = 0; i < n_seqs; i++) {
        if (lens[i] < kmer_size) continue;
        int n = kmer_encode_canonical(seqs[i], lens[i], kmer_size, tmp_kmers);
        int stride_num = (n > 1) ? (lens[i] - kmer_size) : 0;
        int stride_den = (n > 1) ? (n - 1) : 1;
        for (int p = 0; p < n; p++) {
            uint32_t pos = (uint32_t)(((int64_t)p * stride_num + stride_den / 2) / stride_den);
            flat[flat_n].entry = ENTRY(tmp_kmers[p], (uint32_t)i);
            flat[flat_n].pos   = pos;
            flat_n++;
        }
    }
    free(tmp_kmers);

    /* 3. Sort */
    qsort(flat, flat_n, sizeof(KmerEntryPos), cmp_entry_pos);

    /* 4. Scan buckets, emit pairs */
    int *bucket_ids = (int *)malloc((size_t)(max_bucket + 1) * sizeof(int));
    int *bucket_pos = (int *)malloc((size_t)(max_bucket + 1) * sizeof(int));
    if (!bucket_ids || !bucket_pos) {
        free(bucket_ids); free(bucket_pos); free(flat); return 0;
    }

    size_t buf_cap = 65536;
    int *buf_i  = (int *)malloc(buf_cap * sizeof(int));
    int *buf_j  = (int *)malloc(buf_cap * sizeof(int));
    int *buf_pi = (int *)malloc(buf_cap * sizeof(int));
    int *buf_pj = (int *)malloc(buf_cap * sizeof(int));
    if (!buf_i || !buf_j || !buf_pi || !buf_pj) {
        free(buf_i); free(buf_j); free(buf_pi); free(buf_pj);
        free(bucket_ids); free(bucket_pos); free(flat);
        return 0;
    }
    int out_n = 0;

    size_t p = 0;
    while (p < flat_n) {
        uint32_t cur_kmer = ENTRY_KMER(flat[p].entry);
        size_t q = p + 1;
        while (q < flat_n && ENTRY_KMER(flat[q].entry) == cur_kmer) q++;

        int n_ids = 0;
        for (size_t r = p; r < q && n_ids <= max_bucket; r++) {
            int idx = (int)ENTRY_IDX(flat[r].entry);
            int dup = 0;
            for (int s = 0; s < n_ids; s++)
                if (bucket_ids[s] == idx) { dup = 1; break; }
            if (!dup) {
                bucket_ids[n_ids] = idx;
                bucket_pos[n_ids] = (int)flat[r].pos;
                n_ids++;
            }
        }

        if (n_ids >= 2 && n_ids <= max_bucket) {
            for (int a = 0; a < n_ids; a++) {
                for (int b = a + 1; b < n_ids; b++) {
                    int ei = bucket_ids[a], ej = bucket_ids[b];
                    int pi = bucket_pos[a], pj = bucket_pos[b];
                    if (ei > ej) {
                        int t; t = ei; ei = ej; ej = t;
                                t = pi; pi = pj; pj = t;
                    }

                    if ((size_t)out_n >= buf_cap) {
                        buf_cap *= 2;
                        int *ni  = (int *)realloc(buf_i,  buf_cap * sizeof(int));
                        int *nj  = (int *)realloc(buf_j,  buf_cap * sizeof(int));
                        int *npi = (int *)realloc(buf_pi, buf_cap * sizeof(int));
                        int *npj = (int *)realloc(buf_pj, buf_cap * sizeof(int));
                        if (!ni || !nj || !npi || !npj) {
                            free(ni ? ni : buf_i); free(nj ? nj : buf_j);
                            free(npi ? npi : buf_pi); free(npj ? npj : buf_pj);
                            free(bucket_ids); free(bucket_pos); free(flat);
                            return 0;
                        }
                        buf_i = ni; buf_j = nj; buf_pi = npi; buf_pj = npj;
                    }
                    buf_i[out_n]  = ei;
                    buf_j[out_n]  = ej;
                    buf_pi[out_n] = pi;
                    buf_pj[out_n] = pj;
                    out_n++;
                }
            }
        }
        p = q;
    }

    free(flat);
    free(bucket_ids);
    free(bucket_pos);

    if (out_n == 0) {
        free(buf_i); free(buf_j); free(buf_pi); free(buf_pj);
        return 0;
    }
    *out_i = buf_i; *out_j = buf_j; *out_pi = buf_pi; *out_pj = buf_pj;
    return out_n;
}


/* ── kmer_seed_and_anchor_v2_alloc ─────────────────────────────────────── *
 * Same algorithm as kmer_seed_and_anchor_v2 but allocates output arrays   *
 * internally (exact size = alive_count).  */

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
) {
    *out_i = *out_j = *out_pi = *out_pj = *out_rc = NULL;
    if (n_seqs <= 0 || kmer_size < 4 || kmer_size > 15) return 0;

    if (kmer_size > 12) {
        fprintf(stderr, "[kmer_seed_v2_alloc] ERROR: packed mode requires kmer_size<=12 "
                "(got %d)\n", kmer_size);
        return 0;
    }
    if (n_seqs >= (1 << 24)) {
        fprintf(stderr, "[kmer_seed_v2_alloc] ERROR: packed mode requires n_seqs<%d "
                "(got %d)\n", 1 << 24, n_seqs);
        return 0;
    }

    /* ── 1. Compute total k-mer estimate and max sequence length ─────────── */
    size_t total_est = 0;
    int max_len = 0;
    for (int i = 0; i < n_seqs; i++) {
        if (lens[i] >= kmer_size)
            total_est += (size_t)(lens[i] - kmer_size + 1);
        if (lens[i] > max_len) max_len = lens[i];
    }

    uint32_t *tmp_kmers = (uint32_t *)malloc((size_t)max_len * sizeof(uint32_t));
    if (!tmp_kmers) return 0;

    /* ── 2. Prefix chunking ──────────────────────────────────────────────── */
    size_t mem_budget = 256ULL * 1024 * 1024;
    int n_chunks = 1;
    {
        size_t flat_bytes = total_est * sizeof(uint64_t);
        while (flat_bytes / (size_t)n_chunks > mem_budget && n_chunks < 256)
            n_chunks <<= 1;
    }
    int prefix_bits = 0;
    { int tmp = n_chunks; while (tmp > 1) { prefix_bits++; tmp >>= 1; } }
    int prefix_shift = 2 * kmer_size - prefix_bits;

    /* ── 3. Count k-mers per chunk ───────────────────────────────────────── */
    size_t *chunk_counts = (size_t *)calloc((size_t)n_chunks, sizeof(size_t));
    if (!chunk_counts) { free(tmp_kmers); return 0; }

    size_t total_flat_n = 0;
    for (int i = 0; i < n_seqs; i++) {
        if (lens[i] < kmer_size) continue;
        int n = kmer_encode_canonical(seqs[i], lens[i], kmer_size, tmp_kmers);
        for (int p = 0; p < n; p++) {
            int c = (n_chunks > 1) ? (int)(tmp_kmers[p] >> prefix_shift) : 0;
            chunk_counts[c]++;
        }
        total_flat_n += (size_t)n;
    }

    size_t max_chunk = 0;
    for (int c = 0; c < n_chunks; c++)
        if (chunk_counts[c] > max_chunk) max_chunk = chunk_counts[c];
    free(chunk_counts);

    uint64_t *chunk_buf = (uint64_t *)malloc(max_chunk * sizeof(uint64_t));
    if (!chunk_buf) { free(tmp_kmers); return 0; }

    KSLOG("[kmer_seed_v2_alloc:mem] chunk_buf: %zu MB (%zu entries x 8B)  "
          "VmRSS: %ld MB\n",
          max_chunk * sizeof(uint64_t) / (1024*1024), max_chunk, _vmrss_mb());

    /* ── 4. Bucket scratch arrays ────────────────────────────────────────── */
    int *bucket_ids = (int *)malloc((size_t)(max_bucket + 1) * sizeof(int));
    int *bucket_pos = (int *)malloc((size_t)(max_bucket + 1) * sizeof(int));
    if (!bucket_ids || !bucket_pos) {
        free(bucket_ids); free(bucket_pos);
        free(chunk_buf); free(tmp_kmers);
        return 0;
    }

    /* ── 5. Compact hash table for Phase 1a  */
    size_t cht_init = total_flat_n / 128;
    if (cht_init < 1024) cht_init = 1024;
    size_t cht_cap = 1;
    int cht_shift = 64;
    while (cht_cap < cht_init) { cht_cap <<= 1; cht_shift--; }
    cht_cap <<= 1; cht_shift--;

    CompactPair *cht = (CompactPair *)calloc(cht_cap, sizeof(CompactPair));
    if (!cht) {
        free(bucket_ids); free(bucket_pos);
        free(chunk_buf); free(tmp_kmers);
        return 0;
    }
    size_t cht_count = 0;

    KSLOG("[kmer_seed_v2_alloc:mem] compact_ht: %zu MB (%zu slots x %zu B)  "
          "VmRSS: %ld MB\n",
          cht_cap * sizeof(CompactPair) / (1024*1024), cht_cap,
          sizeof(CompactPair), _vmrss_mb());

    /* ── K-mer spectrum histogram (TODO (anant): study this distribution to tune thresholds) ──── */
    int *bucket_hist = (int *)calloc((size_t)(max_bucket + 2), sizeof(int));
    int loc_hist_overflow = 0;

    KSLOG("[kmer_seed_v2_alloc] %zu total k-mers, %d prefix chunks "
          "(budget %zu MB, max_chunk %zu)\n",
          total_flat_n, n_chunks,
          (size_t)(mem_budget / (1024 * 1024)), max_chunk);

    /* ── 6. PHASE 1a: Count shared kmers per pair (compact HT) ─────────── */
    long raw_events = 0;

    for (int c = 0; c < n_chunks; c++) {
        size_t chunk_n = 0;
        for (int i = 0; i < n_seqs; i++) {
            if (lens[i] < kmer_size) continue;
            int n = kmer_encode_canonical(seqs[i], lens[i], kmer_size, tmp_kmers);
            for (int p = 0; p < n; p++) {
                uint32_t kmer = tmp_kmers[p];
                if (n_chunks > 1 && (int)(kmer >> prefix_shift) != c) continue;
                /* Position not needed for counting — store 0 */
                chunk_buf[chunk_n++] = PACK_KEP(kmer, (uint32_t)i, 0);
            }
        }

        qsort(chunk_buf, chunk_n, sizeof(uint64_t), cmp_uint64);

        size_t pp = 0;
        while (pp < chunk_n) {
            uint32_t cur_kmer = KEP_KMER(chunk_buf[pp]);
            size_t qq = pp + 1;
            while (qq < chunk_n && KEP_KMER(chunk_buf[qq]) == cur_kmer) qq++;

            int n_ids = 0;
            for (size_t r = pp; r < qq && n_ids <= max_bucket; r++) {
                int idx = (int)KEP_IDX(chunk_buf[r]);
                int dup = 0;
                for (int s = 0; s < n_ids; s++)
                    if (bucket_ids[s] == idx) { dup = 1; break; }
                if (!dup) {
                    bucket_ids[n_ids] = idx;
                    n_ids++;
                }
            }

            if (bucket_hist) {
                if (n_ids >= 2 && n_ids <= max_bucket)
                    bucket_hist[n_ids]++;
                else if (n_ids > max_bucket)
                    loc_hist_overflow++;
            }

            if (n_ids >= 2 && n_ids <= max_bucket) {
                for (int a = 0; a < n_ids; a++) {
                    for (int b = a + 1; b < n_ids; b++) {
                        int ei = bucket_ids[a], ej = bucket_ids[b];
                        if (ei > ej) { int t = ei; ei = ej; ej = t; }

                        if (filter_n > 0) {
                            int i_set = (ei < filter_n) ? 0 : 1;
                            int j_set = (ej < filter_n) ? 0 : 1;
                            switch (filter_mode) {
                            case 1:
                                if (i_set != j_set &&
                                    (ei % filter_n) == (ej % filter_n))
                                    continue;
                                break;
                            case 2:
                                if (i_set != j_set) continue;
                                break;
                            case 3:
                                if (i_set != 0 || j_set != 0) continue;
                                break;
                            }
                        }

                        if (max_len_ratio > 0.0f) {
                            int li = lens[ei], lj = lens[ej];
                            float ratio = (li >= lj) ? (float)li / (float)lj
                                                     : (float)lj / (float)li;
                            if (ratio > max_len_ratio) continue;
                        }

                        raw_events++;

                        uint64_t key = _ht_pair_key(ei, ej);

                        /* Grow compact HT if needed */
                        if (cht_count * 10 >= cht_cap * 7) {
                            size_t new_cap = cht_cap << 1;
                            int new_shift = cht_shift - 1;
                            CompactPair *new_cht = (CompactPair *)calloc(
                                new_cap, sizeof(CompactPair));
                            if (!new_cht) goto alloc_phase1a_done;

                            for (size_t ri = 0; ri < cht_cap; ri++) {
                                if (cht[ri].key != HT_EMPTY_KEY) {
                                    int dummy;
                                    CompactPair *dst = _cht_find_or_insert(
                                        new_cht, new_cap, new_shift,
                                        cht[ri].key, &dummy);
                                    *dst = cht[ri];
                                }
                            }
                            free(cht);
                            cht = new_cht;
                            cht_cap = new_cap;
                            cht_shift = new_shift;
                            KSLOG("[kmer_seed_v2_alloc:mem] compact_ht rehash → "
                                  "%zu MB (%zu slots x %zu B)  "
                                  "VmRSS: %ld MB\n",
                                  new_cap * sizeof(CompactPair) / (1024*1024),
                                  new_cap, sizeof(CompactPair), _vmrss_mb());
                        }

                        int did_insert;
                        CompactPair *cp = _cht_find_or_insert(
                            cht, cht_cap, cht_shift, key, &did_insert);
                        if (did_insert) {
                            cp->count = 0;
                            cht_count++;
                        }
                        cp->count++;
                    }
                }
            }
            pp = qq;
        }
    }
alloc_phase1a_done:

    if (out_n_raw_events)   *out_n_raw_events   = raw_events;
    if (out_n_unique_pairs) *out_n_unique_pairs  = (long)cht_count;

    if (bucket_hist) {
        if (out_bucket_hist)
            memcpy(out_bucket_hist, bucket_hist,
                   (size_t)(max_bucket + 2) * sizeof(int));
        if (out_hist_overflow)
            *out_hist_overflow = loc_hist_overflow;
        free(bucket_hist);
        bucket_hist = NULL;
    }

    /* ── 6b. Count survivors in compact HT ─────────────────────────────── */
    size_t alive_count = 0;
    for (size_t ri = 0; ri < cht_cap; ri++) {
        if (cht[ri].key != HT_EMPTY_KEY && (int)cht[ri].count >= min_shared)
            alive_count++;
    }

    KSLOG("[kmer_seed_v2_alloc] %ld raw pair events, %zu unique pairs, "
          "%zu survive min_shared=%d\n",
          raw_events, cht_count, alive_count, min_shared);

    if (alive_count == 0) {
        free(cht);
        free(chunk_buf); free(tmp_kmers);
        free(bucket_ids); free(bucket_pos);
        return 0;
    }

    /* ── 6c. Allocate full PairStats HT sized for survivors only ─────── */
    size_t ht_cap = 1;
    int ht_shift = 64;
    {
        size_t ht_target = alive_count * 10 / 7;  /* ~0.7 load factor */
        if (ht_target < 1024) ht_target = 1024;
        while (ht_cap < ht_target) { ht_cap <<= 1; ht_shift--; }
    }

    PairStats *ht = (PairStats *)calloc(ht_cap, sizeof(PairStats));
    if (!ht) {
        free(cht);
        free(bucket_ids); free(bucket_pos);
        free(chunk_buf); free(tmp_kmers);
        return 0;
    }
    KSLOG("[kmer_seed_v2_alloc:mem] full_ht: %zu MB (%zu slots x %zu B, "
          "for %zu survivors)  VmRSS: %ld MB\n",
          ht_cap * sizeof(PairStats) / (1024*1024), ht_cap,
          sizeof(PairStats), alive_count, _vmrss_mb());

    /* ── 7. PHASE 1b: Welford accumulation for survivors only ──────────── *
     * For each pair if it survived min_shared, insert into full HT */
    for (int c = 0; c < n_chunks; c++) {
        size_t chunk_n = 0;
        for (int i = 0; i < n_seqs; i++) {
            if (lens[i] < kmer_size) continue;
            int n = kmer_encode_canonical(seqs[i], lens[i], kmer_size, tmp_kmers);
            int stride_num = (n > 1) ? (lens[i] - kmer_size) : 0;
            int stride_den = (n > 1) ? (n - 1) : 1;
            for (int p = 0; p < n; p++) {
                uint32_t kmer = tmp_kmers[p];
                if (n_chunks > 1 && (int)(kmer >> prefix_shift) != c) continue;
                uint32_t pos = (uint32_t)(((int64_t)p * stride_num
                                           + stride_den / 2) / stride_den);
                if (pos > 0xFFFFu) pos = 0xFFFFu;
                chunk_buf[chunk_n++] = PACK_KEP(kmer, (uint32_t)i, pos);
            }
        }

        qsort(chunk_buf, chunk_n, sizeof(uint64_t), cmp_uint64);

        size_t pp = 0;
        while (pp < chunk_n) {
            uint32_t cur_kmer = KEP_KMER(chunk_buf[pp]);
            size_t qq = pp + 1;
            while (qq < chunk_n && KEP_KMER(chunk_buf[qq]) == cur_kmer) qq++;

            int n_ids = 0;
            for (size_t r = pp; r < qq && n_ids <= max_bucket; r++) {
                int idx = (int)KEP_IDX(chunk_buf[r]);
                int dup = 0;
                for (int s = 0; s < n_ids; s++)
                    if (bucket_ids[s] == idx) { dup = 1; break; }
                if (!dup) {
                    bucket_ids[n_ids] = idx;
                    bucket_pos[n_ids] = (int)KEP_POS(chunk_buf[r]);
                    n_ids++;
                }
            }

            if (n_ids >= 2 && n_ids <= max_bucket) {
                for (int a = 0; a < n_ids; a++) {
                    for (int b = a + 1; b < n_ids; b++) {
                        int ei = bucket_ids[a], ej = bucket_ids[b];
                        int pi_v = bucket_pos[a], pj_v = bucket_pos[b];
                        if (ei > ej) {
                            int t; t = ei; ei = ej; ej = t;
                                    t = pi_v; pi_v = pj_v; pj_v = t;
                        }

                        /* Check compact HT: skip if pair didn't survive */
                        uint64_t key = _ht_pair_key(ei, ej);
                        CompactPair *cp = _cht_find(
                            cht, cht_cap, cht_shift, key);
                        if (!cp || (int)cp->count < min_shared) continue;

                        /* Insert into full HT and accumulate Welford */
                        int did_insert;
                        PairStats *ps = _ht_find_or_insert(
                            ht, ht_cap, ht_shift, key, &did_insert);
                        if (did_insert) {
                            ps->count = 0;
                            ps->flags = 0;
                            ps->p1.mean_d = ps->p1.mean_a = 0.0;
                            ps->p1.M2_d = ps->p1.M2_a = 0.0;
                        }
                        ps->count++;
                        double d = (double)pi_v - (double)pj_v;
                        double a_val = (double)pi_v + (double)pj_v;
                        double delta_d = d - ps->p1.mean_d;
                        ps->p1.mean_d += delta_d / (double)ps->count;
                        ps->p1.M2_d  += delta_d * (d - ps->p1.mean_d);
                        double delta_a = a_val - ps->p1.mean_a;
                        ps->p1.mean_a += delta_a / (double)ps->count;
                        ps->p1.M2_a  += delta_a * (a_val - ps->p1.mean_a);
                    }
                }
            }
            pp = qq;
        }
    }

    /* Free compact HT, since at this point we made the full HT which is all we need now */
    free(cht);
    cht = NULL;

    KSLOG("[kmer_seed_v2_alloc:mem] freed compact_ht  VmRSS: %ld MB\n",
          _vmrss_mb());

    /* ── 7b. Between passes: compute variance, RC, target ────────────────  */
    size_t alive_verify = 0;
    for (size_t ri = 0; ri < ht_cap; ri++) {
        PairStats *ps = &ht[ri];
        if (ps->key == HT_EMPTY_KEY) continue;

        double var_d  = ps->p1.M2_d / (double)ps->count;
        double var_ad = ps->p1.M2_a / (double)ps->count;
        int is_rc = (var_ad < var_d) ? 1 : 0;
        double target = is_rc ? ps->p1.mean_a : ps->p1.mean_d;

        ps->flags = 1u | (is_rc ? 2u : 0u);
        ps->p2.target    = target;
        ps->p2.best_diff = 1e30;
        ps->p2.best_pi   = -1;
        ps->p2.best_pj   = -1;
        alive_verify++;
    }

    /* ── 8. PASS 2: anchor selection ─────────────────────────────────────── */
    for (int c = 0; c < n_chunks; c++) {
        size_t chunk_n = 0;
        for (int i = 0; i < n_seqs; i++) {
            if (lens[i] < kmer_size) continue;
            int n = kmer_encode_canonical(seqs[i], lens[i], kmer_size, tmp_kmers);
            int stride_num = (n > 1) ? (lens[i] - kmer_size) : 0;
            int stride_den = (n > 1) ? (n - 1) : 1;
            for (int p = 0; p < n; p++) {
                uint32_t kmer = tmp_kmers[p];
                if (n_chunks > 1 && (int)(kmer >> prefix_shift) != c) continue;
                uint32_t pos = (uint32_t)(((int64_t)p * stride_num
                                           + stride_den / 2) / stride_den);
                if (pos > 0xFFFFu) pos = 0xFFFFu;
                chunk_buf[chunk_n++] = PACK_KEP(kmer, (uint32_t)i, pos);
            }
        }

        qsort(chunk_buf, chunk_n, sizeof(uint64_t), cmp_uint64);

        size_t pp = 0;
        while (pp < chunk_n) {
            uint32_t cur_kmer = KEP_KMER(chunk_buf[pp]);
            size_t qq = pp + 1;
            while (qq < chunk_n && KEP_KMER(chunk_buf[qq]) == cur_kmer) qq++;

            int n_ids = 0;
            for (size_t r = pp; r < qq && n_ids <= max_bucket; r++) {
                int idx = (int)KEP_IDX(chunk_buf[r]);
                int dup = 0;
                for (int s = 0; s < n_ids; s++)
                    if (bucket_ids[s] == idx) { dup = 1; break; }
                if (!dup) {
                    bucket_ids[n_ids] = idx;
                    bucket_pos[n_ids] = (int)KEP_POS(chunk_buf[r]);
                    n_ids++;
                }
            }

            if (n_ids >= 2 && n_ids <= max_bucket) {
                for (int a = 0; a < n_ids; a++) {
                    for (int b = a + 1; b < n_ids; b++) {
                        int ei = bucket_ids[a], ej = bucket_ids[b];
                        int pi_v = bucket_pos[a], pj_v = bucket_pos[b];
                        if (ei > ej) {
                            int t; t = ei; ei = ej; ej = t;
                                    t = pi_v; pi_v = pj_v; pj_v = t;
                        }

                        uint64_t key = _ht_pair_key(ei, ej);
                        PairStats *ps = _ht_find(ht, ht_cap, ht_shift, key);
                        if (!ps || !(ps->flags & 1u)) continue;

                        double metric = (ps->flags & 2u)
                            ? ((double)pi_v + (double)pj_v)
                            : ((double)pi_v - (double)pj_v);
                        double diff = metric - ps->p2.target;
                        if (diff < 0) diff = -diff;
                        if (diff < ps->p2.best_diff ||
                            (diff == ps->p2.best_diff &&
                             pi_v < ps->p2.best_pi)) {
                            ps->p2.best_diff = diff;
                            ps->p2.best_pi   = pi_v;
                            ps->p2.best_pj   = pj_v;
                        }
                    }
                }
            }
            pp = qq;
        }
    }

    free(chunk_buf);
    free(tmp_kmers);
    free(bucket_ids);
    free(bucket_pos);

    KSLOG("[kmer_seed_v2_alloc:mem] freed chunk_buf+scratch  VmRSS: %ld MB\n",
          _vmrss_mb());

    /* ── 9. collect output from hash table ────────── */
    int *buf_i  = (int *)malloc(alive_verify * sizeof(int));
    int *buf_j  = (int *)malloc(alive_verify * sizeof(int));
    int *buf_pi = (int *)malloc(alive_verify * sizeof(int));
    int *buf_pj = (int *)malloc(alive_verify * sizeof(int));
    int *buf_rc = (int *)malloc(alive_verify * sizeof(int));
    if (!buf_i || !buf_j || !buf_pi || !buf_pj || !buf_rc) {
        free(buf_i); free(buf_j); free(buf_pi); free(buf_pj); free(buf_rc);
        free(ht);
        return 0;
    }

    int out_n = 0;
    for (size_t ri = 0; ri < ht_cap; ri++) {
        PairStats *ps = &ht[ri];
        if (ps->key == HT_EMPTY_KEY || !(ps->flags & 1u)) continue;
        buf_i[out_n]  = (int)((ps->key >> 32) - 1);
        buf_j[out_n]  = (int)((ps->key & 0xFFFFFFFFULL) - 1);
        buf_pi[out_n] = ps->p2.best_pi;
        buf_pj[out_n] = ps->p2.best_pj;
        buf_rc[out_n] = (ps->flags & 2u) ? 1 : 0;
        out_n++;
    }

    free(ht);

    KSLOG("[kmer_seed_v2_alloc:mem] freed hash_table  VmRSS: %ld MB\n",
          _vmrss_mb());

    /* ── 10. Sort output by (i, j) for deterministic order ───────────────── */
    {
        OutAnchor *tmp = (OutAnchor *)malloc((size_t)out_n * sizeof(OutAnchor));
        if (tmp) {
            for (int k = 0; k < out_n; k++)
                tmp[k] = (OutAnchor){buf_i[k], buf_j[k],
                                     buf_pi[k], buf_pj[k], buf_rc[k]};
            qsort(tmp, (size_t)out_n, sizeof(OutAnchor), cmp_out_anchor);
            for (int k = 0; k < out_n; k++) {
                buf_i[k]  = tmp[k].i;
                buf_j[k]  = tmp[k].j;
                buf_pi[k] = tmp[k].pi;
                buf_pj[k] = tmp[k].pj;
                buf_rc[k] = tmp[k].rc;
            }
            free(tmp);
        }
    }

    *out_i = buf_i; *out_j = buf_j;
    *out_pi = buf_pi; *out_pj = buf_pj; *out_rc = buf_rc;
    return out_n;
}
