"""
Feature-aware LTR clustering following Wicker et al.'s 80/80/80 rule.

Pipeline
--------
1. Parse FASTA headers (produced by seq/extract_ltrs) to extract
   per-element sub-sequences: ltr5, internal, ltr3.

2. K-mer seeding (algorithms/kmer_seed.[h,c]):
   Build canonical 12-mer indices over LTR5 sequences. 
   
3. WFA validation (algorithms/wfa_align.[h,c]):
   Run WFA on seeded candidate LTR pairs, those passing min_identity and min_coverage are edges.
   A second WFA pass over the internal region is then run on pairs that pass LTR identity.

4. Graph clustering (networkx):
   Connected components of the similarity graph define TE families.
   Components with < min_cluster_size (default 3) members are discarded.

5. Consensus generation (algorithms/consensus.py):
   Progressive, UPGMA-style profile merging per cluster, with outgroup correction and subfamily splitting.

Output: consensi.fasta + cluster_membership.tsv.
"""

import gc
import os
import re
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import networkx as nx
from Bio import SeqIO

import algorithms.consensus as consensus
from algorithms.c_libs import KmerLib, WFALib
from engine import constants
from utilities import log
from utilities.timing import timed

# ---------------------------------------------------------------------------
# RC helpers
# ---------------------------------------------------------------------------

_COMP = str.maketrans('ACGTacgt', 'TGCAtgca')


def _revcomp(s: str) -> str:
    """Reverse complement over ACGT only; every other symbol is left as-is.
    """
    return s.translate(_COMP)[::-1]


def _subfamily_tag(idx: int) -> str:
    """Our subfamily notation, although indeed according to 80/80/80 these 
    really are "families" and we are treating them as Wicker families. 
    Convert 0-based index to subfamily tag: a..z, aa..zz, aaa..zzz, ..."""
    repeat = idx // 26 + 1        # 1 for a-z, 2 for aa-zz, 3 for aaa-zzz, ...
    letter = chr(ord('a') + idx % 26)
    return letter * repeat

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class FeatureClusterConfig:
    """Configuration for k-mer seeded, WFA-validated LTR clustering.
    """

    # ── Wicker 80/80/80 thresholds ──────────────────────────────────────────
    min_identity: float = constants.MIN_IDENTITY
    min_coverage: float = constants.MIN_COVERAGE

    # ── K-mer seeding ───────────────────────────────────────────────────────
    kmer_size:           int = constants.KMER_SIZE
    min_seed_shared:     int = constants.MIN_SEED_SHARED
    max_bucket:          int = constants.MAX_BUCKET
    max_bucket_internal: int = constants.MAX_BUCKET_INTERNAL

    # ── WFA ─────────────────────────────────────────────────────────────────
    max_edit_frac: float = constants.MAX_EDIT_FRAC

    # Minimum feature length to include in a comparison.
    min_ltr_len:      int = constants.MIN_LTR_LEN
    min_internal_len: int = constants.MIN_INTERNAL_LEN

    # Minimum cluster size retained in the output.
    min_cluster_size: int = constants.MIN_CLUSTER_SIZE

    # ── Consensus generation ────────────────────────────────────────────────
    cons_threshold:       float = constants.CONS_THRESHOLD
    saturation_threshold: float = constants.SATURATION_THRESHOLD
    ambiguity_threshold:  float = constants.AMBIGUITY_THRESHOLD
    use_outgroup:         bool  = constants.USE_OUTGROUP
    split_subfamilies:    bool  = constants.SPLIT_SUBFAMILIES
    min_merge_identity:   float = constants.MIN_MERGE_IDENTITY
    cut_threshold:        float = constants.CUT_THRESHOLD

    # Flank filter: remove edges where flanking sequences match, indicating
    # orthology (same insertion in different genomes) or tandem duplication.
    # Checks all 4 flank pairs (A_L, A_R) × (B_L, B_R) in both orientations.
    min_flank_len: int = constants.MIN_FLANK_LEN

    threads:          int  = 1
    keep_alignments:  bool = False
    debug:            bool = False
    use_flank_filter: bool = True


# ---------------------------------------------------------------------------
# Per-element record with feature sub-sequences
# ---------------------------------------------------------------------------

@dataclass
class LTRRecord:
    """Parsed LTR element with per-feature sub-sequences."""
    name:     str
    full_seq: str
    ltr5:     str   # left LTR  (repeat1 coords)
    internal: str   # between repeat1 and repeat2
    ltr3:     str   # right LTR (repeat2 coords)
    flank_l:  str = ''  # ~100bp upstream of left LTR
    flank_r:  str = ''  # ~100bp downstream of right LTR
    genome:   str = ''  # GML-safe genome label (set when running with --genome-list)
    ltr_divergence: float = -1.0  # intra-element LTR divergence (1 - identity); age proxy

    @property
    def ltr5_len(self)     -> int: return len(self.ltr5)
    @property
    def internal_len(self) -> int: return len(self.internal)


# ---------------------------------------------------------------------------
# FASTA header parsing
# ---------------------------------------------------------------------------

def _parse_coord(tag: str, header: str) -> Optional[Tuple[int, int]]:
    """Extract start-end from a 'tag=start-end' field in a FASTA description."""
    m = re.search(rf'{re.escape(tag)}=(\d+)-(\d+)', header)
    return (int(m.group(1)), int(m.group(2))) if m else None


def _load_records(fasta_path: str) -> List[LTRRecord]:
    """
    Parse the FASTA produced by seq/extract_ltrs.

    Headers contain scaffold-absolute coordinates for flank_l, repeat1,
    repeat2, flank_r.  The extracted sequence starts at flank_l[0], so
    feature offsets in the sequence are coord - flank_l[0].
    """
    records = []
    skipped = 0

    for rec in SeqIO.parse(fasta_path, 'fasta'):
        name   = rec.id
        header = rec.description
        seq    = str(rec.seq).upper()
        n      = len(seq)

        flank_l = _parse_coord('flank_l', header)
        repeat1 = _parse_coord('repeat1', header)
        repeat2 = _parse_coord('repeat2', header)

        if not (flank_l and repeat1 and repeat2):
            skipped += 1
            continue

        extract_start = flank_l[0]

        # Convert scaffold-absolute coords to sequence-relative offsets.
        r1s = max(0, min(repeat1[0] - extract_start, n))
        r1e = max(r1s, min(repeat1[1] - extract_start, n))
        r2s = max(r1e, min(repeat2[0] - extract_start, n))
        r2e = max(r2s, min(repeat2[1] - extract_start, n))

        genome_m = re.search(r'genome=(\S+)', header)
        records.append(LTRRecord(
            name     = name,
            full_seq = seq,
            ltr5     = seq[r1s:r1e],
            internal = seq[r1e:r2s],
            ltr3     = seq[r2s:r2e],
            flank_l  = seq[0:r1s],      # bases upstream of left LTR
            flank_r  = seq[r2e:n],       # bases downstream of right LTR
            genome   = genome_m.group(1) if genome_m else '',
        ))

    log.info(f"  Loaded {len(records)} elements with structural coordinates "
             f"({skipped} skipped)")
    return records


# ---------------------------------------------------------------------------
# C library paths
# ---------------------------------------------------------------------------

_ALGO_DIR = Path(__file__).parent.parent / 'algorithms'


# ---------------------------------------------------------------------------
# WFA validation: filter candidate pairs
# ---------------------------------------------------------------------------

def _wfa_validate(
    candidates:   List[Tuple[int, ...]],  # (i,j) or (i,j,ap_a,ap_b) or (i,j,ap_a,ap_b,is_rc)
    records:      List[LTRRecord],
    attr:         str,              # 'ltr5', 'ltr3', or 'internal'
    min_len:      int,
    config:       FeatureClusterConfig,
    wfa_lib:      WFALib,
    use_anchored: bool = False,
) -> Dict[Tuple[int, int], float]:
    """
    Run WFA on candidate pairs for the given feature region.
    Returns a dict mapping canonical (i, j) → wfa_identity for pairs that
    pass min_identity and min_coverage.
    """
    if isinstance(candidates, np.ndarray):
        if len(candidates) == 0:
            return {}
    elif not candidates:
        return {}

    if isinstance(candidates, np.ndarray) and candidates.ndim == 2:
        needed_idx = sorted(np.unique(candidates[:, :2]).tolist())
    else:
        needed_idx = sorted({idx for pair in candidates for idx in pair[:2]})
    idx_map    = {orig: local for local, orig in enumerate(needed_idx)}
    n_needed   = len(needed_idx)

    fwd_seqs  = []
    lens_list = []
    for orig_idx in needed_idx:
        seq = getattr(records[orig_idx], attr)
        if len(seq) < min_len:
            seq = ''   # will be filtered by coverage check inside WFA
        fwd_seqs.append(seq)
        lens_list.append(len(seq))

    # Append RC sequences so RC pairs can use index (local_j + n_needed)
    if isinstance(candidates, np.ndarray) and candidates.ndim == 2:
        has_rc = candidates.shape[1] >= 5 and candidates[:, 4].any()
    else:
        has_rc = any(len(p) == 5 and p[4] for p in candidates)
    if has_rc:
        rc_seqs  = [_revcomp(s) if s else '' for s in fwd_seqs]
        seqs_bytes = [s.encode() for s in fwd_seqs] + [s.encode() for s in rc_seqs]
        all_lens   = lens_list + lens_list  # same lengths for fwd and rc
    else:
        seqs_bytes = [s.encode() for s in fwd_seqs]
        all_lens   = lens_list

    m = len(candidates)

    # Build pair index arrays — numpy path for anchor arrays, list path for tuples
    if isinstance(candidates, np.ndarray) and candidates.ndim == 2:
        cand_arr = candidates
        remap = np.empty(max(idx_map.keys()) + 1, dtype=np.int32)
        for orig, local in idx_map.items():
            remap[orig] = local
        pi_arr = remap[cand_arr[:, 0]]
        pj_arr = remap[cand_arr[:, 1]]
        api_arr = cand_arr[:, 2].copy()
        apj_arr = cand_arr[:, 3].copy()
        if has_rc:
            rc_mask = cand_arr[:, 4].astype(bool)
            pj_arr[rc_mask] += n_needed
            # Flip anchor j for RC pairs: len_j - apj - kmer_size
            lens_np = np.array(lens_list, dtype=np.int32)
            rc_with_anchor = rc_mask & (apj_arr >= 0)
            if rc_with_anchor.any():
                local_b_rc = remap[cand_arr[rc_with_anchor, 1]]
                apj_arr[rc_with_anchor] = lens_np[local_b_rc] - apj_arr[rc_with_anchor] - config.kmer_size
    else:
        pi_list  = []
        pj_list  = []
        api_list = []
        apj_list = []
        for p in candidates:
            a, b  = p[0], p[1]
            is_rc = has_rc and len(p) == 5 and bool(p[4])
            local_a = idx_map[a]
            local_b = idx_map[b]
            pi_list.append(local_a)
            pj_list.append(local_b + n_needed if is_rc else local_b)

            api = int(p[2]) if len(p) > 2 else -1
            apj = int(p[3]) if len(p) > 3 else -1
            api_list.append(api)
            if is_rc and apj >= 0:
                len_j = lens_list[local_b]
                apj   = len_j - apj - config.kmer_size
            apj_list.append(apj)
        pi_arr = pi_list
        pj_arr = pj_list
        api_arr = api_list
        apj_arr = apj_list

    # Auto-chunk large pair counts to bound peak memory. 
    chunk_size = 10_000 if m > 20_000 else 0
    use_chunking = chunk_size > 0 and m > chunk_size

    if not use_chunking:
        with timed(f"WFA validation ({attr}, {m:,} pairs)"):
            out_id, out_al = wfa_lib.align_batch_anchored(
                seqs_bytes, all_lens, pi_arr, pj_arr, api_arr, apj_arr,
                config.max_edit_frac, config.min_coverage,
            )
    else:
        n_chunks = (m + chunk_size - 1) // chunk_size
        out_id = np.empty(m, dtype=np.float32)
        out_al = np.empty(m, dtype=np.int32)
        import time as _time
        _t0 = _time.monotonic()
        with timed(f"WFA validation ({attr}, {m:,} pairs)"):
            for ci in range(n_chunks):
                lo = ci * chunk_size
                hi = min(lo + chunk_size, m)
                chunk_id, chunk_al = wfa_lib.align_batch_anchored(
                    seqs_bytes, all_lens,
                    pi_arr[lo:hi], pj_arr[lo:hi],
                    api_arr[lo:hi], apj_arr[lo:hi],
                    config.max_edit_frac, config.min_coverage,
                )
                out_id[lo:hi] = chunk_id
                out_al[lo:hi] = chunk_al
                del chunk_id, chunk_al
                _el = _time.monotonic() - _t0
                _eta = _el / (ci + 1) * (n_chunks - ci - 1)
                log.progress(f"    [{attr}] chunk {ci+1}/{n_chunks} "
                             f"({hi:,}/{m:,} pairs, {_el:.0f}s elapsed, "
                             f"~{_eta:.0f}s remaining)")
            log.progress_end()

    passing: Dict[Tuple[int, int], float] = {}
    # Cast threshold to float32 to match out_id dtype (TODO anant bugfix)
    min_id_f32 = np.float32(config.min_identity)
    if isinstance(candidates, np.ndarray) and candidates.ndim == 2:
        pass_mask = out_id >= min_id_f32
        pass_idx = np.where(pass_mask)[0]
        for k in pass_idx:
            a, b = int(candidates[k, 0]), int(candidates[k, 1])
            key = (min(a, b), max(a, b))
            v = float(out_id[k])
            if key not in passing or v > passing[key]:
                passing[key] = v
    else:
        for k, pair in enumerate(candidates):
            a, b = pair[0], pair[1]
            if out_id[k] >= min_id_f32:
                key = (min(a, b), max(a, b))
                if key not in passing or out_id[k] > passing[key]:
                    passing[key] = float(out_id[k])

    return passing


# ---------------------------------------------------------------------------
# Flank-based orthology / tandem filter
# ---------------------------------------------------------------------------

def _flank_filter(
    edges:   Dict[Tuple[int, int], float],
    records: List[LTRRecord],
    config:  FeatureClusterConfig,
    wfa_lib: WFALib,
) -> Dict[Tuple[int, int], float]:
    """
    Remove edges where flanking sequences match, indicating orthology
    (same insertion site in different genomes) or tandem duplication.
    """
    if not edges:
        return edges

    min_fl = config.min_flank_len
    min_id = float(np.float32(config.min_identity))
    min_cov = config.min_coverage

    needed = sorted({idx for pair in edges for idx in pair})
    idx_to_local = {orig: loc for loc, orig in enumerate(needed)}
    n_needed = len(needed)

    seqs_bytes: List[bytes] = []
    lens_list: List[int] = []
    for idx in needed:
        r = records[idx]
        fl, fr = r.flank_l, r.flank_r
        rcl, rcr = _revcomp(fl), _revcomp(fr)
        for s in (fl, fr, rcl, rcr):
            seqs_bytes.append(s.encode())
            lens_list.append(len(s))

    pi_list: List[int] = []
    pj_list: List[int] = []
    edge_idx: List[int] = []  # which edge each pair belongs to
    edge_list = list(edges.keys())

    for ei, (a, b) in enumerate(edge_list):
        la = idx_to_local[a]
        lb = idx_to_local[b]

        if lens_list[la * 4] < min_fl or lens_list[la * 4 + 1] < min_fl:
            continue
        if lens_list[lb * 4] < min_fl or lens_list[lb * 4 + 1] < min_fl:
            continue

        for a_off in (0, 1):      # A_L, A_R
            for b_off in (0, 1, 2, 3):  # B_L, B_R, RC(B_L), RC(B_R)
                si = la * 4 + a_off
                sj = lb * 4 + b_off
                la_len, lb_len = lens_list[si], lens_list[sj]
                if la_len == 0 or lb_len == 0:
                    continue
                if min(la_len, lb_len) / max(la_len, lb_len) < min_cov:
                    continue
                pi_list.append(si)
                pj_list.append(sj)
                edge_idx.append(ei)

    if not pi_list:
        log.info("  Flank filter: no edges removed")
        return edges

    log.info(f"    flank filter: {len(pi_list):,} WFA calls for "
             f"{len(edge_list):,} edges ...")
    out_id, _ = wfa_lib.align_batch(
        seqs_bytes, lens_list, pi_list, pj_list,
        config.max_edit_frac, min_cov,
    )

    edges_to_remove: Set[int] = set()
    for k in range(len(pi_list)):
        if out_id[k] >= min_id:
            edges_to_remove.add(edge_idx[k])

    if edges_to_remove:
        remove_keys = {edge_list[ei] for ei in edges_to_remove}
        filtered = {k: v for k, v in edges.items() if k not in remove_keys}
        log.info(f"  Flank filter: removed {len(edges_to_remove):,} edges "
                 f"(ortholog/tandem, matching flanks at {min_id:.0%} identity)")
        return filtered

    log.info("  Flank filter: no edges removed")
    return edges


# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------

def _add_reuse_edges(
    G:           nx.Graph,
    reuse_edges: Optional[Dict[Tuple[str, str], dict]],
    name_to_idx: Dict[str, int],
) -> None:
    """Merge name-keyed reused edges into the positional graph G.

    Only edges whose *both* endpoints exist in the current record set are added
    Reused edges never collide with fresh new-incident edges: a reused pair has 
    two non-new endpoints, while every fresh pair is incident to a new element.
    """
    if not reuse_edges:
        return
    n_added = 0
    for (na, nb), attrs in reuse_edges.items():
        ia = name_to_idx.get(na)
        ib = name_to_idx.get(nb)
        if ia is None or ib is None:
            continue
        G.add_edge(ia, ib, **attrs)
        n_added += 1
    log.info(f"  [incremental] reused {n_added:,} stable edges")


def _build_graph(
    records:  List[LTRRecord],
    config:   FeatureClusterConfig,
    kmer_lib: KmerLib,
    wfa_lib:  WFALib,
    new_names:    Optional[Set[str]] = None,
    reuse_edges:  Optional[Dict[Tuple[str, str], dict]] = None,
    frozen_params: Optional[dict]    = None,
) -> nx.Graph:
    """
    Build an element similarity graph using k-mer seeding + WFA validation.
    """
    # Deterministic clustering. The graph is sensitive to input record order:
    # anchor-selection tie-breaks in the kmer seeding step
    records.sort(key=lambda r: (r.full_seq, r.name))

    G = nx.Graph()

    # ── Add nodes with structural annotations ───────────────────────────
    has_internal_set = {
        i for i, r in enumerate(records)
        if r.internal_len >= config.min_internal_len
    }
    for i, r in enumerate(records):
        G.add_node(i,
                   label=r.name,           # used as node label in Graphia 
                   name=r.name,
                   ltrlength=r.ltr5_len,
                   internallength=r.internal_len,
                   hasinternal=int(i in has_internal_set),
                   ltrageselfpct=round(r.ltr_divergence * 100.0, 1) if r.ltr_divergence >= 0 else -1.0,
                   # Populated after consensus generation:
                   iscentroid=0,
                   centroiddivergencepct=-1.0,
                   centroidwfafailed=0,
                   genome=r.genome)

    # ── Incremental support ──────────────────────────────────────────────
    name_to_idx = {r.name: i for i, r in enumerate(records)}
    if new_names is not None:
        new_mask = np.fromiter((r.name in new_names for r in records),
                               dtype=bool, count=len(records))
    else:
        new_mask = None

    # ── LTR seeding + WFA ─────────────────────────────────────────────────
    N = len(records)
    virtual: List[LTRRecord] = list(records)   # [0..N-1] use ltr5
    for r in records:
        virtual.append(LTRRecord(               # [N..2N-1] use ltr3
            name=r.name + '_L3', full_seq='',
            ltr5=r.ltr3, internal='', ltr3='',
            genome=r.genome,
        ))

    # same canonical k-mer hits both ltr5 and ltr3
    combined_max_bucket = config.max_bucket * 2

    # ──  kmer_seed_and_anchor ──────────────────────
    all_seqs = []
    all_lens = []
    for rec in virtual:
        seq = rec.ltr5
        all_seqs.append(seq.encode())
        all_lens.append(len(seq))

    # Length-ratio pre-filter 
    max_len_ratio = (1.0 / config.min_coverage
                    if config.min_coverage > 0 else 0.0)

    c_filter_mode = 1

    log.info(f"  LTR seeding: {2*N:,} sequences (ltr5 + ltr3), "
             f"max_bucket={combined_max_bucket}")

    with timed("k-mer seeding + anchoring (LTR)"):
        anchor_arr, n_raw, n_unique, _hist, _hist_overflow = \
            kmer_lib.seed_and_anchor_v2(
                all_seqs, all_lens,
                kmer_size=config.kmer_size,
                max_bucket=combined_max_bucket,
                min_shared=config.min_seed_shared,
                max_len_ratio=max_len_ratio,
                filter_n=N,
                filter_mode=c_filter_mode,
            )
    del all_seqs, all_lens

    log.info(f"    {n_raw:,} raw events, {n_unique:,} unique pairs "
             f"→ {len(anchor_arr):,} anchored pairs")

    anchored = anchor_arr
    del anchor_arr

    # Incremental: keep only anchored pairs incident to a new element.
    if new_mask is not None and isinstance(anchored, np.ndarray) and len(anchored):
        _ri = (anchored[:, 0] % N).astype(np.intp)
        _rj = (anchored[:, 1] % N).astype(np.intp)
        _keep = new_mask[_ri] | new_mask[_rj]
        _n_before = len(anchored)
        anchored = anchored[_keep]
        log.info(f"  [incremental] LTR anchored {_n_before:,} → "
                 f"{len(anchored):,} (new-incident only)")
        del _ri, _rj, _keep

    if isinstance(anchored, np.ndarray) and len(anchored) > 0:
        col0, col1, col4 = anchored[:, 0], anchored[:, 1], anchored[:, 4]
        n_rc    = int(col4.astype(bool).sum())
        both_lt = (col0 < N) & (col1 < N)
        both_ge = (col0 >= N) & (col1 >= N)
        n_55    = int(both_lt.sum())
        n_33    = int(both_ge.sum())
        n_cross = len(anchored) - n_55 - n_33
    else:
        n_55 = n_33 = n_cross = n_rc = 0
        for p in anchored:
            i, j = p[0], p[1]
            if p[4]:
                n_rc += 1
            if i < N and j < N:
                n_55 += 1
            elif i >= N and j >= N:
                n_33 += 1
            else:
                n_cross += 1
    log.step(f"  [kmer_seed_ltr] {2*N:,} seqs → {len(anchored):,} anchored pairs")
    log.info(f"  → {len(anchored):,} anchored pairs "
             f"(5↔5: {n_55:,}, 3↔3: {n_33:,}, cross: {n_cross:,}, "
             f"RC: {n_rc:,})")

    all_edges = _wfa_validate(
        anchored, virtual, 'ltr5', config.min_ltr_len, config, wfa_lib,
        use_anchored=True,
    )

    # Map combined indices → element pairs, best identity wins
    ltr_passing: Dict[Tuple[int, int], float] = {}
    n_55_pass = n_33_pass = n_cross_pass = 0
    for (ci, cj), identity in all_edges.items():
        ri = ci % N
        rj = cj % N
        if ci < N and cj < N:
            n_55_pass += 1
        elif ci >= N and cj >= N:
            n_33_pass += 1
        else:
            n_cross_pass += 1
        key = (min(ri, rj), max(ri, rj))
        if identity > ltr_passing.get(key, 0.0):
            ltr_passing[key] = identity

    log.info(f"  Combined LTR edges: {len(all_edges):,} "
             f"(5↔5: {n_55_pass:,}, 3↔3: {n_33_pass:,}, "
             f"cross: {n_cross_pass:,})")
    log.step(f"  [clustering_ltr] {len(anchored):,} pairs → "
             f"{len(ltr_passing):,} element pairs")

    # Free large LTR-phase intermediates before internal seeding (TODO (anant): mem optimization)
    del virtual, anchored, all_edges

    # ── Flank filter: remove orthologous / tandem edges ──────────────────
    if config.use_flank_filter:
        _n_pre_flank = len(ltr_passing)
        with timed("flank filter"):
            ltr_passing = _flank_filter(ltr_passing, records, config, wfa_lib)
        log.step(f"  Flank filter: {_n_pre_flank:,} → {len(ltr_passing):,} edges")

    log.info("  Internal region check enabled (AND logic)")

    # ── Free LTR-phase intermediates before internal seeding ─────────────
    gc.collect()

    # ── Internal region WFA on LTR-passing pairs only ─────────────────────
    # Wicker 80/80/80: LTR passes AND internal passes.
    ltr_passing_indices = sorted({idx for pair in ltr_passing for idx in pair})
    sub_records = [records[i] for i in ltr_passing_indices]
    global_idx  = {local: g for local, g in enumerate(ltr_passing_indices)}
    ltr_key_set = set(ltr_passing.keys())

    # ── Fused C path for internal seeding ───────────────────────
    int_eligible_idx  = []
    int_eligible_seqs = []
    int_eligible_lens = []
    for ri, rec in enumerate(sub_records):
        seq = rec.internal
        if len(seq) >= config.min_internal_len:
            int_eligible_idx.append(ri)
            int_eligible_seqs.append(seq.encode())
            int_eligible_lens.append(len(seq))

    max_len_ratio_int = (1.0 / config.min_coverage
                        if config.min_coverage > 0 else 0.0)

    # Adaptive min_shared for internal: scale with avg sequence length.
    # Longer sequences produce more coincidental kmer matches, so we
    # need a higher threshold. 
    if frozen_params is not None and 'int_min_shared' in frozen_params:
        int_min_shared   = frozen_params['int_min_shared']
        avg_internal_len = frozen_params.get('avg_internal_len', 0)
    else:
        avg_internal_len = (sum(int_eligible_lens) / len(int_eligible_lens)
                            if int_eligible_lens else 0)
        int_min_shared = max(config.min_seed_shared,
                             int(avg_internal_len // 1000))
        if frozen_params is not None:
            frozen_params['int_min_shared']   = int_min_shared
            frozen_params['avg_internal_len'] = avg_internal_len

    _n_int_eligible = len(int_eligible_seqs)
    log.info(f"  Internal seeding: {len(int_eligible_seqs):,} eligible of "
             f"{len(ltr_passing_indices):,} LTR-passing elements, "
             f"max_bucket={config.max_bucket_internal}, "
             f"min_shared={int_min_shared} (avg_len={avg_internal_len:.0f})")

    with timed("k-mer seeding + anchoring (internal)"):
        int_anchor_arr, int_n_raw, _int_n_unique, _hist, _hist_overflow = \
            kmer_lib.seed_and_anchor_v2(
                int_eligible_seqs, int_eligible_lens,
                kmer_size=config.kmer_size,
                max_bucket=config.max_bucket_internal,
                min_shared=int_min_shared,
                max_len_ratio=max_len_ratio_int,
                filter_n=0,
                filter_mode=0,
            )
    del int_eligible_seqs, int_eligible_lens

    if len(int_anchor_arr) > 0:
        local_remap = np.array(int_eligible_idx, dtype=np.int32)
        int_anchor_arr[:, 0] = local_remap[int_anchor_arr[:, 0]]
        int_anchor_arr[:, 1] = local_remap[int_anchor_arr[:, 1]]
        del local_remap

        global_remap = np.array(ltr_passing_indices, dtype=np.int32)
        int_anchor_arr[:, 0] = global_remap[int_anchor_arr[:, 0]]
        int_anchor_arr[:, 1] = global_remap[int_anchor_arr[:, 1]]
        del global_remap
    del int_eligible_idx

    # Filter to LTR-passing pairs 
    if len(int_anchor_arr) > 0:
        col0 = int_anchor_arr[:, 0].astype(np.int64)
        col1 = int_anchor_arr[:, 1].astype(np.int64)
        mn = np.minimum(col0, col1)
        mx = np.maximum(col0, col1)
        max_idx = int(mx.max()) + 1
        anchor_keys = mn * max_idx + mx
        del col0, col1, mn, mx
        ltr_keys_packed = np.array(
            [a * max_idx + b for a, b in ltr_key_set], dtype=np.int64)
        ltr_keys_packed.sort()
        keep = np.searchsorted(ltr_keys_packed, anchor_keys)
        keep = (keep < len(ltr_keys_packed)) & \
               (ltr_keys_packed[np.minimum(keep, len(ltr_keys_packed) - 1)] == anchor_keys)
        del anchor_keys, ltr_keys_packed
        n_before = len(int_anchor_arr)
        int_anchor_arr = int_anchor_arr[keep]
        del keep
        log.info(f"    {int_n_raw:,} raw events → {n_before:,} anchors "
                 f"→ {len(int_anchor_arr):,} after LTR-passing filter")

    int_anchored = int_anchor_arr
    del int_anchor_arr
    gc.collect()

    if isinstance(int_anchored, np.ndarray) and int_anchored.ndim == 2 and len(int_anchored) > 0:
        n_rc_int = int(int_anchored[:, 4].astype(bool).sum())
    else:
        n_rc_int = sum(1 for p in int_anchored if len(p) == 5 and p[4])
    log.step(f"  [kmer_seed_internal] {_n_int_eligible:,} eligible → "
             f"{len(int_anchored):,} anchored pairs")
    log.info(f"  → {len(int_anchored):,} internal anchored pairs "
             f"({n_rc_int:,} RC)")
    internal_edges: Dict[Tuple[int, int], float] = _wfa_validate(
        int_anchored, records, 'internal',
        config.min_internal_len, config, wfa_lib, use_anchored=True,
    )

    log.step(f"  [clustering_internal] {len(int_anchored):,} pairs → "
             f"{len(internal_edges):,} edges")

    # Add edges: LTR passed AND internal passed (strict intersection).
    # If either element lacks a usable internal region the pair is skipped.
    with timed("edge construction"):
        n_edges = 0
        for (a, b), ltr_id in ltr_passing.items():
            key = (min(a, b), max(a, b))
            ltr_div = round((1.0 - ltr_id) * 100.0, 1)
            if key in internal_edges:
                int_id  = internal_edges[key]
                int_div = round((1.0 - int_id) * 100.0, 1)
                G.add_edge(a, b, ltrpassed=1, internalpassed=1,
                           ltrdivergencepct=ltr_div, internaldivergencepct=int_div)
                n_edges += 1

    _add_reuse_edges(G, reuse_edges, name_to_idx)
    log.step(f"  Graph: {G.number_of_nodes():,} nodes, "
             f"{G.number_of_edges():,} edges")
    log.info(f"    {n_edges:,} of those edges are new-incident this round")
    return G


def _write_gml(G: nx.Graph, output_dir: Optional[Path]) -> None:
    """Write the similarity graph as a GML file with full node/edge attributes."""
    if output_dir is None:
        return
    gml_path = Path(output_dir) / 'similarity_graph.gml'
    with open(gml_path, 'w') as f:
        f.write('graph [\n')
        f.write('  directed 0\n')
        for node_id, attrs in G.nodes(data=True):
            f.write('  node [\n')
            f.write(f'    id {node_id}\n')
            f.write(f'    label "{attrs.get("label", attrs.get("name", str(node_id)))}"\n')
            f.write(f'    ltrlength {attrs.get("ltrlength", 0)}\n')
            f.write(f'    internallength {attrs.get("internallength", 0)}\n')
            f.write(f'    hasinternal {attrs.get("hasinternal", 0)}\n')
            f.write(f'    ltrageselfpct {attrs.get("ltrageselfpct", -1.0):.1f}\n')
            f.write(f'    iscentroid {attrs.get("iscentroid", 0)}\n')
            f.write(f'    centroiddivergencepct {attrs.get("centroiddivergencepct", -1.0):.1f}\n')
            f.write(f'    centroidwfafailed {attrs.get("centroidwfafailed", 0)}\n')
            f.write(f'    genome "{attrs.get("genome", "")}"\n')
            f.write('  ]\n')
        for src, dst, eattrs in G.edges(data=True):
            f.write('  edge [\n')
            f.write(f'    source {src}\n')
            f.write(f'    target {dst}\n')
            f.write(f'    ltrpassed {eattrs.get("ltrpassed", 1)}\n')
            f.write(f'    internalpassed {eattrs.get("internalpassed", -1)}\n')
            if 'ltrdivergencepct' in eattrs:
                f.write(f'    ltrdivergencepct {eattrs["ltrdivergencepct"]:.1f}\n')
            if 'internaldivergencepct' in eattrs:
                f.write(f'    internaldivergencepct {eattrs["internaldivergencepct"]:.1f}\n')
            f.write('  ]\n')
        f.write(']\n')
    log.info(f"  Similarity graph written to {gml_path}")


def write_classified_gml(
    gml_path:         str,
    membership_tsv:   str,
    classified_fasta: str,
    output_path:      str,
    current_round:    int = 1,
) -> None:
    """
    Write a copy of similarity_graph.gml annotated with per-node classifier labels.

    Reads the existing GML, cross-references cluster_membership.tsv and
    classified_fasta to assign a classifier_class to each node, then writes
    similarity_graph_classified.gml.

    Parameters
    ----------
    gml_path          Path to similarity_graph.gml produced by _write_gml().
    membership_tsv    Path to cluster_membership.tsv (consensus → members).
    classified_fasta  Path to consensi_classified.fasta (headers contain class=X).
    output_path       Destination for the annotated GML.
    current_round     Round number for elements that have no r{N}_ prefix.
    """
    G = nx.read_gml(gml_path, label='id')

    # element_name → consensus_name
    element_to_consensus: Dict[str, str] = {}
    with open(membership_tsv) as fh:
        next(fh)  # skip header
        for line in fh:
            cons_name, members_str = line.rstrip('\n').split('\t', 1)
            for elem in members_str.split(','):
                element_to_consensus[elem.strip()] = cons_name

    # consensus_name → class_name  (parsed from FASTA headers: "class=X")
    consensus_to_class: Dict[str, str] = {}
    with open(classified_fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                header = line[1:].strip()
                name   = header.split()[0]
                m      = re.search(r'class=(\S+)', header)
                if m:
                    consensus_to_class[name] = m.group(1)

    for node_id, attrs in G.nodes(data=True):
        elem_name = attrs.get('label', attrs.get('name', ''))
        cons      = element_to_consensus.get(elem_name, '')
        cls       = consensus_to_class.get(cons, 'unassigned')
        G.nodes[node_id]['classifier_class'] = cls
        # Parse round from r{N}_ prefix; plain names belong to current_round
        rnd = current_round
        if elem_name and '_' in elem_name and elem_name[0] == 'r':
            prefix = elem_name[1:elem_name.index('_')]
            if prefix.isdigit():
                rnd = int(prefix)
        G.nodes[node_id]['round'] = rnd

    # Hand-roll GML output matching _write_gml format, plus classifierclass
    with open(output_path, 'w') as f:
        f.write('graph [\n')
        f.write('  directed 0\n')
        for node_id, attrs in G.nodes(data=True):
            f.write('  node [\n')
            f.write(f'    id {node_id}\n')
            f.write(f'    label "{attrs.get("label", str(node_id))}"\n')
            f.write(f'    ltrlength {attrs.get("ltrlength", 0)}\n')
            f.write(f'    internallength {attrs.get("internallength", 0)}\n')
            f.write(f'    hasinternal {attrs.get("hasinternal", 0)}\n')
            f.write(f'    ltrageselfpct {attrs.get("ltrageselfpct", -1.0):.1f}\n')
            f.write(f'    iscentroid {attrs.get("iscentroid", 0)}\n')
            f.write(f'    centroiddivergencepct {attrs.get("centroiddivergencepct", -1.0):.1f}\n')
            f.write(f'    centroidwfafailed {attrs.get("centroidwfafailed", 0)}\n')
            f.write(f'    genome "{attrs.get("genome", "")}"\n')
            f.write(f'    classifierclass "{attrs.get("classifier_class", "unassigned")}"\n')
            f.write(f'    round {attrs.get("round", 1)}\n')
            f.write('  ]\n')
        for src, dst, eattrs in G.edges(data=True):
            f.write('  edge [\n')
            f.write(f'    source {src}\n')
            f.write(f'    target {dst}\n')
            f.write(f'    ltrpassed {eattrs.get("ltrpassed", 1)}\n')
            f.write(f'    internalpassed {eattrs.get("internalpassed", -1)}\n')
            if 'ltrdivergencepct' in eattrs:
                f.write(f'    ltrdivergencepct {eattrs["ltrdivergencepct"]:.1f}\n')
            if 'internaldivergencepct' in eattrs:
                f.write(f'    internaldivergencepct {eattrs["internaldivergencepct"]:.1f}\n')
            f.write('  ]\n')
        f.write(']\n')
    log.info(f"  Classifier-annotated graph written to {output_path}")

    # ── Filtered GML: only components with no "unassigned" or "non-LTR" nodes ─
    _EXCLUDE = {'unassigned', 'non-LTR'}
    keep_nodes = {
        n for n, attrs in G.nodes(data=True)
        if attrs.get('classifier_class', 'unassigned') not in _EXCLUDE
    }
    for comp in nx.connected_components(G):
        if comp & (set(G.nodes) - keep_nodes):
            keep_nodes -= comp
    G_ltr = G.subgraph(keep_nodes)

    ltr_path = str(Path(output_path).with_suffix('')) + '_ltr_only.gml'
    with open(ltr_path, 'w') as f:
        f.write('graph [\n')
        f.write('  directed 0\n')
        for node_id, attrs in G_ltr.nodes(data=True):
            f.write('  node [\n')
            f.write(f'    id {node_id}\n')
            f.write(f'    label "{attrs.get("label", str(node_id))}"\n')
            f.write(f'    ltrlength {attrs.get("ltrlength", 0)}\n')
            f.write(f'    internallength {attrs.get("internallength", 0)}\n')
            f.write(f'    hasinternal {attrs.get("hasinternal", 0)}\n')
            f.write(f'    ltrageselfpct {attrs.get("ltrageselfpct", -1.0):.1f}\n')
            f.write(f'    iscentroid {attrs.get("iscentroid", 0)}\n')
            f.write(f'    centroiddivergencepct {attrs.get("centroiddivergencepct", -1.0):.1f}\n')
            f.write(f'    centroidwfafailed {attrs.get("centroidwfafailed", 0)}\n')
            f.write(f'    genome "{attrs.get("genome", "")}"\n')
            f.write(f'    classifierclass "{attrs.get("classifier_class", "unassigned")}"\n')
            f.write(f'    round {attrs.get("round", 1)}\n')
            f.write('  ]\n')
        for src, dst, eattrs in G_ltr.edges(data=True):
            f.write('  edge [\n')
            f.write(f'    source {src}\n')
            f.write(f'    target {dst}\n')
            f.write(f'    ltrpassed {eattrs.get("ltrpassed", 1)}\n')
            f.write(f'    internalpassed {eattrs.get("internalpassed", -1)}\n')
            if 'ltrdivergencepct' in eattrs:
                f.write(f'    ltrdivergencepct {eattrs["ltrdivergencepct"]:.1f}\n')
            if 'internaldivergencepct' in eattrs:
                f.write(f'    internaldivergencepct {eattrs["internaldivergencepct"]:.1f}\n')
            f.write('  ]\n')
        f.write(']\n')
    log.info(f"  LTR-only graph ({G_ltr.number_of_nodes()} nodes) "
             f"written to {ltr_path}")


# ---------------------------------------------------------------------------
# Consensus generation for one cluster
# ---------------------------------------------------------------------------

def _centroid_idx(comp: List[int], G: nx.Graph) -> int:
    """Pick the highest-degree node as centroid; break ties by record index."""
    return max(range(len(comp)), key=lambda li: G.degree(comp[li]))


def _make_consensus_progressive(
    cluster_records:    List[LTRRecord],
    pairwise_identity:  Dict[Tuple[int, int], float],
    config:             FeatureClusterConfig,
    wfa_lib:            WFALib,
    kmer_lib:           KmerLib,
    debug_plot_path:    Optional[str] = None,
    label:              Optional[str] = None,
) -> Optional[List[Tuple[consensus.ConsensusResult, List[str], Optional[Dict]]]]:
    """Progressive consensus via UPGMA guide tree + profile merging.

    Returns a list of (ConsensusResult, member_names, member_stats) — one per
    subfamily.  A cluster splits into multiple subfamilies when progressive
    merges exceed the edit budget.
    """
    from algorithms.consensus import build_progressive_consensus
    try:
        sequences = [r.full_seq for r in cluster_records]
        all_names = [r.name for r in cluster_records]
        subfamilies = build_progressive_consensus(
            sequences,
            pairwise_identity=pairwise_identity,
            wfa_lib=wfa_lib,
            kmer_lib=kmer_lib,
            cons_threshold=config.cons_threshold,
            saturation_threshold=config.saturation_threshold,
            label=label,
            debug_plot_path=debug_plot_path,
            ltr_divergences=[r.ltr_divergence for r in cluster_records] if debug_plot_path else None,
            element_names=all_names if debug_plot_path else None,
            ambiguity_threshold=config.ambiguity_threshold,
            use_outgroup=config.use_outgroup,
            split_subfamilies=config.split_subfamilies,
            cut_threshold=config.cut_threshold,
            min_merge_identity=config.min_merge_identity,
            keep_alignments=config.keep_alignments,
        )
        if subfamilies is not None:
            results = []
            for cr, member_indices, member_stats in subfamilies:
                member_names = [all_names[i] for i in member_indices]
                if cr.alignment is not None:
                    # The alignment rows are keyed by index; this is the only
                    # place that knows the element names they belong to.
                    cr.alignment.row_names = [
                        all_names[i] for i in cr.alignment.row_indices
                    ]
                results.append((cr, member_names, member_stats))
            return results
    except Exception as e:
        log.warn(f"progressive consensus failed for cluster of "
                 f"{len(cluster_records)} elements: {e}")
    return None


# ---------------------------------------------------------------------------
# Process-pool worker helpers (each subprocess needs its own C library handles)
# ---------------------------------------------------------------------------

_worker_wfa_lib:  Optional[WFALib]  = None
_worker_kmer_lib: Optional[KmerLib] = None


def _worker_init(threads: int) -> None:
    """Initializer for ProcessPoolExecutor workers — load C libs once per process."""
    global _worker_wfa_lib, _worker_kmer_lib
    os.environ['OMP_NUM_THREADS'] = '1'  # each worker is one task; avoid oversubscription
    _worker_kmer_lib = KmerLib(_ALGO_DIR / 'kmer_seed.so')
    _worker_wfa_lib  = WFALib(_ALGO_DIR / 'wfa_align.so')


def _worker_make_consensus(
    cluster_records, config, pairwise_identity, debug_plot_path, label,
):
    """_make_consensus_progressive using per-process library handles."""
    return _make_consensus_progressive(
        cluster_records, pairwise_identity, config,
        _worker_wfa_lib, _worker_kmer_lib,
        debug_plot_path=debug_plot_path, label=label,
    )


def _annotate_consensus_nodes(
    G:            nx.Graph,
    comp:         List[int],
    cent_li:      int,
    member_stats: Optional[Dict],
) -> None:
    """
    Annotate graph nodes with per-member WFA alignment stats collected during
    consensus generation.

    Sets three node attributes:
      is_centroid             — 1 for the centroid node, 0 for all others
      centroid_divergence_pct — (1 - identity) * 100; -1 if WFA failed or
                                stats unavailable (MAFFT path)
      centroid_wfa_failed     — 1 if the member's WFA alignment returned rc < 0
    """
    centroid_global = comp[cent_li]
    for local_i, global_i in enumerate(comp):
        G.nodes[global_i]['iscentroid'] = int(global_i == centroid_global)
        if member_stats is not None and local_i in member_stats:
            identity, failed = member_stats[local_i]
            G.nodes[global_i]['centroiddivergencepct'] = (
                round((1.0 - identity) * 100.0, 1) if identity >= 0.0 else -1.0
            )
            G.nodes[global_i]['centroidwfafailed'] = int(failed)


# ---------------------------------------------------------------------------
# Consensus stage
# ---------------------------------------------------------------------------

def _consensus_name(cr, cluster_num: int, sub_tag: str) -> str:
    """The FASTA name for one consensus.

    Single source of the naming convention: `consensi.fasta`,
    `cluster_membership.tsv` and the per-consensus alignment files all take
    their name from here, so a reader can join them by name.
    """
    return (
        f"CONS_{cluster_num}{sub_tag}-{len(cr.sequence)}"
        f"_clus{cr.num_sequences}"
        f"_tsdl{cr.tsd_length}"
        f"_tsdc{cr.tsd_confidence}"
        f"_tsdm{cr.tsd_motif if cr.tsd_motif else 'none'}"
    )


def _run_consensus_and_write(
    components: List[List[int]],
    records:    List[LTRRecord],
    G:          nx.Graph,
    config:     FeatureClusterConfig,
    wfa_lib:    WFALib,
    out:        Path,
    kmer_lib:   Optional[KmerLib] = None,
) -> Tuple[str, str, list]:
    """
    Build one consensus per cluster, write consensi.fasta and
    cluster_membership.tsv, and return (consensi_path, membership_path,
    all_results).
    """
    record_by_idx = {i: r for i, r in enumerate(records)}
    all_results: List[Tuple[consensus.ConsensusResult, List[str], int, str]] = []

    tasks = []
    for comp in components:
        cluster_recs = [record_by_idx[idx] for idx in comp if idx in record_by_idx]
        if len(cluster_recs) < config.min_cluster_size:
            continue
        cent_li = _centroid_idx(comp, G)   # local index within comp

        # For progressive consensus, extract local pairwise identities from
        # graph edges. 
        pw_id: Optional[Dict[Tuple[int, int], float]] = None
        g2l = {g: l for l, g in enumerate(comp)}
        pw_id = {}
        sub = G.subgraph(comp)
        for u, v, data in sub.edges(data=True):
            li, lj = g2l[u], g2l[v]
            key = (min(li, lj), max(li, lj))
            ltr_div = data.get('ltrdivergencepct', 20.0)
            int_div = data.get('internaldivergencepct', ltr_div)
            ident = 1.0 - (ltr_div + int_div) / 200.0
            pw_id[key] = ident
        tasks.append((cluster_recs, cent_li, comp, pw_id))

    # Set up debug tree plot directory (progressive + debug only)
    tree_dir: Optional[Path] = None
    if config.debug:
        tree_dir = out / 'debug' / 'cluster_trees'
        tree_dir.mkdir(parents=True, exist_ok=True)

    n_tasks    = len(tasks)
    method_tag = 'PROGRESSIVE'

    n_workers    = max(1, config.threads)
    use_parallel = n_workers > 1
    parallel_tag = f"{n_workers} workers" if use_parallel else "serial"
    log.step(f"  Building {n_tasks} consensi "
             f"(method={method_tag}, {parallel_tag}) ...")

    t_consensus = time.perf_counter()
    completed   = 0

    def _collect_subfamilies(subfamilies, task_i, comp, cent_li):
        """Append subfamily results with cluster/subfamily labeling."""
        for sub_i, (cr, member_ids, member_stats) in enumerate(subfamilies):
            # Tag subfamily index so FASTA naming can distinguish them
            sub_tag = _subfamily_tag(sub_i) if len(subfamilies) > 1 else ''
            all_results.append((cr, member_ids, task_i + 1, sub_tag))
            _annotate_consensus_nodes(G, comp, cent_li, member_stats)

    if use_parallel:
        with ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_worker_init,
            initargs=(1,),
        ) as pool:
            futures = {}
            for task_i, (recs, cent, _comp, pw) in enumerate(tasks):
                plot_path = str(tree_dir / f'cluster_{task_i+1}.png') if tree_dir else None
                label = str(task_i + 1)
                futures[pool.submit(
                    _worker_make_consensus, recs, config, pw,
                    plot_path, label,
                )] = task_i
            # process in task order so consensus IDs are deterministic.
            results_by_task = [None] * n_tasks
            for fut in as_completed(futures):
                completed += 1
                if n_tasks > 0 and (completed == n_tasks or
                                    completed * 10 // n_tasks > (completed - 1) * 10 // n_tasks):
                    log.progress(f"  {completed}/{n_tasks} consensi done")
                task_i = futures[fut]
                results_by_task[task_i] = fut.result()
            for task_i, subfamilies in enumerate(results_by_task):
                if subfamilies is not None:
                    _recs, cent_li, comp, _pw = tasks[task_i]
                    _collect_subfamilies(subfamilies, task_i, comp, cent_li)
    else:
        for task_i, (recs, cent, comp, pw) in enumerate(tasks):
            completed += 1
            if n_tasks > 0 and (completed == n_tasks or
                                completed * 10 // n_tasks > (completed - 1) * 10 // n_tasks):
                log.progress(f"  {completed}/{n_tasks} consensi done")
            plot_path = str(tree_dir / f'cluster_{task_i+1}.png') if tree_dir else None
            label = str(task_i + 1)
            subfamilies = _make_consensus_progressive(
                recs, pw, config, wfa_lib, kmer_lib,
                debug_plot_path=plot_path, label=label)
            if subfamilies is not None:
                _collect_subfamilies(subfamilies, task_i, comp, cent)

    log.progress_end()
    log.step(f"  Generated {len(all_results)} consensus sequences")
    log.info(f"  [timing] build_{method_tag.lower()}_consensus ({n_tasks} clusters, "
             f"{parallel_tag}): "
             f"{time.perf_counter() - t_consensus:.2f}s")

    consensi_path   = str(out / 'consensi.fasta')
    membership_path = str(out / 'cluster_membership.tsv')

    aln_dir = None
    if config.keep_alignments:
        from utilities.alignment_output import write_cluster_alignment
        aln_dir = out / 'alignments'
        aln_dir.mkdir(parents=True, exist_ok=True)
    n_aln = n_aln_missing = 0

    with open(consensi_path, 'w') as fa, \
            open(membership_path, 'w') as mem:
        mem.write("consensus\tmembers\n")
        for cr, member_ids, cluster_num, sub_tag in all_results:
            cons_name = _consensus_name(cr, cluster_num, sub_tag)
            fa.write(f">{cons_name}\n{cr.sequence}\n")
            mem.write(f"{cons_name}\t{','.join(member_ids)}\n")
            if aln_dir is not None and cr.alignment is not None:
                write_cluster_alignment(cr.alignment, cons_name, aln_dir)
                n_aln += 1
                n_aln_missing += len(member_ids) - len(cr.alignment.rows)

    log.info(f"  Cluster membership: {membership_path}")
    log.info(f"  Consensus library: {consensi_path}")
    if aln_dir is not None:
        log.step(f"  Member alignments: {n_aln} written to {aln_dir}")
        if n_aln_missing:
            log.info(f"    {n_aln_missing} member(s) omitted — no alignment "
                     f"to their consensus")
    return consensi_path, membership_path, all_results


def cluster_records_incremental(
    records:   List[LTRRecord],
    config:    FeatureClusterConfig,
    kmer_lib:  KmerLib,
    wfa_lib:   WFALib,
    out:       Path,
    *,
    new_names:     Optional[Set[str]] = None,
    reuse_edges:   Optional[Dict[Tuple[str, str], dict]] = None,
    frozen_params: Optional[dict]     = None,
) -> Tuple[str, str, list, nx.Graph]:
    """Cluster a set of LTRRecords into families; write consensi + membership. """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)

    G = _build_graph(records, config, kmer_lib, wfa_lib,
                     new_names=new_names, reuse_edges=reuse_edges,
                     frozen_params=frozen_params)

    with timed("connected components"):
        components = sorted(
            [sorted(comp)
             for comp in nx.connected_components(G)
             if len(comp) >= config.min_cluster_size],
            key=lambda c: c[0],
        )
    log.info(f"  Found {len(components)} families "
             f"(>= {config.min_cluster_size} elements)")

    consensi_path, membership_path, all_results = _run_consensus_and_write(
        components, records, G, config, wfa_lib, out, kmer_lib=kmer_lib,
    )
    return consensi_path, membership_path, all_results, G
