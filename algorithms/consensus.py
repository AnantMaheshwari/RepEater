"""
Consensus generation by progressive profile merging 
(TSD/conservation logic is inspired by pantera).

`build_progressive_consensus` builds a UPGMA guide tree
over the cluster, merges pairs of column profiles along it (optionally splitting
the tree into subfamilies and refining ambiguous columns against an outgroup),
and then for each finished profile:
  1. Generate the majority-rule consensus string
  2. Compute conservation and saturation per column
  3. Find conserved boundaries (adaptive window size)
  4. Trim and filter low-occupancy columns
  5. Detect TSD from the alignment edges
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

import numpy as np
from collections import defaultdict

from engine.constants import (
    AMBIGUITY_THRESHOLD, CONS_THRESHOLD, CONSENSUS_ANCHOR_MAX_BUCKET,
    CUT_THRESHOLD, KMER_SIZE, MERGE_MAX_EDIT_FRAC, MIN_MERGE_IDENTITY,
    SATURATION_THRESHOLD,
)
from utilities import log

if TYPE_CHECKING:
    from algorithms.c_libs import WFALib

from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


# IUPAC ambiguity codes for pairs of bases
_IUPAC = {
    frozenset('AG'): 'R', frozenset('CT'): 'Y', frozenset('GC'): 'S',
    frozenset('AT'): 'W', frozenset('GT'): 'K', frozenset('AC'): 'M',
    frozenset('CGT'): 'B', frozenset('AGT'): 'D', frozenset('ACT'): 'H',
    frozenset('ACG'): 'V', frozenset('ACGT'): 'N',
}

# Row indices in the count matrix
_BASE_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
_IDX_BASE = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}


@dataclass
class ClusterAlignment:
    """The alignment of a subfamily's members against its consensus.
    """
    row_indices:   List[int]
    rows:          List[str]
    consensus_row: str
    row_names:     List[str] = None


@dataclass
class ConsensusResult:
    """Result of consensus generation for one cluster."""
    sequence: str
    num_sequences: int
    tsd_length: int
    tsd_confidence: float
    tsd_motif: str
    raw_alignment_length: int
    alignment: Optional[ClusterAlignment] = None


def consensus_string(
    count_matrix: np.ndarray,
    n_seqs: int,
    threshold: float = CONS_THRESHOLD,
    ambiguity_threshold: float = AMBIGUITY_THRESHOLD,
) -> str:
    """
    Generate a consensus string from a count matrix.
    
    """
    aln_len = count_matrix.shape[1]
    bases_matrix = count_matrix[:4]  # A, C, G, T rows
    result = []

    for j in range(aln_len):
        total = bases_matrix[:, j].sum()
        if total == 0:
            result.append('-')
            continue
        freqs = bases_matrix[:, j] / total
        above = np.where(freqs >= threshold)[0]
        if len(above) == 1:
            result.append(_IDX_BASE[above[0]])
        elif len(above) > 1:
            # Multiple bases above threshold — use IUPAC ambiguity
            base_set = frozenset(_IDX_BASE[i] for i in above)
            result.append(_IUPAC.get(base_set, 'N'))
        elif ambiguity_threshold > 0:
            # No single base dominates — try the lower ambiguity bar
            above_amb = np.where(freqs >= ambiguity_threshold)[0]
            if len(above_amb) > 1:
                base_set = frozenset(_IDX_BASE[i] for i in above_amb)
                result.append(_IUPAC.get(base_set, 'N'))
            elif len(above_amb) == 1:
                result.append(_IDX_BASE[above_amb[0]])
            else:
                result.append('N')
        else:
            result.append('N')

    return ''.join(result)


def compute_conservation(
    count_matrix: np.ndarray, n_seqs: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute per-column conservation and saturation.
    """
    bases = count_matrix[:4]  # A, C, G, T
    conservation = bases.max(axis=0) / n_seqs
    saturation = bases.sum(axis=0) / n_seqs
    return conservation, saturation


def find_conserved_boundaries(
    conservation: np.ndarray,
    threshold: float = 0.8,
    n_seqs: int = 3,
) -> Optional[Tuple[int, int]]:
    """
    Find the first and last conserved runs in the conservation array.
    Returns (start, end) indices (inclusive), or None if no conserved run found.
    """
    min_run = max(2, 12 - math.floor(math.sqrt(max(0, n_seqs - 2))))
    binary = (conservation >= threshold).astype(np.int8)
    binary_str = ''.join(str(b) for b in binary)
    target = '1' * min_run

    first = binary_str.find(target)
    if first == -1:
        return None
    last = binary_str.rfind(target)
    # end is the last position covered by the rightmost run
    end = last + min_run - 1

    return first, end


def trim_and_filter(
    cons: str,
    saturation: np.ndarray,
    start: int,
    end: int,
    saturation_threshold: float = SATURATION_THRESHOLD,
) -> str:
    """
    Trim consensus to conserved boundaries, remove low-occupancy columns,
    and replace gaps with N.

    Returns the final consensus sequence.
    """
    trimmed = cons[start:end + 1]
    sat_slice = saturation[start:end + 1]
    filtered = ''.join(
        base for base, sat in zip(trimmed, sat_slice)
        if sat >= saturation_threshold
    )
    return filtered.replace('-', 'N')


def detect_tsd_from_alignment(
    aligned_seqs: List[str],
    cons_start: int,
    cons_end: int,
    max_tsd_len: int = 13,
) -> Tuple[int, float, str]:
    """
    Detect TSD (Target Site Duplication) from alignment edges.
    For each candidate TSD length k (2..max_tsd_len), count how many
    sequences have right[:k] == reverse(left[:k]).
    """
    edges = []
    for seq in aligned_seqs:
        left = seq[:cons_start].replace('-', '')[::-1]  # reversed, no gaps
        right = seq[cons_end + 1:].replace('-', '')       # no gaps
        if left and right:
            edges.append((left.upper(), right.upper()))

    if not edges:
        return 0, 0.0, ''

    n = len(edges)
    counts = []
    for k in range(2, max_tsd_len + 1):
        matches = sum(
            1 for left, right in edges
            if len(left) >= k and len(right) >= k
            and right[:k] == left[:k][::-1]
        )
        counts.append(matches)

    if not counts or max(counts) == 0:
        return 0, 0.0, ''

    max_count = max(counts)
    # Smallest k with maximum count
    tsd_idx = counts.index(max_count)
    tsd_len = tsd_idx + 2  # offset because range starts at k=2

    tsd_conf = round(max_count / n, 2)

    # Build consensus motif from matching TSDs (IUPAC)
    motifs = []
    for left, right in edges:
        if len(right) >= tsd_len and len(left) >= tsd_len:
            if right[:tsd_len] == left[:tsd_len][::-1]:
                motifs.append(right[:tsd_len])
    tsd_motif = _iupac_consensus(motifs) if motifs else ''

    return tsd_len, tsd_conf, tsd_motif


def _iupac_consensus(motifs: List[str]) -> str:
    """Build an IUPAC consensus string from a list of equal-length motifs."""
    if not motifs:
        return ''
    length = len(motifs[0])
    result = []
    for i in range(length):
        bases = set(m[i] for m in motifs if i < len(m))
        bases.discard('N')
        if not bases:
            result.append('N')
        elif len(bases) == 1:
            result.append(bases.pop())
        else:
            code = _IUPAC.get(frozenset(bases), 'N')
            result.append(code)
    return ''.join(result)


_COMP = str.maketrans('ACGTacgtNnRrYyKkMmSsWwBbDdHhVv',
                      'TGCAtgcaNnYyRrMmKkSsWwVvHhDdBb')


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


_KMER_SIZE = KMER_SIZE


# ------------------------------------------------------------------------------------
# Some ideas around centering the anchor for consensus seeding to get cleaner consensi
# TODO (anant): quantify how much this happens, we see anchors chosen near endpoints
# resulting in consensus artefacts, so this is an approach to try to choose an anchor
# during consensus alignments that is close to the middle of sequences being compared 
# ------------------------------------------------------------------------------------

def _lis_chain(pos_i: np.ndarray, pos_j: np.ndarray) -> np.ndarray:
    import bisect

    n = len(pos_i)
    if n == 0:
        return np.array([], dtype=np.intp)

    order = np.lexsort((pos_j, pos_i))
    pj_sorted = pos_j[order]

    tails = []
    tail_idx = []
    parent = [-1] * n

    for k in range(n):
        val = int(pj_sorted[k])
        pos = bisect.bisect_left(tails, val)
        if pos == len(tails):
            tails.append(val)
            tail_idx.append(k)
        else:
            tails[pos] = val
            tail_idx[pos] = k
        parent[k] = tail_idx[pos - 1] if pos > 0 else -1

    if len(tails) == 0:
        return np.array([], dtype=np.intp)

    chain_sorted_idx = []
    idx = tail_idx[-1]
    while idx >= 0:
        chain_sorted_idx.append(idx)
        idx = parent[idx]
    chain_sorted_idx.reverse()

    return order[np.array(chain_sorted_idx, dtype=np.intp)]


def _anchor_chain_and_orient(
    cons_a: str, cons_b: str, kmer_lib,
) -> Tuple[List[int], List[int], bool]:
    """
    Get a collinear chain of k-mer anchors between two consensus sequences.
    """
    la, lb = len(cons_a), len(cons_b)
    seqs = [cons_a.upper().encode(), cons_b.upper().encode()]
    lens = [la, lb]

    def _ret(anch_a, anch_b, rc):
        return anch_a, anch_b, rc

    events = kmer_lib.emit_pairs_pos(
        seqs, lens,
        kmer_size=_KMER_SIZE, max_bucket=CONSENSUS_ANCHOR_MAX_BUCKET,
    )
    if events.shape[0] == 0:
        return _ret([], [], False)

    mask = (events[:, 0] == 0) & (events[:, 1] == 1)
    ev = events[mask]
    if len(ev) == 0:
        return _ret([], [], False)

    pi = ev[:, 2].astype(np.int64)
    pj = ev[:, 3].astype(np.int64)

    diags = pi - pj
    anti_diags = pi + pj
    var_d = float(np.var(diags)) if len(diags) > 1 else 0.0
    var_ad = float(np.var(anti_diags)) if len(anti_diags) > 1 else 0.0
    is_rc = var_ad < var_d

    if is_rc:
        pj_oriented = lb - pj - _KMER_SIZE
        lis_idx = _lis_chain(pi, pj_oriented)
        if len(lis_idx) == 0:
            return _ret([], [], True)
        chain_pi = pi[lis_idx]
        chain_pj = pj_oriented[lis_idx]
    else:
        lis_idx = _lis_chain(pi, pj)
        if len(lis_idx) == 0:
            return _ret([], [], False)
        chain_pi = pi[lis_idx]
        chain_pj = pj[lis_idx]


    diag_vals = chain_pi.astype(np.int64) - chain_pj.astype(np.int64)
    min_diag_shift = 20 

    keep = []
    last_diag = None
    for k in range(len(lis_idx)):
        if last_diag is None or abs(int(diag_vals[k]) - last_diag) >= min_diag_shift:
            keep.append(k)
            last_diag = int(diag_vals[k])

    def _pick_central():
        best_k, best_score = 0, -1.0
        for k in range(len(lis_idx)):
            pa, pb = int(chain_pi[k]), int(chain_pj[k])
            frac_a = min(pa, la - pa) / la if la > 0 else 0
            frac_b = min(pb, lb - pb) / lb if lb > 0 else 0
            score = min(frac_a, frac_b)
            if score > best_score:
                best_score = score
                best_k = k
        return [int(chain_pi[best_k])], [int(chain_pj[best_k])]

    if len(keep) <= 1:
        anch_a, anch_b = _pick_central()
        return _ret(anch_a, anch_b, is_rc)

    min_end_bp = max(int(0.10 * min(la, lb)), 200)

    while len(keep) > 1 and (int(chain_pi[keep[0]]) < min_end_bp
                              or int(chain_pj[keep[0]]) < min_end_bp):
        keep.pop(0)
    while len(keep) > 1 and (la - int(chain_pi[keep[-1]]) < min_end_bp
                              or lb - int(chain_pj[keep[-1]]) < min_end_bp):
        keep.pop()

    if len(keep) <= 1:
        # After trimming endpoints, only 0-1 anchors remain — use central
        anch_a, anch_b = _pick_central()
        return _ret(anch_a, anch_b, is_rc)

    central_a, central_b = _pick_central()
    anchors_a = central_a + [int(chain_pi[k]) for k in keep]
    anchors_b = central_b + [int(chain_pj[k]) for k in keep]

    return _ret(anchors_a, anchors_b, is_rc)


def _anchor_and_orient(cons_a: str, cons_b: str, kmer_lib,
                       prefer_central: bool = False) -> tuple:
    """
    Return (anchor_a, anchor_b, is_rc) for two sequences using k-mer seeding.
    """
    la, lb = len(cons_a), len(cons_b)
    seqs = [cons_a.upper().encode(), cons_b.upper().encode()]
    lens = [la, lb]

    n_filter = 10 if prefer_central else 1
    anchors, _, _, _, _ = kmer_lib.seed_and_anchor_v2(
        seqs, lens,
        kmer_size=_KMER_SIZE, max_bucket=CONSENSUS_ANCHOR_MAX_BUCKET,
        min_shared=1, filter_n=n_filter,
    )
    if anchors.shape[0] == 0:
        return -1, -1, False

    if prefer_central and anchors.shape[0] > 1:
        mid_a, mid_b = la / 2.0, lb / 2.0
        best_idx = 0
        best_score = -1.0
        for ri in range(anchors.shape[0]):
            pa, pb = int(anchors[ri, 2]), int(anchors[ri, 3])
            is_rc_i = bool(anchors[ri, 4])
            if is_rc_i:
                pb = lb - pb - _KMER_SIZE
            if pb < 0 or pa < 0:
                continue
            frac_a = min(pa, la - pa) / la if la > 0 else 0
            frac_b = min(pb, lb - pb) / lb if lb > 0 else 0
            score = min(frac_a, frac_b)
            if score > best_score:
                best_score = score
                best_idx = ri
        row = anchors[best_idx]
    else:
        row = anchors[0]

    ap_a, ap_b, is_rc = int(row[2]), int(row[3]), bool(row[4])
    if is_rc:
        ap_b = lb - ap_b - _KMER_SIZE
        if ap_b < 0:
            return ap_a, -1, is_rc
    return ap_a, ap_b, is_rc


# ---------------------------------------------------------------------------
# CIGAR helpers
# ---------------------------------------------------------------------------

def _iter_cigar(cigar: str):
    """Yield (count, op) pairs from an RLE CIGAR string like '502M1D198M'."""
    i, n = 0, len(cigar)
    while i < n:
        j = i
        while j < n and cigar[j].isdigit():
            j += 1
        yield int(cigar[i:j]), cigar[j]
        i = j + 1


# Lookup table: ASCII byte → row index (A=0 C=1 G=2 T=3, else 4)
_BYTE_TO_ROW = np.full(256, 4, dtype=np.int8)
for _b, _r in ((ord('A'), 0), (ord('C'), 1), (ord('G'), 2), (ord('T'), 3),
               (ord('a'), 0), (ord('c'), 1), (ord('g'), 2), (ord('t'), 3)):
    _BYTE_TO_ROW[_b] = _r


def _project_cigar(member: str, cigar: str, n_cols: int) -> str:
    """
    Build a projected aligned sequence of length n_cols from the CIGAR.
    """
    member_b = member.encode('ascii')
    member_len = len(member_b)
    projected = bytearray(b'-') * n_cols
    c_pos = m_pos = 0
    for n, op in _iter_cigar(cigar):
        if op == 'M':
            actual = min(n, n_cols - c_pos, member_len - m_pos)
            if actual > 0:
                projected[c_pos:c_pos + actual] = member_b[m_pos:m_pos + actual]
            c_pos += n
            m_pos += n
        elif op == 'D':
            c_pos += n
        elif op == 'I':
            m_pos += n
    return projected.decode('ascii')


def _insert_cap(ins: List[int], budget: int) -> int:
    lo, hi = 0, max(ins)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if sum(v if v < mid else mid for v in ins) <= budget:
            lo = mid
        else:
            hi = mid - 1
    return lo


def _align_members_to_consensus(
    cons_seq:      str,
    member_indices: List[int],
    sequences:     List[str],
    wfa_lib:       'WFALib',
    kmer_lib:      'KmerLib',
    budget_frac:   float = 0.50,
) -> List[Tuple[int, str, str]]:
    """Align each member to the consensus, for --keep-alignments.
    """
    ref_b   = cons_seq.encode()
    ref_len = len(cons_seq)
    out: List[Tuple[int, str, str]] = []
    for idx in member_indices:
        seq = sequences[idx]
        ap_r, ap_m, is_rc = _anchor_and_orient(cons_seq, seq.upper(), kmer_lib)
        member = reverse_complement(seq.upper()) if is_rc else seq.upper()
        budget = int(budget_frac * max(ref_len, len(member))) + 1
        if ap_r >= 0 and ap_m >= 0:
            rc, _id, cigar = wfa_lib.align_cigar_anchored(
                ref_b, member.encode(), ap_r, ap_m,
                max_edit=budget, min_coverage=0.0,
            )
        else:
            rc, _id, cigar = wfa_lib.align_cigar(
                ref_b, member.encode(), max_edit=budget, min_coverage=0.0,
            )
        if rc >= 0 and cigar:
            out.append((idx, member, cigar))
    return out


def _expand_alignment(
    cons_seq: str,
    captured: List[Tuple[int, str, str]],
) -> 'ClusterAlignment':
    """Build a true MSA of cluster members from their CIGARs against *cons_seq*.
    """
    ref_len = len(cons_seq)

    ins = [0] * (ref_len + 1)
    for _idx, member, cigar in captured:
        c_pos = 0
        for n, op in _iter_cigar(cigar):
            if op == 'M':
                c_pos += n
            elif op == 'D':
                c_pos += n
            elif op == 'I':
                p = c_pos if c_pos <= ref_len else ref_len
                if n > ins[p]:
                    ins[p] = n

    total_ins = sum(ins)
    if total_ins > ref_len:
        cap = _insert_cap(ins, ref_len)
        ins = [v if v < cap else cap for v in ins]
        log.detail(f"    alignment: truncated insertions to {cap} col(s) "
                   f"({total_ins} -> {sum(ins)} insert columns, ref {ref_len})")


    width = ref_len + sum(ins)
    ref_col_of   = [0] * ref_len
    ins_start_of = [0] * (ref_len + 1)
    col = 0
    for p in range(ref_len):
        ins_start_of[p] = col
        col += ins[p]
        ref_col_of[p] = col
        col += 1
    ins_start_of[ref_len] = col

    def _blank() -> bytearray:
        return bytearray(b'-') * width

    # Pass 2: place each member's bases.
    rows: List[str] = []
    row_indices: List[int] = []
    for idx, member, cigar in captured:
        member_b   = member.encode('ascii')
        member_len = len(member_b)
        row = _blank()
        c_pos = m_pos = 0
        for n, op in _iter_cigar(cigar):
            if op == 'M':
                # Same defensive clamp as _project_cigar: CIGARs can overrun.
                actual = min(n, ref_len - c_pos, member_len - m_pos)
                for k in range(actual):
                    row[ref_col_of[c_pos + k]] = member_b[m_pos + k]
                c_pos += n
                m_pos += n
            elif op == 'D':
                c_pos += n          # gap in the member; row already '-'
            elif op == 'I':
                p     = c_pos if c_pos <= ref_len else ref_len
                start = ins_start_of[p]
                actual = min(n, ins[p], member_len - m_pos)
                if actual > 0:
                    row[start:start + actual] = member_b[m_pos:m_pos + actual]
                m_pos += n
        rows.append(row.decode('ascii'))
        row_indices.append(idx)

    cons_row = _blank()
    cons_b   = cons_seq.encode('ascii')
    for p in range(ref_len):
        cons_row[ref_col_of[p]] = cons_b[p]

    return ClusterAlignment(
        row_indices=row_indices,
        rows=rows,
        consensus_row=cons_row.decode('ascii'),
    )


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------

def _rc_profile(profile: np.ndarray) -> np.ndarray:
    """Reverse-complement a count matrix: reverse columns, swap A<->T, C<->G."""
    rc = np.empty_like(profile)
    rc[0] = profile[3, ::-1]  # A <- T reversed
    rc[1] = profile[2, ::-1]  # C <- G reversed
    rc[2] = profile[1, ::-1]  # G <- C reversed
    rc[3] = profile[0, ::-1]  # T <- A reversed
    rc[4] = profile[4, ::-1]  # gap row reversed
    return rc


_IDX_BASE_ARR = np.array([_IDX_BASE[i] for i in range(4)] + ['N'], dtype='U1')


def _profile_consensus(profile: np.ndarray,
                       threshold: float = CONS_THRESHOLD) -> str:
    """Quick majority-vote consensus from a profile (no IUPAC, just best base)."""
    n_cols = profile.shape[1]
    if n_cols == 0:
        return ''
    bases = profile[:4]
    totals = bases.sum(axis=0)
    best = bases.argmax(axis=0)       # index of highest-count base per column
    best[totals == 0] = 4             # N for zero-coverage columns
    return ''.join(_IDX_BASE_ARR[best])


def _merge_profiles(
    prof_a:  np.ndarray,
    n_a:     int,
    prof_b:  np.ndarray,
    n_b:     int,
    cigar:   str,
) -> np.ndarray:
    """
    Merge two count-matrix profiles through a CIGAR alignment.

    prof_a : (5, La) — profile for subtree A (n_a sequences)
    prof_b : (5, Lb) — profile for subtree B (n_b sequences)
    cigar  : alignment of consensus_a vs consensus_b
    """
    la = prof_a.shape[1]
    lb = prof_b.shape[1]

    # First pass: compute output length from CIGAR
    out_len = 0
    a_pos = b_pos = 0
    ops = list(_iter_cigar(cigar))
    for count, op in ops:
        if op == 'M':
            n = min(count, la - a_pos, lb - b_pos)
            out_len += max(n, 0)
            a_pos += count
            b_pos += count
        elif op == 'D':
            n = min(a_pos + count, la) - a_pos
            out_len += max(n, 0)
            a_pos += count
        elif op == 'I':
            n = min(b_pos + count, lb) - b_pos
            out_len += max(n, 0)
            b_pos += count

    if out_len == 0:
        return np.zeros((5, 0), dtype=np.int32)

    # Second pass: fill pre-allocated output
    result = np.empty((5, out_len), dtype=np.int32)
    a_pos = b_pos = o_pos = 0
    for count, op in ops:
        if op == 'M':
            n = min(count, la - a_pos, lb - b_pos)
            if n > 0:
                np.add(prof_a[:, a_pos:a_pos + n],
                       prof_b[:, b_pos:b_pos + n],
                       out=result[:, o_pos:o_pos + n])
                o_pos += n
            a_pos += count
            b_pos += count
        elif op == 'D':
            end_a = min(a_pos + count, la)
            n = end_a - a_pos
            if n > 0:
                result[:, o_pos:o_pos + n] = prof_a[:, a_pos:end_a]
                result[4, o_pos:o_pos + n] += n_b
                o_pos += n
            a_pos += count
        elif op == 'I':
            end_b = min(b_pos + count, lb)
            n = end_b - b_pos
            if n > 0:
                result[:, o_pos:o_pos + n] = prof_b[:, b_pos:end_b]
                result[4, o_pos:o_pos + n] += n_b
                o_pos += n
            b_pos += count
    return result[:, :o_pos]


def _outgroup_refine(
    merged_profile: np.ndarray,
    outgroup_cons:  str,
    merged_cons:    str,
    wfa_lib,
    kmer_lib,
    n_merged:       int,
    boost:          float = 1.0,
) -> np.ndarray:
    """
    Refine ambiguous columns in *merged_profile* using an outgroup consensus.

    For each M-column in the alignment of merged_cons vs outgroup_cons,
    if the top-2 bases in merged_profile are within 1.5x of each other
    (genuinely ambiguous), and the outgroup base matches one of the two
    candidates, add *boost* votes for that base.
    """
    # Align outgroup to the merged consensus
    ap_a, ap_b, is_rc = _anchor_and_orient(merged_cons, outgroup_cons, kmer_lib)
    og_oriented = reverse_complement(outgroup_cons) if is_rc else outgroup_cons

    budget = int(0.40 * max(len(merged_cons), len(og_oriented))) + 1
    if ap_a >= 0 and ap_b >= 0:
        rc, _id, cigar = wfa_lib.align_cigar_anchored(
            merged_cons.encode(), og_oriented.encode(), ap_a, ap_b,
            max_edit=budget, min_coverage=0.0,
        )
    else:
        rc, _id, cigar = wfa_lib.align_cigar(
            merged_cons.encode(), og_oriented.encode(),
            max_edit=budget, min_coverage=0.0,
        )
    if rc < 0:
        return merged_profile 

    _BASE_ORD = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3,
                 ord('a'): 0, ord('c'): 1, ord('g'): 2, ord('t'): 3}
    result = merged_profile.copy()
    m_pos = 0
    og_pos = 0
    n_cols = result.shape[1]
    n_refined = 0

    for count, op in _iter_cigar(cigar):
        if op == 'M':
            for _ in range(count):
                if m_pos >= n_cols:
                    break
                col = result[:4, m_pos]
                total = col.sum()
                if total < 2:
                    m_pos += 1
                    og_pos += 1
                    continue
                sorted_idx = np.argsort(col)[::-1]
                top_count = col[sorted_idx[0]]
                second_count = col[sorted_idx[1]]
                if second_count > 0 and top_count < 1.5 * second_count:
                    # Ambiguous — consult outgroup
                    og_byte = ord(og_oriented[og_pos]) if og_pos < len(og_oriented) else -1
                    og_base = _BASE_ORD.get(og_byte, -1)
                    if og_base >= 0 and og_base in (sorted_idx[0], sorted_idx[1]):
                        result[og_base, m_pos] += int(boost)
                        n_refined += 1
                m_pos += 1
                og_pos += 1
        elif op == 'D':
            m_pos += count
        elif op == 'I':
            og_pos += count

    return result


def _finalize_subfamily(
    profile:             np.ndarray,
    n_members:           int,
    cons_str:            str,
    member_indices:      List[int],
    sequences:           List[str],
    wfa_lib:             'WFALib',
    kmer_lib:            'KmerLib',
    cons_threshold:      float,
    saturation_threshold: float,
    ambiguity_threshold: float,
    keep_alignments:     bool = False,
) -> Optional[Tuple[ConsensusResult, List[int], Dict[int, Tuple[float, bool]]]]:
    """Post-process a subfamily profile into a ConsensusResult.
    """
    if len(member_indices) < 2:
        return None

    cons = consensus_string(profile, n_members, threshold=cons_threshold,
                            ambiguity_threshold=ambiguity_threshold)
    conservation, saturation = compute_conservation(profile, n_members)

    boundaries = find_conserved_boundaries(
        conservation,
        threshold=cons_threshold + 0.4,
        n_seqs=n_members,
    )
    if boundaries is None:
        return None

    cons_start, cons_end = boundaries

    # TSD detection: project each member against this subfamily's consensus.
    ref_str = cons_str
    ref_bytes = ref_str.encode()
    ref_len = len(ref_str)
    tsd_budget_frac = 0.25
    aligned_seqs: List[str] = [ref_str]
    member_stats: Dict[int, Tuple[float, bool]] = {}
    for idx in member_indices:
        seq = sequences[idx]
        ap_r, ap_m, is_rc = _anchor_and_orient(ref_str, seq.upper(), kmer_lib)
        member = reverse_complement(seq.upper()) if is_rc else seq.upper()
        tsd_budget = int(tsd_budget_frac * max(ref_len, len(member))) + 1
        if ap_r >= 0 and ap_m >= 0:
            rc_m, id_m, cigar_m = wfa_lib.align_cigar_anchored(
                ref_bytes, member.encode(), ap_r, ap_m,
                max_edit=tsd_budget, min_coverage=0.0,
            )
        else:
            rc_m, id_m, cigar_m = wfa_lib.align_cigar(
                ref_bytes, member.encode(),
                max_edit=tsd_budget, min_coverage=0.0,
            )
        if rc_m < 0:
            member_stats[idx] = (-1.0, True)
        else:
            member_stats[idx] = (id_m, False)
            aligned_seqs.append(_project_cigar(member, cigar_m, ref_len))

    tsd_len, tsd_conf, tsd_motif = detect_tsd_from_alignment(
        aligned_seqs, cons_start, cons_end)

    sequence = trim_and_filter(cons, saturation, cons_start, cons_end,
                               saturation_threshold)
    if not sequence:
        return None

    alignment = None
    if keep_alignments:
        try:
            captured = _align_members_to_consensus(
                sequence, member_indices, sequences, wfa_lib, kmer_lib,
            )
            if captured:
                alignment = _expand_alignment(sequence, captured)
        except Exception as e:
            log.warn(f"keeping alignment failed for a cluster of "
                     f"{n_members} ({e}); consensus kept, alignment dropped")

    return ConsensusResult(
        sequence=sequence,
        num_sequences=n_members,
        tsd_length=tsd_len,
        tsd_confidence=tsd_conf,
        tsd_motif=tsd_motif,
        raw_alignment_length=profile.shape[1],
        alignment=alignment,
    ), member_indices, member_stats


# ---------------------------------------------------------------------------
# Dendrogram-cut helpers for subfamily splitting
# ---------------------------------------------------------------------------

def cut_dendrogram(
    Z: np.ndarray,
    n_leaves: int,
    cut_threshold: float = CUT_THRESHOLD,
) -> Dict[int, List[int]]:
    """Walk UPGMA linkage *Z*; union children at steps with distance <= threshold.
    """
    parent = list(range(n_leaves))
    rank = [0] * n_leaves

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        if rank[ra] == rank[rb]:
            rank[ra] += 1

    rep = list(range(n_leaves)) + [0] * (n_leaves - 1)
    for step in range(len(Z)):
        left, right = int(Z[step, 0]), int(Z[step, 1])
        node = n_leaves + step
        rep[node] = rep[left]
        if Z[step, 2] <= cut_threshold:
            union(rep[left], rep[right])

    groups: Dict[int, List[int]] = defaultdict(list)
    for leaf in range(n_leaves):
        groups[find(leaf)].append(leaf)
    return dict(groups)


def _extract_sub_linkage(
    Z: np.ndarray,
    n_leaves: int,
    subtree_leaves: List[int],
) -> Tuple[np.ndarray, Dict[int, int]]:
    """Extract a sub-linkage matrix for *subtree_leaves* from the full *Z*.
    """
    k = len(subtree_leaves)
    global_to_local = {g: l for l, g in enumerate(subtree_leaves)}
    if k <= 1:
        return np.empty((0, 4)), global_to_local

    node_map: Dict[int, int] = dict(global_to_local)
    next_local = k
    rows: List[List[float]] = []
    for step in range(len(Z)):
        left, right = int(Z[step, 0]), int(Z[step, 1])
        global_node = n_leaves + step
        if left in node_map and right in node_map:
            rows.append([node_map[left], node_map[right], Z[step, 2], 0])
            node_map[global_node] = next_local
            next_local += 1

    sub_Z = np.array(rows, dtype=np.float64) if rows else np.empty((0, 4))
    # Fix count column
    counts: Dict[int, int] = {}
    for i in range(k):
        counts[i] = 1
    for ri, row in enumerate(sub_Z):
        local_node = k + ri
        counts[local_node] = counts[int(row[0])] + counts[int(row[1])]
        sub_Z[ri, 3] = counts[local_node]
    return sub_Z, global_to_local


def _run_subtree_consensus(
    sub_Z: np.ndarray,
    subtree_leaves: List[int],
    sequences: List[str],
    wfa_lib: 'WFALib',
    kmer_lib: 'KmerLib',
    cons_threshold: float,
    merge_max_edit_frac: float,
    use_outgroup: bool,
    min_merge_identity: float,
    label: Optional[str] = None,
) -> Tuple[np.ndarray, int, str, List[int]]:
    """Run bottom-up progressive merging for a single dendrogram subtree.
    """
    k = len(subtree_leaves)

    # ── Initialize leaf profiles in local ID space ───────────────────────
    profiles: Dict[int, Tuple[np.ndarray, int, str]] = {}
    members: Dict[int, List[int]] = {}
    for local_id in range(k):
        global_id = subtree_leaves[local_id]
        seq = sequences[global_id].upper()
        prof = np.zeros((5, len(seq)), dtype=np.int32)
        seq_bytes = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
        rows = _BYTE_TO_ROW[seq_bytes]
        prof[rows, np.arange(len(seq))] = 1
        profiles[local_id] = (prof, 1, seq)
        members[local_id] = [global_id]

    # ── Pre-compute outgroup (sibling) map from sub-linkage ──────────────
    sibling_of: Dict[int, int] = {}
    if use_outgroup:
        for s in range(len(sub_Z)):
            c1, c2 = int(sub_Z[s, 0]), int(sub_Z[s, 1])
            sibling_of[c1] = c2
            sibling_of[c2] = c1

    n_merge_failed = 0

    # ── Bottom-up traversal ──────────────────────────────────────────────
    for step in range(len(sub_Z)):
        idx1, idx2 = int(sub_Z[step, 0]), int(sub_Z[step, 1])
        node_id = k + step

        prof_a, n_a, cons_a = profiles[idx1]
        prof_b, n_b, cons_b = profiles[idx2]

        lmax = max(len(cons_a), len(cons_b))
        merge_budget = int(merge_max_edit_frac * lmax) + 1
        chain_a, chain_b, is_rc = _anchor_chain_and_orient(
            cons_a, cons_b, kmer_lib,
        )

        cons_b_oriented = reverse_complement(cons_b) if is_rc else cons_b
        if is_rc:
            prof_b = _rc_profile(prof_b)

        if len(chain_a) >= 1:
            rc, identity_val, cigar = wfa_lib.align_cigar_anchored(
                cons_a.encode(), cons_b_oriented.encode(),
                chain_a[0], chain_b[0],
                max_edit=merge_budget, min_coverage=0.0,
            )
        else:
            rc, identity_val, cigar = wfa_lib.align_cigar(
                cons_a.encode(), cons_b_oriented.encode(),
                max_edit=merge_budget, min_coverage=0.0,
            )

        do_split = (rc < 0 or not cigar or
                    identity_val < min_merge_identity)

        if do_split:
            # Within a dendrogram subtree: discard smaller side (no orphaning).
            n_merge_failed += 1
            if n_a >= n_b:
                profiles[node_id] = (prof_a, n_a, cons_a)
                members[node_id] = members[idx1]
            else:
                profiles[node_id] = (prof_b, n_b, cons_b_oriented)
                members[node_id] = members[idx2]
        else:
            merged = _merge_profiles(prof_a, n_a, prof_b, n_b, cigar)
            merged_cons = _profile_consensus(merged, threshold=cons_threshold)

            if merged.shape[1] == 0:
                log.warn(f"subtree merge produced empty profile at step={step} "
                         f"A(n={n_a},len={len(cons_a)},cols={prof_a.shape[1]}) "
                         f"B(n={n_b},len={len(cons_b_oriented)},cols={prof_b.shape[1]}) "
                         f"cigar={cigar[:80]}{'...' if len(cigar)>80 else ''} "
                         f"identity={identity_val:.3f}")

            # Outgroup refinement
            if use_outgroup and node_id in sibling_of:
                sib_id = sibling_of[node_id]
                if sib_id in profiles:
                    _sib_prof, _sib_n, sib_cons = profiles[sib_id]
                    merged = _outgroup_refine(
                        merged, sib_cons, merged_cons,
                        wfa_lib, kmer_lib, n_a + n_b,
                    )
                    merged_cons = _profile_consensus(merged,
                                                     threshold=cons_threshold)

            profiles[node_id] = (merged, n_a + n_b, merged_cons)
            members[node_id] = members[idx1] + members[idx2]

        del profiles[idx1]
        del profiles[idx2]
        del members[idx1]
        del members[idx2]

    root_id = 2 * k - 2
    root_prof, root_n, root_cons = profiles[root_id]
    root_members = members[root_id]

    return root_prof, root_n, root_cons, root_members


def build_progressive_consensus(
    sequences:           List[str],
    pairwise_identity:   Dict[Tuple[int, int], float],
    wfa_lib:             'WFALib',
    kmer_lib:            'KmerLib',
    cons_threshold:      float = CONS_THRESHOLD,
    saturation_threshold: float = SATURATION_THRESHOLD,

    merge_max_edit_frac: float = MERGE_MAX_EDIT_FRAC,
    min_merge_identity:  float = MIN_MERGE_IDENTITY,

    split_subfamilies:   bool = False,
    cut_threshold:       float = CUT_THRESHOLD,
    label:               Optional[str] = None,
    debug_plot_path:     Optional[str] = None,
    ltr_divergences:     Optional[List[float]] = None,
    element_names:       Optional[List[str]] = None,
    ambiguity_threshold: float = AMBIGUITY_THRESHOLD,
    use_outgroup:        bool = False,
    keep_alignments:     bool = False,
) -> Optional[List[Tuple[ConsensusResult, List[int], Dict[int, Tuple[float, bool]]]]]:
    """
    Build a consensus using a UPGMA guide tree and progressive profile merging.

    Closer sequences are merged first so each WFA alignment bridges a small
    distance — avoiding the centroid-divergence problem of star consensus.
    kmer_lib is used at each merge step for fast C-speed anchor detection,
    enabling align_cigar_anchored (~2-4x speedup over global WFA for long seqs).

    Parameters
    ----------
    sequences : list of str
        Unaligned cluster member sequences (may be on either strand).
    pairwise_identity : dict
        {(i, j): identity} for graph edges within the cluster (i < j).
        Missing pairs are assigned distance 0.5.
    wfa_lib, kmer_lib : WFALib, KmerLib
        Loaded wfa_align.so and kmer_seed.so wrappers.
    cons_threshold, saturation_threshold : float
        Passed to the downstream consensus pipeline.

    Returns
    -------
    list of (ConsensusResult, member_indices, member_stats) or None
        When a merge exceeds the edit budget the two subtrees become separate
        subfamilies, so a single cluster may yield multiple consensus sequences.
    """
    n_seqs = len(sequences)
    if n_seqs < 2:
        return None  # singletons lack multi-sequence validation

    # ── Build UPGMA guide tree ────────────────────────────────────────────
    # Seed distance matrix with known pairwise identities; unknown pairs
    # start at infinity, then fill them via shortest paths.
    dist_matrix = np.full((n_seqs, n_seqs), np.inf)
    np.fill_diagonal(dist_matrix, 0.0)
    for (i, j), ident in pairwise_identity.items():
        d = 1.0 - ident
        dist_matrix[i, j] = d
        dist_matrix[j, i] = d

    # fill missing pairs with transitive shortest-path
    # distances through known edges.
    for k in range(n_seqs):
        np.minimum(dist_matrix, dist_matrix[:, k:k+1] + dist_matrix[k:k+1, :],
                   out=dist_matrix)

    # Any pairs still unreachable get the 0.5 distance, but note
    # that since we operate over connected components this should never
    # happen (TODO (anant): we should assert here that the graph is connected))
    dist_matrix[np.isinf(dist_matrix)] = 0.5

    condensed = squareform(dist_matrix, checks=False)
    Z = linkage(condensed, method='average')  # UPGMA

    # ── Debug: plot the guide tree if requested ───────────────────────────
    if debug_plot_path is not None and ltr_divergences is not None:
        names = element_names if element_names is not None else [str(i) for i in range(n_seqs)]
        try:
            _plot_cluster_tree(
                Z, dist_matrix, names, ltr_divergences,
                debug_plot_path, label=label or '',
            )
        except Exception as e:
            log.warn(f"debug tree plot failed: {e}")

    # ── Dendrogram-cut: cut tree true to the 80/80/80 rule ────────────────────────
    if split_subfamilies:
        subtree_groups = cut_dendrogram(Z, n_seqs, cut_threshold)

        results: List[Tuple[ConsensusResult, List[int], Dict[int, Tuple[float, bool]]]] = []
        for group_leaves in sorted(subtree_groups.values(), key=lambda g: min(g)):
            group_sorted = sorted(group_leaves)
            k_sub = len(group_sorted)

            if k_sub == 1:
                continue

            sub_Z, _g2l = _extract_sub_linkage(Z, n_seqs, group_sorted)
            root_prof, root_n, root_cons, root_members = _run_subtree_consensus(
                sub_Z, group_sorted, sequences,
                wfa_lib, kmer_lib, cons_threshold,
                merge_max_edit_frac, use_outgroup,
                min_merge_identity, label,
            )
            ret = _finalize_subfamily(
                root_prof, root_n, root_cons, root_members, sequences,
                wfa_lib, kmer_lib, cons_threshold, saturation_threshold,
                ambiguity_threshold, keep_alignments=keep_alignments,
            )
            if ret is not None:
                results.append(ret)

        if not results:
            return None
        return results

    profiles: Dict[int, Tuple[np.ndarray, int, str]] = {}
    members: Dict[int, List[int]] = {}
    for i in range(n_seqs):
        seq = sequences[i].upper()
        prof = np.zeros((5, len(seq)), dtype=np.int32)
        seq_bytes = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
        rows = _BYTE_TO_ROW[seq_bytes]
        prof[rows, np.arange(len(seq))] = 1
        profiles[i] = (prof, 1, seq)
        members[i] = [i]

    n_merge_failed = 0

    sibling_of: Dict[int, int] = {}
    if use_outgroup:
        for s in range(len(Z)):
            c1, c2 = int(Z[s, 0]), int(Z[s, 1])
            node = n_seqs + s
            # c1 and c2 are siblings of each other at this merge
            sibling_of[c1] = c2
            sibling_of[c2] = c1

    # Bottom-up traversal
    for step in range(len(Z)):
        idx1, idx2 = int(Z[step, 0]), int(Z[step, 1])
        node_id = n_seqs + step

        prof_a, n_a, cons_a = profiles[idx1]
        prof_b, n_b, cons_b = profiles[idx2]

        lmax = max(len(cons_a), len(cons_b))
        merge_budget = int(merge_max_edit_frac * lmax) + 1
        chain_a, chain_b, is_rc = _anchor_chain_and_orient(
            cons_a, cons_b, kmer_lib,
        )

        cons_b_oriented = reverse_complement(cons_b) if is_rc else cons_b
        if is_rc:
            prof_b = _rc_profile(prof_b)

        if len(chain_a) >= 1:
            rc, identity_val, cigar = wfa_lib.align_cigar_anchored(
                cons_a.encode(), cons_b_oriented.encode(),
                chain_a[0], chain_b[0],
                max_edit=merge_budget, min_coverage=0.0,
            )
        else:
            rc, identity_val, cigar = wfa_lib.align_cigar(
                cons_a.encode(), cons_b_oriented.encode(),
                max_edit=merge_budget, min_coverage=0.0,
            )

        do_split = (rc < 0 or not cigar)

        if do_split:
            n_merge_failed += 1
            if n_a >= n_b:
                profiles[node_id] = (prof_a, n_a, cons_a)
                members[node_id] = members[idx1]
            else:
                profiles[node_id] = (prof_b, n_b, cons_b_oriented)
                members[node_id] = members[idx2]
        else:
            merged = _merge_profiles(prof_a, n_a, prof_b, n_b, cigar)
            merged_cons = _profile_consensus(merged, threshold=cons_threshold)

            if merged.shape[1] == 0:
                log.warn(f"merge produced empty profile at step={step} "
                         f"A(n={n_a},len={len(cons_a)},cols={prof_a.shape[1]}) "
                         f"B(n={n_b},len={len(cons_b_oriented)},cols={prof_b.shape[1]}) "
                         f"cigar={cigar[:80]}{'...' if len(cigar)>80 else ''} "
                         f"identity={identity_val:.3f}")

            # Outgroup refinement
            if use_outgroup and node_id in sibling_of:
                sib_id = sibling_of[node_id]
                if sib_id in profiles:
                    _sib_prof, _sib_n, sib_cons = profiles[sib_id]
                    merged = _outgroup_refine(
                        merged, sib_cons, merged_cons,
                        wfa_lib, kmer_lib, n_a + n_b,
                    )
                    merged_cons = _profile_consensus(merged,
                                                     threshold=cons_threshold)

            profiles[node_id] = (merged, n_a + n_b, merged_cons)
            members[node_id] = members[idx1] + members[idx2]
            if identity_val < 0.65:
                anchor_tag = (f"anchor=({chain_a[0]},{chain_b[0]})"
                              if chain_a else "no_anchor")
                log.detail(f"    merge_ok   step={step}/{len(Z)-1} "
                           f"A(n={n_a},len={len(cons_a)}) vs "
                           f"B(n={n_b},len={len(cons_b_oriented)}) "
                           f"identity={identity_val:.3f} "
                           f"is_rc={is_rc} {anchor_tag}")

        del profiles[idx1]
        del profiles[idx2]
        del members[idx1]
        del members[idx2]

    root_id = 2 * n_seqs - 2
    root_profile, root_n, _root_cons = profiles[root_id]
    root_members = members[root_id]

    if n_merge_failed > 0:
        tag = f"cluster {label}" if label is not None else "cluster"
        log.warn(f"{tag}: {n_merge_failed} progressive merge(s) "
                 f"failed WFA alignment → discarded smaller subtrees")

    results: List[Tuple[ConsensusResult, List[int], Dict[int, Tuple[float, bool]]]] = []
    root_ret = _finalize_subfamily(
        root_profile, root_n, _root_cons, root_members, sequences,
        wfa_lib, kmer_lib, cons_threshold, saturation_threshold,
        ambiguity_threshold, keep_alignments=keep_alignments,
    )
    if root_ret is not None:
        results.append(root_ret)

    if not results:
        return None
    return results


# ---------------------------------------------------------------------------
# Debug: cluster tree + LTR age visualization
# ---------------------------------------------------------------------------

def _plot_cluster_tree(
    Z:               np.ndarray,
    dist_matrix:     np.ndarray,
    names:           List[str],
    ltr_divergences: List[float],
    output_path:     str,
    label:           str = '',
) -> None:
    """
    Plot a UPGMA dendrogram for a cluster.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram

    n = len(names)
    fig_h = max(4, n * 0.25)
    fig, ax_tree = plt.subplots(1, 1, figsize=(10, fig_h))

    # Truncate long element names for readability
    short_names = [nm[:30] for nm in names]

    # Dendrogram (horizontal, leaves on the left axis)
    dendrogram(
        Z, orientation='left', labels=short_names, ax=ax_tree,
        leaf_font_size=max(5, min(8, 200 // max(n, 1))),
        color_threshold=0.20,  # colour branches at ~80% identity
    )
    ax_tree.set_xlabel('WFA distance (1 - identity)')
    ax_tree.set_title(f'Cluster {label} — UPGMA guide tree ({n} members)')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


