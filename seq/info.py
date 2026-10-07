"""Readers for the sequence and alignment files the pipeline passes around.
"""

import gzip
import subprocess
from pathlib import Path

_REPO = Path(__file__).parent.parent
TACO_BIN     = str(_REPO / 'external' / 'alntools' / 'taco')
ALN_INFO_BIN = str(_REPO / 'seq' / 'aln_info')


def _taco_info(seq_path):
    """Run `taco info` on the .1taco companion of a .1seq and return its stdout."""
    taco_path = str(seq_path)[:-5] + '.1taco'
    return subprocess.run([TACO_BIN, 'info', taco_path],
                          capture_output=True, text=True, check=True).stdout


def aln_info(aln_path):
    """Summarise a .1aln file in one pass of the `aln_info` binary."""
    info = {'alignments': 0, 'total_bp': 0, 'active_bp': 0,
            'active_scaffolds': set()}
    p = Path(aln_path)
    if not p.exists() or p.stat().st_size == 0:
        return info
    out = subprocess.run([ALN_INFO_BIN, str(aln_path)],
                         capture_output=True, text=True, check=True).stdout
    in_list = False
    for line in out.splitlines():
        if in_list:
            if line:
                info['active_scaffolds'].add(line)
        elif line == 'active_scaffolds:':
            in_list = True
        elif line.startswith('alignments='):
            info['alignments'] = int(line.split('=', 1)[1])
        elif line.startswith('total_scaffold_bp='):
            info['total_bp'] = int(line.split('=', 1)[1])
        elif line.startswith('active_scaffold_bp='):
            info['active_bp'] = int(line.split('=', 1)[1])
    return info


def fasta_size(fasta_path):
    """Count total non-header bases in a FASTA file (handles .gz)."""
    opener = gzip.open if str(fasta_path).endswith('.gz') else open
    total = 0
    with opener(fasta_path, 'rt') as f:
        for line in f:
            if not line.startswith('>'):
                total += len(line.strip())
    return total


def seq_size(seq_path):
    """Total scaffold bp for a FASTA or scaffold-aware .1seq."""
    p = str(seq_path)
    if not p.endswith('.1seq'):
        return fasta_size(seq_path)
    for line in _taco_info(p).splitlines():
        if line.startswith('Total:'):
            for tok in line.split():
                if tok.startswith('tacoLen='):
                    return int(tok.split('=', 1)[1])
    raise RuntimeError(f"could not parse tacoLen from `taco info` for {p}")


def scaffold_order(seq_path):
    """Return the ordered list of scaffold IDs in a FASTA or .1seq."""
    p = str(seq_path)
    names = []
    if p.endswith('.1seq'):
        for line in _taco_info(p).splitlines():
            s = line.lstrip()           # "seq <idx> (<name>): origLen=..."
            if not s.startswith('seq '):
                continue
            l = s.find('(')
            r = s.find(')', l + 1)
            if l > 0 and r > l:
                names.append(s[l + 1:r])
        return names
    opener = gzip.open if p.endswith('.gz') else open
    with opener(seq_path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                names.append(line[1:].split()[0])
    return names
