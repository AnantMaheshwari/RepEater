"""Write the alignment behind each consensus, under --keep-alignments.
"""

from pathlib import Path

_WRAP = 60

def _row_names(aln):
    """Row labels, falling back to member indices if names were never set."""
    if aln.row_names:
        return aln.row_names
    return [str(i) for i in aln.row_indices]


def _wrapped(seq, width=_WRAP):
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


def write_alignment_fasta(aln, cons_name, path):
    """Aligned FASTA: the consensus, then one record per member."""
    names = _row_names(aln)
    with open(path, 'w') as fh:
        fh.write(f">{cons_name}\n{_wrapped(aln.consensus_row)}\n")
        for name, row in zip(names, aln.rows):
            fh.write(f">{name}\n{_wrapped(row)}\n")


def write_cluster_alignment(aln, cons_name, out_dir):
    """Write one consensus's member alignment into *out_dir*."""
    write_alignment_fasta(aln, cons_name, Path(out_dir) / f'{cons_name}.fa')
