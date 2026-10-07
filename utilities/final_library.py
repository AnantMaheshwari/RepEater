import re
from pathlib import Path

_LTR_RE = re.compile(r'^LTR(/|$)', re.IGNORECASE)


def is_ltr(class_label):
    return bool(_LTR_RE.match((class_label or '').strip().rstrip('?')))


def parse_header(header):
    header = header.lstrip('>').rstrip('\n')
    name = header.split()[0] if header.split() else ''

    m = re.search(r'\bclass=(\S+)', header)
    if m:
        return name, m.group(1)
    if '#' in name:
        name, cls = name.split('#', 1)
        return name, cls
    return name, 'Unknown'


def write_final_library(classified_fasta, output_dir, base):
    """Write `final/` """
    final_dir = Path(output_dir) / 'final'
    final_dir.mkdir(parents=True, exist_ok=True)

    ltr_fasta = final_dir / f'{base}_all_rounds_classified.fasta'
    rm_fasta = final_dir / f'{base}_all_rounds_classified_rm.fasta'

    src = Path(classified_fasta)
    if not src.exists():
        ltr_fasta.write_text('')
        rm_fasta.write_text('')
        return str(ltr_fasta), str(rm_fasta), 0, 0

    kept = dropped = 0
    seen = {}
    with open(src) as fh, open(ltr_fasta, 'w') as out_ltr, \
            open(rm_fasta, 'w') as out_rm:
        keep = False
        for line in fh:
            if not line.startswith('>'):
                if keep:
                    out_ltr.write(line)
                    out_rm.write(line)
                continue

            name, cls = parse_header(line)
            keep = is_ltr(cls)
            if not keep:
                dropped += 1
                continue

            seen[name] = seen.get(name, 0) + 1
            rm_name = name if seen[name] == 1 else f'{name}_{seen[name]}'

            out_ltr.write(line if line.startswith('>') else f'>{line}')
            out_rm.write(f'>{rm_name}#{cls}\n')
            kept += 1

    return str(ltr_fasta), str(rm_fasta), kept, dropped
