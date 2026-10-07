#!/usr/bin/env bash
# Build seq/aln_info against the SEQUENCE_UTILITIES copies of
# ONElib / utils / array / dict.
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SU="${ONECODE_SEQ_UTILS:-$DIR/../external/ONEcode/SEQUENCE_UTILITIES}"

if [[ ! -d "$SU" ]]; then
  echo "SEQUENCE_UTILITIES dir not found at $SU" >&2
  echo "Set ONECODE_SEQ_UTILS to override." >&2
  exit 1
fi

for obj in ONElib.o utils.o array.o dict.o hash.o; do
  if [[ ! -f "$SU/$obj" ]]; then
    echo "Missing $SU/$obj — run 'make' inside $SU first." >&2
    exit 1
  fi
done

gcc -O3 -Wall -I"$SU" -DONEIO -o "$DIR/aln_info" \
    "$DIR/aln_info.c" \
    "$SU/ONElib.o" \
    "$SU/utils.o" "$SU/array.o" "$SU/dict.o" "$SU/hash.o" \
    -lm -lz -lpthread

echo "Built $DIR/aln_info"
