#!/usr/bin/env bash
# Build algorithms/extract_ltrs against the SEQUENCE_UTILITIES copies of
# ONElib / seqio / utils / array / dict.
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SU="${ONECODE_SEQ_UTILS:-$DIR/../external/ONEcode/SEQUENCE_UTILITIES}"

if [[ ! -d "$SU" ]]; then
  echo "SEQUENCE_UTILITIES dir not found at $SU" >&2
  echo "Set ONECODE_SEQ_UTILS to override." >&2
  exit 1
fi

# Make sure the helper object files exist (the SEQUENCE_UTILITIES Makefile
# builds them when you run `make` in that directory).
for obj in ONElib.o seqio.o utils.o array.o dict.o hash.o; do
  if [[ ! -f "$SU/$obj" ]]; then
    echo "Missing $SU/$obj — run 'make' inside $SU first." >&2
    exit 1
  fi
done

gcc -O3 -Wall -I"$SU" -DONEIO -o "$DIR/extract_ltrs" \
    "$DIR/extract_ltrs.c" \
    "$SU/seqio.o" "$SU/ONElib.o" \
    "$SU/utils.o" "$SU/array.o" "$SU/dict.o" "$SU/hash.o" \
    -lm -lz -lpthread

echo "Built $DIR/extract_ltrs"
