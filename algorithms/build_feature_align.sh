#!/usr/bin/env bash
# Build kmer_seed.so and wfa_align.so in-place.
# Called automatically by c_libs.py when a .so is missing.
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CC="${CC:-gcc}"

# macOS: Apple clang needs libomp from Homebrew for OpenMP support
OPENMP_FLAGS="-fopenmp"
if [[ "$(uname)" == "Darwin" ]]; then
    LIBOMP_PREFIX="$(brew --prefix libomp 2>/dev/null || true)"
    if [[ -n "$LIBOMP_PREFIX" ]]; then
        OPENMP_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include -L${LIBOMP_PREFIX}/lib -lomp"
    else
        echo "WARNING: libomp not found. Install with: brew install libomp" >&2
        echo "         Building wfa_align.so WITHOUT OpenMP (single-threaded)." >&2
        OPENMP_FLAGS=""
    fi
fi

$CC -O3 -Wall -shared -fPIC -o "$DIR/kmer_seed.so" "$DIR/kmer_seed.c"
$CC -O3 -Wall $OPENMP_FLAGS -shared -fPIC -o "$DIR/wfa_align.so" "$DIR/wfa_align.c"
