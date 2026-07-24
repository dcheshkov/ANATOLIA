#!/usr/bin/env bash
# drive_fN.sh N [JOBS] [GAP_MEM]
#
# Compute f(N) -- the number of 2^-closed permutation groups of degree N with no
# fixed point (the "essentially N-spin" realizable symmetry groups) -- by the
# orbit-partition decomposition, without ever building the S_N subgroup lattice.
# Runs one GAP job per orbit partition of N with all parts >= 2 (emit_partition_
# colourings.g), up to JOBS in parallel, then pipes every colouring through
# filter_closed once (its canonical-form dedup merges the block-permutation
# conjugates across partitions).  Then  a(N) = a(N-1) + f(N).
#
#   JOBS     max concurrent GAP jobs      (default: min(nproc-1, 8))
#   GAP_MEM  -o workspace cap per job     (default: 8g; raise for N>=14, the
#            [N-2,2] partition -> S_{N-2} x S_2 is the heaviest)
#
# Needs: GAP with the transitive-groups library (gap-transgrp), and filter_closed
# built here in gap/ (g++ -O3 -march=native -fopenmp -std=c++17 -I../engine \
#  -o filter_closed filter_closed.cpp; canon_ir.h lives in ../engine).
# Validation: f(4..13) = 5,3,16,9,54,41,151,112,554,368.
set -euo pipefail

N="${1:?usage: drive_fN.sh N [JOBS] [GAP_MEM]}"
JOBS="${2:-$(( $(nproc) > 1 ? $(nproc) - 1 : 1 ))}"; [ "$JOBS" -gt 8 ] && [ -z "${2:-}" ] && JOBS=8
GAP_MEM="${3:-8g}"
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$(mktemp -d)"
trap 'rm -rf "$OUT"' EXIT

# partitions of N with min part >= 2, one per line as "12,2"
mapfile -t PARTS < <(gap -q -b -c \
  "for p in Filtered(Partitions($N),x->Minimum(x)>=2) do \
     Print(JoinStringsWithSeparator(List(p,String),\",\"),\"\n\"); od;" 2>/dev/null)

echo "f($N): ${#PARTS[@]} orbit partitions (parts>=2), up to $JOBS parallel jobs, -o $GAP_MEM/job" >&2
i=0
for p in "${PARTS[@]}"; do
  glist="[${p}]"
  gap -o "$GAP_MEM" -q -b -c "NN:=$N;; PART:=$glist;;" \
      "$HERE/emit_partition_colourings.g" > "$OUT/p_${i}.txt" 2>"$OUT/e_${i}.txt" &
  i=$((i+1))
  while [ "$(jobs -r | wc -l)" -ge "$JOBS" ]; do wait -n; done
done
wait

# surface any GAP job that failed (e.g. hit its -o cap)
for e in "$OUT"/e_*.txt; do
  if grep -q 'memory\|Error' "$e" 2>/dev/null; then
    echo "WARNING: a GAP job reported an error ($e):" >&2; tail -2 "$e" >&2
  fi
done

grep -h '^#' "$OUT"/p_*.txt >&2 || true
echo "=== f($N) ===" >&2
cat "$OUT"/p_*.txt | grep -v '^#' | "$HERE/filter_closed" -
