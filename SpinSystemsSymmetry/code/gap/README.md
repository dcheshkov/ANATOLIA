# gap/ — orbit decomposition, closedness filter, catalogue assembly

The route by which `N=14` was obtained, plus an independent re-check of the
lower levels. GAP enumerates the groups, the C++ filter keeps the 2*-closed
ones. All C++ tools include `canon_ir.h` from `../engine`, so they are built
with `-I../engine`.

Requires GAP with the transitive-groups library (`gap-core` + `gap-transgrp`).

## `filter_closed.cpp` — the 2*-closedness filter

Reads lines `"N |G| c_1 … c_E"` (the edge orbit colouring of a group), keeps
the closed ones (`|Aut(colouring)| = |G|`) and deduplicates by canonical form
→ the class count `a(N)`. `DUMP=<file>` additionally writes the canonical
colourings of the closed groups (catalogue assembly). The test is
parallelized over lines (OpenMP).

```
g++ -O3 -march=native -fopenmp -std=c++17 -I../engine -o filter_closed filter_closed.cpp
… | ./filter_closed -          # reads stdin when "-" or no argument
```

## Enumerate-then-filter (re-check of `a(3..13)`)

`emit_subgroup_colourings.g` — enumerates ALL subgroup classes of `S_N`
(GAP, Hulpke's algorithm = OEIS A000638), one representative per line.

```
gap -q -b -c 'NMIN:=3;; NMAX:=13;;' emit_subgroup_colourings.g > col.txt
grep -v '^#' col.txt | OMP_NUM_THREADS=96 ./filter_closed -
```
Check: `N=3..13 → 3,8,11,27,36,90,131,282,394,948,1316`.
The bottleneck is the GAP enumeration (`N=11` ~50 s, `N=12` ~4 min, `N=13`
~8 min); at `N=14` the full subgroup lattice does not fit into GAP's memory
(the `Zuppos` phase, >13 GB) — hence the decomposition below.

## Orbit-partition decomposition (how `a(14)` was obtained)

`a(N) = a(N-1) + f(N)`; `f(N)` (the closed classes with no fixed point) is
assembled over the partitions of `N` with parts `≥ 2`, bypassing the `S_N`
lattice. The driver `drive_fN.sh` sweeps all partitions (per partition —
`emit_partition_colourings.g`; for the balanced two-block partitions `[a,b]`
— `emit_subdirect2.g`, subdirect products, much faster than the full
`S_a × S_b` lattice), then `filter_closed`.

```
./drive_fN.sh 14 [JOBS] [GAP_MEM]              # prints f(N)
# a single partition:
gap -q -b -c 'NN:=14;; PART:=[8,6];;' emit_subdirect2.g > p.txt
```
Check: `f(4..14) = 5,3,16,9,54,41,151,112,554,368,1550`; `a(14)=1316+1550=2866`.
Cross-validation: for every two-block partition the full Young lattice and the
subdirect products give byte-identical sets of closed colourings.

`emit_transitive_colourings.g` — the transitive part only (for when the `S_N`
lattice is out of reach): `degree 14 → 63 groups, 10 closed (single-orbit)`.

## Assembling the full catalogue

`extend_tower.cpp` — the tower extension of a degree-`M` catalogue to `M+1`
(the fixed-point part): append an isolated vertex, colouring its edges with
fresh labels per vertex orbit; prints raw (non-canonicalized) colourings.

```
g++ -O3 -march=native -std=c++17 -I../engine -o extend_tower extend_tower.cpp
./extend_tower 13 ../../materials/13Spins.txt > fixedpoint14.txt   # 1316 entries
```

`catalog_format.cpp` — assembles a catalogue in `materials/*Spins.txt` format
from a stream of colourings: IR canonicalization (`bestNorm`), dedup,
grouping by `|Aut|`, header + `--- order K ---` sections.

```
g++ -O3 -march=native -fopenmp -std=c++17 -I../engine -o catalog_format catalog_format.cpp
cat fixedpoint14.txt f14_ff.txt | ./catalog_format 14 "note" > 14Spins.txt   # 2866 = 1316 + 1550
```
`N=14` uses the IR-canonical representative (`bestNorm`), not the lex-min form
of `N ≤ 13`: lex-min canonicalization degenerates on the nearly asymmetric
colourings of `N=14` (see `../engine/canon_ir.h`). The representative is
isomorphism-invariant, and no consumer of the catalogue depends on the choice.

Method and completeness proof of the decomposition —
`../../docs/symmetry_groups_of_spin_systems_{ru,en}.md` §8.5.
