# engine/ — the core and the BFS enumerators

The core machinery for the symmetries of an edge-weighted `K_N`, plus two
enumerators of all 2*-closed subgroups of `S_N` (the realizable symmetry
groups of spin systems) by the sequential-symmetrization method. The headers
here are also included by the tools in `../gap` (via `-I../engine`).

Conventions: the edges of `K_N` are numbered in upper-triangle column order
(`J12 J13 J23 J14 …`, as in `../../materials/8Spins.txt`); a colouring = an
RGS of length `E = N(N-1)/2`; C++17, `N ≤ 16`.

## Headers (header-only)

- **`spin_graph.h`** — `K_N` edge indexing, RGS encoding, the coarsening step
  `coarsen()` (merging edge classes along the orbits of a permutation), a
  partition generator, a timer. Shared by both drivers.
- **`canon_ir.h`** — nauty-style canonicalization: the canonical form (lex-min
  RGS **and** the IR certificate `bestNorm`), `Aut` generators, vertex orbits,
  `|Aut|` via Schreier–Sims (without enumerating group elements). The core of
  every tool that needs graph symmetry (including those in `../gap`).

## `closed_groups.cpp` — single-node enumerator (OpenMP)

BFS symmetrization starting from the discrete colouring; prints the number of
classes, the order histogram and (with `VERBOSE=1`) the canonical colourings
in `8Spins.txt` format.

```
g++ -O3 -march=native -fopenmp -std=c++17 -o closed_groups closed_groups.cpp
./closed_groups [Nmin [Nmax]]                 # default 3..9
OMP_NUM_THREADS=96 ./closed_groups 12 12
VERBOSE=1 ./closed_groups 10 10 > cat_N10.txt
```
Check: `N=3..12 → 3,8,11,27,36,90,131,282,394,948`.

## `closed_groups_mpi.cpp` — distributed (MPI + OpenMP)

The same algorithm plus batched exchange of raw colourings, checkpoints,
`RESUME`. For `N=13` on several nodes.

```
mpicxx -O3 -march=native -fopenmp -std=c++17 -o closed_groups_mpi closed_groups_mpi.cpp
mpirun -np 4 -H srv2:2,srv3:2 --map-by ppr:2:node --bind-to socket \
  -x OMP_NUM_THREADS=24 -x OMP_PLACES=cores -x OMP_PROC_BIND=close \
  -x VERBOSE=1 ./closed_groups_mpi 13 13 > 13Spins.txt
```
Flags: `RESUME=1` continue from a checkpoint; `NOCKPT=1` write no checkpoints;
`BATCH_MB=<MB>` exchange batch size (default 256).
Check: `mpirun -np 2 ./closed_groups_mpi 3 10 → …,282`.

## `canon_ir_test.cpp` — verifier for `canon_ir.h`

Cross-checks the canonical form, `|Aut|` and the orbits against the legacy
exhaustive search over a catalogue; the `--sstest` mode checks Schreier–Sims
against brute force on 48000 random groups.

```
g++ -O3 -march=native -std=c++17 -o canon_ir_test canon_ir_test.cpp
./canon_ir_test N catalogue.txt        # catalogue cross-check
./canon_ir_test --sstest               # Schreier–Sims regression
```

Method and completeness proof — `../../docs/symmetry_groups_of_spin_systems_{ru,en}.md` §8.
