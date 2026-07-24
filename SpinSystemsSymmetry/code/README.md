# SpinSymmetry / code

The toolkit is organized into three folders, each with its own `README.md`
giving a full description and build/run lines.

Common conventions:
- the spin-system graph = the edge-weighted complete graph `K_N`; edges are
  numbered in upper-triangle column order (`J12 J13 J23 J14 …`), as in
  `../materials/8Spins.txt`;
- an edge colouring = a restricted growth string (RGS) of length
  `E = N(N-1)/2`;
- the C++ tools target C++17, `N` up to 16.

| folder | contents |
|---|---|
| [`engine/`](engine/README.md) | the core: graph indexing (`spin_graph.h`) + nauty-style canonicalizer (`canon_ir.h`) + BFS enumerators of the 2*-closed groups (single-node and MPI) + verifier |
| [`gap/`](gap/README.md) | the route to `N=14`: GAP orbit-partition decomposition, 2*-closedness filter, catalogue assembly. The C++ files include `canon_ir.h` from `engine/` via `-I../engine` |
| [`analysis/`](analysis/README.md) | analysis of the finished catalogues: the machine taxonomy (`classify_catalog`) and the symmetry-exact 3D-embedding tool (`topology2xyz`) |

## Dependency graph

```
engine/spin_graph.h ─┐
                     ├─> engine/closed_groups.cpp, engine/closed_groups_mpi.cpp
engine/canon_ir.h ───┼─> engine/canon_ir_test.cpp
                     └─> gap/filter_closed.cpp, gap/extend_tower.cpp,
                         gap/catalog_format.cpp        (built with -I../engine)

analysis/classify_catalog.cpp — self-contained (does not use engine/).
```

## How to obtain `a(N)`

- `N ≤ 13`: BFS symmetrization — `engine/closed_groups.cpp` (up to 12) and
  `engine/closed_groups_mpi.cpp` (13); or, independently, enumerate-then-filter —
  `gap/emit_subgroup_colourings.g | gap/filter_closed`.
- `N = 14`: the orbit-partition decomposition — `gap/drive_fN.sh` (yields
  `f(14)`, then `a(14) = a(13) + f(14) = 2866`); the full catalogue is
  assembled by `gap/extend_tower.cpp` + `gap/catalog_format.cpp`.

The sequence: `a(N) = 1, 1, 3, 8, 11, 27, 36, 90, 131, 282, 394, 948, 1316, 2866`.

Method, completeness proof and tables — `../docs/symmetry_groups_of_spin_systems_{ru,en}.md`
(§8 — the BFS and its completeness, §8.5 — enumerate-then-filter and the
orbit-partition decomposition). Catalogues and data — `../materials`.