# analysis/ — analysis of the finished catalogues

Tools that work on the already-built catalogues (`../../materials/*Spins.txt`):
the machine taxonomy and the symmetry-exact 3D-embedding tool. Both are
self-contained (no dependence on `../engine`).

## `classify_catalog.cpp` — the machine taxonomy

Reads a catalogue (`8Spins.txt` format) and emits one TSV line per entry:
`id, |Aut|, twins (n_g), |Γq|, orbits, fixed point, cyclicity,
max element order, abelianity, ambivalence`.

```
g++ -O3 -march=native -std=c++17 -o classify_catalog classify_catalog.cpp
./classify_catalog N catalogue.txt >> taxonomy.tsv
```
This is how `../../materials/catalog_taxonomy.tsv` was built (6112 entries,
`N=3..14`). Note: the automorphisms are enumerated by a vertex-ordering scan,
so the rare groups of very large order at `N=14` (e.g. the tower extension of
`S_13`, `|Aut|=13!`) are slow but do finish.

## Python tools

- **`topology2xyz.py`** — symmetry-exact 3D embedding of a topology
  (eigenspaces of the class-weighted adjacency matrix); `.xyz/.mol2/` legend.
  Requires python3 + numpy. Full documentation:
  `../../docs/spin_system_visualization_en.md` (русская версия — `..._ru.md`).
  ```
  python3 topology2xyz.py --file ../../materials/8Spins.txt --order 6 --index 5 --out c6
  ```
