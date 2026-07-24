# SpinSymmetry - realizable symmetry groups of NMR spin systems

## Catalogues `NSpins.txt` (`N = 3..14`)

One file per level (`3Spins.txt` … `14Spins.txt`), all `a(N)` realizable
symmetry classes of `N`-spin systems. The format:

```
N=7: closed subgroup classes = 36 (…)          <- header: level and a(N)
  orders: 1(x1) 2(x3) 3(x1) …                  <- header: |Aut| histogram
--- order K (xM) ---                           <- section: M entries with |Aut| = K
c1 c2 c3 … cE                                  <- one entry = one edge colouring
…
```

Each entry line is the canonical edge colouring of the complete graph `K_N`:
`E = N(N-1)/2` integer class labels in upper-triangle **column order**

```
J12  J13 J23  J14 J24 J34  J15 J25 J35 J45  …
```

(for `N = 4` the six positions are exactly `J12 J13 J23 J14 J24 J34`). Two
couplings belong to the same class iff their labels are equal; the label
values themselves are arbitrary - only the equality pattern matters. Entry
`j` (1-based) of section `K` is referred to throughout the project and the
article as `N/oK.j`.

Two file-specific notes. `8Spins.txt` additionally carries the classical
group names in its section titles (e.g. `--- order 2 (C2,x4) ---`) and a
leading `*` on every entry line — a relic of the cross-check against the
classical catalogue; parsers skip the marker. `14Spins.txt` was assembled by
the orbit-partition decomposition, and its canonical representative is the
IR certificate (`bestNorm`) rather than the lex-min form used for
`N ≤ 13` (see `code/engine/canon_ir.h`); the representative is
isomorphism-invariant, and no consumer depends on the choice.

Control totals: `a(3..14) = 3, 8, 11, 27, 36, 90, 131, 282, 394, 948, 1316,
2866`; `Σ = 6112`.

## `catalog_taxonomy.tsv`

One row per catalogue entry (6112 rows). Columns: `id` (`N/oK.j`), `|Aut|`,
twin sizes `n_g` (comma-separated), `|Γq|` (residual factor-group order,
`|Aut| = |Γq|·∏ n_g!`), vertex orbits (`+`-separated sizes), fixed point
(`y/n`), cyclic (`y/n`), maximal element order, abelian (`y/n`), ambivalent
(`y/n` — whether every element is conjugate to its inverse, i.e.\ whether
all characters are real).

## `8Spins_crossref.txt`

Machine-verified cross-reference of the classical 8-spin catalogue: for
every row of the classical table, the closure of its claimed generators is
computed, `|Aut|` compared with the claimed order, the canonical form
matched into `8Spins.txt`, and the twin groups, `|Γq|`, vertex orbits and
the structural name extracted. (The generator indices in the file header
refer to the permutation list of the classical catalogue; neither that list
nor the scanned PDF is redistributed in this deposit.)