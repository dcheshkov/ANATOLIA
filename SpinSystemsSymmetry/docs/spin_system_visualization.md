# Visualizing spin systems: symmetry-exact 3D embeddings

*Documentation for the tool `code/analysis/topology2xyz.py`. Theory: `symmetry_groups_of_spin_systems_en.md` §7
(Euclidean realization) and Section 6 of the article (spectral embedding,
Proposition 3).*

## 1. The task

A catalogue entry is an abstract edge colouring of the complete graph `K_N`.
To see its symmetry with one's own eyes we need an embedding into
three-dimensional space in which **every graph automorphism acts as a rigid
motion** (a rotation or a rotation-reflection). The tool `topology2xyz`
constructs such an embedding and exports it as a pseudo-molecule
(`.xyz`/`.mol2`): spins are atoms, the element encodes the vertex orbit,
bonds are selected J classes. The files open in any chemistry viewer
(Avogadro, VMD, Jmol, ChemCraft).

## 2. Theoretical backbone: the distorted simplex

The intuition "a spin system is a distorted regular simplex" is literally
exact: the realizable (2*-closed) groups are precisely the isometry groups of
full-dimensional `N`-point configurations in `R^(N−1)`, viewed as permutation
groups of the points (the Euclidean realization theorem; the squared-distance
matrix of the regular simplex is an interior point of the cone of Euclidean
distance matrices, so small per-edge-class perturbations remain legal
configurations). Chirality is included: a twisted antiprism is a lawful point
configuration with the purely rotational group `C₃`. The three-dimensional
picture is a symmetry-consistent "shadow" of that `R^(N−1)` configuration.

## 3. Method: diagonalizing the topology matrix

Each edge class `c` receives its own generic weight (in the code
`w_c = 1 + 0.7318·√(c+2)` — pairwise distinct irrationals) and the
**topology matrix** `W` is assembled: `W_ij = w_c(ij)`, zero diagonal
(vertex classes are redundant for `N ≥ 3` — any vertex partition is induced
by the edges).

Two key facts.

- **`Aut(W) = Aut(G)` exactly.** An automorphism preserves the classes, so
  its permutation matrix commutes with `W`; conversely, a permutation
  commuting with `W` must preserve every class (the weights are pairwise
  distinct). The spectrum of `W` sees the entry's true symmetry — and nothing
  more.
- **The eigenspaces of `W` are `Aut`-invariant.** This is the same
  commutation `[W, P_σ] = 0` that block-diagonalizes the Hamiltonian (the
  L0–L2 methodology), just in vertex space: one group, two matrices, two
  scales (`N` versus `2^N`).

The embedding: take an orthonormal eigenbasis `U`, select a union of
**whole** eigenspaces of total dimension 3, and use the rows of the
corresponding block `U_S` as vertex coordinates. Then for every `σ ∈ Aut`
one has `P_σ·X = X·Q_σ` with orthogonal `Q_σ` — all automorphisms act as
exact isometries, and the Gramian `X·Xᵀ` (the invariant projector onto `S`)
fixes the picture up to a global rotation. Among all whole-eigenspace
selections of total dimension 3 the tool picks the one **maximizing the
minimal inter-vertex distance** (readability + guaranteed vertex
separation), then centres and rescales (`--scale`, default 1.5
"angstrom-like" units).

## 4. Verifying exactness

For `N ≤ 9` the tool builds `Aut` by exhaustive search, computes the vertex
orbits, and for each `σ` recovers the nearest orthogonal matrix `Q_σ`
(orthogonal Procrustes via SVD). The maximal residual `‖P_σX − XQ_σ‖` is
printed in the report; on the demo entries it is of order `1e-14` — machine
precision. For `N ≥ 10` (`--no-verify`) the check is skipped and the orbits
are approximated by coupling-profile classes (the multiset of incident
classes).

## 5. Two by-design "failures" — both informative

1. **Isometric but not faithful.** Distinct automorphisms may land on the
   same isometry: the group may simply not fit into `O(3)`. Example:
   `(C₂)⁴ = D₂×D₂` of order 16 on two tetrads (entry `8/o16.8` — also the
   unique round-4 entry of the symmetrization at `N = 8`) exceeds the largest
   elementary abelian 2-subgroup of `O(3)`, the diagonal sign matrices
   `(C₂)³` of order 8. The report counts the distinct isometries realized.
2. **No separating 3-dimensional invariant subspace.** The group is not a 3D
   point group on the given vertices at all: every 3-dimensional invariant
   choice collapses vertices. The tool then finds the smallest separating
   symmetric embedding of dimension `D > 3` and "folds" the surplus
   coordinates into the third axis with small deterministic coefficients:
   the vertices separate, and the picture stays exact for the subgroup fixing
   the folded directions. The classic example beyond `N = 8` is `C₉` with
   orbits `9+3` (entry `12/o9.3`): the generator advances the 9-ring by one
   step while turning the triangle by 120°; no rigid rotation of `R³` can do
   that.

The flag "does a separating 3-dimensional invariant subspace exist" is also
used in the classification documentation as the 3D-realizability marker
(`classification_and_nomenclature_en.md`).

## 6. Usage

```
# from a catalogue line (edge labels of K_N in column order)
python3 topology2xyz.py --labels "0 1 1 1 1 2 ..." --out name

# from a catalogue file: section "--- order 6 ---", 5th line
python3 topology2xyz.py --file ../../materials/8Spins.txt --order 6 --index 5 --out c6
```

Options: `--bonds auto|all|k1,k2,...` — which classes become bonds
(`auto` = classes with ≤ `N` edges); `--scale d` — minimal inter-vertex
distance in output units; `--no-verify` — skip the `Aut` search. Output
files: `name.xyz` (element = vertex orbit: C, N, O, S, P, ...), `name.mol2`
(same atoms + bonds of the selected classes), `name.txt` (legend: classes,
sizes, what is drawn).

## 7. Gallery: six showcase entries

Six `N = 8` entries with the actual tool reports; each picture is reproduced by
`python3 topology2xyz.py --file ../../materials/8Spins.txt --order K --index j --out name`
with `K.j` taken from the entry identifier (the "entry" column):

| name | entry | \|Aut\| | eigenspaces | report |
|---|---|---|---|---|
| `8-o48-1_cube` | 8/o48.1 | 48 | [3] | exact, faithful (the cube; residual 5·10⁻¹⁵) |
| `8-o8-3_prism` | 8/o8.3 | 8 | [1, 2] | exact, faithful (prism) |
| `8-o8-4_antiprism` | 8/o8.4 | 8 | [1, 2] | exact, faithful (antiprism) |
| `8-o6-5_C6` | 8/o6.5 | 6 | [2, 1] | exact, faithful (chiral twisted triangles + magnetic pair; Aut ≅ C₆) |
| `8-o4-10_C4tw` | 8/o4.10 | 4 | [1, 2] | exact, faithful (twisted square stack, `C₄`) |
| `8-o16-8_D2xD2` | 8/o16.8 | 16 | [1, 1, 1, 1] → fold | no separating 3-dim subspace; 8 of 16 symmetries exact |

The cube is the model case: a single whole 3-dimensional eigenspace yields
the classical solid with the full group of order 48 realized by rigid
motions. `D₂×D₂` is the model counterexample: both caveats of §5 fire on one
entry.

## 8. Cross-references

- Article, Section 6 "Symmetry-Exact Visualization: Diagonalizing the
  Topology Matrix" (Proposition 3 — Euclidean realization; spectral
  embedding).
- `symmetry_groups_of_spin_systems_en.md`, §7 "Euclidean realization: the
  distorted simplex — literally".
- `classification_and_nomenclature_en.md` — the 3D-realizability flag.
- `code/analysis/README.md` — summary of the catalogue-analysis tools.
