# Symmetry groups of spin systems: realizability criterion, classification, and enumeration

*A theory document: which subgroups of the symmetric group can occur as symmetry groups of spin systems with scalar couplings, why, and how to enumerate all of them with a proof of completeness. Companion to the methodological specification of the block-diagonalization pipeline (not included in this deposit; that document explains how to *use* the symmetry for block diagonalization — this one explains which symmetries *exist* in the first place).*

---

## 0. The object and the main question

A system of `N` spin-½ nuclei with scalar couplings is described by the Hamiltonian

```
H = Σ_i ω_i · Iz_i  +  Σ_{i<j} J_ij · (I_i · I_j)
```

and is a weighted undirected complete graph `G`: vertex `i` carries the weight `ω_i` (chemical shift), edge `{i,j}` carries the weight `J_ij`. A spin permutation `σ ∈ S_N` acts by a unitary operator `U_σ`, and

```
[H, U_σ] = 0   ⟺   ω_σ(i) = ω_i  and  J_σ(i)σ(j) = J_ij  ∀ i,j   ⟺   σ ∈ Aut(G).
```

The symmetry group of a spin system is precisely the **automorphism group of the weighted graph**. As in the methodological specification, symmetry is a discrete property: only the equality classes of the weights (their indices) matter, not the numerical values.

> **Main question.** Let `H ≤ S_N` be a subgroup. Does there exist an `N`-spin system whose `Aut(G) = H` — *equal*, not merely containing `H`? Such `H` will be called **realizable**.

The answer turns out to be classical (Wielandt's theory of 2-closures), admits a short proof, explains the entire structure of the catalogue of symmetric topologies, and upgrades the empirical completeness criterion ("the next symmetrization round produces nothing new") to a theorem.

---

## 1. The most symmetric system and the geometric picture

The most symmetric system consists of `N` magnetically equivalent spins: all shifts equal, all `J` equal, `Aut(G) = S_N`. Geometrically this is the **regular simplex**: `N` vertices in a space of dimension `N−1`, with all pairwise "edge lengths" (coupling constants) equal.

Every other symmetric system is an **edge distortion of the simplex that preserves part of the symmetry**: a descent from `S_N` to its subgroups by splitting the equality classes of the weights. This picture is not a metaphor but an exact model (Theorem 5, §7): the realizable groups are precisely the isometry groups of `N`-point configurations in Euclidean space `R^{N−1}`.

The key restriction is visible immediately: not every subgroup is "reachable by distortion". Lowering the symmetry by continuously varying the weights, one can only land on subgroups that are *stabilizers of points* of the weight space — and on no others.

---

## 2. The Galois correspondence and the realizability criterion

**Orbit colouring.** To a subgroup `H ≤ S_N` assign the colouring `c_H`: two pairs `{i,j}` and `{k,l}` (and two vertices `i`, `k`) receive the same label if and only if they lie in the same orbit of `H`. This is the finest labelling of the weights for which every element of `H` is an automorphism.

**Closure.** Set `H^(2*) = Aut(c_H)` — the group of all permutations preserving every orbit of `H` on pairs. Clearly `H ≤ H^(2*)`; the map `H ↦ H^(2*)` is monotone and idempotent — a **closure operator** (the symmetrized, unordered-pairs variant of Wielandt's 2-closure). Groups with `H = H^(2*)` are called **2*-closed** (or simply closed).

> **Theorem 1 (realizability criterion).** `H ≤ S_N` is realizable as the full symmetry group of a spin system if and only if `H = H^(2*)`. Moreover, `H^(2*)` is the smallest realizable group containing `H`.
>
> *Proof.* (⇐) Assign pairwise distinct real weights to the classes of `c_H`. An automorphism of the weighted graph preserves each weight value, i.e. each colour class, hence `Aut(G) = Aut(c_H) = H^(2*) = H`. (⇒) Let `H = Aut(G)` for some weight assignment `w`. The equality classes of `w` are unions of `H`-orbits. The group `H^(2*)` preserves every `H`-orbit, hence preserves `w` pointwise, hence `H^(2*) ≤ Aut(G) = H`. Together with `H ≤ H^(2*)` this gives `H = H^(2*)`. ∎

Structural consequences:

- **Lattice.** An intersection of closed groups is closed; the join of two closed groups is the closure of the subgroup they generate. The realizable groups form a lattice — which is exactly the lattice of symmetry strata of parameter space (§5).
- **Terminology for literature search.** Orbit colourings are the "Schurian" symmetric coherent configurations (D. G. Higman); the closure operator is the 2-closure (H. Wielandt, 1969); algorithmic colour stabilization is Weisfeiler–Leman.
- **Contrast with Frucht's theorem.** Every finite group is abstractly isomorphic to `Aut` of some graph (Frucht, 1939), but our question concerns a concrete *action* on the same `N` vertices, and the answer is restrictive: only closed actions occur.

---

## 3. Vertex weights do not extend the class of realizable groups (N ≥ 3)

It is convenient to keep the topology catalogue by edge labels only (the coupling matrix), treating the shift index as derived — the vertex orbit. No generality is lost:

> **Proposition 2.** For `N ≥ 3` the set of realizable groups is the same whether vertices may be distinguished by shifts or only edges by `J` weights. The single exception is `N = 2`: the trivial group is realizable only by distinct shifts (the transposition of a single edge cannot be removed).
>
> *Proof sketch.* Let `H_e^(2*)` be the group of the purely edge orbit colouring of `H`. We show `H_e^(2*)` preserves every vertex orbit of `H` — then adding vertex labels does not change the closure. (i) An orbit `O` with `|O| ≥ 2`: the classes of internal edges of `O` are `H`-orbits of pairs with support exactly `O`; each class maps to itself (weights differ across classes), hence `O` maps to itself. (ii) Fixed points `u ≠ v`: every class of edges at `u` consists solely of edges incident to `u`; its image must stay in the same class, whence `ĝ(u) = u` — except in the degenerate case `N = 2`. ∎

Remark: for a *specific physical system* the shifts do of course participate in determining its group (a system with "all J equal, shifts split 4+4" has group `S₄ × S₄`), but every group arising this way is also realizable purely by edges — by a different system (for `S₄ × S₄`: two distinct intra-group `J` indices plus an inter-group one).

---

## 4. Structural consequences of the criterion

### 4.1. The chirality dichotomy: a single ring is always dihedral

> **Proposition 3.** Let `C_n` (`n ≥ 3`) act regularly on its `n`-element orbit (a single ring). The orbits on pairs inside the ring are the distance classes `{±d}`, and these are preserved by the reflection `i ↦ −i`. Hence the closure contains the dihedral group `D_n`: **purely rotational symmetry of a single ring is not realizable**.

Physical reading: a scalar coupling is undirected — an edge, not an arrow. A "chiral" permutation symmetry cannot be sustained by a scalar Hamiltonian: the reflection comes for free. This is why the `A′₈` ring in the catalogue yields `D₈` (order 16), not `C₈`. Directed spin interactions — e.g. the antisymmetric Dzyaloshinskii–Moriya exchange `D·(S_i × S_j)` in solid-state systems — change the object to a colouring of *ordered* pairs, and the class of realizable groups widens to Wielandt's 2-closed groups in the original (unsymmetrized) sense.

### 4.2. Yet chiral configurations exist: non-regular Cₙ

A cyclic group is realizable when it acts by **several** cycles with a "twisted" inter-cycle structure. The minimal example is `C₃` on 6 spins: two triangles with weights `a ≠ b` and three classes of matchings between them (the Latin-square structure of `Z₃`):

```
inter-triangle classes:  {14, 25, 36} → c₀,  {15, 26, 34} → c₁,  {16, 24, 35} → c₂
```

With distinct `c₀, c₁, c₂` the only automorphisms are simultaneous rotations of both triangles: `Aut = ⟨(123)(456)⟩ ≅ C₃`. A reflection would have to swap the classes `c₁ ↔ c₂` — but their weights differ. Geometrically these are two parallel triangular "decks" rotated by a generic angle: a **twisted antiprism** — a chiral point configuration (a chiral polygon is impossible; a chiral configuration is perfectly possible).

> **The p-cycle rule.** A closed group of prime order `p` requires `N ≥ 2p` (two `p`-cycles). Check against the tables of §8: order 3 first appears at `N = 6`; orders 5 and 7 first appear exactly at `N = 10` and `N = 14` (both `= 2p`, confirmed by the enumerator), while order 11 is still absent (it requires `N ≥ 22`), whereas `D₅` (order 10) is already present at `N = 5`, `D₇` (order 14) at `N = 7`, and `D₁₁` (order 22) appears exactly at `N = 11`.

> **Minimal realizations of cyclic groups.** The canonical realization of `C_n` is a pair of twisted `n`-gons (two layers): the intra-layer classes (distance classes) are automatically reflection-symmetric, so reflections can only be killed by the inter-layer "skew matchings" — the `n` difference classes `d ∈ Z_n`; pairwise distinct labels always suffice (labelings of `Z_n` with no axial `d ↦ c−d` and no shift `d ↦ d+r` symmetry), and for `n = 3` they are also necessary. Geometrically this is an antiprism at a generic twist angle. But two equal rings are the minimum only for prime `n` and for `C₄`; the general mechanism is a **family of quotient rings**: the generator rotates by `+1` several rings of sizes `m_t | n` with `lcm{m_t} = n` (faithfulness). There are exactly `gcd(m_s, m_t)` skew classes between rings of sizes `m_s` and `m_t`; any reflection acts on them as `k ↦ −k + c`, which for `gcd ≥ 3` cannot preserve the classes, while a ring of size `≥ 3` that belongs to no pair with `gcd ≥ 3` can be flipped independently (which is why `{4,3,3}` does not realize `C₁₂` — the square flips on its own). A family realizes `C_n` ⟺ every ring of size `≥ 3` belongs to a pair with `gcd ≥ 3` and the group of synchronized rotations `{(r_t) : r_s ≡ r_t mod gcd(m_s,m_t)}` has order exactly `n`. Hence for prime powers `n = p^e`: `μ*(C_{p^e}) = p^e + min{p^f : 3 ≤ p^f ≤ p^e}` — `C₃ → 6`, `C₄ → 8`, `C₅ → 10`, `C₇ → 14`, but `C₈ → 8+4 = 12` and `C₉ → 9+3 = 12`. The `N = 12` catalogue confirms this: `C₈` first appears at `N = 12` (orbits `8+4`, entry `12/o8.44`), `C₉` at `N = 12` (orbits `9+3`, `12/o9.3`; element orders `1,3,3,9,…` verified by an independent enumeration), and both are absent for `N ≤ 11` — the earlier `2p^e` estimate (16 and 18) for them is **refuted**. `C₆` additionally has a "ring + quotient triangle" realization `6+3` already at `N = 9` (`9/o6.2`). For **composite** `n` the minimum is an additive assembly of the prime-power components on disjoint supports, `μ*(C_n) = Σ_p μ*(C_{p^{e_p}})`: `C₆ = C₂×C₃` is closed already at `N = 8` (the action with orbits `2+3+3`: `cl⟨(12)(345)(678)⟩ = C₆`, verified; in the `N = 8` catalogue this is the order-6 entry with vertex orbits `3+3+2` — in 8Spins.txt the section header has been corrected to `order 6 (D3x4,C6)`, with the `C₆` line moved to the end of the section following the order-4 convention; in the 8Spins.pdf provenance table this is row 86, named `C3-(A'M')3X2`, generator #5188 = `(12)(345)(678)` — the fragment name with the C3 prefix is correct there, the full group being `C₃ × S₂ ≅ C₆`); `C₁₀ → 12` is confirmed by the catalogue (`12/o10.5`, orbits `5+5+2`); `C₁₂ → 14` (`{4,4,3,3}`) and `C₁₈ = C₂ × C₉ → 2 + 12 = 14`. Parity corollary: `μ*` is always even (except `C₁`), so new cyclic groups debut only at even `N` — the first appearances of `C₂,C₃,C₄,C₆,C₅,C₈,C₉,C₁₀,C₇,C₁₂,C₁₈` are `2,6,8,8,10,12,12,12,14,14,14`. The `N = 14` catalogue confirms the rule and its three predicted debuts `C₇, C₁₂, C₁₈` (all `μ* = 14`): a group of order 7 appears for the first time, exactly at `N = 2·7`.

### 4.3. Overly transitive groups collapse to S_N

> **Proposition 4.** If `H` is transitive on unordered pairs (2-homogeneous), its orbit colouring is trivial — all `J` equal — and `H^(2*) = S_N`.

Casualties include: `A_N` (all `N ≥ 3`); on 8 points — `PSL(2,7)` (order 168), `PGL(2,7)` (336), `AΓL(1,8)` (168), `AGL(3,2)` (1344), `A₈` (20160); on 5 points — `F₂₀` and `A₅`. Hence the characteristic "cliff" at the top of the `N = 8` catalogue: after orders 1152, 1440, 5040 the next entry is 40320.

### 4.4. Products on disjoint supports

If the coloured components are pairwise inequivalent (e.g. of different sizes), the group of the union is the direct product of the component groups: `D₅ × S₂` (order 20), `D₅ × S₃` (order 60 in the `N = 8` catalogue — a pentagon and a triangle), `S₇ × S₁` (order 5040). Identical components with identical labels give a wreath product: two indistinguishable triangles — `S₃ ≀ C₂` of order 72.

### 4.5. Twins: magnetic equivalence as a graph structure

A magnetically equivalent group of spins is a class of **twin vertices** (equal shifts, equal couplings to every external spin, equal mutual couplings; in graph theory — twin vertices, modular decomposition, Gallai 1967). The full automorphism group factorizes:

```
|Aut(G)| = (∏_g n_g!) · |Γ_q|
```

— the internal `S_{n_g}` over twin classes, times the factor group `Γ_q` of automorphisms of the quotient graph. This is exactly the two-level L1/L2 structure of the methodological specification: the Schur–Weyl contraction consumes `∏ n_g!`, orbit weighting and isotypic projection consume `Γ_q`. The computation of automorphism groups of edge-weighted (NMR) graphs and the observation `Γ ⊆ S_{n₁}×…×S_{n_m}` (the group lies in the Young subgroup of the equivalent-vertex classes) go back to Balasubramanian (1994), who solved the **forward** problem — `Aut` of a given graph by scanning `∏ nᵢ!` permutations under the test `A_{p_i p_j} = A_{ij}`, the same method as the author's `WeightedGraphSymmetry` engine; we answer the **inverse** question — which groups are realizable at all.

---

## 5. The lattice of closed groups and the stratification of parameter space

The space of all `N`-spin systems is `w ∈ R^N ⊕ R^{C(N,2)}` (shifts ⊕ couplings). For a subgroup `H`:

- `V_H = { w : H ≤ Aut(w) }` is a **linear subspace** (weights constant on orbits), with `dim V_H` = the number of orbits of `H` (vertex + edge). This is the number of independent parameters of a system with symmetry at least `H`.
- A point `w ∈ V_H` in **general position** has `Aut(w) = H^(2*)` — the closure, not `H`. A non-closed symmetry cannot be "seen" at any parameter values: the stratum of exact symmetry `H` is non-empty ⟺ `H` is closed.
- The strata are ordered by inclusion of closures; their lattice (meet — intersection, join — closure of the generated subgroup) is the lattice of realizable groups. Lowering symmetry by a small perturbation = moving to a lower stratum; a topology catalogue is a census of strata.

---

## 6. Decomposition of distortions into irreducible representations of S_N

The space of edge weights, as a representation of `S_N`, decomposes (for `N ≥ 4`) into three irreducible components:

```
R^{C(N,2)}  ≅  [N] ⊕ [N−1,1] ⊕ [N−2,2],
dim:            1     N−1       N(N−3)/2
```

- `[N]` — a uniform change of all `J`: does not lower the symmetry.
- `[N−1,1]` — "vertex-like" distortions `J_ij = a_i + a_j`: the stabilizer of a generic point is the Young subgroup `∏ S_{n_g}` given by the equality classes of `a_i`. Such distortions generate **only pure magnetic equivalence** (level L1).
- `[N−2,2]` — genuinely pairwise distortions: only these create a non-trivial factor group `Γ_q` (level L2). Example: the square `D₄` is not a Young subgroup, and its colouring (sides/diagonals) is not expressible as `a_i + a_j`.

In other words, the L1/L2 hierarchy of the methodological specification has an exact representation-theoretic meaning: L1 — symmetries already visible in the `[N−1,1]` component, L2 — those requiring `[N−2,2]`.

---

## 7. Euclidean realization: the distorted simplex, literally

> **Theorem 5.** The realizable symmetry groups of `N`-spin systems are precisely the isometry groups of full-dimensional `N`-point configurations in `R^{N−1}`, viewed as permutation groups of the points.
>
> *Proof idea.* The squared-distance matrix of the regular simplex is an interior point of the cone of Euclidean distance matrices (the centred Gram matrix is proportional to `I − J/N` and positive definite on `1⊥`), so all sufficiently small perturbations of the squared distances remain distance matrices of configurations of full affine rank `N−1`. Choose perturbations that separate the classes of the orbit colouring of a closed group `H` — the isometries of the configuration then coincide with `Aut` of the colouring, i.e. with `H`. Conversely, a distance-preserving permutation of the points of a full-dimensional configuration extends uniquely to an isometry of all of `R^{N−1}`. ∎

The intuition "a spin system = a distorted regular simplex" is thereby exact, including the chiral cases (§4.2): the twisted antiprism is a legitimate point configuration with group `C₃`.

*(Caveat: physical `J` values need not satisfy the metric axioms — positivity, the triangle inequality. This is immaterial for symmetry theory: only the combinatorics of equality classes matters, and the Euclidean model reproduces it in full.)*

---

## 8. Enumeration: the algorithm, the completeness proof, the tables

### 8.1. The symmetrization step

> **Lemma 6.** Let `c` be the orbit colouring of a subgroup `K`, and `σ ∈ S_N`. Coarsening `c` to `σ`-invariance (union-find: merge the classes of edges `e` and `σ(e)` until stable) yields the orbit colouring of `⟨K, σ⟩`.
>
> *Proof.* The classes of the result are the connected components of the relation generated by "same colour in `c`" (= `K`-orbits) and "`e ↦ σ(e)`", i.e. the orbits of `⟨K, σ⟩`. ∎

### 8.2. Completeness of sequential symmetrization

> **Theorem 7 (completeness of the catalogue).** BFS: start from the discrete colouring (all weights distinct, trivial group); in each round apply Lemma 6 to all pairs (reached colouring, `σ ∈ S_N`) and canonicalize the result by vertex relabelling. After `k` productive rounds one has obtained the orbit colourings of **all** subgroups generated by `≤ k` elements. Since every subgroup of `S_N` is generated by at most `max(2, ⌊N/2⌋)` elements (McIver–Neumann, 1987), `⌊N/2⌋` rounds exhaust all closed groups.
>
> *Attribution remark.* The `⌊N/2⌋` bound is due to McIver and Neumann; Jerrum's result (Jerrum's filter, 1986) gives the weaker `N − 1`, which is nevertheless also sufficient for finite completeness.

For `N = 8` this means 4 rounds — and the enumerator (`code/engine/closed_groups.cpp`) confirms the prediction in detail:

```
round 1: +21 (total 22)   — closures of cyclic groups
round 2: +59 (total 81)
round 3: +8  (total 89)
round 4: +1  (total 90)   — exactly one group: (C₂)⁴ in its D₂×D₂ action
round 5: +0               — stop
```

The order histogram reproduces the catalogue `materials/8Spins.txt` exactly (31 orders, including the splitting of order 4 into `C₂×C₂ ×9` and `C₄ ×1`). For `N = 9` (`⌊9/2⌋ = 4` rounds) the pattern repeats: `30 → 121 → 130 → 131`, with the fourth round again contributing exactly one group; for `N = 10`: `42 → 243 → 278 → 282` (four groups in round 4); for `N = 11`: `56 → 345 → 389 → 394` (3.5 min on 48 cores) — the bound admits a fifth round for the first time (`⌊11/2⌋ = 5`), but it was not needed; for `N = 12`: `77 → 755 → 921 → 947 → 948` (101 min on 48 cores) — the fifth round is productive for the first time (+1); for `N = 13`: `101 → 1073 → 1284 → 1315 → 1316` (23.4 h on 96 cores of two nodes, MPI) — again a productive fifth round (+1), the sixth empty. The round-4 "straggler" at `N = 8, 9` is `(C₂)⁴` in its `D₂ × D₂` action: two regular Klein four-groups on two tetrads with equal inter-tetrad couplings (the `(D2-A'4)(D2-X'4)` row of the 8Spins.pdf provenance table). Curiously, the Young `(C₂)⁴` of four transpositions, although 4-generated as a group, is reached already in round 3: its orbits are reproduced by the even-weight subgroup of rank 3.

> **The partition skeleton: a bridge between round one and the decomposition of §8.5.** Partitions `λ ⊢ N` play a double role in the enumeration: in the first symmetrization round they index the conjugacy classes of `S_N` — the cycle types of the generators — while in the decomposition of §8.5 they index the vertex-orbit structures ("decks") of the target groups. The link is exact: the closure preserves every vertex orbit, and the orbits of a cyclic group are the cycles of its generator, so the orbit partition of `cl⟨σ⟩` equals the cycle type of `σ`. Hence the catalogue after round one is in **bijection with the partitions** of `N` (conjugate generators give one entry, distinct `λ` give distinct orbit structures), and the round-one total is exactly `p(N)`: `5, 7, 11, 15, 22, 30, 42, 56, 77, 101` for `N = 4…13` — verified by the enumerator at every level. Moreover, `cl⟨σ_λ⟩` is block-transitive on its decks, i.e. it is the cyclic **seed** of stratum `λ` of the §8.5 decomposition: the later BFS rounds grow each stratum from below (the deck structure of a join is the join of the generators' cycle partitions as set partitions), while the decomposition enumerates it from above, inside the Young subgroup `S_λ`. Finally, the tower `a(N) = a(N−1) + f(N)` mirrors, at skeleton level, the elementary partition identity `p(N) = p(N−1) + p_{≥2}(N)` (partitions with a unit part ↔ systems with a fixed spin, inherited from level `N−1`; partitions without units ↔ fixed-point-free classes): at `N = 14` this reads `135 = 101 + 34`, and exactly 34 strata enter the computation of `f(14)`.

### 8.3. How many subgroups survive

| `N` | subgroup classes of `S_N` (OEIS A000638) | closed (realizable) |
|---|---|---|
| 3 | 4   | **3**  |
| 4 | 11  | **8**  |
| 5 | 19  | **11** |
| 6 | 56  | **27** |
| 7 | 96  | **36** |
| 8 | 296 | **90** |
| 9 | 554 | **131** |
| 10 | 1593 | **282** |
| 11 | 3094 | **394** |
| 12 | 10723 | **948** |
| 13 | 20832 | **1316** |
| 14 | 75154 | **2866** |

With the initial terms `a(1) = a(2) = 1` the sequence of realizable classes reads `1, 1, 3, 8, 11, 27, 36, 90, 131, 282, 394, 948, 1316, 2866`; as of July 2026 it is **absent** from the OEIS and has been submitted there.

> **Two more enumeration routes.** Besides the constructive BFS, the catalogue can be reached by *enumerate-then-filter* (list all subgroup classes of `S_N`, keep the closed ones) and — the route that reaches `N = 14` — by the *orbit-partition decomposition*. Both, with a completeness proof of the decomposition, are in §8.5: the enumerate-then-filter route independently reproduces `a(3..13)`, and the decomposition delivers `a(14) = 2866` (`= 1316 + f(14)`, `f(14) = 1550`).

> **The tower of catalogues.** A closed group of degree `N−1`, extended by a fixed vertex (in the catalogue format — an appended column with one fresh label per vertex orbit), is closed of degree `N`; conversely, the restriction along a fixed vertex is closed (for `N−1 ≥ 3`, Proposition 2). The level-`N−1` catalogue is exactly the set of level-`N` closed groups with a fixed point, whence `a(N) = a(N−1) + f(N)`, where `f(N)` is the number of closed classes **without** fixed points (the "essentially N-spin" ones). Verified element by element for `N ≤ 8`; at `N = 3` the identity fails precisely because of the `N = 2` exception. Values: `f(4..14) = 5, 3, 16, 9, 54, 41, 151, 112, 554, 368, 1550` — even `N` are systematically richer (matching-based structures).

> **The removable-vertex criterion (the limit of vertex induction).** Call a vertex `v` of a level-`N` entry `G` *removable* if the colouring `G − v` is itself a canonical level-`(N−1)` entry (the orbit colouring of a closed group). (i) If `G` has a removable vertex, then `G` is built by attaching a "wheel" column to the entry `G − v` of the level-`(N−1)` catalogue (the column labels are `G`'s labels on the edges at `v`, reusing the labels of old classes). (ii) Fixed points of `Aut(G)` are always removable (the tower); the converse fails — in the `S_N` entry every vertex is removable. (iii) There exist entries with no removable vertex at all, not constructible by a wheel over any canonical entry (verified computationally): the `D₆` ring at `N = 6` — the residual after deleting any vertex has three classes, its group is `C₂`, but the orbit colouring of `C₂` has six classes, and any column over a coarsened base leads to `S₆`; and `cl⟨(12)(34)(56)(78)⟩` at `N = 8` — every residual is fully asymmetric, and symmetric columns over a coarsened base impose the extra symmetry `(78)`. **Corollary:** the induction "wheel over the symmetric level-`(N−1)` graphs" is incomplete — it is guaranteed to produce only the inherited part `a(N−1)`, whereas the genuinely new `f(N)` entries (those without fixed points) in general require a full level-`N` search.

Group orders for small `N` (multiplicities in parentheses):

```
N=3:  1, 2, 6
N=4:  1, 2(×2), 4(×2), 6, 8, 24
N=5:  1, 2(×2), 4(×2), 6, 8, 10, 12, 24, 120
N=6:  1, 2(×3), 3, 4(×4), 6(×3), 8(×4), 10, 12(×2), 16, 24, 36, 48(×2), 72, 120, 720
N=7:  1, 2(×3), 3, 4(×4), 6(×3), 8(×4), 10, 12(×3), 14, 16, 20, 24(×3), 36, 48(×3),
      72, 120, 144, 240, 720, 5040
N=8:  catalogue materials/8Spins.txt — 90 topologies, 31 orders from 1 to 40320
N=9:  131 classes, 45 orders from 1 to 362880 (incl. S₃≀S₃ of order 1296; no orders 5 or 7)
N=10: 282 classes, 64 orders from 1 to 3628800; order 5 appears for the first time (×1) —
      exactly at N = 2p, as the p-cycle rule predicts
N=11: 394 classes, 83 orders from 1 to 39916800; order 22 appears for the first time
      (D₁₁ — a single 11-ring); no orders 7 or 11
N=12: 948 classes, 124 orders from 1 to 479001600; order 9 appears for the first time
      (×3): two C₃×C₃ (orbits 3+3+3+3 and the regular action on 9+3) and one C₉
      (orbits 9+3 — refuting the 2p^e estimate); the cyclic C₈ (8+4) and C₁₀
      (5+5+2) also debut; no orders 7 or 11
N=13: 1316 classes, 152 orders from 1 to 6227020800; order 26 debuts (D13 - a
      single 13-ring, exactly at N = 13); no orders 7, 11 or 13 (they need
      N >= 14, 22, 26); no new cyclic group - the parity corollary (new cyclic
      groups debut only at even N) is confirmed
N=14: 2866 classes, 204 orders from 1 to 87178291200 (= 14!); order 7 debuts
      (C7, the first group of order 7 at all, exactly at N = 2*7); three new
      cyclic groups C7, C12, C18 (all mu* = 14) reconfirm the parity
      corollary (new cyclic groups only at even N); orders 11 and 13 still
      absent (need N >= 22, 26)
```

All structural predictions of §4 can be read off the tables directly: order 3 appears from `N = 6`, order 5 from `N = 10` and order 7 from `N = 14` (each at `N = 2p`), the dihedral 10 and 14 are present, and the "cliff" before `N!`.

### 8.4. The case N = 4 worked out: who is closed, who is not

`S₄` has eleven subgroup classes; eight are closed. Edge order — column by column of the upper triangle: `(J12 J13 J23 J14 J24 J34)`.

| Subgroup | Order | Closed? | Realizing colouring |
|---|---|---|---|
| `S₄` | 24 | yes | `a a a a a a` |
| `A₄` | 12 | **no** → `S₄` | (edge-transitive) |
| `D₄` (square 1-2-3-4) | 8 | yes | `a b a a b a` (sides `a`, diagonals `b`) |
| `S₃` (vertex 4 apart) | 6 | yes | `a a a b b b` |
| `C₄` | 4 | **no** → `D₄` | (same orbits as `D₄`) |
| `V₄` regular `⟨(12)(34),(13)(24)⟩` | 4 | yes | `a b c c b a` (three matchings) |
| `V₄' = ⟨(12),(34)⟩` | 4 | yes | `a c c c c b` |
| `C₃` | 3 | **no** → `S₃` | (same orbits as `S₃`) |
| `C₂ = ⟨(12)⟩` | 2 | yes | `c a a b b d` |
| `C₂' = ⟨(12)(34)⟩` | 2 | yes | `c a b b a d` |
| `1` | 1 | yes | all six distinct |

### 8.5. Two more routes: enumerate-then-filter and the orbit-partition decomposition

The sequential symmetrization of §8.1–8.2 *builds* the closed groups from below. Two further routes enumerate them differently; where they overlap they agree exactly, and one of them is what reaches `N = 14`.

**(A) Enumerate-then-filter.** The most direct route is also a complete one: list *every* conjugacy class of subgroups of `S_N` and keep the closed ones. GAP's `ConjugacyClassesSubgroups(SymmetricGroup(N))` (Hulpke's cyclic-extension algorithm; the class counts are the left column of §8.3, OEIS A000638) gives one representative `H` per class; for each we form its edge-orbit colouring `c_H` and test closure directly, `|Aut(c_H)| = |H|` (Theorem 1), with the individualization–refinement canonicalizer (§8.1). The closed count reproduces `a(N)` for **every** `N = 3..13`, and — a closed group *being* the automorphism group of its colouring — the survivors are automatically distinct. Two algorithms of opposite character, *symmetrize-then-BFS* and *enumerate-then-filter*, thus agree on the entire sequence. The route is limited by the lattice, not the filter: at `N = 14` GAP's cyclic-extension lattice exhausts memory in its `Zuppos` phase (prime-order cyclic subgroups) long before reaching the 75154 classes of `S₁₄`. Tools: `code/gap/emit_subgroup_colourings.g`, `code/gap/filter_closed.cpp`.

**(B) Orbit-partition decomposition.** This is how `a(14)` was obtained, and it never materializes the `S_N` lattice. By the tower (§8.3) it suffices to build `f(N)`, the fixed-point-free closed classes; the rest is the level-`(N−1)` catalogue with an isolated vertex appended. The structural fact that makes the assembly complete:

> **Proposition (completeness of the decomposition).** Fix, for each partition `λ = (λ₁, …, λ_k) ⊢ N` with every part `≥ 2`, a splitting of `{1, …, N}` into blocks of those sizes, and let `S_λ = S_{λ₁} × … × S_{λ_k}` be the corresponding Young subgroup; call `H ≤ S_λ` **block-transitive** if it is transitive on each block. Then the fixed-point-free 2*-closed groups of degree `N`, up to conjugacy, are exactly the 2*-closed block-transitive subgroups of the `S_λ` over all such `λ`. Hence `f(N)` is the number of these (deduplicated across conjugacy), and `a(N) = a(N−1) + f(N)`.
>
> *Proof.* Let `G` be fixed-point-free and 2*-closed of degree `N`. Its orbits `O₁, …, O_k` partition `{1, …, N}`; having no fixed point means every `|O_i| ≥ 2`, so `λ := (|O₁|, …, |O_k|)` is a partition of `N` with parts `≥ 2`. `G` stabilizes each orbit setwise, hence `G ≤ Sym(O₁) × … × Sym(O_k) ≅ S_λ`, and by the definition of orbit it is transitive on each `O_i` — i.e. `G` is block-transitive in `S_λ`. Conversely, a block-transitive `H ≤ S_λ` has orbits exactly the blocks, so its orbit-partition is `λ` (parts `≥ 2`, hence no fixed point); if `H` is 2*-closed it is one of the groups sought. The two families therefore coincide. For the count up to `S_N`-conjugacy: two block-transitive subgroups of `S_λ` that are conjugate in `S_N` are conjugate by an element of the normalizer `N_{S_N}(S_λ) = S_λ ⋊ (\text{permutations of equal-size blocks})`; deduplication by canonical colouring — a complete isomorphism/conjugacy invariant — identifies precisely these, counting each class once. Adding the `a(N−1)` fixed-point classes (the tower, §8.3) gives `a(N)`. ∎

It remains to enumerate the block-transitive subgroups of each `S_λ`. Every block-transitive `H` projects onto a transitive group on each block and is a **subdirect product** (Goursat) of those projections, which gives a cheap, library-driven enumeration:

- one block `λ = [N]` → the transitive groups of degree `N` from the tabulated library (for `N = 14`, 63 groups, of which **10** are 2*-closed — the single-orbit, fully coupled systems);
- two blocks `[a, b]` → the subdirect products of a transitive group on each block, read off the small transitive-group libraries — far cheaper than the full lattice of `S_a × S_b`, and decisive for the balanced partitions that otherwise dominate the cost;
- three or more blocks → `S_λ` is a product of small symmetric groups whose subgroup lattice is itself tractable.

> **Remark (projections do not inherit closedness; a two-tier organization of the search).** The components of the subdirect products must range over **all** transitive groups of the respective degrees: the projection of a 2*-closed group onto one of its orbits is transitive but in general **not** 2*-closed. Example — entry `6/o3.1` (the chiral `C₃` on 6 spins): both of its projections are `C₃` on three points, and `cl(C₃) = S₃`; with components restricted to 2*-closed groups this entry cannot be assembled at all — subdirectness demands surjective projections, and the diagonal `C₃` does not project onto `S₃`. Closure and the splitting into orbits do not commute, and the 2*-test is applied exactly once — to the assembled group; it is indispensable in its own right, since a subdirect product of 2*-closed components may fail to be closed (e.g. `S₃ ×_{C₂} S₃` of order 18 at `N = 6` — an order absent from the catalogue). At the same time the induced colouring of a single orbit sees only the **closure** of the projection: `T` and `cl(T)` have identical orbit colourings (idempotence), so what resides inside an orbit is always a 2*-closed transitive subgroup `D = cl(T)`, while the projection `T` itself is pinned down by the cross-orbit classes. Hence an equivalent **two-tier organization** of the search: for each orbit run over **all** 2*-closed transitive subgroups `D_i ≤ S_{λ_i}` (up to conjugacy, including `S_{λ_i}` itself), then a component `T_i` from `Surv(D_i) = { T ≤ D_i : T transitive, cl(T) = D_i }` (e.g. `Surv(S₃) = {S₃, C₃}`); the subdirect products of the tuples `(T₁, …, T_k)` then pass the final 2*-test and the dedup. Since every transitive group has exactly one closure, this search coincides with the search over all transitive components, grouped by closure. For `k = 1` (a single orbit) the scheme degenerates: there are no cross-orbit classes, a non-closed component `T ⊊ D` has nothing to pin it, the final test leaves only `T = D` — and the enumeration reduces to the list of the 2*-closed transitive groups themselves (at `N = 14`: 10 of the 63 transitive groups of degree 14); the rule "a single ring is always dihedral" is precisely this degenerate case.

Each candidate goes through the closure test (§8.1) and is deduplicated by canonical colouring. For `N = 14` this gives `f(14) = 1550` and `a(14) = 1316 + 1550 = 2866`. The decomposition is cross-checked two ways: it reproduces `f(4..13) = 5, 3, 16, 9, 54, 41, 151, 112, 554, 368` exactly, and for **every** two-block partition the full-Young-lattice and subdirect-product enumerations give byte-identical closed-colouring sets. Tools: `code/gap/emit_fixedfree_colourings.g`, `code/gap/emit_partition_colourings.g`, `code/gap/emit_subdirect2.g`, `code/gap/drive_fN.sh`; catalogue assembly `code/gap/extend_tower.cpp` + `code/gap/catalog_format.cpp`.

> **Remark (the fate of the alternating components: the revival of A₄).** The two-tier scheme suggests a natural question: which non-closed components `T ⊊ D` actually "come alive" in the assembled catalogue as projections? For the cyclic `C_n ∈ Surv(D_n)` the answer is the chiral family of §4.2; the next candidates are the natural alternating groups `A_n ∈ Surv(S_n)` (`A_n` is transitive on unordered pairs, its colouring of `K_n` is monochrome, hence `cl(A_n) = S_n`). A direct scan of all 6112 entries (for every twin-free orbit of size `n ≥ 3` the restrictions of `Aut` to the orbit are enumerated explicitly; a subgroup of order `n!/2` in `S_n` is automatically `A_n`, the unique index-2 subgroup, so a twin-carrying orbit is disqualified outright — a twin swap restricts to a transposition) gives a sharp answer. Projections `A₃ = C₃` are plentiful — 561 orbit instances, beginning with `6/o3.1` (the chiral étagères). Projections `A_n` with `n ≥ 5` **never** occur. And `A₄` revives exactly twice: entries `13/o12.10` (orbits `4+6+3`) and `14/o12.32` (the same plus an isolated spin), both with `Aut ≅ A₄` of order 12, faithful. The mechanism is the chiral twist lifted to a quotient group. All three orbits are built from one 4-element set: the 4 points themselves, the 6 pairs, the 3 partitions into two pairs. The point↔pair coupling is the incidence relation (two classes, locking the decks together); the point↔partition coupling is forced monochrome (`A₄` is transitive on such incidences); and the entire parity selection resides in the pair↔partition coupling, which carries **three** classes: "pair inside its partition" plus two distinguished transversal classes. The transposition `(ab)` fixes the partition `ab|cd` but exchanges `ac|bd ↔ ad|bc` — and with them these two classes; distinguishing them kills every odd element: `S₄` acts on the partition triple through `S₄/V₄ ≅ S₃`, the colouring pins that image down to `C₃`, and the preimage of `C₃` in `S₄` is exactly `A₄`. The trick is available to `A₄` alone — the only non-simple `A_n`: for `n ≥ 5` there is no small quotient, and the `A_n`-orbits on couplings between point- and pair-type orbits coincide with the `S_n`-orbits (a transposition of two points outside the configuration always exists; `n = 4` is precisely the tight case in which no free points remain and the class splits; larger `A_n`-actions — on `k`-subsets or coset spaces — no longer fit the `N ≤ 14` budget). Each deck of `13/o12.10` taken separately is a phantom: `cl(A₄ on points) = S₄`, `cl(C₃) = S₃`, `cl(A₄ on pairs) = C₂ ≀ S₃` of order 48 (the octahedral colouring of `K₆`); all three come alive only in the ensemble — the sharpest illustration in the catalogue of the remark above. Incidentally: the **regular** action of `A₄` is closed already at `N = 12` (the single-orbit entry `12/o12.7`; its neighbour `12/o12.9` is the regular `D₆`), but a natural 4-point orbit debuts only at `N = 13`. And since `A₄` is non-ambivalent (its 3-cycle classes split into two mutually inverse ones), both entries carry genuinely complex character tables — a ready-made showcase for the complex-character branch of L2b.

---

## 9. Consequences for spectroscopy

- **Input of the factorization pipeline.** The block-diagonalization algorithm (levels L0–L3 of the methodological specification) consumes `Aut(G)`; the present document describes *which* `Aut(G)` can occur. A catalogue of symmetric topologies is the lattice of closed groups, and its completeness is a theorem (7), not an empirical observation.
- **Complex characters are unavoidable.** Closed `C₃, C₄, C₆` exist (non-regular actions, §4.2), so the isotypic-projection branch L2b with complex characters is a mandatory part of a correct engine, not an exotic corner.
- **Dihedrality of rings** (§4.1) explains the composition of the catalogue: single rings contribute `D_n`, not `C_n`; purely rotational groups arise only from "twisted" multi-orbit actions.
- **The number of independent parameters** of a system with a given symmetry = the number of orbits of the group (the stratum dimension, §5) — a useful invariant when setting up a spectral fit.

---

## 10. Literature and terminology

- H. Wielandt. *Permutation groups through invariant relations and invariant functions.* Lecture notes, Ohio State University, 1969. — `k`-closures; our operator is the symmetrized variant of the 2-closure.
- P. McIver, P. M. Neumann. *Enumerating finite groups.* Quart. J. Math. Oxford (2) **38** (1987) 473–488. — `d(G) ≤ max(2, ⌊n/2⌋)` for `G ≤ S_n`.
- M. Jerrum. *A compact representation for permutation groups.* J. Algorithms **7** (1986) 60–78. — the `n − 1` bound (Jerrum's filter).
- R. Frucht. *Herstellung von Graphen mit vorgegebener abstrakter Gruppe.* Compositio Math. **6** (1939) 239–250. — every finite group is abstractly `Aut` of some graph.
- D. G. Higman. *Coherent configurations.* Geom. Dedicata **4** (1975) 1–32; B. Weisfeiler, A. Leman (1968). — coherent configurations, colour stabilization, "Schurian-ness".
- T. Gallai. *Transitiv orientierbare Graphen.* Acta Math. Acad. Sci. Hung. **18** (1967) 25–66. — modular decomposition (twins).
- C. N. Banwell, H. Primas. Mol. Phys. **6** (1963) 225. — composite particles in NMR; P. L. Corio. *Structure of High-Resolution NMR Spectra.* Academic Press, 1966.
- D. A. Cheshkov, K. F. Sheberstov, D. O. Sinitsyn, V. A. Chertkov. *ANATOLIA: NMR software for spectral analysis of total lineshape.* Magn. Reson. Chem. **56** (2018) 449. — the production context.
- OEIS A000638 — number of conjugacy classes of subgroups of `S_n`; D. F. Holt. *Enumerating subgroups of the symmetric group.* (2010) — enumeration up to `n = 18`.
- The GAP Group. *GAP – Groups, Algorithms, and Programming.* — the subgroup-class enumeration used for the independent check; A. Hulpke. *Constructing transitive permutation groups.* J. Symbolic Comput. **39** (2005) 1–30. — the cyclic-extension algorithm behind `ConjugacyClassesSubgroups`.
- K. Balasubramanian. *Computer generation of automorphism groups of weighted graphs.* J. Chem. Inf. Comput. Sci. **34** (1994) 1146–1150. — direct computation of `Aut` of edge-weighted (NMR) graphs; `Γ ⊆ S_{n₁}×…×S_{n_m}`, the `A_{p_i p_j}=A_{ij}` test.

Computational artifacts in this repository: the closed-group enumerator — `code/engine/closed_groups.cpp`; the `N = 8` catalogue — `materials/8Spins.txt` (the machine-verified cross-reference against the classical catalogue — `materials/8Spins_crossref.txt`).
