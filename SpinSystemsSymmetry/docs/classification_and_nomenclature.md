# Classification and nomenclature of symmetric spin systems

*A classification scheme, canonical identifiers, and a grammar of structural formulas extending the Pople notation. Companion to `symmetry_groups_of_spin_systems_en.md` and to the machine verification `materials/8Spins_crossref.txt`.*

---

## 1. Three layers of classification

**Layer 0 — the canonical identifier.** The complete invariant of a symmetric spin system is its closed group as a permutation action, equivalently the canonical orbit colouring (a catalogue line). The identifier:

```
N/o<group order>.<line number within the section of NSpins.txt>
```

Example: `8/o6.5` — the cyclic C₆ eight-spin system. The ID is uninformative but unambiguous and stable (tied to the published catalogue files) — the CAS-number counterpart of a chemical name.

**Layer 1 — the physical structure (L1/L2).** Every system splits uniquely into magnetic composites (twin classes) and the factor graph with its factor group:

```
|Aut| = (prod_g n_g!) · |Gamma_q|
```

The structure is recursive: the factor graph is itself a smaller symmetric weighted graph. The layer-1 fingerprint (compositions n_g, |Gamma_q|, vertex orbits) is computed for the whole N = 8 catalogue in `8Spins_crossref.txt`.

**Layer 2 — taxonomic axes** (all machine-computable):
- group order and type (cyclic / dihedral / abelian / product / wreath / polyhedral);
- position in the tower: inherited entry (has a fixed point) or essentially new (`f(N)` of them);
- cyclic rank — the round of sequential symmetrization (minimal length of a chain of joins of cyclic closures);
- reality of the characters of Gamma_q — decides whether the engine needs the complex branch of the isotypic projection L2b;
- 3D realizability — whether a separating three-dimensional invariant subspace exists (the criterion of the `topology2xyz` tool).

## 2. The grammar of structural formulas

| construction | meaning | order factor |
|---|---|---|
| `An` | magnetic composite (n equivalent spins) | `n!` |
| `A'`, `AA'`, `A2A2'` | chemical orbits, primes as in Pople | 1 |
| `F*G` | independent assembly (generic couplings between parts) | product |
| `(F)n` | n interchangeable copies (wreath by Sₙ) | `|F|^n · n!` |
| `(F)_G` | copies of a fragment at positions with symmetry G | `|F|^k · |G|` |
| `G<F1, F2, ...>` | group G acting on the listed fragments jointly | `|G|` |
| `@u, @v` | attachment-stabilizer tags (formalizing `(-)`, `(\|)`) | 1 |
| `twn(L1,...,Lk)` | twisted stack of layers — chiral, group Cₙ | `n` |
| `stn^0 / stn^45` | stack at a special angle: eclipsed (□) / staggered (◊) | `2n` |
| `cube(A'8)` | polyhedral primitives (Oh, ...) | `|G|` |

**Checking rule:** the product of the order factors of all constructions must equal the entry's `|Aut|`. **Canonicity:** the formula is human-readable but not unique; only the pair "ID + canonical colouring" is canonical.

## 3. The renaming table of the N = 8 catalogue

Sorted by ID; the "PDF row" column refers to the 8Spins.pdf provenance table.

| ID | \|Aut\| | structural formula | 8Spins.pdf name | PDF row |
|---|---|---|---|---|
| `8/o1.1` | 1 | `1 (asymm.)` | A | 1 |
| `8/o2.1` | 2 | `C2<AA',MM',PP',XX'>` | AA'MM'PP'XX' | 79 |
| `8/o2.2` | 2 | `C2<AA',MM',XX'>` | AA'MM'XX' | 65 |
| `8/o2.3` | 2 | `C2<AA',XX'>` | AA'XX' | 35 |
| `8/o2.4` | 2 | `A2` | A2 | 2 |
| `8/o3.1` | 3 | `tw3(A'3,X'3)` | C3-(A'X')3 | 72 |
| `8/o4.1` | 4 | `D2<A'4; MM'@u, XX'@v>` | D2-(A'4,(-)MM',(|)XX') | 67 |
| `8/o4.2` | 4 | `D2<(A'X')4>` | D2/C2v-(A'X')4 | 80 |
| `8/o4.3` | 4 | `D2<A'4; MM'@u, XX'@u>` | (D2-A'4,(-)MM'(-)XX') | 49 |
| `8/o4.4` | 4 | `D2<A'4,XX'>` | D2-(A'4,XX') | 41 |
| `8/o4.5` | 4 | `C2<AA',MM'>*C2<PP',XX'>` | (AA'MM')(PP'XX') | 45 |
| `8/o4.6` | 4 | `D2<A'4>` | D2-A'4 | 36 |
| `8/o4.7` | 4 | `C2<AA',MM',PP'>*X2` | AA'MM'PP'X2 | 20 |
| `8/o4.8` | 4 | `C2<AA',MM'>*X2` | AA'MM'X2 | 6 |
| `8/o4.9` | 4 | `A2*X2` | A2X2 | 3 |
| `8/o4.10` | 4 | `tw4(A'4,X'4)` | C4-(A'X')4 | 89 |
| `8/o6.1` | 6 | `D3<A'6>` | D3/C3v-A'6 | 66 |
| `8/o6.2` | 6 | `C3v<A'6,XX'>` | C3v-(A'6,XX') | 73 |
| `8/o6.3` | 6 | `st3(A'3,X'3)` | D3/C3v-(A'X')3 | 43 |
| `8/o6.4` | 6 | `A3` | A3 | 34 |
| `8/o6.5` | 6 | `tw3(A'3,M'3)*X2` | C3-(A'M')3X2 | 86 |
| `8/o8.1` | 8 | `C2<A2A2',XX'>` | C2-(A2A2'XX') | 71 |
| `8/o8.2` | 8 | `C2<A2A2',MM',XX'>` | C2-(A2A2'MM'XX') | 85 |
| `8/o8.3` | 8 | `st4^0(A'4,X'4)` | [sq] D4/C4v-(A'X')4 | 60 |
| `8/o8.4` | 8 | `st4^45(A'4,X'4)` | [dia] D4/C4v-(A'X')4 | 70 |
| `8/o8.5` | 8 | `D2<A'4>*X2` | (D2-A'4)X2 | 7 |
| `8/o8.6` | 8 | `D2<A'4>*C2<MM',XX'>` | (D2-A'4)(MM'XX') | 37 |
| `8/o8.7` | 8 | `D2h<A'8>` | D2h-A'8 | 81 |
| `8/o8.8` | 8 | `D4<A'8>` | D4/C4v-A'8 | 84 |
| `8/o8.9` | 8 | `D2<A'4,MM'>*X2` | (D2-A'4,MM')X2 | 9 |
| `8/o8.10` | 8 | `D2h<(A'X')4>` | D2h-(A'X')4 | 46 |
| `8/o8.11` | 8 | `C2<AA',MM'>*P2*X2` | AA'MM'P2X2 | 47 |
| `8/o8.12` | 8 | `A2*M2*X2` | A2M2X2 | 40 |
| `8/o8.13` | 8 | `D4<A'4>` | D4-A'4 | 62 |
| `8/o10.1` | 10 | `D5<A'5>` | D5-A'5 | 64 |
| `8/o12.1` | 12 | `C2<AA',MM'>*X3` | AA'MM'X3 | 75 |
| `8/o12.2` | 12 | `st3(A'3,M'3)*X2` | (D3/C3v-(A'M')3)X2 | 51 |
| `8/o12.3` | 12 | `D3<A'6>*X2` | (D3/C3v-A'6)X2 | 68 |
| `8/o12.4` | 12 | `D3h<A'6,XX'>` | D3h/D3d-(A'6XX') | 87 |
| `8/o12.5` | 12 | `D6<A'6>` | D6/D3h/D3d-A'6 | 74 |
| `8/o12.6` | 12 | `A3*X2` | A3X2 | 63 |
| `8/o14.1` | 14 | `D7<A'7>` | D7-A'7 | 78 |
| `8/o16.1` | 16 | `D8<A'8>` | D8-A'8 | 90 |
| `8/o16.2` | 16 | `D2<A'4>*M2*X2` | (D2-A'4)M2X2 | 42 |
| `8/o16.3` | 16 | `D4<A'4>*C2<MM',XX'>` | (D4-A'4)(MM'XX') | 48 |
| `8/o16.4` | 16 | `A2*M2*P2*X2` | A2M2P2X2 | 8 |
| `8/o16.5` | 16 | `D2<A'4, M2, M2'>` | (D2-A'4,M2,M2') | 54 |
| `8/o16.6` | 16 | `D4h<A'8>` | D4h-A'8 | 83 |
| `8/o16.7` | 16 | `C2<A2A2',MM'>*X2` | C2-(A2A2'MM')X2 | 23 |
| `8/o16.8` | 16 | `D2<A'4>*D2<X'4>` | (D2-A'4)(D2-X'4) | 38 |
| `8/o16.9` | 16 | `D4<A'4>*X2` | (D4-A'4)X2 | 12 |
| `8/o20.1` | 20 | `D5<A'5>*X2` | (D5-A'5)X2 | 76 |
| `8/o24.1` | 24 | `D6<A'6>*X2` | (D6-A'6)X2 | 30 |
| `8/o24.2` | 24 | `D2<A'4>*X3` | (D2-A'4)X3 | 44 |
| `8/o24.3` | 24 | `(A'X')4` | (A'X')4 | 56 |
| `8/o24.4` | 24 | `A2*M2*X3` | A2M2X3 | 16 |
| `8/o24.5` | 24 | `A4` | A4 | 4 |
| `8/o32.1` | 32 | `(D2<A'4>)2` | (D2-A'4)2 | 61 |
| `8/o32.2` | 32 | `D4<A'4>*D2<X'4>` | (D4-A'4)(D2-X'4) | 39 |
| `8/o32.3` | 32 | `(A2)2*M2*X2` | (A2)2M2X2 | 29 |
| `8/o32.4` | 32 | `D4<A'4>*M2*X2` | (D4-A'4)M2X2 | 52 |
| `8/o36.1` | 36 | `A3*X3` | A3X3 | 11 |
| `8/o48.1` | 48 | `cube(A'8)` | Oh-A'8 | 59 |
| `8/o48.2` | 48 | `D4<A'4>*X3` | (D4-A'4)X3 | 77 |
| `8/o48.3` | 48 | `C2<AA',MM'>*X4` | AA'MM'X4 | 21 |
| `8/o48.4` | 48 | `(A2)3` | (A2)3 | 15 |
| `8/o48.5` | 48 | `A2*X4` | A2X4 | 10 |
| `8/o60.1` | 60 | `D5<A'5>*X3` | (D5-A'5)X3 | 88 |
| `8/o64.1` | 64 | `(A2)_D2` | D2-(A2)'4 | 82 |
| `8/o64.2` | 64 | `D4<A'4>*D4<X'4>` | (D4-A'4)(D4-X'4) | 53 |
| `8/o72.1` | 72 | `A2*M3*X3` | A2M3X3 | 22 |
| `8/o72.2` | 72 | `(A3)2` | (A3)2 | 14 |
| `8/o72.3` | 72 | `C2<A3A3',XX'>` | A3A3'XX' | 25 |
| `8/o96.1` | 96 | `(A2)3*X2` | (A2)3X2 | 26 |
| `8/o96.2` | 96 | `D2<A'4>*X4` | (D2-A'4)X4 | 69 |
| `8/o96.3` | 96 | `A3*M2*X2` | A3M2X2 | 50 |
| `8/o120.1` | 120 | `A5` | A5 | 5 |
| `8/o128.1` | 128 | `(A2)_D4` | D4-(A2)'4 | 33 |
| `8/o144.1` | 144 | `(A3)2*X2` | (A3)2X2 | 55 |
| `8/o144.2` | 144 | `A3*X4` | A3X4 | 18 |
| `8/o192.1` | 192 | `D4<A'4>*X4` | (D4-A'4)X4 | 28 |
| `8/o240.1` | 240 | `A2*X5` | A2X5 | 17 |
| `8/o384.1` | 384 | `(A2)4` | (A2)4 | 58 |
| `8/o576.1` | 576 | `A4*X4` | A4X4 | 57 |
| `8/o720.1` | 720 | `A3*X5` | A3X5 | 27 |
| `8/o720.2` | 720 | `A6` | A6 | 13 |
| `8/o1152.1` | 1152 | `(A4)2` | (A4)2 | 32 |
| `8/o1440.1` | 1440 | `A2*X6` | A2X6 | 24 |
| `8/o5040.1` | 5040 | `A7` | A7 | 19 |
| `8/o40320.1` | 40320 | `A8` | A8 | 31 |


## 4. Remarks

- The pictograms of the original PDF are translated into formal constructions: `□ → st4^0`, `◊ → st4^45`, `(-)/(|) → @u/@v`; the decodings are documented in `8Spins_crossref.txt`.
- The single data correction found during verification: `8/o6.5 = tw3(A'3,M'3)*X2` is the cyclic C₆ (the D3 label in the txt section header was fixed; the PDF name, `C3-(A'M')3X2`, was already correct).
- Full machine fingerprints of all 90 entries (twins, |Gamma_q|, orbits) are in `materials/8Spins_crossref.txt`; the taxonomy of the N = 3..13 catalogues (3246 entries) is in `materials/catalog_taxonomy.tsv`.
