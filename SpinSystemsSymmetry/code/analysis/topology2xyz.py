#!/usr/bin/env python3
# topology2xyz.py — symmetric 3D visualization of a spin-system topology.
#
# Takes one catalogue line (edge labels of K_N in column order: J12 J13 J23
# J14 J24 J34 ...) and produces 3D coordinates in which the FULL automorphism
# group of the weighted graph acts by exact isometries.  Method: eigenspaces
# of the class-weighted adjacency matrix (generic distinct weight per label
# class) are Aut-invariant, so any sum of whole eigenspaces of total
# dimension 3 yields a symmetry-exact embedding; among all such sums the tool
# picks the one maximizing the minimal inter-vertex distance.  The 3D action
# is always isometric but need not be faithful (e.g. D2xD2 on two tetrads is
# not a subgroup of O(3)); the report states how many distinct isometries are
# realized.
#
# Outputs (basename.<ext>):
#   .xyz     spins as pseudo-atoms, element = vertex orbit under Aut
#            (view in Avogadro / VMD / Jmol / ChemCraft ...)
#   .mol2    same atoms + bonds for the selected edge classes
#   .txt     legend: classes, sizes, which are drawn, orbit->element map
#
# Usage:
#   python3 topology2xyz.py --labels "0 1 1 1 1 2 1 1 2 2 3 3 4 5 6 3 3 6 4 5 7 3 3 5 6 4 7 7" --out c6
#   python3 topology2xyz.py --file ../materials/8Spins.txt --order 6 --index 5 --out c6
# Options:
#   --bonds auto|all|<k1,k2,...>   which label classes become bonds
#                                  (auto = classes with <= N edges; default)
#   --no-verify                    skip Aut computation (N >= 10), orbits
#                                  then fall back to coupling-profile classes
#   --scale <d>                    minimal inter-vertex distance in output
#                                  units (default 1.5, "angstrom-like")
#
# Requires: numpy.

import argparse
import re
import sys
from itertools import permutations, combinations

import numpy as np

ELEMENTS = ["C", "N", "O", "S", "P", "F", "Cl", "B", "Br", "I", "Se", "Si"]


def infer_n(nlabels):
    n = int((1 + (1 + 8 * nlabels) ** 0.5) / 2)
    if n * (n - 1) // 2 != nlabels:
        sys.exit(f"label count {nlabels} is not N(N-1)/2 for integer N")
    return n


def build_edges(N):
    edges = [(i, j) for j in range(1, N) for i in range(j)]
    eidx = {}
    for k, (i, j) in enumerate(edges):
        eidx[(i, j)] = k
        eidx[(j, i)] = k
    return edges, eidx


def rgs(c):
    lab, out = {}, []
    for v in c:
        if v not in lab:
            lab[v] = len(lab)
        out.append(lab[v])
    return tuple(out)


def automorphisms(c, N, edges, eidx):
    E = len(edges)
    out = []
    for p in permutations(range(N)):
        ok = True
        for e in range(E):
            i, j = edges[e]
            if c[eidx[(p[i], p[j])]] != c[e]:
                ok = False
                break
        if ok:
            out.append(p)
    return out


def vertex_orbits(A, N):
    parent = list(range(N))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for p in A:
        for v in range(N):
            a, b = find(v), find(p[v])
            if a != b:
                parent[b] = a
    roots = {}
    orb = [0] * N
    for v in range(N):
        r = find(v)
        if r not in roots:
            roots[r] = len(roots)
        orb[v] = roots[r]
    return orb


def profile_classes(c, N, eidx):
    # fallback vertex classification without Aut: multiset of incident labels
    prof = {}
    orb = [0] * N
    for v in range(N):
        key = tuple(sorted(c[eidx[(v, u)]] for u in range(N) if u != v))
        if key not in prof:
            prof[key] = len(prof)
        orb[v] = prof[key]
    return orb


def embed(c, N, edges):
    E = len(edges)
    m = max(c) + 1
    w = 1.0 + np.sqrt(np.arange(2, m + 2, dtype=float)) * 0.7318
    W = np.zeros((N, N))
    for e, (i, j) in enumerate(edges):
        W[i, j] = W[j, i] = w[c[e]]
    ev, U = np.linalg.eigh(W)
    groups = []
    for i in range(N):
        if groups and abs(ev[i] - ev[groups[-1][-1]]) < 1e-8:
            groups[-1].append(i)
        else:
            groups.append([i])

    def best_of_dim(D):
        best = None
        for r in range(1, D + 1):
            for sel in combinations(range(len(groups)), r):
                cols = [j for s in sel for j in groups[s]]
                if len(cols) != D:
                    continue
                X = U[:, cols]
                dmin = min(np.linalg.norm(X[i] - X[j])
                           for i in range(N) for j in range(i + 1, N))
                if best is None or dmin > best[0]:
                    best = (dmin, sel, X)
        return best

    SEP = 1e-6
    best3 = best_of_dim(3)
    if best3 and best3[0] > SEP:
        dmin, sel, X = best3
        return X, dmin, [len(groups[s]) for s in sel], True

    # No separating 3-dim invariant subspace exists: the group is not a 3D
    # point-group action on these vertices.  Find the smallest separating
    # D-dim symmetric embedding and fold the extra coordinates into the
    # third axis with small deterministic coefficients: vertices separate,
    # the picture stays exact for the subgroup fixing the folded directions.
    for D in range(4, N):
        bD = best_of_dim(D)
        if bD and bD[0] > SEP:
            dmin, sel, X = bD
            Y = X[:, :3].copy()
            for k in range(3, D):
                Y[:, 2] += (0.35 / (k - 1)) * X[:, k]
            dmin = min(np.linalg.norm(Y[i] - Y[j])
                       for i in range(N) for j in range(i + 1, N))
            return Y, dmin, [len(groups[s]) for s in sel], False
    sys.exit("no separating invariant subspace found (unexpected)")


def isometry_report(X, A, tol=1e-6):
    qs = set()
    resmax = 0.0
    nexact = 0
    scale = max(1.0, float(np.abs(X).max()))
    for p in A:
        P = np.zeros((len(p), len(p)))
        for v in range(len(p)):
            P[p[v], v] = 1
        M = X.T @ P @ X
        Uu, _, Vt = np.linalg.svd(M)
        Q = Uu @ Vt
        r = float(np.abs(P @ X - X @ Q).max())
        resmax = max(resmax, r)
        if r < tol * scale:
            nexact += 1
            qs.add(tuple(np.round(Q, 6).ravel()))
    return resmax, nexact, len(qs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", help="edge labels, space separated")
    ap.add_argument("--file", help="catalogue file (8Spins.txt format)")
    ap.add_argument("--order", type=int, help="section order in the catalogue file")
    ap.add_argument("--index", type=int, help="1-based line index inside the section")
    ap.add_argument("--out", default="topology", help="output basename")
    ap.add_argument("--bonds", default="auto")
    ap.add_argument("--no-verify", action="store_true")
    ap.add_argument("--scale", type=float, default=1.5)
    args = ap.parse_args()

    if args.labels:
        c = [int(x) for x in args.labels.split()]
    elif args.file and args.order and args.index:
        cur, pos, c = None, 0, None
        for line in open(args.file):
            if line.startswith("---"):
                cur = int(re.search(r"order (\d+)", line).group(1))
                pos = 0
                continue
            if line.strip() and cur == args.order:
                pos += 1
                if pos == args.index:
                    c = [int(x) for x in line.replace("*", "").split()]
                    break
        if c is None:
            sys.exit("catalogue line not found")
    else:
        sys.exit("give either --labels or --file/--order/--index")

    N = infer_n(len(c))
    edges, eidx = build_edges(N)
    c = list(rgs(tuple(c)))
    E = len(edges)

    verify = not args.no_verify and N <= 9
    if verify:
        A = automorphisms(c, N, edges, eidx)
        orb = vertex_orbits(A, N)
    else:
        A = []
        orb = profile_classes(c, N, eidx)

    X, dmin, chosen, exact = embed(c, N, edges)
    X = X - X.mean(axis=0)
    X = X * (args.scale / dmin)

    print(f"N={N}, classes={max(c)+1}", end="")
    if verify:
        res, nexact, ndistinct = isometry_report(X, A)
        if exact:
            faithful = ("faithful" if ndistinct == len(A)
                        else f"NOT faithful ({ndistinct} of {len(A)} isometries distinct)")
            print(f", |Aut|={len(A)}, eig-dims chosen={chosen}, "
                  f"isometry residual={res:.1e}, 3D action {faithful}")
        else:
            print(f", |Aut|={len(A)}, eig-dims chosen={chosen}: no separating "
                  f"3-dim invariant subspace (group is not a 3D point-group "
                  f"action on these vertices) - higher-dim symmetric embedding "
                  f"folded to 3D, {nexact} of {len(A)} symmetries remain exact")
    else:
        print(f", eig-dims chosen={chosen} (Aut verification skipped)"
              + ("" if exact else "; folded higher-dim embedding"))

    # bond class selection
    from collections import defaultdict
    cls = defaultdict(list)
    for e in range(E):
        cls[c[e]].append(edges[e])
    if args.bonds == "all":
        drawn = sorted(cls)
    elif args.bonds == "auto":
        drawn = sorted(k for k in cls if len(cls[k]) <= N)
    else:
        drawn = sorted(int(x) for x in args.bonds.split(","))

    elems = [ELEMENTS[orb[v] % len(ELEMENTS)] for v in range(N)]

    with open(args.out + ".xyz", "w") as f:
        f.write(f"{N}\n{args.out}: spin-system topology, element = vertex orbit\n")
        for v in range(N):
            f.write(f"{elems[v]:<3} {X[v,0]:12.6f} {X[v,1]:12.6f} {X[v,2]:12.6f}\n")

    bonds = [(i, j, k) for k in drawn for (i, j) in cls[k]]
    with open(args.out + ".mol2", "w") as f:
        f.write("@<TRIPOS>MOLECULE\n" + args.out + f"\n{N} {len(bonds)} 0 0 0\nSMALL\nNO_CHARGES\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for v in range(N):
            f.write(f"{v+1:>4} {elems[v]}{v+1:<3} {X[v,0]:10.4f} {X[v,1]:10.4f} "
                    f"{X[v,2]:10.4f} {elems[v]:<3} 1 SPIN 0.0\n")
        f.write("@<TRIPOS>BOND\n")
        for bi, (i, j, k) in enumerate(bonds, 1):
            f.write(f"{bi:>4} {i+1:>4} {j+1:>4} 1\n")

    with open(args.out + ".txt", "w") as f:
        f.write(f"{args.out}: legend\n")
        f.write(f"vertex orbits -> elements: " +
                ", ".join(f"orbit {o}={ELEMENTS[o % len(ELEMENTS)]}"
                          for o in sorted(set(orb))) + "\n")
        for k in sorted(cls):
            mark = "drawn" if k in drawn else "omitted"
            f.write(f"class {k}: {len(cls[k])} edges [{mark}]  "
                    + " ".join(f"({i+1},{j+1})" for i, j in cls[k]) + "\n")

    print(f"written: {args.out}.xyz, {args.out}.mol2, {args.out}.txt "
          f"(bond classes: {drawn})")


if __name__ == "__main__":
    main()
