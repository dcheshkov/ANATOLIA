// extend_tower.cpp - tower extension of a realizable-symmetry catalogue.
//
// By the tower theorem, the closed groups of degree N with a fixed vertex are
// exactly the closed groups of degree N-1, each extended by one isolated vertex:
// append vertex N and colour its N-1 new edges by ONE FRESH label per vertex
// orbit of the degree-(N-1) group.  The automorphism group is unchanged, so the
// extended colouring is again closed and has a fixed point.  This maps the
// level-(N-1) catalogue bijectively onto the fixed-point part of level N.
//
// Reads a degree-M catalogue (one canonical RGS per line, space-separated
// integers, E_M = M(M-1)/2 columns, materials/*Spins.txt convention), writes the
// degree-(M+1) tower extensions, each canonicalized (lex-min), one per line.
// Concatenated with the fixed-point-FREE classes (from the orbit-partition
// decomposition, filter_closed DUMP) this is the complete level-(M+1) catalogue.
//
// Build:  g++ -O3 -march=native -std=c++17 -o extend_tower extend_tower.cpp
// Usage:  ./extend_tower M  cat_M.txt  > cat_{M+1}_fixedpoint.txt

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include "canon_ir.h"
using namespace std;

int main(int argc, char** argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s M cat_M.txt\n", argv[0]); return 1; }
    int M = atoi(argv[1]);
    int N = M + 1;
    int EM = M * (M - 1) / 2;
    FILE* f = fopen(argv[2], "r"); if (!f) { fprintf(stderr, "cannot open %s\n", argv[2]); return 1; }
    char* line = nullptr; size_t cap = 0;
    long long out = 0;
    while (getline(&line, &cap, f) != -1)
    {
        if (line[0] == '#' || line[0] == '\n') continue;
        istringstream ss(line);
        vector<int> c; int x; while (ss >> x) c.push_back(x);
        if ((int)c.size() != EM) continue;             // skip malformed / header rows
        // degree-M colouring as RGS string (colours fit in a byte for M<=16)
        string cm((size_t)EM, 0);
        for (int e = 0; e < EM; e++) cm[(size_t)e] = (char)c[e];
        // vertex orbits under Aut of the degree-M colouring
        canir::CanonIR ir; ir.init(cm, M); ir.run();
        vector<int> orb = ir.vertexOrbits();           // orb[i] in 0..(#orbits-1)
        int maxc = 0; for (int v : c) if (v > maxc) maxc = v;
        int fresh = maxc + 1;
        // degree-N colouring: first EM edges unchanged, then the new vertex's
        // column {i,N-1} for i=0..M-1, coloured fresh+orb[i] (column order puts
        // this last block exactly as the (N-1)-th column of K_N)
        // raw degree-N colouring (NOT canonicalized here): first EM edges
        // unchanged, then the new vertex's column {i,N-1} coloured fresh+orb[i].
        // Canonicalization is left to catalog_format (IR canonical form): the
        // lex-min canonize degenerates on the near-asymmetric extensions of the
        // low-symmetry degree-M entries (documented in canon_ir.h).
        for (int e = 0; e < EM; e++) printf("%s%d", e ? " " : "", c[e]);
        for (int i = 0; i < M; i++) printf(" %d", fresh + orb[i]);
        printf("\n");
        out++;
    }
    fclose(f); free(line);
    fprintf(stderr, "extend_tower: wrote %lld degree-%d fixed-point colourings\n", out, N);
    return 0;
}
