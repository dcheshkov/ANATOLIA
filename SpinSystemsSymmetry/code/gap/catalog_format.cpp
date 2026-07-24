// catalog_format.cpp - assemble a realizable-symmetry catalogue in the
// materials/*Spins.txt format from a stream of colourings.
//
// Reads degree-N edge colourings (space-separated integers, E = N(N-1)/2 per
// line; lines that do not parse to E integers are skipped, so headers and
// "--- order ---" separators of an existing catalogue pass through harmlessly).
// Each colouring is put into the IR canonical form (canon_ir.h, bestNorm) and
// tagged with |Aut|; duplicates (same canonical form) are removed; the distinct
// classes are grouped by |Aut|, orders ascending, and printed in the catalogue
// layout used by closed_groups.cpp -VERBOSE:
//
//   N=<N>: closed subgroup classes = <count> (<note>)
//     orders: o1(xk1) o2(xk2) ...
//   --- order o (xk) ---
//   <canonical colouring>            (space-separated RGS integers)
//   ...
//
// The IR canonical form is used (not the lex-min of the other catalogues)
// because lex-min canonize is impractical on the near-asymmetric N=14 colourings
// (canon_ir.h, lines ~494); it is an isomorphism-invariant canonical
// representative all the same, so the catalogue is well defined and every
// consumer (classify_catalog, mirror_catalog, topology2xyz) is
// representative-independent.
//
// Build:  g++ -O3 -march=native -fopenmp -std=c++17 -o catalog_format catalog_format.cpp
// Usage:  cat parts...  | ./catalog_format N "note"  > NSpins.txt

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <sstream>
#include "canon_ir.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

int main(int argc, char** argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s N [note]\n", argv[0]); return 1; }
    int N = atoi(argv[1]);
    const char* note = (argc > 2) ? argv[2] : "orbit-partition decomposition";
    int E = N * (N - 1) / 2;

    // ---- read all colourings, RGS-normalized ----
    vector<string> raw;
    char* line = nullptr; size_t cap = 0;
    while (getline(&line, &cap, stdin) != -1)
    {
        if (line[0] == '#' || line[0] == '\n') continue;
        istringstream ss(line);
        vector<int> v; int x; while (ss >> x) v.push_back(x);
        if ((int)v.size() != E) continue;
        int hi = 0; for (int c : v) if (c > hi) hi = c;
        vector<int> m(hi + 1, -1); int nx = 0;
        string c((size_t)E, 0);
        for (int e = 0; e < E; e++) { if (m[v[e]] < 0) m[v[e]] = nx++; c[e] = (char)m[v[e]]; }
        raw.push_back(std::move(c));
    }
    free(line);

    // ---- IR-canonicalize + |Aut| (parallel) ----
    long long M = (long long)raw.size();
    vector<string> canon((size_t)M);
    vector<long long> ord((size_t)M);
    #pragma omp parallel for schedule(dynamic, 64)
    for (long long i = 0; i < M; i++)
    {
        canir::CanonIR ir; ir.init(raw[(size_t)i], N); ir.run();
        canon[(size_t)i] = ir.bestNorm;
        ord[(size_t)i] = ir.autOrder();
    }

    // ---- dedup by canonical form, keep its order ----
    unordered_map<string, long long> seen;   // canon -> order
    for (long long i = 0; i < M; i++) seen.emplace(canon[(size_t)i], ord[(size_t)i]);

    // ---- group by order ----
    map<long long, vector<const string*>> byOrder;
    for (auto& kv : seen) byOrder[kv.second].push_back(&kv.first);
    for (auto& kv : byOrder) sort(kv.second.begin(), kv.second.end(),
                                  [](const string* a, const string* b){ return *a < *b; });

    size_t total = seen.size();
    printf("N=%d: closed subgroup classes = %zu (%s)\n", N, total, note);
    printf("  orders:");
    for (auto& kv : byOrder) printf(" %lld(x%zu)", kv.first, kv.second.size());
    printf("\n");
    for (auto& kv : byOrder)
    {
        printf("--- order %lld (x%zu) ---\n", kv.first, kv.second.size());
        for (const string* c : kv.second)
        {
            for (int e = 0; e < E; e++) printf("%s%d", e ? " " : "", (int)(*c)[(size_t)e]);
            printf("\n");
        }
    }
    fprintf(stderr, "catalog_format: %lld inputs -> %zu distinct classes\n", M, total);
    return 0;
}
