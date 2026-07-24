// filter_closed.cpp - filter a list of subgroup edge-colourings by 2^-closure.
//
// Reads lines "N |G| c_1 c_2 ... c_E" (E = N(N-1)/2, column order) produced by
// a subgroup enumerator (e.g. GAP): c_k is the pair-orbit index of edge k under
// the group G.  A group is 2^-closed (realizable as a spin-system symmetry) iff
// |Aut| of its edge orbit-colouring equals |G| (for N>=3 the vertex orbits are
// induced by the edges, so the edge colouring alone decides it).  Deduplicates
// the closed ones by canonical form and prints the count per N.
//
// Build:  g++ -O3 -march=native -fopenmp -std=c++17 -o filter_closed filter_closed.cpp
// Usage:  ./filter_closed colourings.txt          (reads stdin if - or omitted)
//
// The 2^-closure test (canon_ir per group) is OpenMP-parallel over the input
// lines; a single line's canonicalization is not the bottleneck, the count of
// groups is (75154 at N=14), so the fan-out is over lines.

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <sstream>
#include "canon_ir.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

int main(int argc, char** argv)
{
    const char* path = (argc >= 2) ? argv[1] : "-";
    FILE* f = (path[0] == '-' && path[1] == 0) ? stdin : fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return 1; }

    // ---- read + RGS-normalize all colourings (single pass, cheap) ----
    vector<int> ns;            // degree N of each entry
    vector<long long> gs;      // |G| of each entry
    vector<string> cols;       // RGS-normalized colouring, length E(N)
    char* line = nullptr; size_t cap = 0;
    while (getline(&line, &cap, f) != -1)
    {
        if (line[0] == '#' || line[0] == '\n') continue;
        istringstream ss(line);
        int n; long long g; if (!(ss >> n >> g)) continue;
        int E = n * (n - 1) / 2;
        vector<int> lab; int x; while (ss >> x) lab.push_back(x);
        if ((int)lab.size() != E) continue;
        // RGS-normalize (GAP colour labels run 1..#orbits, up to E, so the
        // relabel map must cover the full label range)
        string c((size_t)E, 0);
        int hi = 0; for (int v : lab) if (v > hi) hi = v;
        vector<int> m(hi + 1, -1); int nx = 0;
        for (int e = 0; e < E; e++) { if (m[lab[e]] < 0) m[lab[e]] = nx++; c[e] = (char)m[lab[e]]; }
        ns.push_back(n); gs.push_back(g); cols.push_back(std::move(c));
    }
    if (f != stdin) fclose(f);
    free(line);

    // ---- 2^-closure test per entry (parallel) ----
    long long M = (long long)ns.size();
    vector<char> isClosed((size_t)M, 0);
    vector<string> canon((size_t)M);       // canonical form of the closed ones
    #pragma omp parallel for schedule(dynamic, 1)
    for (long long i = 0; i < M; i++)
    {
        canir::CanonIR ir; ir.init(cols[(size_t)i], ns[(size_t)i]); ir.run();
        if (ir.autOrder() == gs[(size_t)i])
        {
            isClosed[(size_t)i] = 1;
            // dedup key: the TRUE lex-min canonical form (CanonLex, reusing the
            // automorphisms just found for pruning).  ir.bestNorm is only a
            // strong invariant, not a complete canonical form, so it can split
            // S_N-conjugate colourings apart -- harmless when the input has one
            // rep per class (the full-lattice route), but the orbit-partition
            // route emits Y_lambda-classes, several of which may be S_N-conjugate
            // (equal blocks), and those MUST be merged here.  CanonLex runs only
            // on the closed groups, so the cost is small.
            canir::CanonLex cl; cl.init(cols[(size_t)i], ns[(size_t)i], ir.gens);
            canon[(size_t)i] = cl.run();
        }
    }

    // ---- aggregate per degree N ----
    map<int, long long> tested, closed;
    map<int, unordered_set<string>> closedCanon;
    for (long long i = 0; i < M; i++)
    {
        int n = ns[(size_t)i];
        tested[n]++;
        if (isClosed[(size_t)i]) { closed[n]++; closedCanon[n].insert(canon[(size_t)i]); }
    }
    for (auto& kv : tested)
    {
        int n = kv.first;
        printf("n=%d: tested %lld, 2^-closed %lld, distinct closed classes %zu\n",
               n, kv.second, closed[n], closedCanon[n].size());
    }

    // ---- optional: dump the distinct closed canonical colourings ----
    // DUMP=<file> writes one line per distinct closed class: the canonical RGS
    // as space-separated integers (the materials/*Spins.txt convention), grouped
    // by degree N.  This is the realizable-symmetry catalogue for that N (the
    // fixed-point-free part when the input came from emit_fixedfree/partition).
    const char* dump = getenv("DUMP");
    if (dump)
    {
        FILE* o = fopen(dump, "w");
        if (!o) { fprintf(stderr, "cannot write DUMP=%s\n", dump); return 1; }
        for (auto& kv : closedCanon)
        {
            for (const string& c : kv.second)
            {
                for (size_t e = 0; e < c.size(); e++)
                    fprintf(o, "%s%d", e ? " " : "", (int)c[e]);
                fprintf(o, "\n");
            }
        }
        fclose(o);
        fprintf(stderr, "DUMP: wrote distinct closed colourings to %s\n", dump);
    }
    return 0;
}
