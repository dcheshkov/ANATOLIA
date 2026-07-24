// canon_ir_test.cpp - exhaustive correctness + speed check of canon_ir.h
// against the reference canon()/autOrder()/orbit computation used by the
// existing tools, over a catalogue file (VERBOSE format).
//
// For every entry: IR canonical form must equal the reference lex-min RGS,
// IR |Aut| must equal the reference automorphism count, and IR vertex orbits
// must equal the reference orbits (as a partition).  Also times both.
//
// Build:  g++ -O3 -march=native -std=c++17 -o canon_ir_test canon_ir_test.cpp
// Usage:  ./canon_ir_test N catalogue.txt

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <unordered_set>
#include <random>
#include "canon_ir.h"

using namespace std;

static int N, E;
static int eidxA[16][16];

static double now()
{ using namespace std::chrono; return duration<double>(steady_clock::now().time_since_epoch()).count(); }

static void buildEdges()
{ E = 0; for (int j = 1; j < N; j++) for (int i = 0; i < j; i++) { eidxA[i][j] = eidxA[j][i] = E; E++; } }

static string rgsOf(const int* lab)
{
    string out((size_t)E, (char)0);
    int m[512]; for (int i = 0; i < E; i++) m[i] = -1; int nx = 0;
    for (int i = 0; i < E; i++) { int v = lab[i]; if (m[v] < 0) m[v] = nx++; out[(size_t)i] = (char)m[v]; }
    return out;
}

// --- reference canon() (old prefix-pruned lex-min), verbatim in spirit ---
static string refCanon(const string& c)
{
    int v[16]; for (int i = 0; i < N; i++) v[i] = i;
    string best = c, cur((size_t)E, (char)0); int map_[512];
    for (;;)
    {
        if (!next_permutation(v, v + N)) break;
        int nextLab = 0; for (int i = 0; i < E; i++) map_[i] = -1;
        int pos = 0, cmp = 0, failCol = -1;
        for (int j = 1; j < N && failCol < 0; j++)
            for (int i = 0; i < j; i++)
            {
                int col = (unsigned char)c[(size_t)eidxA[v[i]][v[j]]];
                if (map_[col] < 0) map_[col] = nextLab++;
                char lb = (char)map_[col]; cur[(size_t)pos] = lb;
                if (cmp == 0) { char bb = best[(size_t)pos]; if (lb > bb) { failCol = j; break; } if (lb < bb) cmp = -1; }
                pos++;
            }
        if (failCol >= 0) { sort(v + failCol + 1, v + N, greater<int>()); continue; }
        if (cmp < 0) best = cur;
    }
    return best;
}

// --- reference |Aut| and orbits (old prefix-pruned search), with S_N shortcut ---
static void refAutOrbits(const string& c, long long& order, vector<int>& orbSizes)
{
    bool oneColour = true;
    for (int e = 1; e < E && oneColour; e++) if (c[(size_t)e] != c[0]) oneColour = false;
    if (oneColour && N >= 1)
    {
        order = 1; for (int k = 2; k <= N; k++) order *= k;
        orbSizes.assign(1, N);
        return;
    }
    int vp[16]; for (int i = 0; i < N; i++) vp[i] = i;
    int par[16]; for (int i = 0; i < N; i++) par[i] = i;
    auto find = [&](int x) { while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; } return x; };
    order = 0;
    for (;;)
    {
        int failCol = -1;
        for (int j = 1; j < N && failCol < 0; j++)
            for (int i = 0; i < j; i++)
                if (c[(size_t)eidxA[vp[i]][vp[j]]] != c[(size_t)eidxA[i][j]]) { failCol = j; break; }
        if (failCol < 0) { order++; for (int v = 0; v < N; v++) { int a = find(v), b = find(vp[v]); if (a != b) par[b] = a; } }
        else sort(vp + failCol + 1, vp + N, greater<int>());
        if (!next_permutation(vp, vp + N)) break;
    }
    int cnt[16] = { 0 }; for (int v = 0; v < N; v++) cnt[find(v)]++;
    orbSizes.clear(); for (int v = 0; v < N; v++) if (cnt[v]) orbSizes.push_back(cnt[v]);
    sort(orbSizes.rbegin(), orbSizes.rend());
}

// Schreier-Sims regression guard: |group| from canon_ir vs brute-force
// closure, over random generator sets (the 2^-closed catalogue groups do not
// exercise the nested-stabilizer case that broke the first SS implementation).
static long long bruteOrder(int n, const std::vector<canir::Perm>& gens)
{
    auto key = [&](const canir::Perm& p) { std::string s; for (int i = 0; i < n; i++) s += (char)p[i]; return s; };
    canir::Perm id; for (int i = 0; i < 16; i++) id[i] = i;
    unordered_set<string> seen{ key(id) }; vector<canir::Perm> q{ id };
    for (size_t i = 0; i < q.size(); i++)
        for (const auto& g : gens)
        {
            canir::Perm h; for (int k = 0; k < n; k++) h[k] = g[q[i][k]]; for (int k = n; k < 16; k++) h[k] = k;
            if (seen.insert(key(h)).second) q.push_back(h);
            if (q.size() > 3000000) return -1;
        }
    return (long long)q.size();
}
static int ssSelfTest()
{
    // the reviewer witness: <(1 2)(3 4), (0 3)> on 5 points has order 12
    canir::Perm w1, w2; for (int i = 0; i < 16; i++) { w1[i] = i; w2[i] = i; }
    w1[1] = 2; w1[2] = 1; w1[3] = 4; w1[4] = 3; w2[0] = 3; w2[3] = 0;
    if (canir::groupOrder(5, { w1, w2 }) != 12) { fprintf(stderr, "SS witness FAILED\n"); return 1; }
    std::mt19937 rng(20260713);
    int fails = 0, total = 0;
    for (int n = 3; n <= 8; n++)
        for (int t = 0; t < 8000; t++)
        {
            int ng = 1 + rng() % 3; vector<canir::Perm> gens;
            for (int k = 0; k < ng; k++)
            {
                canir::Perm p; for (int i = 0; i < 16; i++) p[i] = i;
                for (int i = n - 1; i > 0; i--) std::swap(p[i], p[(int)(rng() % (i + 1))]);
                gens.push_back(p);
            }
            long long b = bruteOrder(n, gens); if (b < 0) continue;
            long long s = canir::groupOrder(n, gens); total++;
            if (s != b) { if (++fails <= 3) fprintf(stderr, "SS FAIL n=%d ss=%lld brute=%lld\n", n, s, b); }
        }
    fprintf(stderr, "SS self-test: %d/%d failures\n", fails, total);
    return fails ? 1 : 0;
}

int main(int argc, char** argv)
{
    if (argc == 2 && strcmp(argv[1], "--sstest") == 0) return ssSelfTest();
    if (argc < 3) { fprintf(stderr, "usage: %s N catalogue.txt  |  %s --sstest\n", argv[0], argv[0]); return 1; }
    N = atoi(argv[1]); buildEdges();
    FILE* f = fopen(argv[2], "r"); if (!f) { fprintf(stderr, "cannot open\n"); return 1; }
    char buf[8192];
    long long nEntries = 0, canonFail = 0, autFail = 0, orbFail = 0;
    double tRef = 0, tIR = 0;
    while (fgets(buf, sizeof(buf), f))
    {
        if (strncmp(buf, "---", 3) == 0) continue;
        string line(buf), cleaned; for (char ch : line) if (ch != '*') cleaned += ch;
        vector<int> lab; { const char* p = cleaned.c_str(); char* end; for (;;) { long v = strtol(p, &end, 10); if (end == p) break; lab.push_back((int)v); p = end; } }
        if ((int)lab.size() != E) continue;
        nEntries++;
        string c = rgsOf(lab.data());

        double a0 = now();
        string rc = refCanon(c);
        long long ro; vector<int> rorb; refAutOrbits(c, ro, rorb);
        tRef += now() - a0;

        double b0 = now();
        canir::CanonIR ir; ir.init(c, N); ir.run();
        long long io = ir.autOrder();
        vector<int> iorbLbl = ir.vertexOrbits();
        canir::CanonLex cl; cl.init(c, N, ir.gens);
        string ic = cl.run();               // exact lex-min, auto-accelerated
        tIR += now() - b0;

        int cnt[16] = { 0 }; for (int v = 0; v < N; v++) cnt[iorbLbl[v]]++;
        vector<int> iorb; for (int v = 0; v < N; v++) if (cnt[v]) iorb.push_back(cnt[v]);
        sort(iorb.rbegin(), iorb.rend());

        if (ic != rc) { canonFail++; if (canonFail <= 3) fprintf(stderr, "CANON MISMATCH entry %lld\n", nEntries); }
        if (io != ro) { autFail++; if (autFail <= 5) fprintf(stderr, "AUT MISMATCH entry %lld: ref=%lld ir=%lld\n", nEntries, ro, io); }
        if (iorb != rorb) { orbFail++; if (orbFail <= 5) fprintf(stderr, "ORBIT MISMATCH entry %lld\n", nEntries); }
    }
    fclose(f);
    printf("N=%d entries=%lld  canonFail=%lld autFail=%lld orbFail=%lld  refTime=%.2fs irTime=%.2fs speedup=%.2fx\n",
           N, nEntries, canonFail, autFail, orbFail, tRef, tIR, tRef / (tIR > 0 ? tIR : 1e-9));
    return (canonFail || autFail || orbFail) ? 2 : 0;
}
