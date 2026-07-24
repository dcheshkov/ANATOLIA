// classify_catalog.cpp — taxonomy classifier for closed-group catalogues.
//
// Reads a catalogue in the 8Spins.txt / VERBOSE format of closed_groups
// (sections "--- order K ... ---" followed by lines of E = N(N-1)/2 edge
// labels in column order) and emits one TSV row per entry:
//
//   id            N/o<order>.<index>              (the canonical identifier)
//   order         computed |Aut| (verified against the section header)
//   twins         magnetic composite sizes, comma separated (descending)
//   gq            |Gamma_q| = |Aut| / prod(n_g!)
//   orbits        vertex orbit sizes under Aut, plus separated (descending)
//   fixedpt       y = has a fixed vertex (inherited from level N-1 by the
//                 tower theorem), n = fixed-point-free (essentially new)
//   cyclic        y/n (an element of order |Aut| exists) - always computed
//   maxelt        maximal element order
//   abelian       y/n/-   (computed for |Aut| <= 2000, else '-')
//   ambivalent    y/n/-   every element conjugate to its inverse
//                 (real character table <=> y); |Aut| <= 2000, else '-'
//
// Automorphisms are enumerated by the same prefix-pruned ordering search as
// in closed_groups.cpp, so huge groups (up to S_N) are handled: the vertex
// orbits, fixed points and element orders are accumulated on the fly, and
// the full element list is stored only when |Aut| <= 2000.
//
// Build:  g++ -O3 -march=native -o classify_catalog classify_catalog.cpp
// Usage:  ./classify_catalog N catalogue.txt >> taxonomy.tsv

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace std;

static int N, E;
static int eidxA[16][16];

static void buildEdges(void)
{
    E = 0;
    for (int j = 1; j < N; j++)
        for (int i = 0; i < j; i++)
        { eidxA[i][j] = eidxA[j][i] = E; E++; }
}

static long long gcdll(long long a, long long b) { while (b) { long long t = a % b; a = b; b = t; } return a; }

static long long permOrder(const int* p)
{
    bool seen[16] = { false };
    long long o = 1;
    for (int v = 0; v < N; v++)
    {
        if (seen[v]) continue;
        int len = 0, u = v;
        while (!seen[u]) { seen[u] = true; len++; u = p[u]; }
        o = o / gcdll(o, len) * len;
    }
    return o;
}

struct Entry
{
    long long order = 0;
    vector<int> twins, orbits;
    long long gq = 0;
    bool fixedpt = false, cyclic = false;
    long long maxelt = 1;
    int abelian = -1, ambivalent = -1;   // -1 unknown, 0 no, 1 yes
};

static const long long STORE_LIMIT = 2000;

// Landau function g(N): maximal element order in S_N (max lcm over partitions)
static long long landauRec(int rem, int minPart, long long lcmSoFar)
{
    long long best = lcmSoFar;
    for (int p = minPart; p <= rem; p++)
    {
        long long lcmNext = lcmSoFar / gcdll(lcmSoFar, p) * p;
        long long v = landauRec(rem - p, p + 1, lcmNext);
        if (v > best) best = v;
    }
    return best;
}

static Entry classify(const string& c)
{
    Entry R;

    // single-colour pattern: Aut = S_N; the pruned search would walk all N!
    // successes (6.2e9 at N=13), while every column is known in closed form
    {
        bool oneColour = true;
        for (int e = 1; e < E && oneColour; e++) if (c[(size_t)e] != c[0]) oneColour = false;
        if (oneColour && N >= 3)
        {
            R.order = 1; for (int k = 2; k <= N; k++) R.order *= k;
            R.twins.push_back(N);
            R.orbits.push_back(N);
            R.gq = 1;
            R.fixedpt = false;
            R.maxelt = landauRec(N, 1, 1);
            R.cyclic = false;
            R.abelian = R.ambivalent = (R.order <= STORE_LIMIT) ? 0 : -1;
            if (R.order <= STORE_LIMIT) R.ambivalent = 1;   // S_N is ambivalent
            return R;
        }
    }

    // twins: identical coupling rows off the mutual edge
    {
        int parent[16];
        for (int v = 0; v < N; v++) parent[v] = v;
        auto find = [&](int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
        for (int i = 0; i < N; i++)
            for (int j = i + 1; j < N; j++)
            {
                bool same = true;
                for (int k = 0; k < N && same; k++)
                {
                    if (k == i || k == j) continue;
                    if (c[(size_t)eidxA[i][k]] != c[(size_t)eidxA[j][k]]) same = false;
                }
                if (same) { int a = find(i), b = find(j); if (a != b) parent[b] = a; }
            }
        int cnt[16] = { 0 };
        for (int v = 0; v < N; v++) cnt[find(v)]++;
        for (int v = 0; v < N; v++) if (cnt[v]) R.twins.push_back(cnt[v]);
        sort(R.twins.rbegin(), R.twins.rend());
    }

    // automorphism enumeration (prefix-pruned ordering search)
    int vp[16]; for (int i = 0; i < N; i++) vp[i] = i;
    int vparent[16]; for (int i = 0; i < N; i++) vparent[i] = i;
    auto vfind = [&](int x) { while (vparent[x] != x) { vparent[x] = vparent[vparent[x]]; x = vparent[x]; } return x; };
    bool moved[16] = { false };
    vector<vector<int>> elems;

    for (;;)
    {
        int failCol = -1;
        for (int j = 1; j < N && failCol < 0; j++)
            for (int i = 0; i < j; i++)
                if (c[(size_t)eidxA[vp[i]][vp[j]]] != c[(size_t)eidxA[i][j]]) { failCol = j; break; }
        if (failCol < 0)
        {
            R.order++;
            for (int v = 0; v < N; v++)
            {
                if (vp[v] != v) moved[v] = true;
                int a = vfind(v), b = vfind(vp[v]);
                if (a != b) vparent[b] = a;
            }
            long long o = permOrder(vp);
            if (o > R.maxelt) R.maxelt = o;
            if (R.order <= STORE_LIMIT) elems.push_back(vector<int>(vp, vp + N));
        }
        else sort(vp + failCol + 1, vp + N, greater<int>());
        if (!next_permutation(vp, vp + N)) break;
    }

    {
        int cnt[16] = { 0 };
        for (int v = 0; v < N; v++) cnt[vfind(v)]++;
        for (int v = 0; v < N; v++) if (cnt[v]) R.orbits.push_back(cnt[v]);
        sort(R.orbits.rbegin(), R.orbits.rend());
    }
    R.fixedpt = false;
    for (int v = 0; v < N; v++) if (!moved[v]) R.fixedpt = true;
    R.cyclic = (R.maxelt == R.order);
    long long pf = 1;
    for (int t : R.twins) for (int k = 2; k <= t; k++) pf *= k;
    R.gq = (R.order % pf == 0) ? R.order / pf : -1;

    if (R.order <= STORE_LIMIT)
    {
        int ne = (int)elems.size();
        auto comp = [&](const int* a, const int* b, int* out)
        { for (int i = 0; i < N; i++) out[i] = a[b[i]]; };
        // abelian
        R.abelian = 1;
        int ab[16], ba[16];
        for (int x = 0; x < ne && R.abelian; x++)
            for (int y = x + 1; y < ne && R.abelian; y++)
            {
                comp(elems[x].data(), elems[y].data(), ab);
                comp(elems[y].data(), elems[x].data(), ba);
                if (memcmp(ab, ba, N * sizeof(int)) != 0) R.abelian = 0;
            }
        // ambivalence: every g conjugate to g^{-1}
        R.ambivalent = 1;
        int gi[16], t1[16], t2[16];
        for (int x = 0; x < ne && R.ambivalent; x++)
        {
            for (int i = 0; i < N; i++) gi[elems[x][i]] = i;   // inverse of g
            bool found = false;
            for (int h = 0; h < ne && !found; h++)
            {
                // h g h^{-1} == g^{-1} ?
                int hi[16];
                for (int i = 0; i < N; i++) hi[elems[h][i]] = i;
                comp(elems[x].data(), hi, t1);   // g h^{-1}
                comp(elems[h].data(), t1, t2);   // h g h^{-1}
                if (memcmp(t2, gi, N * sizeof(int)) == 0) found = true;
            }
            if (!found) R.ambivalent = 0;
        }
    }
    return R;
}

int main(int argc, char** argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s N catalogue.txt\n", argv[0]); return 1; }
    N = atoi(argv[1]);
    buildEdges();
    FILE* f = fopen(argv[2], "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", argv[2]); return 1; }

    char buf[4096];
    long long secOrder = -1;
    int idx = 0;
    printf("id\torder\ttwins\tgq\torbits\tfixedpt\tcyclic\tmaxelt\tabelian\tambivalent\n");
    while (fgets(buf, sizeof(buf), f))
    {
        if (strncmp(buf, "---", 3) == 0)
        {
            const char* p = strstr(buf, "order ");
            secOrder = p ? atoll(p + 6) : -1;
            idx = 0;
            continue;
        }
        // pattern line: digits (possibly with leading '*')
        string line(buf);
        string cleaned;
        for (char ch : line) if (ch != '*') cleaned += ch;
        vector<int> lab;
        {
            const char* p = cleaned.c_str();
            char* end;
            for (;;)
            {
                long v = strtol(p, &end, 10);
                if (end == p) break;
                lab.push_back((int)v);
                p = end;
            }
        }
        if ((int)lab.size() != E) continue;
        idx++;
        string c((size_t)E, (char)0);
        { int map_[128]; for (int i = 0; i < E; i++) map_[i] = -1; int nx = 0;
          for (int e = 0; e < E; e++) { if (map_[lab[e]] < 0) map_[lab[e]] = nx++; c[(size_t)e] = (char)map_[lab[e]]; } }
        Entry R = classify(c);
        if (secOrder > 0 && R.order != secOrder)
            fprintf(stderr, "WARNING: N=%d order-%lld #%d computed |Aut|=%lld\n", N, secOrder, idx, R.order);
        string tw, ob;
        for (size_t i = 0; i < R.twins.size(); i++) { if (i) tw += ","; tw += to_string(R.twins[i]); }
        for (size_t i = 0; i < R.orbits.size(); i++) { if (i) ob += "+"; ob += to_string(R.orbits[i]); }
        printf("%d/o%lld.%d\t%lld\t%s\t%lld\t%s\t%c\t%c\t%lld\t%s\t%s\n",
               N, R.order, idx, R.order, tw.c_str(), R.gq, ob.c_str(),
               R.fixedpt ? 'y' : 'n', R.cyclic ? 'y' : 'n', R.maxelt,
               R.abelian < 0 ? "-" : (R.abelian ? "y" : "n"),
               R.ambivalent < 0 ? "-" : (R.ambivalent ? "y" : "n"));
    }
    fclose(f);
    return 0;
}
