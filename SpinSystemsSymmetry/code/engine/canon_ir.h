// canon_ir.h - nauty-style individualization-refinement canonicalization of
// edge-coloured complete graphs K_N, specialized for the spin-symmetry code.
//
// Given a symmetric edge colouring col[i][j] (0..K-1) of K_N, computes in one
// pass:
//   * the canonical form: the lexicographically minimal restricted-growth
//     string (RGS) of the colouring over all vertex orderings, read column by
//     column (c(v_i,v_j), j=1..N-1, i<j) - identical definition to the old
//     canon() in closed_groups.cpp, so canonical forms are interchangeable;
//   * generators of the automorphism group Aut (colour-preserving vertex
//     permutations);
//   * the vertex orbits (connected components under the generators);
//   * |Aut|, computed by Schreier-Sims from the generators (no element
//     enumeration - handles Aut up to S_N).
//
// Method (McKay-Piperno individualization-refinement):
//   1. colour refinement (1-WL with edge colours) to an equitable partition;
//   2. individualize one vertex of the first non-singleton cell, refine,
//      recurse - a DFS whose leaves are discrete partitions = vertex orders;
//   3. at each leaf, the RAW colour sequence is the leaf certificate; a leaf
//      whose raw certificate equals the first leaf's yields an automorphism
//      (their orders differ by a colour-preserving permutation); the canonical
//      form is the lex-min of the RGS-NORMALIZED certificate over all leaves;
//   4. automorphism pruning: within a target cell only one vertex per orbit
//      (under the automorphisms found so far that fix the individualization
//      path pointwise) is branched on - so highly symmetric graphs, the worst
//      case for the old prefix search, are the cheapest here.
//
// The two certificates are deliberately different: normalized (RGS) for the
// canonical form / dedup, RAW for automorphism detection.  Using the RGS for
// automorphisms would admit colour-permuting maps that are not automorphisms
// of the fixed colouring.
//
// Correctness is checked exhaustively against the old canon()/autOrder() over
// the whole catalogue N=3..13 in canon_ir_test.cpp.

#pragma once
#include <array>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <cstring>

namespace canir {

static const int MAXN = 16;

// ---------------- permutation helpers ----------------

typedef std::array<int, MAXN> Perm;   // p[i] = image of i (points 0..N-1)

static inline bool isIdentity(const Perm& p, int N)
{
    for (int i = 0; i < N; i++) if (p[i] != i) return false;
    return true;
}

// ---------------- Schreier-Sims group order (no element enumeration) --------
//
// Deterministic BSGS builder (Sims): returns the group order as the product
// of the basic orbit lengths.  The control flow follows Sims exactly - after a
// Schreier generator fails to sift, its residue is added as a strong generator
// and the outer level index jumps to the drop level (i := dropLevel), so all
// affected levels are reprocessed to a fixpoint.  (An earlier hand-rolled
// single-pass re-sift undercounted groups whose base needs extension via a
// deep Schreier generator, e.g. <(1 2)(3 4),(0 3)> on 5 points -> 6 not 12;
// this version is verified against brute-force closure on 120000 random
// groups over N=3..8 with zero mismatches.)

static Perm ssCompose(const Perm& a, const Perm& b, int N)   // (a after b): i -> a[b[i]]
{ Perm r; for (int i = 0; i < N; i++) r[i] = a[b[i]]; return r; }
static Perm ssInverse(const Perm& a, int N)
{ Perm r; for (int i = 0; i < N; i++) r[a[i]] = i; return r; }

struct SchreierSims
{
    int N;
    std::vector<int> base;
    std::vector<std::vector<Perm>> S;        // S[i]: strong gens fixing base[0..i-1]
    std::vector<std::vector<int>> orbit;     // orbit[i]: points, orbit[i][0] = base[i]
    std::vector<std::vector<int>> pos;        // per level: point -> index in orbit (or -1)
    std::vector<std::vector<Perm>> u, uinv;   // transversal reps u[i][k](base[i]) = orbit[i][k]

    void newLevel(int b)
    {
        base.push_back(b); S.push_back({}); orbit.push_back({});
        pos.push_back(std::vector<int>(N, -1)); u.push_back({}); uinv.push_back({});
    }

    void computeOrbit(int i)
    {
        int b = base[i];
        orbit[i].assign(1, b);
        std::fill(pos[i].begin(), pos[i].end(), -1);
        pos[i][b] = 0;
        Perm id; for (int k = 0; k < N; k++) id[k] = k;
        u[i].assign(1, id); uinv[i].assign(1, id);
        for (size_t oi = 0; oi < orbit[i].size(); oi++)
        {
            int p = orbit[i][oi];
            for (const Perm& s : S[i])
            {
                int q = s[p];
                if (pos[i][q] < 0)
                {
                    pos[i][q] = (int)orbit[i].size();
                    orbit[i].push_back(q);
                    Perm rep = ssCompose(s, u[i][oi], N);
                    u[i].push_back(rep); uinv[i].push_back(ssInverse(rep, N));
                }
            }
        }
    }

    bool sift(Perm g, Perm& res, int& dropLevel)
    {
        for (size_t i = 0; i < base.size(); i++)
        {
            int img = g[base[i]];
            if (pos[i][img] < 0) { res = g; dropLevel = (int)i; return false; }
            g = ssCompose(uinv[i][pos[i][img]], g, N);
        }
        res = g; dropLevel = (int)base.size();
        return isIdentity(g, N);
    }

    int firstMoved(const Perm& p) { for (int i = 0; i < N; i++) if (p[i] != i) return i; return -1; }

    long long build(int n, const std::vector<Perm>& generators)
    {
        N = n;
        base.clear(); S.clear(); orbit.clear(); pos.clear(); u.clear(); uinv.clear();
        std::vector<Perm> gg;
        for (const Perm& g : generators) if (!isIdentity(g, N)) gg.push_back(g);
        for (const Perm& g : gg)
        {
            bool movesBase = false;
            for (int b : base) if (g[b] != b) { movesBase = true; break; }
            if (!movesBase) newLevel(firstMoved(g));
        }
        for (const Perm& g : gg)
            for (size_t i = 0; i < base.size(); i++)
            {
                bool fixes = true;
                for (size_t t = 0; t < i; t++) if (g[base[t]] != base[t]) { fixes = false; break; }
                if (fixes) S[i].push_back(g);
            }
        for (size_t i = 0; i < base.size(); i++) computeOrbit((int)i);

        int i = (int)base.size() - 1;
        while (i >= 0)
        {
            bool restart = false;
            std::vector<Perm> Ssnap = S[i];       // snapshot: orbit/pos only grow, safe to read live
            std::vector<int> Osnap = orbit[i];
            for (int beta : Osnap)
            {
                Perm ub = u[i][pos[i][beta]];
                for (const Perm& s : Ssnap)
                {
                    int gamma = s[beta];
                    Perm sg = ssCompose(uinv[i][pos[i][gamma]], ssCompose(s, ub, N), N);
                    Perm res; int dl;
                    if (!sift(sg, res, dl))
                    {
                        if (dl == (int)base.size()) newLevel(firstMoved(res));
                        for (int j = i + 1; j <= dl; j++) S[j].push_back(res);
                        for (int j = i + 1; j <= dl; j++) computeOrbit(j);
                        i = dl; restart = true; break;      // reprocess from the drop level
                    }
                }
                if (restart) break;
            }
            if (!restart) i--;
        }
        long long o = 1;
        for (size_t k = 0; k < base.size(); k++) o *= (long long)orbit[k].size();
        return o;
    }
};

static long long groupOrder(int N, const std::vector<Perm>& generators)
{
    SchreierSims ss; return ss.build(N, generators);
}

// ---------------- the canonicalizer ----------------

struct CanonIR
{
    int N, E, K;
    int col[MAXN][MAXN];                 // symmetric edge colours; diag unused
    std::vector<std::pair<int,int>> edgeV;   // column order
    int eidx[MAXN][MAXN];

    // search state
    std::string bestNorm;                // lex-min normalized certificate
    std::vector<int> firstOrd;           // first leaf vertex order
    std::string firstRaw;                // its raw certificate
    bool haveFirst = false;
    std::vector<Perm> gens;              // automorphism generators

    void init(const std::string& rgs, int n)
    {
        N = n; E = 0; edgeV.clear();
        for (int j = 1; j < N; j++)
            for (int i = 0; i < j; i++) { eidx[i][j] = eidx[j][i] = E; edgeV.push_back({i, j}); E++; }
        K = 0;
        for (int e = 0; e < E; e++)
        {
            int c = (unsigned char)rgs[(size_t)e];
            auto [i, j] = edgeV[e];
            col[i][j] = col[j][i] = c;
            if (c + 1 > K) K = c + 1;
        }
        bestNorm.clear(); haveFirst = false; gens.clear();
    }

    // raw certificate: colours read in column order under vertex order ord
    std::string rawCert(const std::vector<int>& ord)
    {
        std::string s((size_t)E, (char)0);
        int p = 0;
        for (int j = 1; j < N; j++)
            for (int i = 0; i < j; i++)
                s[(size_t)p++] = (char)col[ord[i]][ord[j]];
        return s;
    }

    static std::string normalize(const std::string& raw)
    {
        std::string out(raw.size(), (char)0);
        int map_[256]; std::memset(map_, -1, sizeof(map_));
        int nx = 0;
        for (size_t i = 0; i < raw.size(); i++)
        {
            int v = (unsigned char)raw[i];
            if (map_[v] < 0) map_[v] = nx++;
            out[i] = (char)map_[v];
        }
        return out;
    }

    // colour refinement to an equitable ordered partition
    void refine(std::vector<std::vector<int>>& P)
    {
        bool changed = true;
        while (changed)
        {
            changed = false;
            std::vector<std::vector<int>> NP;
            NP.reserve(P.size());
            for (auto& C : P)
            {
                if (C.size() <= 1) { NP.push_back(C); continue; }
                std::vector<std::pair<std::vector<int>, int>> sv;
                sv.reserve(C.size());
                for (int v : C)
                {
                    std::vector<int> sig;
                    sig.reserve(P.size() * K);
                    for (auto& D : P)
                    {
                        int cnt[128] = { 0 };
                        for (int u : D) if (u != v) cnt[col[v][u]]++;
                        for (int k = 0; k < K; k++) sig.push_back(cnt[k]);
                    }
                    sv.push_back({std::move(sig), v});
                }
                std::sort(sv.begin(), sv.end());
                size_t i = 0;
                bool split = false;
                while (i < sv.size())
                {
                    size_t j = i + 1;
                    while (j < sv.size() && sv[j].first == sv[i].first) j++;
                    std::vector<int> cell;
                    for (size_t t = i; t < j; t++) cell.push_back(sv[t].second);
                    if (cell.size() < C.size()) split = true;
                    NP.push_back(std::move(cell));
                    i = j;
                }
                if (split) changed = true;
            }
            P.swap(NP);
        }
    }

    void addAuto(const std::vector<int>& baseOrd, const std::vector<int>& ord)
    {
        Perm g;
        for (int i = 0; i < N; i++) g[baseOrd[i]] = ord[i];
        if (isIdentity(g, N)) return;
        for (const Perm& h : gens) { bool eq = true; for (int i = 0; i < N; i++) if (h[i] != g[i]) { eq = false; break; } if (eq) return; }
        gens.push_back(g);
    }

    void leaf(const std::vector<int>& ord)
    {
        std::string raw = rawCert(ord);
        std::string norm = normalize(raw);
        if (!haveFirst) { haveFirst = true; firstOrd = ord; firstRaw = raw; bestNorm = norm; return; }
        if (raw == firstRaw) addAuto(firstOrd, ord);
        if (norm < bestNorm) bestNorm = norm;
    }

    // orbits of target cell W under automorphisms fixing the path pointwise
    bool prunedByAuto(const std::vector<int>& W, int w,
                      const std::vector<int>& path)
    {
        int parent[MAXN]; for (int x : W) parent[x] = x;
        auto find = [&](int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
        for (const Perm& g : gens)
        {
            bool fixes = true;
            for (int pv : path) if (g[pv] != pv) { fixes = false; break; }
            if (!fixes) continue;
            for (int x : W) { int gx = g[x]; int a = find(x), b = find(gx); if (a != b) parent[b] = a; }
        }
        int rw = find(w);
        for (int x : W) if (x < w && find(x) == rw) return true;   // smaller rep exists
        return false;
    }

    void search(std::vector<std::vector<int>> P, std::vector<int> path)
    {
        refine(P);
        int t = -1;
        for (size_t i = 0; i < P.size(); i++) if (P[i].size() > 1) { t = (int)i; break; }
        if (t < 0)
        {
            std::vector<int> ord; ord.reserve(N);
            for (auto& C : P) ord.push_back(C[0]);
            leaf(ord);
            return;
        }
        std::vector<int> W = P[t];
        std::sort(W.begin(), W.end());
        for (int w : W)
        {
            if (prunedByAuto(W, w, path)) continue;
            std::vector<std::vector<int>> Q;
            Q.reserve(P.size() + 1);
            for (int i = 0; i < (int)P.size(); i++)
            {
                if (i != t) { Q.push_back(P[i]); continue; }
                std::vector<int> rest;
                for (int u : P[i]) if (u != w) rest.push_back(u);
                Q.push_back(std::vector<int>{w});
                Q.push_back(rest);
            }
            std::vector<int> np = path; np.push_back(w);
            search(std::move(Q), std::move(np));
        }
    }

    // run: fills bestNorm, gens
    void run()
    {
        std::vector<std::vector<int>> P(1);
        for (int i = 0; i < N; i++) P[0].push_back(i);
        std::vector<int> path;
        search(P, path);
    }

    // vertex orbits (component labels) under the generators
    std::vector<int> vertexOrbits()
    {
        int parent[MAXN]; for (int i = 0; i < N; i++) parent[i] = i;
        std::function<int(int)> find = [&](int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
        for (const Perm& g : gens)
            for (int i = 0; i < N; i++) { int a = find(i), b = find(g[i]); if (a != b) parent[b] = a; }
        std::vector<int> lbl(N, -1); int nx = 0;
        for (int i = 0; i < N; i++) { int r = find(i); if (lbl[r] < 0) lbl[r] = nx++; }
        std::vector<int> orb(N);
        for (int i = 0; i < N; i++) orb[i] = lbl[find(i)];
        return orb;
    }

    long long autOrder() { return groupOrder(N, gens); }
};

// ---------------- lex-min canonical string, automorphism-pruned ----------------
//
// Produces the SAME canonical form as the old canon() (lexicographically
// minimal RGS over all vertex orderings), but accelerated: the automorphism
// generators (computed by CanonIR) let us skip, at each ordering position,
// all but one vertex per orbit of the prefix-stabilizer subgroup - the exact
// symmetry that made the old next_permutation walk slow on symmetric graphs.
// Soundness: skipped candidates are images under a genuine automorphism that
// fixes the placed prefix, so their certificate continuations are identical
// and cannot beat the representative already tried.

struct CanonLex
{
    int N, E;
    int col[MAXN][MAXN];
    std::vector<Perm> gens;
    std::string best;                    // length E
    int ord[MAXN];
    bool used[MAXN];

    void init(const std::string& rgs, int n, const std::vector<Perm>& g)
    {
        N = n; E = 0;
        std::vector<std::pair<int,int>> ev;
        for (int j = 1; j < N; j++) for (int i = 0; i < j; i++) { ev.push_back({i, j}); E++; }
        for (int e = 0; e < E; e++) { auto [i, j] = ev[e]; int c = (unsigned char)rgs[(size_t)e]; col[i][j] = col[j][i] = c; }
        gens = g;
        best.assign((size_t)E, (char)127);   // +inf sentinel
        for (int i = 0; i < N; i++) used[i] = false;
    }

    // place vertex at position j; column j = (c(ord[i],ord[j]) : i<j).  The
    // comparison against best is recomputed over the whole written prefix each
    // time, NOT carried down: best changes as leaves are found, so an inherited
    // cmp would go stale (a prefix judged "below best" against an earlier best
    // must be re-judged against the current one, else a larger leaf overwrites
    // a smaller best).
    void rec(int j, int* mp, int nlab, char* cur)
    {
        int base = j * (j - 1) / 2;       // certificate offset of column j
        if (j == N)
        {
            for (int p = 0; p < E; p++)
            {
                unsigned char a = (unsigned char)cur[p], b = (unsigned char)best[(size_t)p];
                if (a < b) { best.assign(cur, (size_t)E); return; }
                if (a > b) return;
            }
            return;                       // equal to best: keep best
        }
        // orbit pruning: skip w if a smaller unused vertex lies in the same
        // orbit under generators fixing ord[0..j-1] pointwise
        int par[MAXN]; for (int i = 0; i < N; i++) par[i] = i;
        std::function<int(int)> find = [&](int x) { while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; } return x; };
        for (const Perm& g : gens)
        {
            bool fixes = true;
            for (int t = 0; t < j; t++) if (g[ord[t]] != ord[t]) { fixes = false; break; }
            if (!fixes) continue;
            for (int x = 0; x < N; x++) if (!used[x]) { int a = find(x), b = find(g[x]); if (a != b) par[b] = a; }
        }
        bool repDone[MAXN] = { false };

        for (int w = 0; w < N; w++)
        {
            if (used[w]) continue;
            int r = find(w);
            if (repDone[r]) continue;     // another vertex of this orbit already tried
            repDone[r] = true;

            int mp2[512]; std::memcpy(mp2, mp, sizeof(int) * 512);
            int nl2 = nlab;
            for (int i = 0; i < j; i++)
            {
                int c = col[ord[i]][w];
                if (mp2[c] < 0) mp2[c] = nl2++;
                cur[base + i] = (char)mp2[c];
            }
            // compare the written prefix [0, base+j) against best; prune if
            // strictly greater (best is a complete string once the first leaf
            // is in, so a larger prefix cannot yield a smaller full string)
            bool worse = false;
            for (int p = 0; p < base + j; p++)
            {
                unsigned char a = (unsigned char)cur[p], b = (unsigned char)best[(size_t)p];
                if (a < b) break;
                if (a > b) { worse = true; break; }
            }
            if (worse) continue;
            ord[j] = w; used[w] = true;
            rec(j + 1, mp2, nl2, cur);
            used[w] = false;
        }
    }

    std::string run()
    {
        char cur[256];
        int mp[512]; std::memset(mp, -1, sizeof(mp));
        rec(0, mp, 0, cur);
        return best;
    }
};

// convenience: canonical RGS string of an RGS colouring.  IR-mode (default)
// returns the IR canonical form (fast, isomorphism-invariant, but a different
// representative than the old canon).  lexMin=true returns the exact old
// lex-min RGS, automorphism-accelerated - use this to keep catalogue
// representatives stable.
//
// Worst case (lexMin): a near-rainbow colouring (many distinct colours) with
// trivial automorphism group makes both prunings inert and the lex-min search
// degenerates to the full N! ordering walk - the SAME worst case as the old
// next_permutation canon (not a regression).  It does not arise in the
// enumerator: the discrete colouring is never canonized (it seeds the BFS
// directly), and BFS colourings are coarsenings with few colours and ample
// symmetry.  Callers that might feed large asymmetric colourings should pass
// lexMin=false (IR form) instead.
static inline std::string canonize(const std::string& rgs, int N, bool lexMin = true)
{
    CanonIR ir; ir.init(rgs, N); ir.run();
    if (!lexMin) return ir.bestNorm;
    CanonLex cl; cl.init(rgs, N, ir.gens); return cl.run();
}

} // namespace canir
