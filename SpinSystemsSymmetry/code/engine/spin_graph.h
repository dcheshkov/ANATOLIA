// spin_graph.h - shared primitives for the closed-group enumerators.
//
// K_N edge indexing in column order (J12 J13 J23 J14 ...), restricted-growth
// strings, and the coarsening step (union edge-colour classes along a
// permutation's orbits so that the permutation becomes an automorphism).
// Included by both closed_groups.cpp (single node) and closed_groups_mpi.cpp
// (distributed); the canonical form and |Aut| live in canon_ir.h.
//
// The graph state (N, E, eidxA, edgeV) is file-local (static) in each
// translation unit and initialized by buildEdges(); the two enumerators are
// separate programs, so there is no cross-TU sharing to worry about.

#pragma once
#include <string>
#include <vector>
#include <utility>
#include <chrono>
#include <algorithm>

static int N, E;                                  // number of spins; number of edges
static int eidxA[16][16];                         // eidxA[i][j] = edge index of {i,j}
static std::vector<std::pair<int, int>> edgeV;    // edges in column order

static double now(void)
{
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

static void buildEdges(void)
{
    E = 0;
    edgeV.clear();
    for (int j = 1; j < N; j++)
        for (int i = 0; i < j; i++)
        {
            eidxA[i][j] = eidxA[j][i] = E;
            edgeV.push_back({i, j});
            E++;
        }
}

// colouring -> restricted-growth string (labels renumbered by first occurrence)
static inline std::string rgsOf(const int* lab)
{
    std::string out((size_t)E, (char)0);
    int map_[128]; for (int i = 0; i < E; i++) map_[i] = -1;
    int next = 0;
    for (int i = 0; i < E; i++)
    {
        int v = lab[i];
        if (map_[v] < 0) map_[v] = next++;
        out[(size_t)i] = (char)map_[v];
    }
    return out;
}

// coarsen colouring c so that vertex permutation sig becomes an automorphism:
// union same-colour edges, then union e with sig(e); classes = edge orbits of
// <colours, sig>.
static std::string coarsen(const std::string& c, const int* sig)
{
    int parent[128];
    for (int e = 0; e < E; e++) parent[e] = e;
    auto find = [&](int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };

    int firstOf[128]; for (int i = 0; i < E; i++) firstOf[i] = -1;
    for (int e = 0; e < E; e++)
    {
        int col = (unsigned char)c[(size_t)e];
        if (firstOf[col] < 0) firstOf[col] = e;
        else { int a = find(firstOf[col]), b = find(e); if (a != b) parent[b] = a; }
    }
    for (int e = 0; e < E; e++)
    {
        int a = find(e), b = find(eidxA[sig[edgeV[e].first]][sig[edgeV[e].second]]);
        if (a != b) parent[b] = a;
    }
    int lab[128];
    for (int e = 0; e < E; e++) lab[e] = find(e);
    return rgsOf(lab);
}

// all integer partitions of n with parts <= maxPart (cycle types of S_N)
static void genPartitions(int n, int maxPart, std::vector<int>& cur, std::vector<std::vector<int>>& out)
{
    if (n == 0) { out.push_back(cur); return; }
    for (int p = std::min(n, maxPart); p >= 1; p--)
    {
        cur.push_back(p);
        genPartitions(n - p, p, cur, out);
        cur.pop_back();
    }
}
