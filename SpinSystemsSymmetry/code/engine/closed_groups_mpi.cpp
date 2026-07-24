// closed_groups_mpi.cpp — distributed (MPI + OpenMP) version of
// closed_groups.cpp: enumerates, up to conjugacy, the subgroups of S_N that
// are automorphism groups of edge-weighted complete graphs (2^-closed groups,
// = symmetry groups of N-spin systems with scalar couplings).
//
// Same algorithm as closed_groups.cpp v2 (sequential-symmetrization BFS,
// cycle-type round 1, prefix-pruned canonical form), distributed as a BSP:
//   per round:  1) SWEEP    — the (frontier entry x permutation-chunk) work
//                             units are striped over ranks, OpenMP inside;
//                             raw colourings are deduplicated locally;
//               2) SHUFFLE  — every raw string is routed to its owner rank
//                             (FNV-1a hash mod nranks, MPI_Alltoallv), where
//                             it is deduplicated against the persistent raw
//                             shard of that owner (global dedup across rounds);
//               3) CANON    — owners canonicalize their fresh raws (OpenMP);
//               4) SYNC     — new canonical entries are allgathered, sorted
//                             (determinism) and appended to the replicated
//                             catalogue; they form the next frontier.
// State is checkpointed after every completed round (NOCKPT=1 disables);
// RESUME=1 continues an interrupted run (the rank count must match the
// checkpoint).  The owner-shuffle is exchanged in size-capped batches
// (BATCH_MB megabytes per batch, default 256) so that every MPI_Alltoallv
// count and displacement stays below 2^31: from N=13 on, a round's raw
// volume exceeds 2 GB and a single exchange would overflow the int
// counts/displacements MPI requires (the segfault of the first version).
//
// Build:  mpicxx -O3 -march=native -fopenmp -std=c++17 -o closed_groups_mpi closed_groups_mpi.cpp
// Run (two dual-socket nodes, one rank per socket, 24 threads per rank):
//   mpirun -np 4 -H srv3:2,srv4:2 --map-by ppr:2:node --bind-to socket \
//     -x OMP_NUM_THREADS=24 -x OMP_PLACES=cores -x OMP_PROC_BIND=close \
//     ./closed_groups_mpi 13 13
// Requirements: passwordless ssh between the nodes, the same binary path on
// both, and a shared working directory (NFS) — checkpoint/seen/frontier files
// are read by all ranks.  Self-check first:
//   mpirun -np 2 ./closed_groups_mpi 3 10   ->  3, 8, 11, 27, 36, 90, 131, 282
// Known values (N=3..12): 3, 8, 11, 27, 36, 90, 131, 282, 394, 948.
// VERBOSE=1 prints the catalogue (rank 0) in the 8Spins.txt format;
// TRACE=1 prints new entries of rounds >= 4 with their group orders.

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "spin_graph.h"     // N, E, eidxA, edgeV, now, buildEdges, rgsOf, coarsen, genPartitions
#include "canon_ir.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

static int mpiRank = 0, mpiSize = 1;

static string canon(const string& c)
{
    // nauty-style individualization-refinement (canon_ir.h): byte-identical
    // to the old lex-min RGS but 100-20000x faster on symmetric colourings,
    // which is exactly where the old N!-ordering walk stalled (verified over
    // N=3..13).  The single-colour S_N shortcut is subsumed (refinement makes
    // it a one-leaf search).
    return canir::canonize(c, N, true);
}

static long long autOrder(const string& c)
{
    // |Aut| by Schreier-Sims from the IR generators - no element enumeration,
    // so the old S_N shortcut and the N!-ordering walk are both retired.
    canir::CanonIR ir; ir.init(c, N); ir.run();
    return ir.autOrder();
}

static long long factorialLL(int n)
{
    long long f = 1;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

// k-th lexicographic permutation of 0..N-1 (factorial number system)
static void unrankPerm(long long k, int* v)
{
    bool used[16] = { false };
    for (int pos = 0; pos < N; pos++)
    {
        long long fp = factorialLL(N - 1 - pos);
        int idx = (int)(k / fp); k %= fp;
        int cnt = 0;
        for (int val = 0; val < N; val++)
            if (!used[val])
            {
                if (cnt == idx) { v[pos] = val; used[val] = true; break; }
                cnt++;
            }
    }
}

static inline uint64_t fnv1a(const string& s)
{
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char ch : s) { h ^= ch; h *= 1099511628211ULL; }
    return h;
}
static inline int ownerOf(const string& s) { return (int)(fnv1a(s) % (uint64_t)mpiSize); }

// ---------------- checkpoint helpers (fixed E-byte records) ----------------

static string ckName(const char* what, int forRank = -1)
{
    char buf[128];
    if (forRank >= 0) snprintf(buf, sizeof(buf), "ck_cg_N%d_%s_r%dof%d.bin", N, what, forRank, mpiSize);
    else              snprintf(buf, sizeof(buf), "ck_cg_N%d_%s.bin", N, what);
    return string(buf);
}

static void writeRecords(const string& fname, const vector<string>& recs)
{
    FILE* f = fopen(fname.c_str(), "wb");
    if (!f) { fprintf(stderr, "cannot write %s\n", fname.c_str()); MPI_Abort(MPI_COMM_WORLD, 3); }
    long long n = (long long)recs.size();
    fwrite(&n, sizeof(n), 1, f);
    for (const string& r : recs) fwrite(r.data(), 1, (size_t)E, f);
    fclose(f);
}

static bool readRecords(const string& fname, vector<string>& recs)
{
    FILE* f = fopen(fname.c_str(), "rb");
    if (!f) return false;
    long long n = 0;
    if (fread(&n, sizeof(n), 1, f) != 1) { fclose(f); return false; }
    recs.clear(); recs.reserve((size_t)n);
    vector<char> buf((size_t)E);
    for (long long i = 0; i < n; i++)
    {
        if (fread(buf.data(), 1, (size_t)E, f) != (size_t)E) { fclose(f); return false; }
        recs.emplace_back(buf.data(), (size_t)E);
    }
    fclose(f);
    return true;
}

// ---------------- distributed round machinery ----------------

// route raw strings to their owners; owners dedup against the persistent
// shard, canonicalize the fresh ones; new canonical entries are allgathered,
// sorted and appended to the replicated catalogue; returns the new frontier.
// The exchange runs in batches of at most gBatchBytes send volume per rank:
// every rank sends at most nPerDest strings to each destination per batch,
// so all int counts and displacements MPI_Alltoallv requires stay < 2^31
// regardless of the round's total raw volume (at N>=13 it exceeds 2 GB).
static long long gBatchBytes = 256LL << 20;

static vector<string> processRaws(vector<string>& local,
                                  unordered_set<string>& rawShard,
                                  unordered_set<string>& seen,
                                  vector<string>& seenList)
{
    // group by owner (moves, no copies)
    vector<vector<string>> byOwner((size_t)mpiSize);
    for (string& s : local) byOwner[(size_t)ownerOf(s)].push_back(std::move(s));
    local.clear(); local.shrink_to_fit();

    const long long nPerDest = (gBatchBytes / mpiSize / E) > 0 ? gBatchBytes / mpiSize / E : 1;

    vector<size_t> pos((size_t)mpiSize, 0);
    vector<string> fresh;
    for (;;)
    {
        long long left = 0;
        for (int r = 0; r < mpiSize; r++) left += (long long)byOwner[(size_t)r].size() - (long long)pos[(size_t)r];
        long long leftAll = 0;
        MPI_Allreduce(&left, &leftAll, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (leftAll == 0) break;

        vector<int> scounts(mpiSize, 0), sdispls(mpiSize, 0), rcounts(mpiSize, 0), rdispls(mpiSize, 0);
        vector<long long> take((size_t)mpiSize, 0);
        long long sendTotal = 0;
        for (int r = 0; r < mpiSize; r++)
        {
            long long avail = (long long)byOwner[(size_t)r].size() - (long long)pos[(size_t)r];
            take[(size_t)r] = avail < nPerDest ? avail : nPerDest;
            scounts[(size_t)r] = (int)(take[(size_t)r] * E);
            sdispls[(size_t)r] = (int)sendTotal;          // <= gBatchBytes, int-safe
            sendTotal += take[(size_t)r] * E;
        }
        vector<char> sendbuf((size_t)sendTotal);
        for (int r = 0; r < mpiSize; r++)
            for (long long i = 0; i < take[(size_t)r]; i++)
            {
                string& s = byOwner[(size_t)r][pos[(size_t)r] + (size_t)i];
                memcpy(sendbuf.data() + (size_t)sdispls[(size_t)r] + (size_t)(i * E), s.data(), (size_t)E);
                string().swap(s);                         // free as we pack
            }
        for (int r = 0; r < mpiSize; r++) pos[(size_t)r] += (size_t)take[(size_t)r];

        MPI_Alltoall(scounts.data(), 1, MPI_INT, rcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
        long long recvTotal = 0;                          // <= mpiSize*nPerDest*E = gBatchBytes
        for (int r = 0; r < mpiSize; r++) { rdispls[(size_t)r] = (int)recvTotal; recvTotal += rcounts[(size_t)r]; }
        vector<char> recvbuf((size_t)recvTotal);
        MPI_Alltoallv(sendbuf.data(), scounts.data(), sdispls.data(), MPI_BYTE,
                      recvbuf.data(), rcounts.data(), rdispls.data(), MPI_BYTE, MPI_COMM_WORLD);
        sendbuf.clear(); sendbuf.shrink_to_fit();

        // owner-side dedup against the persistent shard
        for (long long off = 0; off < recvTotal; off += E)
        {
            string r(recvbuf.data() + off, (size_t)E);
            if (rawShard.insert(r).second) fresh.push_back(std::move(r));
        }
    }
    byOwner.clear(); byOwner.shrink_to_fit();

    // canonicalize fresh raws
    vector<string> cans(fresh.size());
    #pragma omp parallel for schedule(dynamic)
    for (long long i = 0; i < (long long)fresh.size(); i++)
        cans[(size_t)i] = canon(fresh[(size_t)i]);
    sort(cans.begin(), cans.end());
    cans.erase(unique(cans.begin(), cans.end()), cans.end());
    vector<string> mineNew;
    for (const string& c : cans) if (!seen.count(c)) mineNew.push_back(c);

    // allgather the new canonical entries (tiny volumes; guarded regardless)
    long long myBytesLL = (long long)mineNew.size() * E;
    if (myBytesLL > (1LL << 30)) { fprintf(stderr, "allgather volume overflow\n"); MPI_Abort(MPI_COMM_WORLD, 4); }
    int myBytes = (int)myBytesLL;
    vector<int> gcounts(mpiSize, 0), gdispls(mpiSize, 0);
    MPI_Allgather(&myBytes, 1, MPI_INT, gcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    long long gTotal = 0;
    for (int r = 0; r < mpiSize; r++) { gdispls[r] = (int)gTotal; gTotal += gcounts[r]; }
    if (gTotal > (1LL << 31) - 1) { fprintf(stderr, "allgather displacement overflow\n"); MPI_Abort(MPI_COMM_WORLD, 4); }
    vector<char> gsend((size_t)myBytes), grecv((size_t)gTotal);
    for (size_t i = 0; i < mineNew.size(); i++) memcpy(gsend.data() + i * (size_t)E, mineNew[i].data(), (size_t)E);
    MPI_Allgatherv(gsend.data(), myBytes, MPI_BYTE,
                   grecv.data(), gcounts.data(), gdispls.data(), MPI_BYTE, MPI_COMM_WORLD);

    vector<string> allNew;
    for (long long off = 0; off < gTotal; off += E) allNew.emplace_back(grecv.data() + off, (size_t)E);
    sort(allNew.begin(), allNew.end());          // determinism: identical order on all ranks
    allNew.erase(unique(allNew.begin(), allNew.end()), allNew.end());

    vector<string> next;
    for (const string& c : allNew)
        if (seen.insert(c).second) { seenList.push_back(c); next.push_back(c); }
    return next;
}

// stripe the (entry x permutation-chunk) units over ranks; OpenMP inside
static vector<string> sweepLocal(const vector<string>& frontier)
{
    const long long NF = factorialLL(N);
    const long long CHUNK = 4LL << 20;
    const long long cpe = (NF + CHUNK - 1) / CHUNK;
    const long long total = (long long)frontier.size() * cpe;

    unordered_set<string> lset;
    #pragma omp parallel
    {
        unordered_set<string> tset;
        #pragma omp for schedule(dynamic)
        for (long long u = 0; u < total; u++)
        {
            if (u % (long long)mpiSize != (long long)mpiRank) continue;
            const string& c = frontier[(size_t)(u / cpe)];
            long long k0 = (u % cpe) * CHUNK;
            long long k1 = (NF < k0 + CHUNK) ? NF : (k0 + CHUNK);
            int sig[16];
            unrankPerm(k0, sig);
            for (long long k = k0; k < k1; k++)
            {
                tset.insert(coarsen(c, sig));
                if (!next_permutation(sig, sig + N)) break;
            }
            if (tset.size() > 500000)
            {
                #pragma omp critical(mergeraws)
                lset.insert(tset.begin(), tset.end());
                tset.clear();
            }
        }
        #pragma omp critical(mergeraws)
        lset.insert(tset.begin(), tset.end());
    }
    // move the strings out instead of copying (halves the peak footprint)
    vector<string> out;
    out.reserve(lset.size());
    while (!lset.empty()) out.push_back(std::move(lset.extract(lset.begin()).value()));
    return out;
}

int main(int argc, char** argv)
{
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    int Nmin = 3, Nmax = 9;
    if (argc > 1) Nmin = Nmax = atoi(argv[1]);
    if (argc > 2) Nmax = atoi(argv[2]);
    bool verbose = getenv("VERBOSE") != NULL;
    bool trace   = getenv("TRACE") != NULL;
    bool resume  = getenv("RESUME") != NULL;
    bool ckpt    = getenv("NOCKPT") == NULL;
    if (const char* bm = getenv("BATCH_MB"))
    { long long v = atoll(bm); if (v >= 1) gBatchBytes = v << 20; }

    for (N = Nmin; N <= Nmax; N++)
    {
        double t0 = now();
        buildEdges();

        int lab0[128]; for (int e = 0; e < E; e++) lab0[e] = e;
        string discrete = rgsOf(lab0);

        unordered_set<string> seen, rawShard;
        vector<string> seenList, frontier;
        int rounds = 0;
        bool restored = false;

        if (resume)
        {
            // meta: N rounds mpiSize (text); shards and lists in binary files
            FILE* f = fopen(ckName("meta").c_str(), "r");
            if (f)
            {
                int mN = 0, mR = 0, mS = 0;
                if (fscanf(f, "%d %d %d", &mN, &mR, &mS) == 3 && mN == N)
                {
                    if (mS != mpiSize)
                    { if (!mpiRank) fprintf(stderr, "RESUME: checkpoint has %d ranks, running with %d - abort\n", mS, mpiSize); MPI_Abort(MPI_COMM_WORLD, 5); }
                    vector<string> tmp;
                    if (readRecords(ckName("seen"), tmp))
                    {
                        for (const string& s : tmp) { seen.insert(s); seenList.push_back(s); }
                        readRecords(ckName("frontier"), frontier);
                        vector<string> shard;
                        readRecords(ckName("raw", mpiRank), shard);
                        for (string& s : shard) rawShard.insert(std::move(s));
                        rounds = mR;
                        restored = true;
                        if (!mpiRank)
                            fprintf(stderr, "N=%d RESUME at round %d: total %zu, frontier %zu, my shard %zu\n",
                                    N, rounds, seenList.size(), frontier.size(), rawShard.size());
                    }
                }
                fclose(f);
            }
        }

        auto saveCheckpoint = [&](void)
        {
            if (ckpt)
            {
                writeRecords(ckName("raw", mpiRank), vector<string>(rawShard.begin(), rawShard.end()));
                if (!mpiRank)
                {
                    writeRecords(ckName("seen"), seenList);
                    writeRecords(ckName("frontier"), frontier);
                    FILE* f = fopen(ckName("meta").c_str(), "w");
                    fprintf(f, "%d %d %d\n", N, rounds, mpiSize);
                    fclose(f);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        };

        if (!restored)
        {
            seen.insert(discrete); seenList.push_back(discrete);

            // ---- round 1: one candidate permutation per cycle type (rank 0) ----
            rounds = 1;
            vector<string> local;
            if (!mpiRank)
            {
                vector<vector<int>> parts; vector<int> tmp;
                genPartitions(N, N, tmp, parts);
                unordered_set<string> r1;
                for (auto& pt : parts)
                {
                    bool ident = true;
                    for (int L : pt) if (L > 1) { ident = false; break; }
                    if (ident) continue;
                    int sig[16]; int o = 0;
                    for (int L : pt) { for (int k = 0; k < L; k++) sig[o + k] = o + (k + 1) % L; o += L; }
                    r1.insert(coarsen(discrete, sig));
                }
                local.assign(r1.begin(), r1.end());
            }
            frontier = processRaws(local, rawShard, seen, seenList);
            if (!mpiRank)
                fprintf(stderr, "N=%d round 1 [np=%d]: cycle types, new %zu, total %zu (%.1fs)\n",
                        N, mpiSize, frontier.size(), seenList.size(), now() - t0);
            saveCheckpoint();
        }

        // ---- rounds >= 2 ----
        while (!frontier.empty())
        {
            rounds++;
            double tr = now();
            vector<string> local = sweepLocal(frontier);
            long long nloc = (long long)local.size(), nlocSum = 0;
            MPI_Reduce(&nloc, &nlocSum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            vector<string> next = processRaws(local, rawShard, seen, seenList);
            if (!mpiRank)
                fprintf(stderr, "N=%d round %d [np=%d]: frontier %zu, raw(local-unique) %lld, new %zu, total %zu (%.1fs)\n",
                        N, rounds, mpiSize, frontier.size(), nlocSum, next.size(), seenList.size(), now() - tr);
            if (trace && rounds >= 4 && !mpiRank)
                for (const string& cc : next)
                {
                    fprintf(stderr, "TRACE round %d order %lld:", rounds, autOrder(cc));
                    for (int e = 0; e < E; e++) fprintf(stderr, " %d", (int)cc[(size_t)e]);
                    fprintf(stderr, "\n");
                }
            frontier.swap(next);
            saveCheckpoint();
        }

        // ---- order histogram: stripe entries over ranks ----
        long long nSeen = (long long)seenList.size();
        vector<long long> ords((size_t)nSeen, 0), ordsAll((size_t)nSeen, 0);
        #pragma omp parallel for schedule(dynamic)
        for (long long i = 0; i < nSeen; i++)
            if (i % (long long)mpiSize == (long long)mpiRank)
                ords[(size_t)i] = autOrder(seenList[(size_t)i]);
        MPI_Reduce(ords.data(), ordsAll.data(), (int)nSeen, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

        if (!mpiRank)
        {
            map<long long, vector<const string*>> byOrder;
            for (long long i = 0; i < nSeen; i++) byOrder[ordsAll[(size_t)i]].push_back(&seenList[(size_t)i]);
            printf("N=%d: closed subgroup classes = %lld (productive rounds: %d, %.1fs, np=%d)\n",
                   N, nSeen, rounds - 1, now() - t0, mpiSize);
            printf("  orders:");
            for (auto& kv : byOrder) printf(" %lld(x%zu)", kv.first, kv.second.size());
            printf("\n");
            if (verbose)
                for (auto& kv : byOrder)
                {
                    printf("--- order %lld (x%zu) ---\n", kv.first, kv.second.size());
                    for (const string* c : kv.second)
                    {
                        for (int e = 0; e < E; e++) printf("%d ", (int)(*c)[(size_t)e]);
                        printf("\n");
                    }
                }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
