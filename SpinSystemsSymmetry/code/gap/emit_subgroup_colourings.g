# emit_subgroup_colourings.g  —  GAP driver for the direct 2^-closure filter.
#
# For each degree n in a range, enumerate ALL conjugacy classes of subgroups of
# the symmetric group S_n (Hulpke's cyclic-extension algorithm, = OEIS A000638),
# and print one representative per class as an edge-orbit colouring of K_n:
#
#     n  |G|  c_1 c_2 ... c_E        (E = n(n-1)/2)
#
# where c_k is the orbit index (under G, action OnSets on unordered pairs) of the
# k-th edge in upper-triangle column order  J12 J13 J23 J14 J24 J34 ...  — the
# same edge convention as materials/8Spins.txt and the C++ enumerator.
#
# The companion filter_closed.cpp reads these lines and keeps the 2^-closed ones
# (those realizable as spin-system symmetry groups): a group is closed iff the
# automorphism group of this colouring has order |G|.  The count of closed
# classes per n is a(n).  This reproduces the BFS catalogue (closed_groups.cpp)
# by an independent route — enumerate-then-filter instead of symmetrize-then-BFS.
#
# Two tweaks make the output usable:
#   * AtlasRep remote access is disabled (subgroup enumeration for n>=6 would
#     otherwise try to fetch data over the network and abort);
#   * automatic 80-column line wrapping is turned off (otherwise long colouring
#     lines are split with "\" continuations and the parser mis-reads them).
#
# Run:   gap -q -b -c 'NMIN:=3;; NMAX:=13;;' emit_subgroup_colourings.g > col.txt
#        grep -v '^#' col.txt | ./filter_closed /dev/stdin
# (NMIN/NMAX default to 3..13 if not set on the command line.)

SetUserPreference("AtlasRep", "AtlasRepAccessRemoteFiles", false);
SetPrintFormattingStatus("*stdout*", false);

emitN := function(n)
  local cc, G, orb, M, line, i, j, k, pr;
  cc := List(ConjugacyClassesSubgroups(SymmetricGroup(n)), Representative);
  Print("# n=", n, " classes=", Length(cc), "\n");
  for G in cc do
    orb := Orbits(G, Combinations([1..n], 2), OnSets);
    M := NullMat(n, n);                       # M[i][j] = orbit index of pair {i,j}
    for k in [1..Length(orb)] do
      for pr in orb[k] do
        M[pr[1]][pr[2]] := k;
        M[pr[2]][pr[1]] := k;
      od;
    od;
    line := Concatenation(String(n), " ", String(Size(G)));
    for j in [2..n] do for i in [1..j-1] do    # upper triangle, column order
      line := Concatenation(line, " ", String(M[i][j]));
    od; od;
    Print(line, "\n");
  od;
end;;

if not IsBound(NMIN) then NMIN := 3; fi;
if not IsBound(NMAX) then NMAX := 13; fi;
for n in [NMIN..NMAX] do
  t := Runtime();;
  emitN(n);
  Print("# n=", n, " done in ", Int((Runtime() - t) / 1000), "s\n");
od;
QUIT;
