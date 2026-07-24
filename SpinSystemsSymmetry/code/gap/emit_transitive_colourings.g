# emit_transitive_colourings.g  —  transitive part of a(N).
#
# When the full subgroup lattice of S_N is out of reach (GAP's cyclic-extension
# ConjugacyClassesSubgroups exhausts memory in its Zuppos phase around N=14),
# the *transitive* closed groups are still immediate: they are among the
# tabulated transitive groups of degree N (GAP transitive-groups library).
# This driver emits one edge-orbit colouring per transitive group of degree N,
# in the same "N |G| c_1 ... c_E" format as emit_subgroup_colourings.g, so
# filter_closed keeps the 2*-closed ones — the single-orbit (fully coupled)
# N-spin systems, one part of a(N) = a(N-1) + f(N).
#
# Run:   gap -q -b -c 'NN:=14;;' emit_transitive_colourings.g > trans.txt
#        grep -v '^#' trans.txt | ./filter_closed -
# Known: degree 14 has 63 transitive groups, 10 of them 2*-closed.

SetPrintFormattingStatus("*stdout*", false);
if not IsBound(NN) then NN := 14; fi;
n := NN;;
G := AllTransitiveGroups(NrMovedPoints, n);;
Print("# transitive degree ", n, " groups=", Length(G), "\n");
for H in G do
  orb := Orbits(H, Combinations([1..n], 2), OnSets);;
  M := NullMat(n, n);;
  for k in [1..Length(orb)] do
    for pr in orb[k] do M[pr[1]][pr[2]] := k; M[pr[2]][pr[1]] := k; od;
  od;
  line := Concatenation(String(n), " ", String(Size(H)));
  for j in [2..n] do for i in [1..j-1] do
    line := Concatenation(line, " ", String(M[i][j]));
  od; od;
  Print(line, "\n");
od;
QUIT;
