# emit_subdirect2.g  —  fast block-transitive emitter for a 2-block partition.
#
# For a two-block orbit partition [a,b] of N=a+b, the block-transitive subgroups
# of S_a x S_b are exactly the subdirect products of a transitive group on the
# a-block and a transitive group on the b-block (Goursat): every such subgroup
# projects onto a transitive group on each block, and is determined by a common
# quotient of the two.  Iterating over the (small) transitive-groups libraries
# and calling GAP's SubdirectProducts is far cheaper than building the entire
# subgroup lattice of the Young subgroup S_a x S_b (the balanced partitions
# [8,6], [9,5], [7,7] are the bottleneck of the full-lattice route).
#
# DirectProduct(A,B) places the transitive A on [1..a] and B on [a+1..a+b], so
# each subdirect product already acts on [1..N] with the right block placement.
# Output format and downstream use are identical to emit_partition_colourings.g;
# the filter_closed canonical-form dedup merges the conjugates.
#
# Run:  gap -q -b -c 'NN:=14;; PART:=[8,6];;' emit_subdirect2.g > p_8_6.txt

SetPrintFormattingStatus("*stdout*", false);
SetUserPreference("AtlasRep", "AtlasRepAccessRemoteFiles", false);

emitH := function(H, n)
  local orb, M, k, pr, i, j, line;
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
end;;

if not IsBound(NN) then Error("set NN"); fi;
if not IsBound(PART) or Length(PART) <> 2 then Error("PART must be a 2-block partition [a,b]"); fi;
n := NN;; a := PART[1];; b := PART[2];;
if a + b <> n then Error("PART must sum to NN"); fi;

t := Runtime();;
Alist := AllTransitiveGroups(NrMovedPoints, a);;
Blist := AllTransitiveGroups(NrMovedPoints, b);;
cnt := 0;;
for A in Alist do
  for B in Blist do
    for H in SubdirectProducts(A, B) do
      emitH(H, n);
      cnt := cnt + 1;
    od;
  od;
od;
Print("# N=", n, " PART=", PART, " subdirect(", Length(Alist), "x", Length(Blist),
      ")=", cnt, " (", Int((Runtime()-t)/1000), "s)\n");
QUIT;
