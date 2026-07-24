# emit_partition_colourings.g  —  one orbit-partition's block-transitive subgroups.
#
# The per-partition worker behind the f(N) decomposition (see
# emit_fixedfree_colourings.g).  Given a single orbit partition PART of N (a
# GAP list of block sizes, each >= 2) it emits the edge colourings of the
# block-transitive subgroups of the Young subgroup prod S_{PART_i}; for the
# single-block partition [N] it uses the transitive-groups library.  Splitting
# f(N) this way lets the partitions run as independent parallel jobs (the
# heaviest for N=14 is [12,2] -> S_12 x S_2); concatenate all outputs and pipe
# through filter_closed once to get f(N).
#
# Run one partition:
#   gap -o 32g -q -b -c 'NN:=14;; PART:=[12,2];;' emit_partition_colourings.g > p_12_2.txt
# Drive all partitions in parallel, then filter (bash):
#   N=14
#   gap -q -b -c "N:=$N;; Print(JoinStringsWithSeparator(List(Filtered(Partitions(N),
#       p->Minimum(p)>=2), p->JoinStringsWithSeparator(List(p,String),\",\")),\"\\n\"),\"\\n\");" \
#     | while read p; do echo "$p"; done            # list of partitions "12,2" etc.
#   # (or just: for each partition run the gap line above with PART:=[...])
#   cat p_*.txt | grep -v '^#' | ./filter_closed -  # -> f(N)

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

if not IsBound(NN) then Error("set NN (degree)"); fi;
if not IsBound(PART) then Error("set PART (partition, list of block sizes)"); fi;
n := NN;;
if Sum(PART) <> n then Error("PART must sum to NN"); fi;

t := Runtime();;
if Length(PART) = 1 then
  reps := AllTransitiveGroups(NrMovedPoints, n);;
  for H in reps do emitH(H, n); od;
  Print("# N=", n, " PART=[", n, "] transitive=", Length(reps),
        " (", Int((Runtime()-t)/1000), "s)\n");
else
  Y := DirectProduct(List(PART, SymmetricGroup));;
  reps := List(ConjugacyClassesSubgroups(Y), Representative);;
  reps := Filtered(reps, H -> Length(Orbits(H, [1..n])) = Length(PART));;
  for H in reps do emitH(H, n); od;
  Print("# N=", n, " PART=", PART, " block-transitive=", Length(reps),
        " (", Int((Runtime()-t)/1000), "s)\n");
fi;
QUIT;
