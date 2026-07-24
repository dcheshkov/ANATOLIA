# emit_fixedfree_colourings.g  —  orbit-partition decomposition for f(N).
#
# Computes the "essentially N-spin" closed groups without a full S_N lattice.
# By the tower theorem  a(N) = a(N-1) + f(N),  where f(N) counts the 2^-closed
# groups of degree N with NO fixed point, i.e. whose orbit partition has every
# part >= 2.  A subgroup with orbit partition lambda lives in the Young subgroup
# Y_lambda = prod S_{lambda_i} and is transitive on each block; so for each
# lambda of N with min part >= 2 we enumerate the block-transitive subgroups of
# Y_lambda and emit their edge colourings.  filter_closed then keeps the
# 2^-closed ones and deduplicates by canonical colouring (which also merges the
# block-permutation conjugates that Y_lambda-classes do not).
#
# This never touches the lattice of S_N itself: the single-block case lambda=[N]
# uses the transitive-groups library, and every multi-block Y_lambda has largest
# block <= N-2 (for N=14, <= 12), all tractable where the S_14 lattice is not.
#
# Run:   gap -o 16g -q -b -c 'NN:=14;;' emit_fixedfree_colourings.g > ff14.txt
#        grep -v '^#' ff14.txt | ./filter_closed -           # -> f(N)
# Then   a(N) = a(N-1) + f(N).
# Validation: f(4..13) = 5, 3, 16, 9, 54, 41, 151, 112, 554, 368.

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

emitFN := function(n)
  local parts, lam, k, Y, reps, H, t;
  parts := Filtered(Partitions(n), lam -> Minimum(lam) >= 2);;
  Print("# f-decomposition n=", n, " partitions(parts>=2)=", Length(parts), "\n");
  for lam in parts do
    t := Runtime();;
    k := Length(lam);
    if k = 1 then
      reps := AllTransitiveGroups(NrMovedPoints, n);;
      for H in reps do emitH(H, n); od;
      Print("# lambda=[", n, "] transitive groups=", Length(reps),
            " (", Int((Runtime()-t)/1000), "s)\n");
    else
      Y := DirectProduct(List(lam, SymmetricGroup));;
      reps := List(ConjugacyClassesSubgroups(Y), Representative);;
      reps := Filtered(reps, H -> Length(Orbits(H, [1..n])) = k);;  # block-transitive
      for H in reps do emitH(H, n); od;
      Print("# lambda=", lam, " Y-classes->block-transitive=", Length(reps),
            " (", Int((Runtime()-t)/1000), "s)\n");
    fi;
  od;
end;;

if not IsBound(NN) then NN := 6; fi;
emitFN(NN);
QUIT;
