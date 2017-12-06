Program for automated NMR spectra total lineshape analysis.
MRC...
Working with Bruker NMR format (dataset) only.
Can be run in standalone mode or via AU programm
'anatolia' from Bruker TopSpin.

When using with TopSpin it should be located
in TopSpinHome/prog/anatolia directory.

All source code currently located in single c++
file 'anatolia.cpp'.
For matrix diagonalization and multidimentional
function minimization it utilizes procedures from
GNU scientific library (https://www.gnu.org/software/gsl/).

Versions for ALL platforms compiled and linked 
in 'static' mode, so they are a distributed as
single executable files.
