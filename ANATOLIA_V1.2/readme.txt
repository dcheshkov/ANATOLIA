*************************************************************************
***                        ANATOLIA V1.2                              ***
***              Program for total lineshape analysis                 ***
***                         of NMR spectra                            ***
*************************************************************************
***         (C) 2021 Dmitry Cheshkov, Dmitry Sinitsyn,                ***
***                     Kirill Sheberstov                             ***
*************************************************************************
***                     dcheshkov@gmail.com                           ***
***                  http://anatolia.nmrclub.ru                       ***
***             https://github.com/dcheshkov/ANATOLIA                 ***
*************************************************************************
***   D.A. Cheshkov, K.F. Sheberstov, D.O. Sinitsyn, V.A. Chertkov,   ***
***  ANATOLIA: NMR software for spectral analysis of total lineshape, ***
***      Magn. Reson. Chem., 2018, 56, 449, DOI: 10.1002/mrc.4689.    ***
*************************************************************************

      V1.2 vs. V1.0 & V1.1 contains some improvements with respect
         to performance and efficiency of calculating algorithm.
           No changes were made to the format of input files.

            Working with Bruker NMR format (dataset) only.

           Can be run in standalone mode or via AU programm
                   'anatolia' from Bruker TopSpin.

              When using with TopSpin it should be located
                in TopSpinHome/prog/anatolia directory.
 
            Versions for ALL platforms compiled and linked 
            in 'static' mode, so they are a distributed as
                        single executable files.
 
            All source code currently located in single c++
                          file 'anatolia.cpp'.
 
        For matrix diagonalization it utilizes procedures from
       GNU scientific library (https://www.gnu.org/software/gsl/).
          Powell's BOBYQA algoritm used for multidimentional
                        function minimization.
 
         ANATOLIA reads the input files from current directory,
           when called without arguments, or can take path
            to working directory as command line argument.