*************************************************************************
***                             ANATOLIA                              ***
***              Program for total lineshape analysis                 ***
***                         of NMR spectra                            ***
*************************************************************************
***           (C) 2017 Dmitry Cheshkov, Dmitry Sinitsyn,              ***
***                        Kirill Sheberstov                          ***
*************************************************************************
***                       dcheshkov@gmail.com                         ***
***                    http://anatolia.nmrclub.ru                     ***
***               https://github.com/dcheshkov/ANATOLIA               ***
*************************************************************************
***   D.A. Cheshkov, K.F. Sheberstov, D.O. Sinitsyn, V.A. Chertkov,   ***
***  ANATOLIA: NMR software for spectral analysis of total lineshape, ***
***      Magn. Reson. Chem., 2018, 56, 449, DOI: 10.1002/mrc.4689.    ***
*************************************************************************

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

Our recent papers:
Total lineshape analysis of a-tetrahydrofuroic acid 1H NMR spectra - http://doi.org/10.48550/arXiv.2209.03708
Total line shape analysis of high-resolution NMR spectra -  http://doi.org/10.1016/bs.arnmr.2019.11.001