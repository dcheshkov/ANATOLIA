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


Our recent papers and preprints:
Facilitating Total Lineshape Analysis of Complex NMR Spectra With FOMA and ANATOLIA-X Multiplet Fitting Tools,
Exemplified by the Vinyl Norbornene Case - http://doi.org/10.1002/mrc.70111
The Origin of Mirror Symmetry in High-Resolution Nuclear Magnetic Resonance Spectra - http://doi.org/10.5194/mr-7-15-2026
Mirror Symmetry of the NMR Spectrum and the Connection with the Structure of Spin Hamiltonian Matrix Representations - http://doi.org/10.48550/arXiv.2602.03871
Empirical Reevaluation of Computational Limits in Exact NMR Spectral Simulation: 16 Spin-½ Nuclei on Standard Hardware - http://doi.org/10.13140/RG.2.2.21355.20008
Total lineshape analysis of a-tetrahydrofuroic acid 1H NMR spectra - http://doi.org/10.48550/arXiv.2209.03708
Total line shape analysis of high-resolution NMR spectra -  http://doi.org/10.1016/bs.arnmr.2019.11.001