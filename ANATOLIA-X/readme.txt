*************************************************************************
***                     ANATOLIA-X based on V1.2                      ***
***              Program for total lineshape analysis                 ***
***                         of NMR spectra                            ***
*************************************************************************
***         (C) 2025 Dmitry Cheshkov, Dmitry Sinitsyn,                ***
***                     Kirill Sheberstov                             ***
*************************************************************************
***                     dcheshkov@gmail.com                           ***
***                  http://anatolia.nmrclub.ru                       ***
***             https://github.com/dcheshkov/ANATOLIA                 ***
*************************************************************************

             ANATOLIA-X derived from ANATOLIA V1.2 code.
The spectrum calculation module uses the weak-coupling (X) approximation.

    The algorythm automatically detects strongly coupled spin groups
    based on the criterion |deltaResFreq/J| < 10 and performs explicit
         strong-coupling treatment within these spin groups.

        When specifying a spin system in input files, strongly
              coupled spins must be listed sequentially. 
       For convenience and correct detection, it is recommended
           to specify the spin system by listing spins in
          order of their resonant frequencies (ascending
           or descending), ensuring that strongly coupled
                    spins appear consecutively.

           No changes were made to the format of input files.

            Working with Bruker NMR format (dataset) only.

            Versions for ALL platforms compiled and linked 
            in 'static' mode, so they are a distributed as
                        single executable files.
 
            All source code currently located in single c++
                          file 'anatolia_x.cpp'.
 
        For matrix diagonalization it utilizes procedures from
       GNU scientific library (https://www.gnu.org/software/gsl/).
          Powell's BOBYQA algoritm used for multidimentional
                        function minimization.
 
         ANATOLIA reads the input files from current directory,
           when called without arguments, or can take path
            to working directory as command line argument.