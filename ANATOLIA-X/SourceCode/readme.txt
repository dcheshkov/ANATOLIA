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

    The algorithm automatically detects strongly coupled spin groups
    based on the criterion |deltaResFreq/J| < 10 and performs explicit
         strong-coupling treatment within these spin groups.

        When specifying a spin system in input files, strongly
              coupled spins must be listed sequentially. 
       For convenience and correct detection, it is recommended
           to specify the spin system by listing spins in
          order of their resonant frequencies (ascending
           or descending), ensuring that strongly coupled
                    spins appear consecutively.

           No changes were made to the format of input
                   files compared to ANATOLIA V1.0.

            Working with Bruker NMR format (dataset) only.

        For matrix diagonalization it utilizes procedures from
       GNU scientific library (https://www.gnu.org/software/gsl/).
          Powell's BOBYQA algorithm used for multidimensional
                        function minimization.

          All source code is currently located in a single c++
                        file 'anatolia_x.cpp'.

      ANATOLIA-X binary versions for ALL platforms are compiled
       and linked in 'static' mode, so they are distributed as
          single executable files and can be downloaded from
 https://github.com/dcheshkov/ANATOLIA/tree/master/ANATOLIA-X/Binaries

     ANATOLIA-X reads the input files from the current directory,
           when called without arguments, or can take path
            to working directory as command line argument.


ANATOLIA-X (based on ANATOLIA V1.2)
anatolia_x.cpp 96175 bytes,
MD5  (anatolia_x.cpp) = f1c70db953cfacc42a8c639d5962a567
SHA1 (anatolia_x.cpp) = 65af89f01f7c406eb1df5cdb8c22fe412d8259e1


                Instructions for ANATOLIA-X compilation

For compilation on Unix-type systems, the g++ compiler and
GNU scientific library (GSL, for code development) should be installed.

Debian Linux compilation commands:
apt-get install g++ libgsl-dev (with root privileges)
wget http://anatolia.nmrclub.ru/ANATOLIA-X/SourceCode/anatolia_x.cpp
g++ -s -static anatolia_x.cpp -lgsl -lgslcblas -o ANATOLIA-X

FreeBSD compilation commands:
pkg install gsl (with root privileges)
wget http://anatolia.nmrclub.ru/ANATOLIA-X/SourceCode/anatolia_x.cpp
c++ -s -static -I/usr/local/include -L/usr/local/lib anatolia_x.cpp -lgsl -lgslcblas -o ANATOLIA-X

MacOS GSL installation & ANATOLIA compilation commands:
curl -O ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
gzip -dc gsl-latest.tar.gz | tar xvf -
cd gsl-X.Y (for now it's gsl-2.8)
./configure
make
make install (with root privileges)
cd ..
rm -rf gsl-*
curl -O http://anatolia.nmrclub.ru/ANATOLIA-X/SourceCode/anatolia_x.cpp
g++ anatolia_x.cpp /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a -o ANATOLIA-X

Windows
ANATOLIA-X MS Visual Studio Projects with compiled GSL library (v2.8) can be
downloaded from the following links:
http://anatolia.nmrclub.ru/ANATOLIA-X/SourceCode/MSVC2026_Project_x64.7z
https://github.com/dcheshkov/ANATOLIA/blob/master/ANATOLIA-X/SourceCode/MSVC2026_Project_x64.7z
for Microsoft Visual Studio 2026 (x64).

Users are encouraged to use the appropriate compiled binary file, but if there
is no such file for your platform, feel free to contact Dmitry Cheshkov
(dcheshkov@gmail.com) regarding any problems with program compilation
and usage.