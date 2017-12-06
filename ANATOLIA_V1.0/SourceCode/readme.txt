*************************************************************************
***                        ANATOLIA V1.0                              ***
***              Program for total lineshape analysis                 ***
***                         of NMR spectra                            ***
*************************************************************************
***         (C) 2017 Dmitry Cheshkov, Dmitry Sinitsyn,                ***
***                     Kirill Sheberstov                             ***
*************************************************************************
***                     dcheshkov@gmail.com                           ***
***                  http://anatolia.nmrclub.ru                       ***
***             https://github.com/dcheshkov/ANATOLIA                 ***
*************************************************************************
***   D.A. Cheshkov, K.F. Sheberstov, D.O. Sinitsyn, V.A. Chertkov,   ***
***  ANATOLIA: NMR software for spectral analysis of total lineshape, ***
***                    Magn. Reson. Chem., 2017.                      ***
*************************************************************************

                 Instructions for ANATOLIA compilation

Program utilizes GNU scientific library.
All source code is placed in single file anatolia.cpp

ANATOLIA V_1.0
anatolia.cpp 88272 bytes,
MD5  (anatolia.cpp) = 5bfc48f957cdb36681c2f14c6e7a3585
SHA1 (anatolia.cpp) = 711fa924aa7a09332647636648c1b423d0977cf5

For compilation on Unix type systems g++ compiller and
GNU scientific library (GSL, for code development) should be installed.

Debian Linux compilation commands:
apt-get install g++ libgsl0-dev (with root privileges)
wget http://anatolia.nmrclub.ru/ANATOLIA_v1.0/SourceCode/anatolia.cpp
g++ -std=c++11 -static anatolia.cpp -lgsl -lgslcblas -o ANATOLIA

FreeBSD compilation commands:
pkg install gsl (with root privileges)
wget http://anatolia.nmrclub.ru/ANATOLIA_v1.0/SourceCode/anatolia.cpp
c++ -std=c++11 -static -I/usr/local/include -L/usr/local/lib anatolia.cpp -lgsl -lgslcblas -o ANATOLIA

MacOS GSL instalation & ANATOLIA compilation commands:
curl -O ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
gzip -dc gsl-latest.tar.gz | tar xvf -
cd gsl-X.Y (for now it's gsl-2.4)
./configure
make
make install (with root privileges)
cd ..
rm -rf gsl-*
curl -O http://anatolia.nmrclub.ru/ANATOLIA_v1.0/SourceCode/anatolia.cpp
g++ -std=c++11 anatolia.cpp /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a -o ANATOLIA

Windows
ANATOLIA MS Visual Studio Projects with compiled GSL library can be downloaded from
the following links:
http://anatolia.nmrclub.ru/ANATOLIA_v1.0/SourceCode/MSVC17_Project.rar
for Microsoft Visual Studio 2017 (x86 and x64), and
http://anatolia.nmrclub.ru/ANATOLIA_v1.0/SourceCode/MSVC10_Project_WinXP.rar
for Microsoft Visual Studio 2010 (x86, which we use for CLR4.0 WinXP binaries generation).

Users are offered to use appropriate compiled binary file, but if there is no such file
for your platform, feel free to contact with Dmitry Cheshkov (dcheshkov@gmail.com) on
any problems with program compilation.