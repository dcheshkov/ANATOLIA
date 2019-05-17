ANATOLIA Analysis of NMR spectra total lineshape.

Program utilizes GNU scientific library.
All source code is placed in one file anatolia.cpp

ANATOLIA V_1.0
MD5 (anatolia.cpp) = bfe3a0d5bb7445315644e5f4cc7656b6
SHA1 (anatolia.cpp) = dc5cf942fe3e7a9d69cc9f3ba08d4207b76b22c6

For compilation on Unix type systems g++ compiller and
GNU scientific library (GSL, for code development) should be installed.

Debian Linux compilation commands:
apt-get install g++ libgsl0-dev (with root privileges)
wget http://anatolia.nmrclub.ru/ANATOLIA_v1.0/source_code/anatolia.cpp
g++ -std=c++11 -static anatolia.cpp -lgsl -lgslcblas -o ANATOLIA

FreeBSD compilation commands:
pkg install gsl (with root privileges)
wget http://anatolia.nmrclub.ru/ANATOLIA_v1.0/source_code/anatolia.cpp
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
curl -O http://anatolia.nmrclub.ru/ANATOLIA_v1.0/source_code/anatolia.cpp
g++ -std=c++11 anatolia.cpp /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a -o ANATOLIA

Windows
ANATOLIA MS Visual Studio Projects with compiled GSL library can be downloaded from
the following links:
http://anatolia.nmrclub.ru/ANATOLIA_v1.0/source_code/MSVC17_Project.rar
for Microsoft Visual Studio 2017, and
http://anatolia.nmrclub.ru/ANATOLIA_v1.0/source_code/MSVC10_Project_WinXP.rar
for Microsoft Visual Studio 2010 (which we use for CLR4.0 WinXP binaries generation).

Users are offered to use appropriate compiled binary file, but if there is no such file
for your platform, feel free to contact with Dmitry Cheshkov
(dcheshkov at gmail dot com) on any problems with program compilation.