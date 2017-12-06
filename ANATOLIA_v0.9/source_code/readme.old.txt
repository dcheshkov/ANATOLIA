ANATOLIA Analysis of NMR spectra total lineshape.

Program utilizes GNU scientific library.
All source code is placed in one file anatolia.cpp

For V_1.0
MD5 (anatolia.cpp) = bfe3a0d5bb7445315644e5f4cc7656b6
SHA1 (anatolia.cpp) = dc5cf942fe3e7a9d69cc9f3ba08d4207b76b22c6

Linux compilation:
g++ compiller and GNU scientific library (devel) should be installed.
Use the command:
g++ -std=c++11 -static anatolia.cpp -lgsl -lgslcblas -o ANATOLIA
