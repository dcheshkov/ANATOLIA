ANATOLIA V1.2
Statically compiled and linked. Single EXEcutable file
Don't forget to set executable flag (chmod 755 ANATOLIA).

Compilled by the command (with GSL-2.7 library):
c++ -s -static -O2 -ffast-math -fomit-frame-pointer -fstrict-aliasing -funroll-loops anatolia.cpp -lgsl -lgslcblas -o ANATOLIA
