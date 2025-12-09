ANATOLIA-X
Statically compiled and linked. Single EXEcutable file
Don't forget to set executable flag (chmod 755 ANATOLIA-X).

Compilled by the command (with GSL-2.7 library):
c++ -s -static -O2 -DNDEBUG -ffast-math -fomit-frame-pointer -fstrict-aliasing -funroll-loops anatolia_x.cpp -lgsl -lgslcblas -o ANATOLIA-X