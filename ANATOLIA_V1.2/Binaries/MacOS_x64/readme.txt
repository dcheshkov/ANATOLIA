ANATOLIA V1.2
Statically compiled and linked. Single EXEcutable file
Don't forget to set executable flag (chmod 755 ANATOLIA).

Compilled by the command (with statically linked GSL-2.7 library):
c++ -O2 -ffast-math -fomit-frame-pointer -fstrict-aliasing -funroll-loops anatolia.cpp /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a -o ANATOLIA
