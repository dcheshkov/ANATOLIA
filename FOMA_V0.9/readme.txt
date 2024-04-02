*************************************************************************
***                           FOMA V0.9                               ***
***      Program for total lineshape analysis of first-order NMR      ***
***           multiplets (due to couplings with spins-1/2)            ***
***       designed as AU program for Buker TopSpin for Windows        ***
*************************************************************************
***         (C) 2024 Dmitry Cheshkov, Dmitry Sinitsyn                 ***
*************************************************************************
***                     dcheshkov@gmail.com                           ***
***                  http://anatolia.nmrclub.ru                       ***
***             https://github.com/dcheshkov/ANATOLIA                 ***
*************************************************************************

The First Order Multiplet Analyzer (FOMA) is actually a lightweight version
of ANATOLIA, designed for total lineshape analysis of a single multiplet
using an additional broadening approach to eliminate local minima.

FOMA is Bruker TopSpin AU program. File 'foma' should be placed to
'TopSpinHome\exp\stan\nmr\au\src\user\' folder.

To compile, file makeau 'TopSpinHome\exp\stan\nmr\au\makeau' shold modified
in a such a way, that The GDI32 library must be added at link time.

you need to change the file TopSpin Home\exp\stan\nmr\au\makeau. 

In oder to compile, file makeau 'TopSpinHome\exp\stan\nmr\au\makeau' should
modified. GDI32 library should be added on linking stage. 


In section 'my %ldflags' in line containing 'Windows_GCC'
-lgdi32 should be added.

Original section:

my %ldflags = (
	    'LINUX_INTEL' => ' -lm -ldl',
	    'MAC_INTEL'	=> ' -lm -ldl',
	    'Windows_GCC'=> ' -lm -lwinspool -lwsock32',
	    'Windows_NT'=> " /nologo /INCREMENTAL:no"
			  ." kernel32.lib advapi32.lib gdi32.lib comdlg32.lib"
			  ." user32.lib shell32.lib winspool.lib wsock32.lib"
    );

Modified:

my %ldflags = (
	    'LINUX_INTEL' => ' -lm -ldl',
	    'MAC_INTEL'	=> ' -lm -ldl',
	    'Windows_GCC'=> ' -lm -lwinspool -lwsock32 -lgdi32',
	    'Windows_NT'=> " /nologo /INCREMENTAL:no"
			  ." kernel32.lib advapi32.lib gdi32.lib comdlg32.lib"
			  ." user32.lib shell32.lib winspool.lib wsock32.lib"
    );

To run the program, one needs to define a single integral region containing the multiplet of
interest. The program takes one command line parameter corresponding to the number of coupling
constants or multiplet type in notation like am3x2n. After the program starts up, an interactive
window with parameters to optimize is opened. For example, to analyze the 'A' multiplet of
an 'AM3' spin system, 'foma 3' or 'foma am3' should be written. In the interactive window,
it is possible to define starting values for coupling constants, but for many cases,
the default values are suitable.