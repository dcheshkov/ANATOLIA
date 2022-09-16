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

                          Tutorial examples

There are several examples of 1H NMR spectra analysis (ODCB, styrene, etc).
Spectra are stored in the Bruker NMR data format as datasets. So the folder
names 'ODCB' and 'Styrene' corresponds to 'ExpName'. All the necessary input
files for ANATOLIA are located also there.

User should copy the entire dataset to the place of TopSpin spectral data
(e.g. '/u/data/user/nmr/' or 'D:\data\user\nmr\').

Then, put the AU program 'anatolia' (from the TopSpinAU folder) to the
'TopSpinHomeDirectory/exp/stan/nmr/au/src/user', create the folder
'TopSpinHomeDirectory/prog/anatolia' and put there ANATOLIA executable
file. For now, ANATOLIA can be executed from TopSpin command line by
entering the 'anatolia' command.

Alternatively, ANATOLIA executable file can be placed to the
'ExpNo' folder of the given examples and executed there directly.