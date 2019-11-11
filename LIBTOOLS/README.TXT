Release tag : LIBTOOLS-CNRM-4-8-a  on September 24 2009

Welcome to MESONH Libtools...

Branch : LIBTOOLS-CNRM-branch
TAG    : $Name$


Support for gfortran ver > 4.3 and NCL/NCAR 5.1.1
with :

export ARCH=LXgfortran

Documentation can be found in 'readme' directory.

LaTeX sources of documentation can be found in 'readme/LATEX'
directory : type make in the LATEX directory to build the 
postscript documentation file : tools.ps


How to compile ?

 cd lib
 export ARCH=LXgfortran
 make COMPRESS
 make NEWLFI
 make SURCOUCHE
 make MPIvide
 make RAD2
 make vis5d
 cd gribex_1302b
 export ARCH=linux
 export CNAME=_gfortran
 export A64=A64
 make
 
 cd ../../tools
 cd diachro
 export ARCH=LXgfortran
 export MNH_LIBTOOLS= absolute path for libtools directory
 make
# The executables are in the directory LXgfortran_64 conv2dia and LXgfortran_32 for the others


 cd ../fmmore
 make

# The executable is in the directory LXgfortran_64

 cd ../lfiz
 make

# The executables ate in the directory LXgfortran

