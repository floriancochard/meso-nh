#MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
#MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
#MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
#MNH_LIC for details. version 1.
##########################################################
#                                                        #
# Compiler Options                                       #
#                                                        #
##########################################################
#OBJDIR_PATH=${WORKDIR}
# -qsigtrap -qfloat=nans
# -qflttrap=enable:overflow:zerodivide:invalid
# -qextname 
#OPT_BASE  = -qflttrap=enable:overflow:zerodivide:invalid \
#            -qfloat=nans -qarch=450 -qmoddir=$(OBJDIR) \
#            -qautodbl=dbl4 -qzerosize -g -qfullpath -qspillsize=32648 \
#            -qinitauto=0 -qdpc=e -qmaxmem=-1

#OPT_BASE  = -qmoddir=$(OBJDIR) -qautodbl=dbl4 -qzerosize  
#OPT_BASE  = -g -qautodbl=dbl4 -qzerosize -qextname=flush -qnohot -qnoescape 
OPT_BASE  = -g -qautodbl=dbl4 -qdpc=e -qzerosize -qextname -qnohot -qnoescape -qarch=auto -qtune=auto -qspillsize=32648 -qfloat=nomaf
OPT_BASE_R8  = -g -qrealsize=8 -qdpc=e -qzerosize -qextname -qnohot -qnoescape -qarch=auto -qtune=auto -qspillsize=32648 -qfloat=nomaf
#            -qsimd=noauto -qfloat=nomaf:nofold:norsqrt
 
OPT_PERF0   = -O0 -qnooptimize -qkeepparm -qfullpath 
OPT_PERF2   = -O2 -qmaxmem=-1
OPT_CHECK = -C 
OPT_NONAN =  -qsigtrap -qflttrap=qpxstore:overflow:zerodivide:invalid:enable -qfloat=nans
# -qflttrap=qpxstore
OPT_I8      = -qintsize=8 -qxlf77=intarg
OPT_I4      = -qintsize=4 -qxlf77=intarg
#
# Real/Integer 4/8 option
#
MNH_REAL  ?=8
MNH_INT   ?=4
LFI_RECL  ?=512
#
OPT_BASE_I4       := $(OPT_BASE) $(OPT_I4)
ifeq "$(MNH_INT)" "8"
OPT_BASE         += $(OPT_I8)
LFI_INT           ?=8
MNH_MPI_RANK_KIND ?=8
else
OPT_BASE         += $(OPT_I4)
MNH_MPI_RANK_KIND ?=4
LFI_INT           ?=4
endif
#
OPT       = $(OPT_BASE) $(OPT_PERF2) 
OPT0      = $(OPT_BASE) $(OPT_PERF0) 
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF2)
#
ifeq "$(OPTLEVEL)" "DEBUG"
OPT       = $(OPT_BASE) $(OPT_PERF0) $(OPT_CHECK)
OPT0      = $(OPT_BASE) $(OPT_PERF0) $(OPT_CHECK)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF0)
#LIBS     += -L/bglocal/prod/TotalView/8.10.0-0/linux-power/lib/ -ltvheap_bluegene_p
endif
#
ifeq "$(OPTLEVEL)" "O2"
OPT       = $(OPT_BASE) $(OPT_PERF2) $(OPT_NONAN)
OPT0      = $(OPT_BASE) $(OPT_PERF0) $(OPT_NONAN)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF2) $(OPT_NONAN)
endif
#
ifeq "$(OPTLEVEL)" "O2NAN"
OPT       = $(OPT_BASE) $(OPT_PERF2)
OPT0      = $(OPT_BASE) $(OPT_PERF0)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF2)
endif

# Problem xlf , no option to promote only not kinded real
# use -qrealsize=8 inplace
# OBJSR8 = spll_modi_gather_and_write_mpi.o gather_and_write_mpi.o
#$(OBJSR8) :  OPT =  $(OPT_BASE_R8) $(OPT_PERF2) 
#
#
# File with problem OPTLEVEL > 2
#
# Compilation problems
OBJS2 = spll_mode_dustopt.o spll_mode_saltopt.o spll_write_lfifm1_for_diag.o \
        spll_mode_thermo.o spll_asymtx.o
#
# Runtime Problems
OBJS2 += spll_read_sso_canopy_n.o
#
ifeq "$(OPTLEVEL)" "O3"
OPT_PERF3    = -O3 -qstrict -qmaxmem=-1
OPT       = $(OPT_BASE) $(OPT_PERF3) 
OPT0      = $(OPT_BASE) $(OPT_PERF0) 
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF3)
$(OBJS2) :  OPT =  $(OPT_BASE) $(OPT_PERF2) 
endif
#            
ifeq "$(OPTLEVEL)" "O3SMP"
OPT_PERF3 = -O3 -qsmp -qstrict -qmaxmem=-1
OPT       = $(OPT_BASE) $(OPT_PERF3) 
OPT0      = $(OPT_BASE) $(OPT_PERF0)    
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF3)
$(OBJS2) :  OPT =  $(OPT_BASE)  $(OPT_PERF2) 
endif
#   
ifeq "$(OPTLEVEL)" "O2SMP"
OPT_PERF2SMP = -O2 -qsmp -qstrict -qmaxmem=-1
OPT       = $(OPT_BASE) $(OPT_PERF2SMP)
OPT0      = $(OPT_BASE) $(OPT_PERF0)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF2SMP)
$(OBJS2) :  OPT =  $(OPT_BASE) $(OPT_PERF2) 
endif         
#
ifeq "$(OPTLEVEL)" "O4"
OPT_PERF4    = -O4 
OPT       = $(OPT_BASE) $(OPT_PERF4) 
OPT0      = $(OPT_BASE) $(OPT_PERF0) 
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF4)
endif
#
#
F90 = mpixlf95_r
F90FLAGS =       $(OPT) -qfree=f90 -qsuffix=f=f90 
F77 = $(F90)
F77FLAGS      =  $(OPT) -qfixed
FX90 = $(F90)
FX90FLAGS     =  $(OPT) -qfixed
#
# compiler & flags for compilation of grib_api
# for reproductibility need : -qfloat=nomaf
#
FC = mpixlf95_r
FCFLAGS = -qfloat=nomaf
#CC = mpixlc_r
#CFLAGS =  $(FCFLAGS)
CC = powerpc64-bgq-linux-gcc
export CC FCFLAGS CFLAGS
#
LDFLAGS   =  $(OPT) -Wl,--relax   
AR = /bgsys/drivers/ppcfloor/gnu-linux/powerpc64-bgq-linux/bin/ar
#
# preprocessing flags 
#
CPP = cpp -P -traditional -Wcomment

#
CPPFLAGS_SURFEX    =
#CPPFLAGS_SURCOUCHE = -DMNH_MPI_DOUBLE_PRECISION -DMNH_LINUX -DMNH_SP4 -DMNH_MPI_ISEND -DMNH_MPI_RANK_KIND=$(MNH_MPI_RANK_KIND)
CPPFLAGS_SURCOUCHE = -DMNH_MPI_DOUBLE_PRECISION -DMNH_LINUX -DMNH_SP4 -DMNH_MPI_BSEND -DMNH_MPI_RANK_KIND=$(MNH_MPI_RANK_KIND) -DSNGL=REAL
CPPFLAGS_RAD       =
CPPFLAGS_NEWLFI    = -DLINUX  -DLFI_INT=${LFI_INT} -DLFI_RECL=${LFI_RECL}
CPPFLAGS_MNH       = -DAMAX1=MAX -DMNH -DSFX_MNH
#
# Rules for GA = Global Array
#
ifdef VER_GA
CPPFLAGS_SURCOUCHE += -DMNH_GA
INC                += -I${GA_ROOT}/include
LIBS               += -L${GA_ROOT}/lib -larmci -lga -lgfortran
endif
#
# Gribex flags
#
#TARGET_GRIBEX=rs6000
TARGET_GRIBEX=ibm_power4
CNAME_GRIBEX=""
#A64=A64
# Gribapi flags
GRIBAPI_CONF= --host=powerpc64-bgq-linux 
#
# LIBTOOLS flags
#
#if MNH_TOOLS exists => compile the tools
#MNH_TOOLS = no
#
## COMPRESS flag
#
#if MNH_COMPRESS exists => compile the COMPRESS library (for LFI files)
#MNH_COMPRESS=no
#
## S4PY flag
#
#if MNH_S4PY exists => compile the libs4py library (for epygram)
#MNH_S4PY=no
#
##########################################################
#                                                        #
# Source of MESONH PACKAGE  Distribution                 #
#                                                        #
##########################################################
#DIR_SURFEX      += ARCH_SRC/surfex 
#DIR_MNH      += ARCH_SRC/bug_mnh_BG
#
include Makefile.MESONH.mk
#
##########################################################
#                                                        #
# extra VPATH, Compilation flag modification             #
#         systeme module , etc ...                       #
#         external precompiled module librairie          #
#         etc ...                                        #
#                                                        #
##########################################################
OBJS_NOCB += spll_prep_ideal_case.o spll_mesonh.o
$(OBJS_NOCB) : OPT = $(OPT_NOCB)
#
#IGNORE_OBJS += spll_abort.o spll_ch_make_lookup.o \
#spll_compute_ver_grid.o spll_convlfi.o spll_diag.o spll_example_fwd.o spll_latlon_to_xy.o \
#spll_prep_nest_pgd.o spll_prep_pgd.o spll_prep_real_case.o \
#spll_prep_surfex.o spll_rad1driv.o spll_rttov_ascii2bin_coef.o spll_rttovcld_testad.o spll_rttovcld_test.o \
#spll_rttovscatt_test.o spll_spawning.o spll_test_2_coef.o spll_test_coef.o spll_test_errorhandling.o \
#spll_test_q2v.o spll_xy_to_latlon.o spll_zoom_pgd.o 

ifneq "$(findstring 8,$(LFI_INT))" ""
OBJS_I8=spll_NEWLFI_ALL.o
$(OBJS_I8) : OPT = $(OPT_BASE) $(OPT_PERF2) $(OPT_I8)
endif

ifeq "$(MNH_INT)" "8"
OBJS_I4=spll_modd_netcdf.o
$(OBJS_I4) : OPT = $(OPT_BASE_I4)
endif

