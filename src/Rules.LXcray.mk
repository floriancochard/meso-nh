#MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
#MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
#MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
#MNH_LIC for details. version 1.
##########################################################
#                                                        #
# Compiler Options                                       #
#                                                        #
##########################################################
#OBJDIR_PATH=/home/escj/azertyuiopqsdfghjklm/wxcvbn/azertyuiopqsdfghjklmwxcvbn
#
#OPT_BASE   =  -r8 -g -w -assume nosource_include -assume byterecl -fpe0 -ftz -fpic -traceback  -fp-model precise -switch fe_inline_all_arg_copy_inout
OPT_BASE   =  -sreal64 -hpic -em -ef
OPT_PERF0  =  -O0 -g
OPT_PERF2  =  -O2 -hflex_mp=intolerant -Ofp0 -hnoomp
OPT_PERF1  =  -O1 -hflex_mp=intolerant -Ofp0 -hnoomp -hcpu=istanbul -hfp0 -K trap=fp
# –hcpu=Istanbul –hfp0
#OPT_CHECK  =  -CB -ftrapuv 
OPT_CHECK  =  -Rbc
OPT_I8     =  -sdefault64
#
# Integer 4/8 option
#
MNH_INT   ?=4
LFI_RECL  ?=512
#
OPT_BASE_I4       := $(OPT_BASE)
ifeq "$(MNH_INT)" "8"
#OPT_BASE         += $(OPT_I8)
OPT_BASE           = -sdefault64 -hpic -em -ef
LFI_INT           ?=8
MNH_MPI_RANK_KIND ?=8
else
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
CFLAGS   += -g
endif
ifeq "$(OPTLEVEL)" "O2PAR"
PAR= -homp
OPT       = $(OPT_BASE) $(OPT_PERF2) $(PAR)
OPT0      = $(OPT_BASE) $(OPT_PERF0) $(PAR)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF2) $(PAR)
endif
ifeq "$(OPTLEVEL)" "O2NOVEC"
OPT       = $(OPT_BASE) $(OPT_PERF2) -O vector0
OPT0      = $(OPT_BASE) $(OPT_PERF0) -O vector0
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF2) -O vector0
endif
ifeq "$(OPTLEVEL)" "O1"
OPT       = $(OPT_BASE) $(OPT_PERF1) 
OPT0      = $(OPT_BASE) $(OPT_PERF0) 
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF1) 
endif
#
#
FC = ftn -em -ef
FCFLAGS = -em -ef
CC=gcc
export FC CC FCFLAGS
F90 = ftn
F90FLAGS  =  $(OPT)
F77  = $(F90)
F77FLAGS  =  $(OPT) 
# -132
FX90 = $(F90)
FX90FLAGS =  $(OPT)
# -132 
#
#LDFLAGS    =  -Wl,-noinhibit-exec  -Wl,-warn-once $(PAR)
LDFLAGS    =   -Wl,-warn-once $(PAR) $(OPT_BASE)
#
# preprocessing flags 
#
CPP = cpp -P -traditional -Wcomment
#
CPPFLAGS_SURFEX    =
CPPFLAGS_SURCOUCHE = -DMNH_MPI_DOUBLE_PRECISION -DMNH_LINUX -DDEV_NULL -DMNH_MPI_RANK_KIND=$(MNH_MPI_RANK_KIND) 
CPPFLAGS_RAD       =
CPPFLAGS_NEWLFI    = -DSWAPIO -DLINUX -DLFI_INT=${LFI_INT} -DLFI_RECL=${LFI_RECL}
CPPFLAGS_MNH       = -DMNH -DSFX_MNH 
ifdef VER_GA
CPPFLAGS_SURCOUCHE += -DMNH_GA
INC                += -I${GA_ROOT}/include
LIBS               += -L${GA_ROOT}/lib -larmci -lga -lgfortran
endif
#
# Gribex flags
#
TARGET_GRIBEX=linux
CNAME_GRIBEX=_gfortran
#
# GRIB_API
#
GRIBAPI_CONF="FCFLAGS= -em -ef "

#
# LIBTOOLS flags
#
#if MNH_TOOLS exists => compile the tools
MNH_TOOLS = yes
#
## COMPRESS flag
#
#if MNH_COMPRESS exists => compile the COMPRESS library (for LFI files)
MNH_COMPRESS=yes
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
#DIR_SURCOUCHE   += ARCH_SRC/bug_surcouche
#DIR_MNH         += ARCH_SRC/bug_mnh
#DIR_RAD         += ARCH_SRC/bug_rad
#DIR_SURFEX      += ARCH_SRC/surfex
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
# Juan & Maud 20/03/2008 --> Ifort 10.1.008 Bug O2 optimization
#OPT_PERF1  =  -O1
OBJS_O1= spll_schu.o spll_ps2str.o spll_p_abs.o spll_ini_one_way_n.o spll_urban_solar_abs.o
$(OBJS_O1) : OPT = $(OPT_BASE) $(OPT_PERF1)
OBJS_O0= spll_mode_gridproj.o spll_ini_dynamics.o spll_sunpos_n.o spll_ground_param_n.o

$(OBJS_O0) : OPT = $(OPT_BASE) $(OPT_PERF0)

ifneq "$(findstring 8,$(LFI_INT))" ""
OBJS_I8=spll_NEWLFI_ALL.o
$(OBJS_I8) : OPT = $(OPT_BASE) $(OPT_PERF2) $(OPT_I8)
endif

ifeq "$(MNH_INT)" "8"
OBJS_I4=spll_modd_netcdf.o
$(OBJS_I4) : OPT = $(OPT_BASE_I4)
endif

#mpi.mod : 
#	ln -sf /opt/cray/mpt/5.6.3/gni/mpich2-cray/74/include/MPI.mod $(OBJDIR_MASTER)/mpi.mod

