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
OPT_BASE  = -r8 -g -fno-second-underscore -fpic
OPT_PERF0 = -O0
OPT_PERF2 = -O2
OPT_CHECK = -fbounds-check -ftrace=full -freal=nan -fpointer=invalid
OPT_I8     =  -i8
#
# Integer 4/8 option
#
MNH_INT   ?=4
LFI_RECL  ?=512
#
OPT_BASE_I4       := $(OPT_BASE)
ifeq "$(MNH_INT)" "8"
OPT_BASE          += $(OPT_I8)
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
OBJS_NOCB += spll_mode_cover_301_573.o
endif
#
#
FC = g95
ifeq "$(VER_MPI)" "MPIAUTO"
F90 = mpif90   
else     
F90 = g95
endif
#
F90FLAGS      =  $(OPT) 
F77 = $(F90)
#F77FLAGS      =  $(OPT) -ffixed-line-length-132
F77FLAGS      =  $(OPT) 
FX90 = $(F90)
#FX90FLAGS     =  $(OPT) -ffixed-line-length-132
FX90FLAGS     =  $(OPT) 
#
LDFLAGS   =  -Wl,-warn-once
#
# preprocessing flags 
#
CPP = cpp -P -traditional -Wcomment
#
LFI_INT  ?=4
LFI_RECL ?=512
CPPFLAGS_SURFEX    =
CPPFLAGS_SURCOUCHE = -DMNH_MPI_DOUBLE_PRECISION -DMNH_LINUX -DMNH_MPI_BSEND -DDEV_NULL -DMNH_MPI_RANK_KIND=$(MNH_MPI_RANK_KIND)
CPPFLAGS_RAD       =
CPPFLAGS_NEWLFI    = -DSWAPIO -DLINUX  -DLFI_INT=${LFI_INT} -DLFI_RECL=${LFI_RECL}
CPPFLAGS_MNH       = -DAINT=INT -DAMOD=MOD -DMNH -DSFX_MNH
#
# Gribex flags
#
TARGET_GRIBEX=linux
CNAME_GRIBEX=g95
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

ifneq "$(findstring 8,$(LFI_INT))" ""
OBJS_I8=spll_NEWLFI_ALL.o
$(OBJS_I8) : OPT = $(OPT_BASE) $(OPT_PERF2) $(OPT_I8)
endif

ifeq "$(MNH_INT)" "8"
OBJS_I4=spll_modd_netcdf.o
$(OBJS_I4) : OPT = $(OPT_BASE_I4)
endif

