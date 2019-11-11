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
#
OPT_BASE   =  -g -Ktrap=fp -Mbackslash
# -Munixlogical 
# -Mrecursive -mcmodel=medium
OPT_PERF0  =  -O0 -Kieee
OPT_PERF2  =  -O2 -Kieee
#OPT_CUDA  =  -O2 -Mcuda=keepgpu -ta=nvidia,cc20,cuda3.1,host,time -Minfo=accel,intensity,all,ccff  
#OPT_CUDA  =  -O3 -fast -ta=nvidia,cc20,cuda4.2,keepgpu,host -Minfo=accel,all,intensity,ccff 
OPT_CUDA  =  -O2 -Kieee -nofma -ta=host,nvidia,nofma,cc20,cc35,cuda5.5 -Minfo=ccff,all,intensity -Mprof=ccff 
#OPT_CUDA  =  -O2 -Kieee -ta=host,nvidia,cc20,cuda4.2 -Minfo=ccff,all,intensity

OPT_CHECK  =  -C 
OPT_PROF   =  -Mprof=time,ccff
OPT_I8     =  -i8
OPT_R8     =  -r8
#
IGNORE_OBJS += pgprof.o
#
# Real/integer 4/8 option
#
MNH_REAL  ?=8
MNH_INT   ?=4
LFI_RECL  ?=512
#
ifneq "$(MNH_REAL)" "4"
OPT_BASE           += $(OPT_R8)
CPPFLAGS_SURCOUCHE += -DMNH_MPI_DOUBLE_PRECISION
endif
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
#
OPT       = $(OPT_BASE) $(OPT_PERF2)
OPT0      = $(OPT_BASE) $(OPT_PERF0)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF2)
#
ifeq "$(OPTLEVEL)" "O2PROF"
OPT       = $(OPT_BASE) $(OPT_PERF2) $(OPT_PROF)
OPT0      = $(OPT_BASE) $(OPT_PERF0) $(OPT_PROF)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF2) $(OPT_PROF)
endif
ifeq "$(OPTLEVEL)" "DEBUG"
OPT       = $(OPT_BASE) $(OPT_PERF0) $(OPT_CHECK)
OPT0      = $(OPT_BASE) $(OPT_PERF0) $(OPT_CHECK)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF0)
endif

ifeq "$(OPTLEVEL)" "CUDA"
OPT       = $(OPT_BASE) $(OPT_CUDA) 
OPT0      = $(OPT_BASE) $(OPT_CUDA) $(OPT_PERF0)
OPT_NOCB  = $(OPT_BASE) $(OPT_CUDA) 
endif

ifeq "$(OPTLEVEL)" "CUDA_DB"
OPT_CUDA  =  -O0 -Kieee -ta=host,nvidia,cc20,cuda4.2 -Minfo=ccff,all,intensity
OPT       = $(OPT_BASE) $(OPT_CUDA) 
OPT0      = $(OPT_BASE) $(OPT_CUDA)
OPT_NOCB  = $(OPT_BASE) $(OPT_CUDA)
endif

#-Mcuda -ta=nvidia,host,time -Minfo=accel,intensity
#
CC = pgcc
FC = pgf90
ifeq "$(VER_MPI)" "MPIAUTO"
F90 = mpif90
else
F90 = pgf90
endif
#
F77FLAGS  =  $(OPT)
F77 = $(F90)
F90FLAGS  =  $(OPT)
FX90 = $(F90)
FX90FLAGS =  $(OPT)
#
LDFLAGS    =   -Wl,-warn-once $(OPT)
#
# preprocessing flags 
#
CPP = cpp -P -traditional -Wcomment
#

CPPFLAGS_SURFEX    =
CPPFLAGS_SURCOUCHE += -DMNH_LINUX -DMNH_MPI_RANK_KIND=$(MNH_MPI_RANK_KIND)
CPPFLAGS_RAD       =
CPPFLAGS_NEWLFI    = -DSWAPIO -DLINUX -DLFI_INT=${LFI_INT} -DLFI_RECL=${LFI_RECL}
CPPFLAGS_MNH       = -DMNH -DMNH_PGI -DSFX_MNH

#
# Gribex flags
#
TARGET_GRIBEX=linux
CNAME_GRIBEX=_pgf77
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
#DIR_SURFEX   += ARCH_SRC/surfex.MNH-462

OBJS_NOCB +=  spll_isba.o
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
OPT_PERF1  =  -O1 -Kieee -g
OBJS_O1= spll_modd_isba_n.o spll_pack_isba_patch_n.o spll_mode_construct_ll.o \
         spll_init_surf_atm_n.o spll_mode_scatter_ll.o spll_convert_patch_teb.o \
         spll_define_mask_n.o spll_del1dfield_ll.o spll_mode_fm.o spll_mode_gather_ll.o \
	 spll_phys_param_n.o \
	 spll_convect_updraft.o spll_convect_updraft_shal.o
$(OBJS_O1) : OPT = $(OPT_BASE) $(OPT_PERF1)

#
#MODULE_SYSTEM = /opt/F95_42/lib/
#VPATH += $(MODULE_SYSTEM)
#

ifneq "$(findstring 8,$(LFI_INT))" ""
OBJS_I8=spll_NEWLFI_ALL.o
$(OBJS_I8) : OPT = $(OPT_BASE) $(OPT_PERF2) $(OPT_I8)
endif

ifeq "$(MNH_INT)" "8"
OBJS_I4=spll_modd_netcdf.o
$(OBJS_I4) : OPT = $(OPT_BASE_I4)
endif

