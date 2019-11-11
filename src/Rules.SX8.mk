#MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
#MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
#MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
#MNH_LIC for details. version 1.
##########################################################
#                                                        #
# Compiler Options                                       #
#                                                        #
##########################################################
ifneq "$(findstring brodie,$(shell uname -n))" ""
OBJDIR_PATH=${workdir}
endif
#
#OPT_LIST = -R5 -Wf"-pvctl fullmsg -O infomsg -L fmtlist" -Wf"-L transform" -Wf"-L summary"
#OPT_BASE = -dw -Wf,'-A idbl4' -to $(OBJDIR) -P stack -Wf"-pvctl noassume loopcnt=1500000 vwork=stack " $(OPT_LIST)
OPT_BASE = -dw -Wf,-Nesc,'-A idbl4' -P stack -Wf"-pvctl noassume loopcnt=1500000 vwork=stack " $(OPT_LIST)
OPT_PERF0 = -C debug 
OPT_PERF2 = -C vsafe 
OPT_PERF3 = -C vopt
OPT_PERF4 = -C hopt
OPT_CHECK = -ePR -Wf"-init stack=nan" -Wf"-init heap=nan"
OPT_I8    = -ew
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
ifeq "$(OPTLEVEL)" "O2"
OPT       = $(OPT_BASE) $(OPT_PERF2) 
OPT0      = $(OPT_BASE) $(OPT_PERF0) 
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF2)
endif
#
ifeq "$(OPTLEVEL)" "O2PROF"
OPT       = $(OPT_BASE) $(OPT_PERF2) -ftrace
OPT0      = $(OPT_BASE) $(OPT_PERF0) -ftrace
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF2) -ftrace
endif
#
ifeq "$(OPTLEVEL)" "DEBUG"
OPT       = $(OPT_BASE) $(OPT_PERF0) $(OPT_CHECK)
OPT0      = $(OPT_BASE) $(OPT_PERF0) $(OPT_CHECK)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF0)
#LDFLAGS   +=  -Wl"-f nan" 
endif
#
ifeq "$(OPTLEVEL)" "O4"
OPT       = $(OPT_BASE) $(OPT_PERF4) 
OPT0      = $(OPT_BASE) $(OPT_PERF0)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF4)
# routines with probleme with inlining
OBJS_NI += spll_ecmwf_radiation_vers2.o spll_slow_terms.o
$(OBJS_NI) :  OPT =  $(OPT_BASE)  $(OPT_PERF4) -Wf,-Ni
# problem for NPROC >= 4 & reprod_sum
OBJS2 += spll_mode_construct_ll.o \
spll_mode_repro_sum.o spll_dotprod.o \
spll_ini_radiations_ecmwf.o spll_mass_leak.o spll_mode_sum_ll.o spll_modi_sum_ll.o \
spll_p_abs.o spll_reduce_sum_0dd_ll.o spll_reduce_sum_1dd_ll.o spll_relaxdef.o \
spll_set_ref.o spll_surf_solar_sum.o spll_test_double_double.o spll_trid.o
$(OBJS2) :  OPT =  $(OPT_BASE)  $(OPT_PERF2) 
endif

F90 = sxmpif90 -clear
F90FLAGS  =    $(OPT) 
F77 = $(F90)     
F77FLAGS  =    -f0 $(OPT)
FX90 = $(F90)        
FX90FLAGS =    -f0 $(OPT)
# 
#LDFLAGS   +=  $(OPT) -Wl"-h lib_cyclic -h nodefs" 
LDFLAGS   +=  $(OPT) -Wl"-h lib_cyclic " 
#
# preprocessing flags 
#
CPP = cpp -P -traditional -Wcomment
AR=sxar
#
CPPFLAGS_SURFEX    =
CPPFLAGS_SURCOUCHE = -DMNH_MPI_DOUBLE_PRECISION -DMNH_SX5 -DMNH_MPI_BSEND -DMNH_MPI_RANK_KIND=$(MNH_MPI_RANK_KIND)
CPPFLAGS_RAD       =
CPPFLAGS_NEWLFI    = -DMNH_SX5 -DLFI_INT=${LFI_INT} -DLFI_RECL=${LFI_RECL}
CPPFLAGS_MNH       = -DMNH -DSFX_MNH
#
# Gribex flags
#
#ARCH_GRIBEX=NEC
TARGET_GRIBEX=NEC
CNAME_GRIBEX=sxmpif90
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
#
OBJS0 += spll_ch_jac.o spll_ch_terms.o  \
         spll_fm_writ_ll.o spll_NEWLFI_ALL.o spll_fm_read_ll.o
$(OBJS0) : OPT_CHECK = 

ifneq "$(findstring 8,$(LFI_INT))" ""
OBJS_I8=spll_NEWLFI_ALL.o
$(OBJS_I8) : OPT = $(OPT_BASE) $(OPT_PERF2) $(OPT_I8)
endif

ifeq "$(MNH_INT)" "8"
OBJS_I4=spll_modd_netcdf.o
$(OBJS_I4) : OPT = $(OPT_BASE_I4)
endif

