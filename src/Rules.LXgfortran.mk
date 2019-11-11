#MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
#MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
#MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
#MNH_LIC for details. version 1.
##########################################################
#                                                        #
# Compiler Options                                       #
#                                                        #
##########################################################
#
#   Gfortran version
GFV=$(shell  gfortran --version | grep -E -m1 -o ' [[:digit:]\.]{2,}( |$$)'  | sed 's/\.//g' )
#
#OBJDIR_PATH=/home/escj/azertyuiopqsdfghjklm/wxcvbn/azertyuiopqsdfghjklmwxcvbn
#
OPT_BASE  =  -g -fno-second-underscore -fPIC  -ffpe-trap=overflow,zero,invalid  -fbacktrace
#
OPT_PERF0 = -O0
OPT_PERF2 = -O2
OPT_CHECK = -fbounds-check -finit-real=nan
OPT_I8    = -fdefault-integer-8 
OPT_R8    = -fdefault-real-8 -fdefault-double-8
#
#
# Real/Integer 4/8 option
#
MNH_REAL  ?=8
MNH_INT   ?=4
LFI_RECL  ?=512
#
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
ifeq "$(OPTLEVEL)" "DEBUG"
OPT       = $(OPT_BASE) $(OPT_PERF0) $(OPT_CHECK)
OPT0      = $(OPT_BASE) $(OPT_PERF0) $(OPT_CHECK)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF0)
CFLAGS    += -g -O0
endif
#
#  
CC = gcc
CXX = g++
FC = gfortran 
ifeq "$(VER_MPI)" "MPIAUTO"
F90 = mpif90
CPPFLAGS_SURCOUCHE += -DMNH_USE_MPI_STATUSES_IGNORE
else         
F90 = gfortran
endif
#
F90FLAGS      =  $(OPT) 
F77 = $(F90)
F77FLAGS      =  $(OPT) 
FX90 = $(F90)
FX90FLAGS     =  $(OPT) 
#
#LDFLAGS   =   -Wl,-warn-once
#
# preprocessing flags 
#
CPP = cpp -P -traditional -Wcomment
#
CPPFLAGS_SURFEX    =
CPPFLAGS_SURCOUCHE += -DMNH_LINUX -DDEV_NULL  -DMNH_MPI_RANK_KIND=$(MNH_MPI_RANK_KIND)
CPPFLAGS_RAD       =
CPPFLAGS_NEWLFI    = -DSWAPIO -DLINUX -DLFI_INT=${LFI_INT} -DLFI_RECL=${LFI_RECL}
CPPFLAGS_MNH       = -DMNH -DSFX_MNH
ifdef VER_GA
CPPFLAGS_SURCOUCHE += -DMNH_GA
INC                += -I${GA_ROOT}/include
LIBS               += -L${GA_ROOT}/lib -larmci -lga
endif
#
# Gribex flags
#
TARGET_GRIBEX=linux
CNAME_GRIBEX=_gfortran
#
# Netcdf/HDF5 flags
#
HDF_CONF= CFLAGS=-std=c99
HDF_OPT ?= -fPIC
NETCDF_OPT ?= -fPIC
#
## LIBTOOLS flags
#
#if MNH_TOOLS exists => compile the tools if gfortran >= 5.X
ifeq "$(MNH_INT)" "4"
ifeq ($(shell test $(GFV) -ge 500 ; echo $$?),0)
MNH_TOOLS=yes
endif
endif
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

ifeq ($(shell test $(GFV) -le 482 ; echo $$?),0)
ifneq "$(OPTLEVEL)" "DEBUG"
OBJS_O0= spll_lima_phillips_integ.o
$(OBJS_O0) : OPT = $(OPT_BASE) $(OPT_PERF0)
endif
endif

