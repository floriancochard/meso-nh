#MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
#MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
#MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
#MNH_LIC for details. version 1.
#-----------------------------------------------------------------
#--------------- special set of characters for CVS information
#-----------------------------------------------------------------
# $Source$
# $Name$ 
# $Revision$ 
# $Date$
#-----------------------------------------------------------------
#-----------------------------------------------------------------

VERSION=-7.4.2    
VPATH = src:$(DIR_OBJ)
#######################################
LIBSURC = libsurcouche$(VERSION).a

DIR_OBJ = ./$(ARCH)$(VERSION)
INC = -I$(DIR_OBJ)


DIR_CONF:=$(shell pwd|sed -e 's/lib\/.*/conf/')

OBJS =  modd_mpif.o mode_ll.o modi_fm.o \
	modi_fmread.o modi_fmwrit.o modi_io.o modi_scatter.o \
	extern_userio.o fmread_ll.o fmwrit_ll.o modd_argslist_ll.o \
	modd_io.o modd_type_date.o mode_fm.o mode_io.o \
	modi_adddelfield2_ll.o modi_addndfield_ll.o modi_delndfield_ll.o modi_get_ll.o \
	modi_init_ll.o modi_location_ll.o modi_nest_ll.o modi_remap_ll.o \
	modi_sum_ll.o modi_update_ll.o extern_usersurc_ll.o \
	fmreadwrit.o modd_errcodes.o modd_parameters_ll.o modd_structure2_ll.o \
	mode_allocbuff.o mode_argslist2_ll.o mode_argslist_ll.o mode_boundaries_ll.o \
	mode_distriblb.o mode_exchange2_ll.o mode_exchange_ll.o mode_fd.o \
	mode_gather.o mode_init_ll.o mode_lb_ll.o mode_ls_ll.o \
	mode_nest_ll.o mode_scatter.o mode_sum2_ll.o mode_sum_ll.o \
	mode_tools_ll.o modi_gather.o modd_dim_ll.o modd_structure_ll.o \
	modd_var_ll.o mode_construct_ll.o mode_splitting_ll.o 

include $(DIR_CONF)/config.$(ARCH)
include Rules.$(ARCH)


%.o:%.f90
	$(CPP) $(INC) $(CPPFLAGS) $< >  $(DIR_OBJ)/cpp_$(*F).f90
	$(F90) $(INC) -c $(F90FLAGS) $(DIR_OBJ)/cpp_$(*F).f90 -o $(DIR_OBJ)/$(*F).o
	-@mv  *.mod $(DIR_OBJ)/. 2> /dev/null || echo pas de module dans $*.f90

%.o:%.c
	$(CC) $(INC) $(CFLAGS) $(CPPFLAGS) -c $< -o $(DIR_OBJ)/$(*F).o

$(LIBSURC) : $(DIR_OBJ)/.dummy $(OBJS)
	cd $(DIR_OBJ);$(AR) crv $@ $(OBJS)

$(DIR_OBJ)/.dummy :
	mkdir -p $(DIR_OBJ)
	@touch $(DIR_OBJ)/.dummy

modd_argslist_ll.o: modd_argslist_ll.f90 modd_structure2_ll.o 
modd_io.o: modd_io.f90 
modd_type_date.o: modd_type_date.f90 
modd_errcodes.o: modd_errcodes.f90 
modd_parameters_ll.o: modd_parameters_ll.f90 
modd_structure2_ll.o: modd_structure2_ll.f90 
modd_dim_ll.o: modd_dim_ll.f90 
modd_structure_ll.o: modd_structure_ll.f90 modd_argslist_ll.o \
	modd_parameters_ll.o 
modd_var_ll.o: modd_var_ll.f90 modd_structure_ll.o 

modi_fm.o: modi_fm.f90 
modi_fmread.o: modi_fmread.f90 modd_type_date.o 
modi_fmwrit.o: modi_fmwrit.f90 modd_type_date.o 
modi_io.o: modi_io.f90 modd_io.o 
modi_scatter.o: modi_scatter.f90 
modi_adddelfield2_ll.o: modi_adddelfield2_ll.f90 modd_argslist_ll.o 
modi_addndfield_ll.o: modi_addndfield_ll.f90 modd_argslist_ll.o 
modi_delndfield_ll.o: modi_delndfield_ll.f90 modd_argslist_ll.o 
modi_get_ll.o: modi_get_ll.f90 
modi_init_ll.o: modi_init_ll.f90 
modi_location_ll.o: modi_location_ll.f90 
modi_nest_ll.o: modi_nest_ll.f90 
modi_remap_ll.o: modi_remap_ll.f90 
modi_sum_ll.o: modi_sum_ll.f90 
modi_update_ll.o: modi_update_ll.f90 modd_argslist_ll.o 
modi_gather.o: modi_gather.f90 

extern_userio.o: extern_userio.f90 fmread_ll.o \
	fmwrit_ll.o modd_io.o modd_type_date.o \
	mode_fm.o mode_gather.o mode_io.o 
extern_usersurc_ll.o: extern_usersurc_ll.f90 modd_argslist_ll.o \
	mode_argslist2_ll.o mode_argslist_ll.o mode_boundaries_ll.o \
	mode_exchange2_ll.o mode_exchange_ll.o mode_init_ll.o \
	mode_lb_ll.o mode_ls_ll.o mode_nest_ll.o \
	mode_sum2_ll.o mode_sum_ll.o mode_tools_ll.o \
	mode_scatter.o


mode_ll.o: mode_ll.f90 modd_argslist_ll.o \
	modi_adddelfield2_ll.o modi_addndfield_ll.o modi_delndfield_ll.o \
	modi_get_ll.o modi_init_ll.o modi_location_ll.o \
	modi_nest_ll.o modi_remap_ll.o modi_sum_ll.o \
	modi_update_ll.o 


mode_io.o: mode_io.f90 modd_errcodes.o \
	modd_io.o mode_fd.o 
mode_fd.o: mode_fd.f90 modd_io.o 
mode_fm.o: mode_fm.f90 modd_errcodes.o \
	modd_io.o mode_fd.o mode_io.o 
fmread_ll.o: fmread_ll.f90 fmreadwrit.o \
	modd_io.o modd_parameters_ll.o modd_type_date.o \
	mode_allocbuff.o mode_distriblb.o mode_fd.o \
	mode_fm.o mode_scatter.o mode_tools_ll.o 
fmwrit_ll.o: fmwrit_ll.f90 fmreadwrit.o \
	modd_io.o modd_parameters_ll.o modd_type_date.o \
	mode_allocbuff.o mode_distriblb.o mode_fd.o \
	mode_fm.o mode_gather.o mode_tools_ll.o 
fmreadwrit.o: fmreadwrit.f90 
mode_allocbuff.o: mode_allocbuff.f90 modd_io.o \
	modd_parameters_ll.o mode_tools_ll.o 
mode_gather.o: mode_gather.f90 modd_io.o \
	modd_structure_ll.o modd_var_ll.o mode_tools_ll.o 
mode_scatter.o: mode_scatter.f90 modd_io.o \
	modd_structure_ll.o modd_var_ll.o 

mode_init_ll.o: mode_init_ll.f90 modd_dim_ll.o \
	modd_parameters_ll.o modd_structure_ll.o modd_var_ll.o \
	mode_construct_ll.o mode_nest_ll.o mode_splitting_ll.o 
mode_argslist2_ll.o: mode_argslist2_ll.f90 modd_argslist_ll.o 
mode_argslist_ll.o: mode_argslist_ll.f90 modd_argslist_ll.o 
mode_boundaries_ll.o: mode_boundaries_ll.f90 modd_argslist_ll.o \
	modd_var_ll.o mode_exchange_ll.o 
mode_distriblb.o: mode_distriblb.f90 modd_parameters_ll.o \
	modd_structure_ll.o modd_var_ll.o mode_tools_ll.o 
mode_exchange2_ll.o: mode_exchange2_ll.f90 modd_argslist_ll.o \
	modd_structure_ll.o modd_var_ll.o mode_exchange_ll.o 
mode_exchange_ll.o: mode_exchange_ll.f90 modd_argslist_ll.o \
	modd_dim_ll.o modd_parameters_ll.o modd_structure_ll.o \
	modd_var_ll.o mode_argslist_ll.o mode_tools_ll.o 
mode_lb_ll.o: mode_lb_ll.f90 modd_argslist_ll.o \
	modd_dim_ll.o modd_parameters_ll.o modd_structure_ll.o \
	modd_var_ll.o mode_argslist_ll.o mode_construct_ll.o \
	mode_exchange_ll.o mode_nest_ll.o mode_splitting_ll.o \
	mode_tools_ll.o 
mode_ls_ll.o: mode_ls_ll.f90 modd_structure_ll.o \
	modd_var_ll.o mode_argslist_ll.o mode_construct_ll.o \
	mode_exchange_ll.o mode_nest_ll.o 
mode_nest_ll.o: mode_nest_ll.f90 modd_dim_ll.o \
	modd_parameters_ll.o modd_structure_ll.o modd_var_ll.o \
	mode_construct_ll.o mode_splitting_ll.o mode_tools_ll.o 
mode_sum2_ll.o: mode_sum2_ll.f90 modd_structure_ll.o \
	modd_var_ll.o 
mode_sum_ll.o: mode_sum_ll.f90 modd_parameters_ll.o \
	modd_structure_ll.o modd_var_ll.o mode_tools_ll.o 
mode_tools_ll.o: mode_tools_ll.f90 modd_dim_ll.o \
	modd_parameters_ll.o modd_structure_ll.o modd_var_ll.o 
mode_construct_ll.o: mode_construct_ll.f90 modd_dim_ll.o \
	modd_parameters_ll.o modd_structure_ll.o modd_var_ll.o \
	mode_argslist_ll.o mode_tools_ll.o 
mode_splitting_ll.o: mode_splitting_ll.f90 modd_parameters_ll.o \
	modd_structure_ll.o 

clean:
	(if [ -d $(DIR_OBJ) ] ; then cd $(DIR_OBJ); rm -f cpp_*.f90 *.o ; fi)

distclean:
	rm -rf $(DIR_OBJ)
