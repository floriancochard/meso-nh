#MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
#MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
#MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
#MNH_LIC for details. version 1.
VPATH = src:$(DIR_OBJ)
#######################################
LIBLFI = libNEWLFI_ALL.a
DIR_OBJ = ./$(ARCH)

DIR_CONF:=$(shell pwd|sed -e 's/lib\/.*/conf/')

include $(DIR_CONF)/config.$(ARCH)
include Rules.$(ARCH)

%.o : %.f
	$(CPP) $(CPPFLAGS) $< > $(DIR_OBJ)/cpp_$(*F).f
	$(F77) -c $(F77FLAGS) -o $(DIR_OBJ)/$(*F).o $(DIR_OBJ)/cpp_$(*F).f

%.o : %.f90
	$(CPP) $(INC) $(CPPFLAGS) $< > $(DIR_OBJ)/cpp_$(*F).f90
	$(F90) $(INC) -c $(F90FLAGS) -o $(DIR_OBJ)/$(*F).o $(DIR_OBJ)/cpp_$(*F).f90

%.o : %.c
	$(CC) $(INC) $(CFLAGS) $(CPPFLAGS) -c -o $(DIR_OBJ)/$(*F).o  $<



$(LIBLFI) : $(DIR_OBJ)/.dummy $(OBJS)
	cd $(DIR_OBJ); $(AR) crv $@ $(OBJS)

$(DIR_OBJ)/.dummy :
	mkdir -p $(DIR_OBJ)
	@touch $(DIR_OBJ)/.dummy

clean:
	(if [ -d $(DIR_OBJ) ] ; then cd $(DIR_OBJ); rm -f cpp_*.f* *.o ; fi)

distclean:
	rm -rf $(DIR_OBJ)
