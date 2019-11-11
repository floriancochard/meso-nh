#MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
#MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
#MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
#MNH_LIC for details. version 1.
VPATH=$(DIR_OBJ)
#######################################
DIR_OBJ = ./$(ARCH)

DIR_CONF:=$(shell pwd|sed -e 's/lib\/.*/conf/')

include $(DIR_CONF)/config.$(ARCH)
include Rules.$(ARCH)

libmpivide.a : $(DIR_OBJ)/.dummy mpivide.o
	cd $(DIR_OBJ); $(AR) -cr $@ mpivide.o


mpivide.o : mpivide.c
	$(CC) -c $(CFLAGS) -DFUJI -I./include $< -o $(DIR_OBJ)/$(*F).o


$(DIR_OBJ)/.dummy:
	mkdir -p $(DIR_OBJ)
	@touch $(DIR_OBJ)/.dummy

clean:
	(if [ -d $(DIR_OBJ) ] ; then cd $(DIR_OBJ); rm -f cpp_*.f90 *.o ; fi)

distclean:
	rm -rf $(DIR_OBJ)
