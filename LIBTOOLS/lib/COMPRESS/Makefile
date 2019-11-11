LIBCOMP = liblficomp.a
#######################################
DIR_OBJ = ./$(ARCH)

VPATH = src:$(DIR_OBJ)
INC = -I$(DIR_OBJ)

DIR_CONF:=$(shell pwd|sed -e 's/lib\/.*/conf/')

include $(DIR_CONF)/config.$(ARCH)
include Rules.$(ARCH)


%.o:%.f90
	$(CPP) $(INC) $(CPPFLAGS) $< >  $(DIR_OBJ)/cpp_$(*F).f90
	$(F90) $(INC) -c $(F90FLAGS) $(DIR_OBJ)/cpp_$(*F).f90 -o $(DIR_OBJ)/$(*F).o
	-@mv  *.mod $(DIR_OBJ)/. 2> /dev/null || echo pas de module dans $*.f90

%.o:%.c
	$(CC) $(INC) $(CFLAGS) $(CPPFLAGS) -c $< -o $(DIR_OBJ)/$(*F).o

$(LIBCOMP) : $(DIR_OBJ)/.dummy $(OBJS)
	cd $(DIR_OBJ);$(AR) crv $@ $(OBJS)

$(DIR_OBJ)/.dummy :
	mkdir -p $(DIR_OBJ)
	@touch $(DIR_OBJ)/.dummy

compress.o : searchgrp.o comppar.o
decompress.o : searchgrp.o comppar.o

clean:
	(if [ -d $(DIR_OBJ) ] ; then cd $(DIR_OBJ); rm -f cpp_*.f90 *.o ; fi)

distclean:
	rm -rf $(DIR_OBJ)


