VPATH=src:$(DIR_OBJ)
###########################
DIR_OBJ = ./$(ARCH)

include ../where.Libs

include $(DIR_CONF)/config.$(ARCH)
#include Rules.$(ARCH)

INC = -I $(DIR_OBJ)

PROG = vergrid

OBJS = mode_pos.o

%.o:%.f90 $(DIR_OBJ)/.dummy
	$(CPP) $(INC) $(CPPFLAGS) $< >  $(DIR_OBJ)/cpp_$(*F).f90
	$(F90) $(INC) -c $(F90FLAGS) $(DIR_OBJ)/cpp_$(*F).f90 -o $(DIR_OBJ)/$(*F).o
	-@mv  *.mod $(DIR_OBJ)/. 2> /dev/null || echo pas de module dans $*.f90

all: $(PROG)

$(PROG): $(PROG).o $(OBJS)  
	cd $(DIR_OBJ); $(F90) $(LDFLAGS) $(patsubst $(DIR_OBJ)/%,%,$^) -o $@ $(LIBS)

$(DIR_OBJ)/.dummy:
	mkdir $(DIR_OBJ)
	@touch $(DIR_OBJ)/.dummy

clean:
	(if [ -d $(DIR_OBJ) ] ; then cd $(DIR_OBJ); rm -f cpp_*.f90 cpp_*.f *.o *.mod ; fi)

distclean:
	(if [ -d $(DIR_OBJ) ] ; then  rm -rf $(DIR_OBJ) ;fi)


$(PROG).o: $(PROG).f90 mode_pos.o
