B = 64
DIR_OBJ=./$(ARCH)_$(B)

ifeq ($(origin MNH_LIBTOOLS), undefined)
include ../where.Libs
else
include $(MNH_LIBTOOLS)/tools/where.Libs
endif

VPATH=src:$(DIR_DIA)/$(DIR_OBJ)

#INC = -I $(DIR_OBJ)
INC = -I $(DIR_OBJ) -I $(DIR_DIA)/$(DIR_OBJ)

PROG = fmmore

OBJS = readuntouch.o 

include $(DIR_CONF)/config.$(ARCH)
include Rules.$(ARCH)

%.o:%.f90 $(DIR_OBJ)/.dummy
	$(CPP) $(INC) $(CPPFLAGS) $< >  $(DIR_OBJ)/cpp_$(*F).f90
	$(F90) $(INC) -c $(F90FLAGS) $(DIR_OBJ)/cpp_$(*F).f90 -o $(DIR_OBJ)/$(*F).o

$(PROG): $(PROG).o $(OBJS) $(LIBDIA) $(LIBLFI) $(LIBCOMP)
	cd $(DIR_OBJ); $(F90) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo executable $(PROG) available under $(DIR_OBJ)

$(DIR_OBJ)/.dummy :
	mkdir $(DIR_OBJ)
	@touch $(DIR_OBJ)/.dummy

$(LIBLFI):
	$(MAKE) -C $(DIR_LFI)

$(LIBCOMP):
	$(MAKE) -C $(DIR_COMP)

clean:
	(if [ -d $(DIR_OBJ) ] ; then cd $(DIR_OBJ); rm -f cpp_*.f90 *.o ; fi)

distclean:
	(if [ -d $(DIR_OBJ) ] ; then  rm -rf $(DIR_OBJ) ;fi)

