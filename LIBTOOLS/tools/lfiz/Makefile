VPATH=src:$(DIR_OBJ)
###########################
DIR_OBJ = ./$(ARCH)

include ../where.Libs

INC = -I$(DIR_COMP)/$(ARCH)

include $(DIR_CONF)/config.$(ARCH)
include Rules.$(ARCH)

PROGS = lfiz unlfiz testlibcomp

%.o:%.f90 $(DIR_OBJ)/.dummy
	$(CPP) $(INC) $(CPPFLAGS) $< >  $(DIR_OBJ)/cpp_$(*F).f90
	$(F90) $(INC) -c $(F90FLAGS) $(DIR_OBJ)/cpp_$(*F).f90 -o $(DIR_OBJ)/$(*F).o

all: $(PROGS)

$(PROGS): %:%.o  $(LIBCOMP) $(LIBLFI)
	cd $(DIR_OBJ); $(F90) $(LDFLAGS) -o $@ $@.o $(LIBLFI) $(LIBCOMP) $(LIBS)

$(DIR_OBJ)/.dummy:
	mkdir -p $(DIR_OBJ)
	@touch $(DIR_OBJ)/.dummy

$(LIBLFI): $(DIR_LFI)
	$(MAKE) -C $(DIR_LFI)

$(LIBCOMP): $(DIR_COMP)
	$(MAKE) -C $(DIR_COMP)

$(DIR_LFI) $(DIR_COMP):        
	@echo "ERROR : COMPRESS and/or NEWLFI directory can't be found"
	@echo "        from root directory DIR_LIB = $(DIR_LIB)";echo
	@echo "Please check SRC_MESONH or DIR_LIB environment variable"
	@echo "and try again...";exit 1

clean:
	(if [ -d $(DIR_OBJ) ] ; then cd $(DIR_OBJ); rm -f cpp_*.f90 cpp_*.f *.o *.mod ; fi)

distclean:
	rm -rf $(DIR_OBJ)
