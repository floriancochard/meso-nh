VPATH = src:$(DIR_OBJ)
#######################################
DIR_OBJ = ./$(ARCH)

ifeq ($(origin SRC_MESONH), undefined)
SRC_MESONH := $(shell pwd|sed -e 's/\/tools\/.*//')
endif

ifeq ($(origin DIR_LIB), undefined)
DIR_LIB := $(SRC_MESONH)/lib
endif

DIR_LFI = $(DIR_LIB)/NEWLFI
LIBLFI  = $(DIR_LFI)/$(ARCH)/libNEWLFI_ALL.a

DIR_COMP = $(DIR_LIB)/COMPRESS
LIBCOMP  = $(DIR_COMP)/$(ARCH)/liblficomp.a


OBJS = mode_options.o lfi2cdf.o modd_ncparam.o mode_dimlist.o fieldtype.o mode_util.o
PROGS = lfi2cdf 

INC = -I$(DIR_OBJ) -DLFI_INT=$(LFI_INT)

DIR_CONF:=$(SRC_MESONH)/conf

include $(DIR_CONF)/config.$(ARCH)
include Rules.$(ARCH)

%.o:%.f90 $(DIR_OBJ)/.dummy
	$(CPP) $(INC) $(CPPFLAGS) $< >  $(DIR_OBJ)/cpp_$(*F).f90
	$(F90) $(INC) -c $(F90FLAGS) $(DIR_OBJ)/cpp_$(*F).f90 -o $(DIR_OBJ)/$(*F).o
	-@mv  *.mod $(DIR_OBJ)/. 2> /dev/null || echo pas de module dans $*.f90

%.o:%.c $(DIR_OBJ)/.dummy
	$(CC) $(INC) $(CFLAGS) $(CPPFLAGS) -c $< -o $(DIR_OBJ)/$(*F).o

all : $(PROGS) cdf2cdf cdf2lfi

cdf2cdf: $(PROGS) 
	cd $(DIR_OBJ); rm -f cdf2cdf; ln -s $(PROGS) cdf2cdf

cdf2lfi: $(PROGS) 
	cd $(DIR_OBJ); rm -f cdf2lfi; ln -s $(PROGS) cdf2lfi

$(PROGS): $(OBJS) $(LIBLFI) $(LIBCOMP)
	cd $(DIR_OBJ); $(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBLFI) $(LIBCOMP) $(LIBCDF)

#$(OBJS): $(LIBCDF)
$(OBJS): 

$(DIR_OBJ)/.dummy :
	mkdir -p $(DIR_OBJ)
	@touch $(DIR_OBJ)/.dummy

$(LIBLFI): $(DIR_LFI)
	$(MAKE) -C $(DIR_LFI)

$(LIBCOMP): $(DIR_COMP)
	$(MAKE) -C $(DIR_COMP)

$(DIR_LFI):
	@echo "ERROR : NEWLFI directory can't be found"
	@echo "        from root directory DIR_LIB = $(DIR_LIB)";echo
	@echo "please check SRC_MESONH or DIR_LIB (= \$$SRC_MESONH/lib) env. variable"
	@echo "and try again...";exit 1

$(LIBCDF):
	@echo "*************      NETCDF library not found        ***************";\
	echo "Please, give NETCDFHOME variable in Rules.$(ARCH) the right path !";\
	echo "******************************************************************";\
	exit 1

clean:
	(if [ -d $(DIR_OBJ) ] ; then cd $(DIR_OBJ); rm -f cpp_*.f90 cpp_*.f *.o *.mod ; fi)

distclean:
	rm -rf $(DIR_OBJ) 


lfi2cdf.o: lfi2cdf.f90 mode_util.o
modd_ncparam.o: modd_ncparam.f90
mode_dimlist.o: mode_dimlist.f90
mode_util.o: mode_util.f90 modd_ncparam.o fieldtype.o mode_dimlist.o
fieldtype.o: fieldtype.f90 modd_ncparam.o

