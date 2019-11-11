DIR_OBJ= ./$(ARCH)

ifndef ARCH
VALID_ARCH=$(subst Rules.,,$(wildcard Rules.*))
dummy %:
	@echo "Error : ARCH variable is not set ! Valid values are :"
	@echo $(VALID_ARCH)

else
PROGALL = conv2dia diaprog extractdia  compute_r00_pc obs2mesonh mesonh2obs concat_time_diafile
include Rules.$(ARCH)

all : $(PROGALL)

#
# l ordre est a respecter
# 
conv2dia:
	$(MAKE) -f Makefile.conv2dia B=64
	$(MAKE) -f Makefile.conv2dia B=32

diaprog: 
	$(MAKE) -f Makefile.diaprog B=32

extractdia: 
	$(MAKE) -f Makefile.extractdia B=32

	
#exrwdia:
#	$(MAKE) -f Makefile.exrwdia
compute_r00_pc:
	$(MAKE) -f Makefile.exrwdia PROG=$@
obs2mesonh:
	$(MAKE) -f Makefile.exrwdia PROG=$@
mesonh2obs:
	$(MAKE) -f Makefile.exrwdia PROG=$@
concat_time_diafile:
	$(MAKE) -f Makefile.exrwdia PROG=$@


clean :
	@for dir in $(DIR_OBJ)_* ; do\
	(if [ -d $$dir ] ; then cd $$dir; rm -f cpp_*.f90 cpp_*.f *.o *.mod ; fi) \
	done

distclean :
	rm -rf $(DIR_OBJ)_*


endif
