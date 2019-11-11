GRIB_DIR=$(wildcard gribex*)
SUBDIRS = NEWLFI COMPRESS MPIvide RAD2 SURCOUCHE vis5d
.PHONY: subdirs $(SUBDIRS) $(GRIB_DIR)

ifndef ARCH
VALID_ARCH=$(subst ../conf/config.,,$(wildcard ../conf/config.*))
dummy %:
	@echo "ERROR : ARCH variable is not set !";echo
	@echo "Please, choose one of these statements then try again :";echo " "
	@for i in $(VALID_ARCH); do echo export ARCH=$$i; done

else	
subdirs: $(SUBDIRS) $(GRIB_DIR)

$(SUBDIRS):
	$(MAKE) -C $@

$(GRIB_DIR):
	@echo "==========================================================================="
	@echo "GRIB library : please go into $@ directory and see README files"
	@echo "               in order to generate manually the GRIB library."
	@echo "==========================================================================="

clean distclean:
	@for dir in $(SUBDIRS) $(GRIB_DIR); do \
	$(MAKE) -C $$dir $@; \
	done

endif



