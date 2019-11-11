SUBDIRS = lfiz lfi2cdf diachro fmmore vergrid
.PHONY: subdirs $(SUBDIRS)

ifndef ARCH

VALID_ARCH=$(subst ../conf/config.,,$(wildcard ../conf/config.*))
dummy %:
	@echo "ERROR : ARCH variable is not set !";echo
	@echo "Please, choose one of these statements then try again :";echo " "
	@for i in $(VALID_ARCH); do echo export ARCH=$$i; done
else	

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean distclean:
	@for dir in $(SUBDIRS); do \
	$(MAKE) -C $$dir $@; \
	done

endif



