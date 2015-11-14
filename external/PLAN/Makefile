# include Make.inc for machine specific commands and options
include Make.inc
#
.DEFAULT:
	@$(ECHO) "Unknown target $@"
	@$(ECHO) "The following targes are defined: lib plib clean"
#
lib liblanso.a:
	$(CD) lanso1; $(MAKE) lib
#
plib libplan.a:
	$(CD) plan; $(MAKE) lib
#
clean:
	$(CD) lanso1; $(MAKE) clean
	$(CD) plan; $(MAKE) clean
	$(CD) examples; $(MAKE) clean
#
