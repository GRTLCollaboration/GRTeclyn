define SEARCH_FOR_MAKE
$(wildcard $1/*/GNUmakefile) $(wildcard $1/*/makefile) $(wildcard $1/*/Makefile)
endef

GRTECLYN_HOME = $(realpath .)

TestsDir := $(GRTECLYN_HOME)/Tests
ExampleDirsWithGNUmakefile := $(call SEARCH_FOR_MAKE, $(GRTECLYN_HOME)/Examples)
ExampleDirs := $(dir $(ExampleDirsWithGNUmakefile))
CleanExampleDirs := $(ExampleDirs:%=clean-%)
CleanConfigTestsDir := $(TestsDir:%=cleanconfig-%)
CleanConfigExampleDirs := $(ExampleDirs:%=cleanconfig-%)

.PHONY: all examples tests tests-config run clean cleanconfig $(ExampleDirs)

ECHO?=@ # set this to null on the command line to increase verbosity

tests: tests-config
	$(info ################# Making Tests #################)
	$(ECHO)$(MAKE) -C $(TestsDir) --no-print-directory

# Separate this out from the tests target just to avoid the race condition
# when doing make all (see below).
tests-config:
	$(ECHO)$(MAKE) -C $(TestsDir) --no-print-directory AMReX_Config.H
	$(ECHO)$(MAKE) -C $(TestsDir) --no-print-directory AMReX_Version.H

run: tests
	$(info ################# Running Tests #################)
	$(ECHO)$(MAKE) -C $(TestsDir) --no-print-directory run

examples: $(ExampleDirs)

all: tests examples

clean: clean-testsdir $(CleanExampleDirs)

cleanconfig: $(CleanConfigTestsDir) $(CleanConfigExampleDirs)

# We add the tests-config dependency just for the case where the build
# configuration is the same for tests and examples to avoid the race condition
# in generating AMReX_Config.H and AMReX_Version.H when doing make all.
$(ExampleDirs): tests-config
	$(info ################# Making example $@ #################)
	$(ECHO)$(MAKE) -C $@ --no-print-directory

clean-testsdir:
	$(ECHO)$(MAKE) -C $(TestsDir) --no-print-directory clean

$(CleanExampleDirs):
	$(ECHO)$(MAKE) -C $(@:clean-%=%) --no-print-directory clean 

$(CleanConfigTestsDir):
	$(ECHO)$(MAKE) -C $(@:cleanconfig-%=%) --no-print-directory cleanconfig 

$(CleanConfigExampleDirs):
	$(ECHO)$(MAKE) -C $(@:cleanconfig-%=%) --no-print-directory cleanconfig 
