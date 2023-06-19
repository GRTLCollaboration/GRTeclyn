define SEARCH_FOR_MAKE
$(wildcard $1/*/GNUmakefile) $(wildcard $1/*/makefile) $(wildcard $1/*/Makefile)
endef

GRAMREX_HOME = $(realpath .)

TestsDir := $(GRAMREX_HOME)/Tests
ExampleDirsWithGNUmakefile := $(call SEARCH_FOR_MAKE, $(GRAMREX_HOME)/Examples)
ExampleDirs := $(dir $(ExampleDirsWithGNUmakefile))
CleanExampleDirs := $(ExampleDirs:%=clean-%)
CleanConfigTestsDir := $(TestsDir:%=cleanconfig-%)
CleanConfigExampleDirs := $(ExampleDirs:%=cleanconfig-%)

.PHONY: all examples tests run clean cleanconfig $(ExampleDirs)

ECHO?=@ # set this to null on the command line to increase verbosity

tests:
	$(info ################# Making Tests #################)
	$(ECHO)$(MAKE) -C $(TestsDir) --no-print-directory

run: tests
	$(info ################# Running Tests #################)
	$(ECHO)$(MAKE) -C $(TestsDir) --no-print-directory run

examples: $(ExampleDirs)

all: tests examples

clean: clean-testsdir $(CleanExampleDirs)

cleanconfig: $(CleanConfigTestsDir) $(CleanConfigExampleDirs)

$(ExampleDirs):
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
