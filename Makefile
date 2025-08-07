# Makefile for tswHist Project (https://github.com/cyber-g/tswHist)
#
# Description:
#   Automates the build, clean, debug, and test process for the tswHist MATLAB/Octave MEX functions.
#
# Targets:
#   all      : Build all MEX files
#   debug    : Build debug versions of all MEX files
#   clean    : Remove all built MEX files
#   test     : Run all MATLAB test scripts in the test directory
#
# Variables:
#   MEX      : MATLAB/Octave mex compiler (default: /usr/local/bin/mex)
#   MEXEXT   : Extension for MEX files (default: mexa64)
#
# Author: Germain PHAM
# Date: August 2025
# License: GNU General Public License v3.0


MEX:= /usr/local/bin/mex
# Override at command line with:
# make MEX=/path/to/mex

MEXEXT:= mexa64
# Override at command line with:
# make MEXEXT=mexa64

# MEX functions to compile
SRC := $(wildcard *.c)
HDR := $(wildcard *.h)

# MEX objects to create
MEXOBJ := $(patsubst %.c,%.$(MEXEXT),$(SRC))

# Default target to build all MEX files
all: $(MEXOBJ)

# MEX compilation flags
MEXFLAGS := -largeArrayDims 

# Compilation and linking rule for MEX files
%.$(MEXEXT): %.c $(HDR)
	$(MEX) $(MEXFLAGS) $< -output $*

# Debug MEX files
DEBUGOBJ := $(patsubst %.c,%_debug.$(MEXEXT),$(SRC))
# Debug compilation flags
MEXDEBUGFLAGS := -g

debug: $(DEBUGOBJ)

# Debug compilation rule
%_debug.$(MEXEXT): %.c $(HDR)
	$(MEX) $(MEXFLAGS) $(MEXDEBUGFLAGS) $< -output $*

clean:
	rm -f $(MEXOBJ) $(DEBUGOBJ)
	@echo "Cleaned up MEX files."

test: $(MEXOBJ)
	@echo "Running tests..."
	cd test && \
	for file in *.m; do \
		echo "Running $$file..."; \
		matlab -batch "$$(basename $$file .m)"; \
	done

.PHONY: all clean debug test