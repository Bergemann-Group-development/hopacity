# Build HOP as a standalone shared library (libhop.so) with no DISPATCH
# framework dependencies.  The framework singletons (io, os, trace) are
# replaced by the -DHOP_STANDALONE branch of ../hop_env_mod.f90, and every
# source below lives under HOP/ -- nothing is pulled in from ../../SRC, so
# the package is self-contained.  This directory is invisible to the main
# experiment build (its wildcard only scans SRC/*.f90 HOP/*.f90).
#
# Usage:  make [-j N] [OPT=-O2]
BUILD := build
LIB   := $(BUILD)/libhop.so
FC    := gfortran
OPT   ?= -O1
FFLAGS = -x f95-cpp-input -DHOP_STANDALONE -ffree-line-length-0 -w \
         -fno-range-check -fallow-argument-mismatch $(OPT) -g -fbacktrace \
         -fopenmp -fPIC -I$(BUILD) -J$(BUILD)

# Topologically ordered source list (verified USE graph)
SRCS = ../hop_env_mod.f90 ../const_mod.f90 ../prec_mod.f90 ../op_math_mod.f90 \
       ../interpol_mod.f90 ../handyfuncs_mod.f90 ../hvcs_mod.f90 \
       ../level_mod.f90 ../line_types_mod.f90 ../bbtrans_mod.f90 \
       ../ion_mod.f90 ../element_mod.f90 ../molecule_mod.f90 \
       ../hprofl_mod.f90 ../voigt_av_mod.f90 ../linelist_mod.f90 \
       ../bftrans_mod.f90 ../trans_mod.f90 ../spectrum_mod.f90 \
       ../optab_mod.f90 ../osmet_mod.f90 ../composition_mod.f90 hop_c_api.f90

OBJS = $(addprefix $(BUILD)/,$(notdir $(SRCS:.f90=.o)))

vpath %.f90 . ..

all: $(LIB)

# Chain each object to its predecessor: the topological list stands in for
# per-module .mod dependencies, so editing any module rebuilds everything
# after it (keeps incremental and parallel builds consistent).
prev :=
$(foreach o,$(OBJS),$(eval $(o): $(prev))$(eval prev := $(o)))

$(LIB): $(OBJS)
	$(FC) -shared -fopenmp -o $@ $(OBJS)

$(BUILD)/%.o: %.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $@

clean:
	rm -rf $(BUILD)

.PHONY: all clean
