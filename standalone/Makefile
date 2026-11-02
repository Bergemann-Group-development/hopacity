# Build HOP as a standalone shared library (libhop.so), the hop_prune.x
# line-list pruning executable and the atom_convert.x model-atom converter, with
# no DISPATCH framework dependencies.  The
# framework singletons (io, os, trace) are replaced by the -DHOP_STANDALONE
# branch of ../hop_env_mod.f90, and every source below lives under HOP/ --
# nothing is pulled in from ../../SRC, so the package is self-contained.  This
# directory is invisible to the main experiment build (its wildcard only scans
# SRC/*.f90 HOP/*.f90); the top-level Makefile's hop_standalone target invokes it
# with its own COMPILER and OPTS, and the objects here (compiled with
# -DHOP_STANDALONE) never mix with the framework-enabled ones linked into dispatch.x.
#
# The default target also copies libhop.so into ../hopacity, where the Python package
# loads it, and unpacks the ABO table that HOP reads from ../hopacity/data; the wheel
# build (../setup.py) makes $(BUILD)/libhop.so alone instead.
#
# Usage:  make [-j N] [COMPILER=ifort] [OPT=-O2] [BUILD=dir]   # libhop.so, hop_prune.x, atom_convert.x
#         make bindings                                        # rewrite libhop.h after an API change
BUILD ?= build
LIB   := $(BUILD)/libhop.so
PKGLIB := ../hopacity/libhop.so
ANSTEE := ../hopacity/data/Anstee-1802-2014.dat
HDR   := libhop.h
PRUNE := $(BUILD)/hop_prune.x
CONVERT := $(BUILD)/atom_convert.x
DISPATCH_TOP ?= ../../../..
COMPILER ?= gfortran
OPTS ?= hoppe
# Match the application build's optimisation and architecture flags.  Outside a Dispatch
# checkout (the hopacity repository) the fallbacks take its optimisation without
# -ffast-math and no architecture flags, so the library runs on any host of the ISA.
-include $(DISPATCH_TOP)/config/compiler/$(COMPILER)/$(OPTS).mkf
OPT   ?= -O3
HW    ?=
# Hide the compiler lines as the DISPATCH build does; `make E=` shows them.
E = @

ifeq ($(COMPILER),gfortran)
  FC      = gfortran
  DEBUG  ?= -fbacktrace
  FFLAGS  = -x f95-cpp-input -DHOP_STANDALONE -ffree-line-length-0 -w \
            -fno-range-check -fallow-argument-mismatch -fno-realloc-lhs -fopenmp \
            $(OPT) $(HW) $(DEBUG) -fPIC -I$(BUILD) -J$(BUILD)
  LDFLAGS = $(OPT) $(DEBUG) -fopenmp
else ifneq ($(filter ifort ifx,$(COMPILER)),)
  FC      = $(COMPILER)
  DEBUG  ?= -traceback
  FFLAGS  = -fpp -DHOP_STANDALONE -assume byterecl,protect_parens -diag-disable remark \
            -diag-disable=16219,10448 -qopenmp $(OPT) $(HW) $(DEBUG) -fPIC \
            -I$(BUILD) -module $(BUILD)
  LDFLAGS = $(OPT) $(DEBUG) -qopenmp
else
  $(error COMPILER=$(COMPILER) is not supported here; use gfortran, ifort or ifx)
endif

# Topologically ordered source list (verified USE graph)
SRCS = ../console_mod.f90 ../namelist_check_mod.f90 ../qlg_mod.f90 ../hop_env_mod.f90 ../const_mod.f90 ../prec_mod.f90 ../op_math_mod.f90 \
       ../interpol_mod.f90 ../handyfuncs_mod.f90 ../linelist_reader_mod.f90 ../anstee_mod.f90 \
       ../hvcs_mod.f90 \
       ../level_mod.f90 ../line_types_mod.f90 ../bbtrans_mod.f90 \
       ../ion_mod.f90 ../element_mod.f90 ../molecule_mod.f90 \
       ../hprofl_mod.f90 ../linelist_mod.f90 \
       ../bftrans_mod.f90 ../trans_mod.f90 ../spectrum_mod.f90 \
       ../optab_mod.f90 ../osmet_mod.f90 ../composition_mod.f90 ../atmos_reader_mod.f90 \
       ../atom_convert_mod.f90 hop_c_api.f90

OBJS = $(addprefix $(BUILD)/,$(notdir $(SRCS:.f90=.o)))
# The pruning program is the library objects plus its main program.
PRUNE_SRCS = hop_prune.f90
PRUNE_OBJS = $(addprefix $(BUILD)/,$(notdir $(PRUNE_SRCS:.f90=.o)))

vpath %.f90 . ..

all: $(LIB) $(PKGLIB) $(ANSTEE) $(PRUNE) $(CONVERT)
	@:

# The header is the C prototypes of the BIND(C) functions, written by gen_bindings.py
# from the shim source and committed.  It is source, not a build product: `make bindings`
# rewrites it after an interface change, and tests/test_hop_abi.py fails on a stale
# one, so building the library needs neither Python nor a header-writing compiler.
bindings:
	python3 gen_bindings.py hop_c_api.f90 > $(HDR)

# Chain each object to its predecessor: the topological list stands in for
# per-module .mod dependencies, so editing any module rebuilds everything
# after it (keeps incremental and parallel builds consistent).
prev :=
$(foreach o,$(OBJS),$(eval $(o): $(prev))$(eval prev := $(o)))

$(LIB): $(OBJS)
	@if [ "x$(E)" == "x@" ] ; then echo " [ LD ]" $(notdir $@) ; fi
	$(E)$(FC) -shared $(LDFLAGS) -o $@ $(OBJS)

# Copied by rename, so a process that has the old library loaded keeps a valid file
$(PKGLIB): $(LIB)
	$(E)cp $< $@.tmp && mv -f $@.tmp $@

# Tracked compressed
$(ANSTEE): $(ANSTEE).gz
	$(E)gunzip -c $< > $@

# Fortran INCLUDE files are not in the chain above; name their includers
$(BUILD)/level_mod.o $(BUILD)/bftrans_mod.o: ../exp_neg.inc

# hop_prune.x: the library objects plus the program
$(BUILD)/hop_prune.o: $(lastword $(OBJS))
$(PRUNE): $(OBJS) $(PRUNE_OBJS)
	@if [ "x$(E)" == "x@" ] ; then echo " [ LD ]" $(notdir $@) ; fi
	$(E)$(FC) $(LDFLAGS) -o $@ $(OBJS) $(PRUNE_OBJS)

hop_prune.x: $(PRUNE)

# atom_convert.x: the library objects plus the program
$(BUILD)/atom_convert.o: $(lastword $(OBJS))
$(CONVERT): $(OBJS) $(BUILD)/atom_convert.o
	@if [ "x$(E)" == "x@" ] ; then echo " [ LD ]" $(notdir $@) ; fi
	$(E)$(FC) $(LDFLAGS) -o $@ $(OBJS) $(BUILD)/atom_convert.o

atom_convert.x: $(CONVERT)

$(BUILD)/%.o: %.f90 | $(BUILD)
	@if [ "x$(E)" == "x@" ] ; then echo " [ FC ]" $(notdir $@) ; fi
	$(E)$(FC) $(FFLAGS) -c $< -o $@

$(BUILD):
	$(E)mkdir -p $@

# gfortran searches the compilation directory for .mod files ahead of -I and -J, so a
# module file left in this source directory shadows the freshly built one in $(BUILD)
# and makes later compiles fail on stale derived types.  Nothing here writes .mod files
# into the source directory on purpose, so clean removes any that a stray compiler run
# has dropped.
clean:
	rm -rf $(BUILD) $(PKGLIB)
	rm -f *.mod

.PHONY: all clean bindings hop_prune.x atom_convert.x
