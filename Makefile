# dependencies
composition_mod.o 	: bftrans_mod.o osmet_mod.o os_mod.o handyfuncs_mod.o optab_mod.o hprofl_mod.o
element_mod.o		    : op_math_mod.o ion_mod.o
molecule_mod.o		  : op_math_mod.o ion_mod.o
ion_mod.o			      : os_mod.o prec_mod.o level_mod.o bbtrans_mod.o
level_mod.o			    : prec_mod.o
bftrans_mod.o	      : element_mod.o molecule_mod.o os_mod.o
bbtrans_mod.o       : op_math_mod.o level_mod.o
op_math_mod.o		    : prec_mod.o
osmet_mod.o         : prec_mod.o op_math_mod.o
optab_mod.o         : prec_mod.o spectrum_mod.o
prec_mod.o          : const_mod.o trace_mod.o
handyfuncs_mod.o    : os_mod.o op_math_mod.o
hprofl_mod.o        : const_mod.o trans_mod.o
trans_mod.o         : bbtrans_mod.o mmath_mod.o
const_mod.o         :