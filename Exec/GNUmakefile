AMREX_HOME = ../Submodules/amrex

## --- Compiler Options --- ##
DEBUG        = FALSE
COMP         = gnu
USE_OMP      = FALSE
USE_MPI      = TRUE
##########################

## --- User Inputs --- ##
DIM          = 3
EBASE        = gflame
#########################

## --- Include --- ##
include $(AMREX_HOME)/Tools/GNUMake/Make.defs
include ../Src/Make.package
VPATH_LOCATIONS  += ../Src
INCLUDE_LOCATIONS += ../Src
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Tools/GNUMake/Make.rules
######################
