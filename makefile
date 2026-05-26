# ============================================================================
# Makefile for the dyablo FoF tool.
#
# Default build: macOS / Apple Silicon with Homebrew HDF5 + Open MPI.
#   brew install hdf5 open-mpi
# On a cluster, override HDF5_DIR and MPIFC on the command line or via env.
# ============================================================================

HDF5_DIR ?= /opt/homebrew

MPIFC := mpif90
FC    := mpif90

# for flanders

MPIFC := h5pfc
FC    := h5pfc

# -DLONGINT activates 8-byte particle IDs (needed if N_particles > 2^31).
OPTS    ?= -O3 -cpp -DLONGINT \
           -I$(HDF5_DIR)/include \
           -ffree-line-length-none

LIBS    ?= -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5

MODS = 			\
	modconst.mod	\
	modvariable.mod	\
	modtiming.mod	\
	modparam.mod	\
	modmpicom.mod	\
	modio.mod	\
	modsort.mod	\
	modfofpara.mod

OBJS = 			\
	modconst.o	\
	modvariable.o	\
	modtiming.o	\
	modparam.o	\
	modmpicom.o	\
	modio.o 	\
	modsort.o	\
	modfofpara.o	\
	fof.o


.SUFFIXES: .f90

.f90.o:
	$(FC) $(OPTS) -c $<

all : fof

fof : $(OBJS)
	$(MPIFC) $(OBJS) -o fof $(LIBS)

# Module dependency order
modconst.o    : modconst.f90
modvariable.o : modvariable.f90 modconst.o
modparam.o    : modparam.f90    modconst.o
modtiming.o   : modtiming.f90   modconst.o
modmpicom.o   : modmpicom.f90   modconst.o
modio.o       : modio.f90       modconst.o modparam.o modvariable.o modmpicom.o modtiming.o
modsort.o     : modsort.f90     modconst.o modmpicom.o
modfofpara.o  : modfofpara.f90  modconst.o modparam.o modvariable.o modmpicom.o modtiming.o modsort.o modio.o
fof.o         : fof.f90         modconst.o modparam.o modvariable.o modmpicom.o modtiming.o modio.o modfofpara.o

clean :
	rm -f *.mod *.o *~ fof
