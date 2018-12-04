# Master makefile for DL_POLY_2.0 
# Author: W. Smith January 2009
# 
# Modified for pydlpoly
# 
# this is a special Makefile for building the python wrapped
# version of dl_poly
# instead of dlpoly.f we build pydlpoly.f, which is just a module
# with subroutines wrapped via f2py
# The trick is: we need only to make the .o files
# f2py will build the shared object 
#
#=======================================================================
# Define default settings
#=======================================================================

SHELL=/bin/sh
MAKE=make -f Makefile

#=====================================================================
# Define object files

OBJ_MOD = parse_module.o setup_module.o error_module.o \
	config_module.o site_module.o pair_module.o \
	molecule_module.o basic_comms.o \
	utility_module.o solvation_module.o tether_module.o \
	vdw_module.o property_module.o rigid_body_module.o \
	angles_module.o bonds_module.o shake_module.o \
	inversion_module.o dihedral_module.o core_shell_module.o \
	exclude_module.o nlist_builders_module.o \
	ewald_module.o coulomb_module.o\
	external_field_module.o four_body_module.o \
	hkewald_module.o metal_module.o ensemble_tools_module.o \
	temp_scalers_module.o three_body_module.o spme_module.o \
	tersoff_module.o neu_coul_module.o timing_module.o\
	forces_module.o \
	lf_motion_module.o lf_rotation1_module.o \
	lf_rotation2_module.o vv_motion_module.o \
	vv_rotation1_module.o vv_rotation2_module.o \
	pmf_module.o integrator_module.o optimiser_module.o \
	hyper_dynamics_module.o driver_module.o \
	define_system_module.o \
	lbfgs_module.o metamd_module.o

OBJ_SRC = pydlpoly.o

OBJ_PAR = merge_tools.o pass_tools.o

# never forget the PIC -> we build an so!

INTEL = False 
ifneq ($(HOSTNAME),cmcc.cmcc.cluster) 
#================== GNU Fortran, MPI version ==============================
# gfortran (SUSE version , fftw2 is called libdfftw
FC=mpif90
LD=mpif90 -o
LDFLAGS=-O2 -ffast-math
#FFLAGS=-c -O2 -ffast-math -fPIC -cpp -DDEBUG
FFLAGS=-c -O2 -ffast-math -fPIC -cpp
# For FFTW2
# LDLIBS=-lfftw
# For FFTW3
# LDLIBS=-lfftw3 -lblas
LDLIBS=-lblas -lfftw3
F2PY_COMP=gnu95
else
ifeq ($(INTEL),True)
## cmcc version using ifort and fftw and MKL Blas
FC=mpiifort
LD=mpiifort -o
LDFLAGS=-O2 $(LDLIBS)
LDLIBS=-lblas -L${FFTW_PATH}/lib -lfftw3
FFLAGS= -FI -c -O2 -xhost -fPIC -cpp
F2PY_COMP=intelem
else
## cmcc version using gfortran and MKL Blas
FC=mpif90
LD=mpif90 -o
LDFLAGS=-O2 -ffast-math
#FFLAGS=-c -O2 -ffast-math -fPIC -cpp -DDEBUG
FFLAGS=-c -O2 -ffast-math -fPIC -cpp
# For FFTW2
# LDLIBS=-lblas -L${FFTW_PATH}/lib -lfftw
# For FFTW3
LDLIBS=-lblas -L${FFTW_PATH}/lib -lfftw3
F2PY_COMP=gnu95
endif
endif

#=====================================================================
all: _pydlpoly.so

_pydlpoly.so: pydlpoly.o $(OBJ_PAR) _pydlpoly.pyf
	f2py --fcompiler=$(F2PY_COMP) --f90exec=$(FC) _pydlpoly.pyf $(OBJ_MOD) $(OBJ_SRC) $(OBJ_PAR) -c $(LDLIBS) --link-lapack_mkl 
	cp _pydlpoly.so ../so/_pydlpoly.so


#=====================================================================
# Clean up the source directory
clean:
	rm -f $(OBJ_MOD) $(OBJ_PAR) $(OBJ_SRC) *.mod

#=====================================================================
# Declare dependencies
.f.o: 
	$(FC) $(FFLAGS) $*.f
.c.o: 
	$(CC) -c $*.c

#=====================================================================
# Declare dependency on module files

$(OBJ_SRC): $(OBJ_MOD)

