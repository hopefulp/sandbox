
# Commands
FC=gfortran
F2PY=f2py
PYTHON=python

# Library name
NAME=ric

# Name of wrapped module file
MODNAME=$(NAME).f90

# Name of Python module
PYNAME=_$(NAME)

# Name of F2Py signature file
SIGNNAME=$(PYNAME).pyf

# Name of library
LIBNAME=$(PYNAME).so

# LAPACK
LAPACK=-llapack

# F2PY debug flags
F2PY_DEBUG=

# List of all Fortran source files
SRC= ric_common.f90 \
     ric_math.f90 \
     ric_stretches.f90 \
     ric_in_bends.f90 \
     ric_out_bends.f90 \
     ric_lin_bends.f90 \
     ric_torsions.f90 \
     ric_eckart.f90 \
     $(MODNAME)

OBJ=$(SRC:.f90=.o)

all: $(LIBNAME)

debug: F2PY_DEBUG += --debug-capi -DF2PY_REPORT_ATEXIT -DF2PY_REPORT_ON_ARRAY_COPY=1
debug: $(LIBNAME)

$(LIBNAME): $(SIGNNAME) $(OBJ)
	$(F2PY) $(F2PY_DEBUG) -c $(SIGNNAME) $(SRC) $(LAPACK)

$(SIGNNAME): $(MODNAME)
	$(F2PY) -m $(PYNAME) -h $(SIGNNAME) $(MODNAME)

%.o: %.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $<

check:
	cd tests; \
	$(PYTHON) -tt3m unittest discover

doc:
	cd ..; \
	$(PYTHON) -c 'import ric; help(ric)'

clean:
	$(RM) $(SIGNNAME)
	$(RM) *.o
	$(RM) *.mod

realclean: clean
	$(RM) $(PYNAME).so
	$(RM) *.pyc
	$(RM) tests/*.pyc

distclean: realclean
	$(RM) *~
	$(RM) tests/*~

