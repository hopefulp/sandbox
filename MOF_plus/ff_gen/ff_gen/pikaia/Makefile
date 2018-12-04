
ifneq ($(HOSTNAME),cmcc.cmcc.cluster) 
#================== GNU Fortran, MPI version ==============================
F2PY_COMP=gnu95
else
## cmcc version using ifort
F2PY_COMP=intelem
endif


all: pikaia.f
	f2py --fcompiler=$(F2PY_COMP) -c pikaia.pyf pikaia.f 


