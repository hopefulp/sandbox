#!/bin/tcsh
set NPP = '4'
if ($#argv == 1) then
        set NPP = $argv[1]
endif

if ($HOST == kdft) then
	set NPP = 12
	set MPI = "/opt/intel/Compiler/composer_xe_2011_sp1.13.367/mpi/openmpi-1.6.3/bin/mpirun -np $NPP "
	set EXEC = "/opt/applic/vasp/bin/vasp-5.3.3-xe11-static-openmpi-1.6.3-full"
	setenv LD_LIBRARY_PATH /opt/intel/Compiler/composer_xe_2011_sp1.13.367/compiler/lib/intel64:/opt/intel/Compiler/composer_xe_2011_sp1.13.367/mkl/lib/intel64:/opt/intel/Compiler/composer_xe_2011_sp1.13.367/mpi/openmpi-1.6.3/lib:$LD_LIBRARY_PATH

else if ($HOST == psi) then
	set NPP = 16
	set MPI = "/opt/mpi/intel-12.1.6/openmpi-1.6.3/bin/mpirun -np $NPP "
	set EXEC = "/opt/applic/vasp/bin/vasp-5.3.5-xe11-static-openmpi-1.6.3-full"
	setenv LD_LIBRARY_PATH /opt/intel/composer_xe_2011_sp1.13.367/compiler/lib/intel64:/opt/intel/composer_xe_2011_sp1.13.367/mkl/lib/intel64:/opt/mpi/intel-12.1.6/openmpi-1.6.3/lib:$LD_LIBRARY_PATH
endif

else if ($HOST == rho) then
        set NPP = 16
        set MPI = "/opt/mpi/intel-12/openmpi-1.8.1/bin/mpirun -np $NPP "
        set EXEC = "/opt/applic/vasp/bin/vasp-5.3.5-if12-openmpi-1.8.1-full"
        setenv LD_LIBRARY_PATH /opt/intel/composer_xe_2011_sp1.13.367/compiler/lib/intel64:/opt/intel/composer_x
e_2011_sp1.13.367/mkl/lib/intel64:/opt/mpi/intel-12.1.6/openmpi-1.6.3/lib:$LD_LIBRARY_PATH
endif


$MPI $EXEC
