###
login shell: bash

qsub script is written by tcsh
    bash has problem

qchem
    qsub -N optm1 -v job=1-PP-A-opt /gpfs/home/joonho/sandbox_gl/pypbs/sge_qchem.tcsh
    tcsh
        /gpfs/opt/qchem/bin/qchem $job.in $job.out

