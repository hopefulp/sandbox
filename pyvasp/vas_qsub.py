from envvasp import get_hostname
from server_env import nXn
import sys
import re

hostname = get_hostname()

def qsub_command(ndir, X=3, nnode=4, np=40, nmpi=None):
    if hostname == 'mlet':
        s = f"qsub -N {ndir} -pe numa {np} -v np={np} -v dir={ndir} -v vas=gam $SB/pypbs/sge_vasp_exe.csh"
    elif hostname == 'kisti':
        nnode=20
        np=40
        if not nmpi:
            if re.search('bd', ndir) or re.search('band', ndir):
                nmpi = np/2
            else:
                nmpi = np
        if np == nmpi:
            s = f"qsub -N {ndir} $SB/pypbs/pbs_vasp_kisti_skl.sh"
        ### due to memory prob for band calculation of large sc, use nmpi = np/2
        else:
            s = f"qsub -N {ndir} -l select={nnode}:ncpus={np}:mpiprocs={nmpi}:ompthreads=1  $SB/pypbs/pbs_vasp_kisti_skl.sh"
    elif hostname == 'pt':
        nproc = nnode * nXn[X]
        s = f"sbatch -J {ndir} -p X{X} -N {nnode} -n {nproc} /home/joonho/sandbox/pypbs/slurm_sbatch.sh"
    else:
        print(f"No qsub command for {hostname}")
        s=''
        #sys.exit(1)
    return s        


