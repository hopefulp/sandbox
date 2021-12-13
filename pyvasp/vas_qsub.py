from envvasp import get_hostname
from server_env import nXn
import sys

hostname = get_hostname()

def qsub_command(ndir, np=24, qopt=None, X=3, nnode=4):
    if hostname == 'mlet':
        s = f"qsub -N {ndir} -pe numa {np} -v np={np} -v dir={ndir} -v vas=gam $SB/pypbs/sge_vasp_exe.csh"
    elif hostname == 'kisti':
        if qopt:
            s = f"qsub -N {ndir} $SB/pypbs/pbs_vasp_kisti_skl2.sh"
        else:
            s = f"qsub -N {ndir} $SB/pypbs/pbs_vasp_kisti_skl.sh"
    elif hostname == 'pt':
        nproc = nnode * nXn[X]
        s = f"sbatch -J {ndir} -p X{X} -N {nnode} -n {nproc} /home/joonho/sandbox/pypbs/slurm_sbatch.sh"
    else:
        print(f"No qsub command for {hostname}")
        s=''
        #sys.exit(1)
    return s        


