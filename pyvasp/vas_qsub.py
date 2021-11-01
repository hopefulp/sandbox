from myvasp import get_hostname
import sys

hostname = get_hostname()

def qsub_command(ndir, np=24, qopt=None):
    if hostname == 'mlet':
        s = f"qsub -N {ndir} -pe numa {np} -v np={np} -v dir={ndir} -v vas=gam $SB/pypbs/sge_vasp_exe.csh"
    elif hostname == 'kisti':
        if qopt:
            s = f"qsub -N {ndir} $SB/pypbs/pbs_vasp_kisti_skl2.sh"
        else:
            s = f"qsub -N {ndir} $SB/pypbs/pbs_vasp.sh"
    else:
        print(f"No qsub command for {hostname}")
        s=''
        #sys.exit(1)
    return s        


