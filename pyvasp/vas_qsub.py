#!/home/joonho/anaconda3/bin/python

from envvasp import get_hostname
from server_env import nXn
import sys
import re
import subprocess

hostname = get_hostname()

def get_queue_pt(qx=None):
    '''
    obtain empty queue and nodes: convert linux function to python
    '''
    s = 'pestat'
    popen = subprocess.Popen(s, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    #data = popen.read().strip()
    (stdoutdata, stderrdata) = popen.communicate()
    dataline = stdoutdata.decode('utf-8').split("\n")   # split is not working 
    #sline = '\n'.join(dataline)
    #print(f"{sline}") # nlen {len(dataline)}")
    free_node={}
    for line in dataline:
        linestrip = line.strip()
        if re.match('n', linestrip):
        #    print (f"{line}")
            ele = linestrip.split()
            if ele[1][:2] not in free_node.keys():
                free_node[ele[1][:2]] = 0
            if ele[2] == 'idle':
                print (f"{line}")
                free_node[ele[1][:2]] += 1
    print(f"{free_node}")
    if qx:
        qname = 'X' + str(qx)
        return qx, free_node[qname]
    else:
        if free_node['X5'] >= 2:
            return 5, free_node['X5'] 
        elif free_node['X4'] >= 3:
            return 4, free_node['X4']
        elif free_node['X3'] >= 4:
            return 3, free_node['X3']
        elif free_node['X2'] >= 20:
            return 2, free_node['X2']
        else:
            return free_node

def qsub_command(ndir, X=3, nnode=4, np=40, nmpi=None):
    if hostname == 'mlet':
        s = f"qsub -N {ndir} -pe numa {np} -v np={np} -v dir={ndir} -v vas=gam $SB/pypbs/sge_vasp_exe.csh"
    elif hostname == 'kisti':
        nnode=20
        np=40
        if not nmpi:
            if re.search('bd', ndir) or re.search('band', ndir):
                nmpi = int(np/2)
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

if __name__ == '__main__':
    get_queue_pt()

