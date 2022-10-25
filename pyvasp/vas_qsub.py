#!/home/joonho/anaconda3/bin/python

from envvasp import get_hostname
from server_env import nXn
import sys
import re
import os
import subprocess
from common import yes_or_no
hostname = get_hostname()

### Now this is being depricated
def run_vasp(dirname, qx, qN, np, issue=None):
    #print(f'qx {qx} and qN {qN} in run_vasp()')
    if get_hostname()=='pt' and ( not qx or not qN):
        qx, qN = get_queue_pt(qx=qx)

    s = qsub_command(dirname,X=qx,nnode=qN,np=np, issue=issue)
    print(s)
    if yes_or_no("Will you run"):
        os.system(s)
    return 0

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
    print(f"Idle nodes: {free_node}")
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
            return 1, free_node

def qsub_command(ndir, X=3, nnode=4, np=None, issue=None):
    if hostname == 'kisti':
        nnode=20
        np=40
        if issue != 'mem' and ( re.search('bd', ndir) or re.search('band', ndir)):
            print("Use -m for half nproc for memory problem")
        if issue == 'mem':
            hproc = int(np/2)
            ### due to memory prob for band calculation of large sc, use hmem = np/2
            s = f"qsub -N {ndir} -l select={nnode}:ncpus={np}:mpiprocs={hproc}:ompthreads=1  $SB/pypbs/pbs_vasp_kisti_skl.sh"
        elif issue == 'opt':
            s = f"qsub -N {ndir} $SB/pypbs/pbs_vasp_kisti_sklopt.sh"
        else:
            s = f"qsub -N {ndir} $SB/pypbs/pbs_vasp_kisti_skl.sh"
    elif hostname == 'pt':
        if not X or not nnode:
            qx, qN = get_queue_pt(qx=qx)
        if np:
            nproc = np
        else:
            nproc = nnode * nXn[X]
        if issue == 'mem':
            hproc = int(nXn[X]/2)
            s = f"sbatch -J {ndir} -p X{X} -N {nnode} -c {hproc} --export=hmem=1 /home/joonho/sandbox/pypbs/slurm_sbatch.sh"
        elif issue == 'opt':
            s = f"sbatch -J {ndir} -p X{X} -N {nnode} -n {nproc} /home/joonho/sandbox/pypbs/slurm_sbatch_vaspopt.sh"
        elif issue == 'sim':
            s = f"sbatch -J {ndir} -p X{X} -N {nnode} -n {nproc} /home/joonho/sandbox/pypbs/slurm_sbatch_sim.sh"
        else:
            s = f"sbatch -J {ndir} -p X{X} -N {nnode} -n {nproc} /home/joonho/sandbox/pypbs/slurm_sbatch.sh"
        print(s)
        if yes_or_no("Will you run"):
            os.system(s)
    else:
        print(f"No qsub command for {hostname}")
        s=''
        #sys.exit(1)
    return s        

if __name__ == '__main__':
    get_queue_pt()

