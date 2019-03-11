#!/usr/bin/env python

from fireworks.core.firework import FWAction, Firework, FireTaskBase
from fireworks import FileTransferTask, ScriptTask
from fireworks.utilities.fw_utilities import explicit_serialize
import os
import shutil
import glob
import string

__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2015, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Sept 8, 2015'

def check_aoforce(fname):
    nfreq = False
    with open(fname, 'r') as f:
        stop = False
        while not stop:
            line = string.split(f.readline())
            if len(line) > 0 and line[0] == '[FREQ]':
                stop = True
        for i in range(6):
            line = string.split(f.readline())
            freq = float(line[0])
            if freq  < 0.0:
                print 'CAUTION: IMAGINARY FREQUENCY AT %6.2f CM^-1 DETECTED' % freq
                nfreq = True
    if nfreq:
        return False
    else:
        return True




@explicit_serialize
class JobexTask(FireTaskBase):

    _fw_name = "JobexTask"
    required_params = ["wf_name"]

    def run_task(self, fw_spec):
        job_info_array = fw_spec['_job_info']
        prev_job_info  = job_info_array[-1]
        path2setup     = prev_job_info['launch_dir']
        for file in glob.glob(path2setup+'/*'):
            if string.split(file,'.')[-1] == 'json': continue
            shutil.copy(file, '.')
        ### print infos to file 
        ST = ScriptTask.from_str('rm JOB-*;touch JOB-' + self["wf_name"] + '-opt')
        Action = ST.run_task(fw_spec)
        ### run jobex
        ST = ScriptTask.from_str('jobex -ri -c 200 > jobex.out')
        Action = ST.run_task(fw_spec)
        return Action 

@explicit_serialize
class AoforceTask(FireTaskBase):

    _fw_name = "AoforceTask"
    required_params = ["wf_name"]

    def run_task(self, fw_spec):
        job_info_array = fw_spec['_job_info']
        prev_job_info  = job_info_array[-1]
        path2setup     = prev_job_info['launch_dir']
        for file in glob.glob(path2setup+'/*'):
            if string.split(file,'.')[-1] == 'json': continue
            shutil.copy(file, '.')
        ### print infos to file 
        ST = ScriptTask.from_str('rm JOB-*;touch JOB-' + self["wf_name"] + '-hessian')
        Action = ST.run_task(fw_spec)
        ### prepare controle file
        ST = ScriptTask.from_str('prep_aoforce')
        Action = ST.run_task(fw_spec)
        ### run aoforce
        ST = ScriptTask.from_str('aoforce > force.out')
        Action = ST.run_task(fw_spec)
        ### run Merz Kollman fit
        ST = ScriptTask.from_str('dscf > charge.out')
        Action = ST.run_task(fw_spec)
        ### tm2molden ###
        ST = ScriptTask.from_str('tm2molden < tm2molden.prot > /dev/null')
        Action = ST.run_task(fw_spec)
        ### check for negative frequencies
        noneg = check_aoforce('molden.input')
        return Action 


@explicit_serialize
class ExtractTask(FireTaskBase):

    _fw_name = "ExtractTask"
    required_params = ["fxyz", "fref"]

    def run_task(self, fw_spec):
        job_info_array = fw_spec['_job_info']
        prev_job_info  = job_info_array[-1]
        path2setup     = prev_job_info['launch_dir']
        #print path2setup
        shutil.copy(path2setup + '/' + self["fxyz"], ".")
        #print path2setup
        ST = ScriptTask.from_str('create_ref -c ' + self["fxyz"]+ ' -r ' + 
                self["fref"] +' -p ' + path2setup + ' > create_ref.out 2> create_ref.err')
        Action = ST.run_task(fw_spec)
        shutil.copy(path2setup + '/final.xyz', ".")
        return Action 

@explicit_serialize
class CopyTask(FireTaskBase):

    _fw_name = "CopyTask"
    #required_params = ["fxyz", "fref"]

    def run_task(self, fw_spec):
        job_info_array = fw_spec['_job_info']
        prev_job_info  = job_info_array[-1]
        path2setup     = prev_job_info['launch_dir']
        ### get the files
        for file in glob.glob(path2setup+'/*'):
            if string.split(file,'.')[-1] == 'json': continue
            shutil.copy(file, '.')

@explicit_serialize
class RemoteCopyTask(FireTaskBase):

    _fw_name = "CopyTask"
    required_params = ["server", "files"]

    def run_task(self, fw_spec):
        import paramiko
        job_info_array = fw_spec['_job_info']
        prev_job_info  = job_info_array[-1]
        path2setup     = prev_job_info['launch_dir']
        cwd            = os.getcwd()
        ### setup paramiko
        ssh = paramiko.SSHClient()
        ssh.load_host_keys(os.path.join("/home/johannes", ".ssh", "known_hosts"))
        ssh.connect(self['server'], username=self.get('johannes'),
                key_filename=os.path.join("/home/johannes", ".ssh", "id_rsa.pub"))
        sftp = ssh.open_sftp()
        ### get the files
        for f in self["files"]:
            sftp.get(os.path.join(path2setup,f),os.path.join(cwd,f))
 
