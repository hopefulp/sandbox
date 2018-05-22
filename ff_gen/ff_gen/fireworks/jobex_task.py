#!/usr/bin/env python

from fireworks.core.firework import FWAction, Firework, FireTaskBase
from fireworks import FileTransferTask, ScriptTask
from fireworks.utilities.fw_utilities import explicit_serialize
import os
import shutil
import glob

__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2015, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Sept 8, 2015'

@explicit_serialize
class JobexTask(FireTaskBase):

    _fw_name = "JobexTask"

    def run_task(self, fw_spec):
        job_info_array = fw_spec['_job_info']
        prev_job_info  = job_info_array[-1]
        path2setup     = prev_job_info['launch_dir']
        for file in glob.glob(path2setup+'/*'):
            shutil.copy(file, '.')
        #ST = ScriptTask.from_str('jobex -ri -c 200 > jobex.out')
        ST = ScriptTask.from_str('touch foo')
        Action = ST.run_task(fw_spec)
        return Action 
        #shutil.copy(path2setup+'/*', '.')
 #      os.system('echo hallo')



#        print('The name of the previous job was: {}'.format(prev_job_info['name']))
#        print('The id of the previous job was: {}'.format(prev_job_info['fw_id']))
#        print('The location of the previous job was: {}'.format(prev_job_info['launch_dir']))
#       print('The location of the previous job was: {}'.format(dir))

