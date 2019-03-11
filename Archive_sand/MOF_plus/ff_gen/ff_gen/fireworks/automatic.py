#!/usr/bin/env python
# -*- coding: utf-8 -*-


from fireworks import Firework, Workflow, FWorker, LaunchPad, ScriptTask, FileTransferTask
from fireworks.core.rocket_launcher import rapidfire
from ff_gen.fireworks.ffgen_tasks import JobexTask, AoforceTask, ExtractTask, CopyTask, RemoteCopyTask
import string

def info(wfname, fwname):
    return 'rm JOB-*;touch JOB-%s-%s &> /dev/null' %(wfname, fwname)

def define(fprot):
    return "define <" + fprot +" > define.prot 2> define.err"

def x2t(fxyz):
    return "x2t %s > coord 2> x2t.err" % fxyz

class automatic:

    def __init__(self, wfname, fxyz, fprot, reset = False):
        self.launchpad = LaunchPad()
        if reset:
            self.launchpad.reset('', require_password=False)
        self.wfname = wfname
        self.fxyz = fxyz
        self.fprot = fprot
        self.ftxyz = wfname +'.txyz'
        self.fkey = wfname +'.key'
        self.fref = wfname +'.hdf5'
        return

    def copy2cmcc(self):
        # copy to CMCC
        copy1 = FileTransferTask({'files':['../'+self.fxyz, '../'+self.fprot], 'dest':'ffgen/fireworks', 'mode':'rtransfer',
            'server':'cmcc'})
        copy =  Firework([copy1], name = 'copy to cmcc', spec={"pass_job_info":True,})
        return copy

    def TMsetup(self, parent):
        # setup TM
        setup0 = FileTransferTask({'files':['../'+self.fxyz, '../'+self.fprot], 'dest':'.', 'mode':'copy'})
        setup1 = ScriptTask.from_str(x2t(self.fxyz)) # x2t
        setup2 = ScriptTask.from_str(define(self.fprot)) # define
        setup3 = ScriptTask.from_str(info(self.wfname,'TM_setup'))
        setup = Firework([setup0,setup1,setup2,setup3], name = 'TM_setup', 
                spec={"_pass_job_info": True,'_fworker':'cmcc_front'}, parents = [parent])
        return setup

    def opt(self,parent):
        #optimization
        opt1 = JobexTask(wf_name = self.wfname) 
        opt  = Firework([opt1], name = 'opt', parents = [parent], 
                spec={"_pass_job_info":True,'_fworker':'cmcc_queue'})
        return opt

    def hessian(self,parent):
        # aoforce
        force1 = AoforceTask(wf_name = self.wfname)
        force  = Firework([force1], name = 'hessian', parents = [parent], 
                spec={"_pass_job_info":True,'_fworker':'cmcc_queue'})
        return force

    def extract(self,parent):
        # extract infos
        extract1 = ExtractTask(fxyz = self.fxyz, fref = self.wfname)
        extract2 = ScriptTask.from_str(info(self.wfname,'extract'))
        extract  = Firework([extract1, extract2], name = 'extract', parents = [parent], 
                spec={"_pass_job_info":True,'_fworker':'cmcc_front'})
        return extract

    def FFsetup(self,parent):
        # setup FF
        FF1 = CopyTask() 
        FF2 = ScriptTask.from_str(info(self.wfname,'FF_setup'))
        FF3 = ScriptTask.from_str('xyz2txyz -i final.xyz -o %s > /dev/null 2> xyz2txyz.err' % self.ftxyz )
        FF4 = ScriptTask.from_str('atomtyper -c %s -o %s > /dev/null 2> atomtyper.err' %(self.ftxyz, self.ftxyz))
        FF5 = ScriptTask.from_str('create_key -i %s -o %s -r %s > key.out 2> key.err' % (self.ftxyz, self.fkey, self.wfname))
        FF  = Firework([FF1,FF2,FF3,FF4,FF5], name = 'FF_setup', parents = [parent], 
                spec={"_pass_job_info":True,'_fworker':'cmcc_front'})
        return FF

    def FFfit(self,parent):
        # fit FF
        fit0 = ScriptTask.from_str(info(self.wfname,'FF_fit'))
        fit1 = RemoteCopyTask(server = 'cmcc', files = [self.fkey, self.ftxyz, self.fref])
        fit2 = ScriptTask.from_str('fit-cma -c %s -k %s -r %s > fit.out \
                2> fit.err' %(self.ftxyz, self.fkey, self.wfname))
        fit3 = ScriptTask.from_str('compare-freqs -c %s -k %s/opt.key -r %s > \
                /dev/null 2> compare-freqs.err' %(self.ftxyz, self.wfname, self.wfname))
        fit  = Firework([fit0,fit1,fit2,fit3], name = 'FF_fit', parents = [parent], 
                spec={"_pass_job_info":True,'_fworker':'merry'})
        return fit

    def __call__(self):
        copy = self.copy2cmcc()
        setup = self.TMsetup(copy)
        opt = self.opt(setup)
        force = self.hessian(opt)
        extract = self.extract(force)
        FF = self.FFsetup(extract)
        fit = self.FFfit(FF)
        workflow = Workflow([copy, setup, opt, force, extract, FF,fit], name = self.wfname)
        # store workflow and launch it locally
        #launchpad.add_wf(setup)
        self.launchpad.add_wf(workflow)
        rapidfire(self.launchpad, FWorker())
        return

