#!/home/joonho/anaconda3/bin/python
'''
Make amp_run.py job string
Filenames for consistency in several modules
'''
### Filenames used in several module: 

score_file      = 'score.dat'
fitness_onegen  = 'ga_gen_fit.txt'
train_finish    = 'train_finished'

from ase.calculators.calculator import Parameters

import re
import numpy as np
from common import get_digits_4str, whereami


### Imported by amp_wrapper.py, gaamp/__init__.py 
class Ga_string:
    '''
    Make amp_run.py job string
    jsubmit = ['qsub', 'node']
    flist: force list = [force rmse, force coefficient for weight w.r.t. energy], f_maxres is not included
    hack_force: making force rmse per image, force file for a maxres image
    '''
    def __init__(self, jsubmit='qsub', script_f='mllorenz.py', queuename='test', mem='3G', ifname=None, ofname=None, mljob='trga',\
                ncore='1', hl='4', elimit='0.001', train_f='0.12 0.04', sftype='dellog', dtype='int1000', dlist=None,\
                add_ml_kw=None, max_iter=5000, ndtotal=4000, nam=None):
        p = self.parameters = Parameters()
        p.jsubmit       = jsubmit
        p.script_f      = script_f
        p.queuename     = queuename
        p.mem           = mem
        p.ifname        = ifname
        p.ofname        = ofname
        p.mljob         = mljob
        p.ncore         = ncore
        p.hl            = hl
        p.elimit        = elimit
        p.train_f       = train_f
        p.sftype        = sftype
        p.dtype         = dtype
        p.dlist         = dlist
        p.max_iter      = max_iter
        p.nam           = nam
        p.ndtotal       = ndtotal
        #print(f"{p.jsubmit}")
        ### overwrite keyword
        if add_ml_kw:
            for key in add_ml_kw.keys():       # add_amp_string.__dict__.keys()
                if  key in p.keys():                # p.__dict__.keys()
                    p[key] = add_ml_kw[key]    # use var as key in p[key], in p.key, key is new-key with name of 'key'
                    #print(f"p[{key}] = {p[key]}")
        if p.jsubmit == 'node':
            if re.search('tr', p.mljob):
                self.gastring = self.mkstr_tr()
            elif re.search('te', p.mljob):
                self.gastring = self.mkstr_te()

        elif p.jsubmit == 'qsub':
            ### write a script or run script
            Lwrite = 1
            if Lwrite:
                    self.qscript = self.queue_write_script()
                    self.gastring = self.queue_string()
            else:
                if re.search('tr',p.mljob):
                    self.gastring = self.make_quetrainstring()
                elif re.search('te', p.mljob):
                    self.gastring = self.make_queteststring()
    def __call__(self):
        return self.gastring

    def mkstr_tr(self):
        p   = self.parameters
        s   = ' tr '
        #s   += ' -inf ' + p.fname
        #s   += ' -j ' + p.mnjob
        #s   += ' -nc ' + p.ncore
        s   += f" -hl {p.hl}"
        #s   += ' -el '    + p.elimit
        #s   = ' -ms ' + p.ofname
        #if p.train_f:
        #    s   += f" -fl {p.train_f}"
        s   += ' -nd ' + str(p.ndtotal)
        s   += ' -mi ' + str(p.max_iter)
        #s   += self.symmetry_function()
        #s   += self.data_selection()
        comm = p.script_f + s + ' -dbp 2  '
        return comm

    def mkstr_te(self):
        '''
        nc = 1 for test calculation
        data_selection should be same for training and test
        only for test:
            nam = number of atoms for energy scaling
        not used in test:
            gaussian symmetry function
        '''
        p = self.parameters
        s   = ' -inf ' + p.fname
        s   += ' -j ' + p.mljob
        if p.nam:
            s   += ' -nam ' + str(p.nam)
        if p.train_f:
            s   += f" -fl {p.train_f}"
        s   += self.data_selection()
        comm = "amp_run.py " + s + " -g &"
        return comm

def score(cost, ref=1.0e-6):
    logref = np.log10(ref)
    logcost = np.log10(cost)
    score = logcost - logref
    return 1./score

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="make amp jobs string")
    parser.add_argument('-j', '--job', choices=['node','qsub'], help='job in node or qsub')
    parser.add_argument('-a', '--append', help='additional argument')
    args = parser.parse_args()

    ampstr = Amp_string(jobname='chromo00', hl='3 5 7')
    print(ampstr.make_queuestring())
