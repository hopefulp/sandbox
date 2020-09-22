#!/home/joonho/anaconda3/bin/python
'''
Make amp_run.py job string
Filenames for consistency in several modules
'''
### Filenames used in several module: 
ampdb   = ['amp-fingerprint-primes.ampdb', 'amp-neighborlists.ampdb', 'amp-fingerprints.ampdb', 'OUTCAR']
amp_amp = [ "amp.amp", "amp-untrained-parameters.amp" ]
ampout_chromosome_fitness       = "ga_fit.txt"              # used in amp_run.py
ampout_onegeneration_fitness    = 'ga_gen_fit.txt'
ampout_te_e_chk                 = 'test_energy.pkl'         # used in amp_run.py, ampplot_stat_dir.py
ampout_te_f_chk                 = 'test_force.touch'        # used in amp_util.py 
#ampout_images_frmse            = 'f_rmse.dat'
#ampout_image_maxforce          = "forces.dat"
from ase.calculators.calculator import Parameters

import re

### Imported by amp_wrapper.py, gaamp/__init__.py 
class Amp_string:
    '''
    Make amp_run.py job string
    jsubmit = ['qsub', 'node']
    flist: force list = [force rmse, force coefficient for weight w.r.t. energy], f_maxres is not included
    hack_force: making force rmse per image, force file for a maxres image
    '''
    def __init__(self, jsubmit='qsub', queuename='test', mem='3G', fname='OUTCAR', ampjob='trga', ncore='1', hl='4',
                       elimit='0.001', train_f='0.1 0.04', sftype='dellog10', dtype='int1000', dlist=None, add_amp_kw=None, 
                       max_iter=10000, nam=None):
        p = self.parameters = Parameters()
        p.jsubmit       = jsubmit
        p.queuename     = queuename
        p.mem           = mem
        p.fname         = fname
        p.ampjob        = ampjob
        p.ncore         = ncore
        p.hl            = hl
        p.elimit        = elimit
        p.train_f       = train_f
        p.sftype        = sftype
        p.dtype         = dtype
        p.dlist         = dlist
        p.max_iter      = max_iter
        p.nam           = nam
        if add_amp_kw:
            for key in add_amp_kw.keys():       # add_amp_string.__dict__.keys()
                if  key in p.keys():                # p.__dict__.keys()
                    p[key] = add_amp_kw[key]    # use var as key in p[key], in p.key, key is new-key with name of 'key'
                    #print(f"p[{key}] = {p[key]}")
        if p.jsubmit == 'node':
            if re.search('tr',p.ampjob):
                self.ampstring = self.mkstr_tr()
            elif re.search('te', p.ampjob):
                self.ampstring = self.mkstr_te()

        elif p.jsubmit == 'qsub':
            ### write a script or run script
            Lwrite = 1
            if Lwrite:
                if re.search('tr',p.ampjob):
                    self.qscript = self.queue_write_script('tr')
                    self.ampstring = self.queue_string('tr')
                elif re.search('te', p.ampjob):
                    self.qscript = self.queue_write_script('te')
                    self.ampstring = self.queue_string('te')
            else: 
                if re.search('tr',p.ampjob):
                    self.ampstring = self.make_quetrainstring()
                elif re.search('te', p.ampjob):
                    self.ampstring = self.make_queteststring()
    
    def mkstr_tr(self):
        p   = self.parameters
        s   = ' -inf ' + p.fname
        s   += ' -j ' + p.ampjob
        s   += ' -nc ' + p.ncore
        s   += f" -hl {p.hl}"
        s   += ' -el '    + p.elimit
        if p.train_f:
            s   += f" -fl {p.train_f}"
        s   += ' -mi ' + str(p.max_iter)
        s   += self.symmetry_function()
        s   += self.data_selection()
        return s

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
        s   += ' -j ' + p.ampjob
        if p.nam:
            s   += ' -nam ' + str(p.nam)
        if p.train_f:
            s   += f" -fl {p.train_f}"
        s   += self.data_selection()
        return s


    def mlet_script(self, job):
        s   =   '#!/usr/bin/csh\n'
        s   +=  '#$ -cwd\n'
        s   +=  '#$ -N amp\n'
        s   +=  '#$ -V\n\n'

        s   +=  f'set logfile = "$SGE_O_WORKDIR/../sge_{job}$JOB_ID"\n'
        s   +=  'set nproc = $NSLOTS\n'
        s   +=  'set sandbox = "$HOME/sandboxg"\n'
        s   +=  'setenv PYTHONPATH $sandbox/pyamp:$sandbox/pycommon:$sandbox/chem_mod\n'
        s   +=  'set PYTHON = "$HOME/anaconda3/bin/python"\n'
        s   +=  'set EXE = "$sandbox/pyamp/amp_run.py"\n'

        s   +=  'date >> $logfile\n'
        return s

    def queue_write_script(self, job):
        p = self.parameters
        s = self.mlet_script(p.ampjob)

        ### run amp_run.py
        st   = ' -inf ' + p.fname
        st   += ' -j ' + p.ampjob
        st   += ' -nc ' + p.ncore
        if job == 'tr':
            st   += f" -hl {p.hl}"
            st   += ' -el '    + p.elimit
        if p.train_f:
            st   += f" -fl {p.train_f}"
        if job == 'tr':
            st   += ' -mi ' + str(p.max_iter)
            st   += self.symmetry_function()
        st   += self.data_selection()
        st   += " -g "
        
        s   += f'set string = "{st}"\n\n'

        s   +=  '$PYTHON $EXE $string\n\n'

        ### write to logfile
        s   += f'echo "Job options: $string" >> $logfile\n'
        s   +=  'echo "HOSTNAME    JOB_NAME    NSLOTS(nproc)" >> $logfile\n'
        s   +=  'date >> $logfile'
        return s

    def queue_string(self, job):
        p = self.parameters

        s   =   'qsub '
        s   +=  ' -N '          + p.queuename
        s   +=  ' -pe numa '    + p.ncore
        s   +=  ' -l mem='      + p.mem
        if job == 'tr':
            s += ' mlet_tr.csh'
        elif job == 'te':
            s += ' mlet_te.csh'
        return s            

    def make_quetrainstring(self):
            p = self.parameters
            s    = self.__queuestring = ""
            s    = 'qsub '
            s   += ' -N '       + p.jobname
            s   += ' -pe numa ' + p.ncore
            s   += ' -l mem='   + p.mem
            s   += ' -v fname=' + p.fname
            s   += ' -v pyjob=' + p.ampjob
            s   += " -v hl='"   + p.hl  +   "' "
            s   += ' -v el='    + str(p.elimit) 
            if p.train_f:
                s   += " -v flist='" + p.train_f + "' "
            s   += self.symmetry_function()
            s   += self.data_selection()
            return s

    def make_queteststring(self):
        p = self.parameters
        s    = self.__queuestring = ""
        s    = 'qsub '
        s   += ' -N '       + p.jobname
        s   += ' -pe numa ' + str(p.nqueue)
        s   += ' -l mem='   + p.mem
        s   += ' -v fname=' + p.fname
        s   += ' -v pyjob=' + p.ampjob
        s   += " -v hl='"   + p.hl  +   "' "
        s   += ' -v el='    + str(p.elimit) 
        if p.hack_force:
            s   += ' -v hf=any'
        if p.train_f:
            s   += " -v flist='" + p.train_f + "' "
        s   += self.symmetry_function()
        s   += self.data_selection()
        return s

    def __call__(self):
        return self.ampstring

    def data_selection(self):
        p = self.parameters
        ndata = 4000                # now fp was calculated only for 4000 images
        if p.dlist:
            dlist = p.dlist
        else:
            dlist = [1000, 2500, 3500, 3600]
        nd_tr   = dlist[1] - dlist[0]
        str_data     = ""
        #if p.jsubmit == 'qsub':
        #    tmp = f' -v nt={ndata} -v ntr={nd_tr}'
        #else:
        tmp = f' -nt {ndata} -ntr {nd_tr}'
        str_data += tmp
        if 'int' in p.dtype:
            #if p.jsubmit == 'qsub':
            #    tmp = f" -v dtype=int -v dlist='{dlist[0]} {dlist[1]}'"
            #else:
            tmp = f" -dtype int -dl {dlist[0]} {dlist[1]}"
            if len(dlist) >= 3: tmp += f" {dlist[2]}"
            if len(dlist) == 4: tmp += f" {dlist[3]}"
            str_data += tmp
        return str_data

    def symmetry_function(self):
        p = self.parameters
        str_sf=""
        if 'log' in p.sftype:
            #if p.jsubmit == 'qsub':
            #    tmp = ' -v des=gs -v pf=log10'
            #else:
            tmp = ' -des gs -pf log10'
            str_sf += tmp
        if 'del' in p.sftype:
            #if p.jsubmit == 'qsub':
            #    tmp = ' -v pmod=del'
            #else:
            tmp = ' -pmod del'
            str_sf += tmp
        #if p.jsubmit == 'qsub':
        #    tmp = " -v pmm='0.05 200.0' -v pn=10"
        #else:
        tmp = " -pmm 0.05 200.0 -pn 10"
        str_sf += tmp
        return str_sf

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="make amp jobs string")
    parser.add_argument('-j', '--job', choices=['node','qsub'], help='job in node or qsub')
    parser.add_argument('-a', '--append', help='additional argument')
    args = parser.parse_args()

    ampstr = Amp_string(jobname='chromo00', hl='3 5 7')
    print(ampstr.make_queuestring())
