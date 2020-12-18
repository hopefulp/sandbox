#!/home/joonho/anaconda3/bin/python
'''
Make amp_run.py job string
Filenames for consistency in several modules
'''
### Filenames used in several module: 
ampdb   = ['amp-fingerprint-primes.ampdb', 'amp-neighborlists.ampdb', 'amp-fingerprints.ampdb', 'OUTCAR']
amp_amp = [ "amp.amp", "amp-untrained-parameters.amp" ]
amptr_backup = [ "score.dat", "test_energy.dat", "test_energy.pkl", "test_force_stat.dat", "test_fstat_acc.pkl", "test_fstat_acc.txt"]
ampout_score       = 'score.dat'               # "ga_fit.txt" refer to output file from  amp_run.py
ampout_onegeneration_fitness    = 'ga_gen_fit.txt'
ampout_te_e_chk                 = 'test_energy.pkl'         # used in amp_run.py, ampplot_stat_dir.py
ampout_te_f_chk                 = 'test_fstat_acc.pkl'      # used in ampplot_stat_dir.py
amptrain_finish                 = 'train_finished'
#ampout_images_frmse            = 'f_rmse.dat'
#ampout_image_maxforce          = "forces.dat"
from ase.calculators.calculator import Parameters

import re
from common import get_digits_4str, whereami


### Imported by amp_wrapper.py, gaamp/__init__.py 
class Amp_string:
    '''
    Make amp_run.py job string
    jsubmit = ['qsub', 'node']
    flist: force list = [force rmse, force coefficient for weight w.r.t. energy], f_maxres is not included
    hack_force: making force rmse per image, force file for a maxres image
    '''
    def __init__(self, jsubmit='qsub', queuename='test', mem='3G', fname='OUTCAR', ampjob='trga', ncore='1', hl='4',
                       elimit='0.001', train_f='0.12 0.04', sftype='dellog', dtype='int1000', dlist=None, add_amp_kw=None, 
                       max_iter=5000, ndtotal=4000, nam=None):
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
        p.ndtotal       = ndtotal
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
                    self.qscript = self.queue_write_script()
                    self.ampstring = self.queue_string()
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
        comm = "amp_run.py " + s + " -g &"
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
        s   += ' -j ' + p.ampjob
        if p.nam:
            s   += ' -nam ' + str(p.nam)
        if p.train_f:
            s   += f" -fl {p.train_f}"
        s   += self.data_selection()
        comm = "amp_run.py " + s + " -g &"
        return comm


    def mlet_script(self):
        p = self.parameters
        s   =   '#!/usr/bin/csh\n'
        s   +=  '#$ -cwd\n'
        s   +=  '#$ -N amp\n'
        s   +=  '#$ -V\n\n'
        s   +=  'set dirname = `basename $SGE_O_WORKDIR`\n'
        if re.search('ga', p.ampjob):
            s   +=  f'set logfile = "$SGE_O_WORKDIR/../sge_{p.ampjob}$dirname"\n'
        else:
            s   +=  f'set logfile = "$SGE_O_WORKDIR/sge_{p.ampjob}$dirname"\n'
        s   +=  'set nproc = $NSLOTS\n'
        s   +=  'set sandbox = "$HOME/sandboxg"\n'
        s   +=  'setenv PYTHONPATH $sandbox/pyamp:$sandbox/pycommon:$sandbox/chem_mod\n'
        s   +=  'set PYTHON = "$HOME/anaconda3/bin/python"\n'
        s   +=  'set EXE = "$sandbox/pyamp/amp_run.py"\n\n'
        s   +=  'echo "HOSTNAME    JOB_NAME    NSLOTS(nproc)" >> $logfile\n'
        s   +=  'echo "$HOSTNAME    $JOB_NAME    $nproc " >> $logfile\n'

        return s

    def queue_write_script(self):
        p = self.parameters
        s = self.mlet_script()

        ### run amp_run.py
        st   = ' -inf ' + p.fname
        st   += ' -j ' + p.ampjob
        st   += ' -nc ' + p.ncore
        if re.search('tr', p.ampjob):
            st   += f" -hl {p.hl}"
            st   += ' -el '    + p.elimit
        if p.train_f:
            st   += f" -fl {p.train_f}"
        if re.search('tr', p.ampjob):
            st   += ' -mi ' + str(p.max_iter)
            st   += self.symmetry_function()
        st   += self.data_selection()
        st   += " -g "
        
        s   += f'set string = "{st}"\n'
        s   +=  'echo "Job options: $string" >> $logfile\n\n'
        s   +=  'date >> $logfile\n'

        s   +=  '$PYTHON $EXE $string\n\n'

        s   +=  'date >> $logfile'
        return s

    def queue_string(self):
        p = self.parameters

        s   =   'qsub '
        s   +=  ' -N '          + p.queuename
        s   +=  ' -pe numa '    + p.ncore
        s   +=  ' -l mem='      + p.mem
        if re.search('tr', p.ampjob):
            s += ' mlet_tr.csh'
        elif re.search('te', p.ampjob):
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
        ndata = p.ndtotal                # now fp was calculated only for 4000 images
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
        str_sf=" -des gs"
        print(f"sf type in {whereami()}: {p.sftype}")
        if 'log' in p.sftype:
            str_sf += ' -pf log10'
            if 'del' in p.sftype:
                tmp = ' -pmod del'
                str_sf += tmp

            if 'p' in p.sftype:        # number of parameters of eta
                word = 'p'
                if 'pn' in p.sftype:
                    word = 'pn'
                nparam = get_digits_4str(word, p.sftype)
            else:
                nparam = '10'
            if 'm' in p.sftype:
                word = 'm'
                if 'max' in p.sftype:
                    word = 'max'
                pmax = get_digits_4str(word, p.sftype)
            else:
                pmax = '200'
            tmp = f" -pmm 0.05 {float(pmax):5.1f} -pn {nparam}"
            str_sf += tmp
        elif 'NN' in p.sftype:
            str_sf += ' -pf powNN'
            nparam = re.sub(r'\D+', '', p.sftype)   # using substitution of word to null for one pair
            str_sf +=  f' -pn {nparam}'
            
        return str_sf

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="make amp jobs string")
    parser.add_argument('-j', '--job', choices=['node','qsub'], help='job in node or qsub')
    parser.add_argument('-a', '--append', help='additional argument')
    args = parser.parse_args()

    ampstr = Amp_string(jobname='chromo00', hl='3 5 7')
    print(ampstr.make_queuestring())
