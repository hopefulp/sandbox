#!/gpfs/home/joonho/anaconda3/bin/python

import argparse
import os
import re
from common import dir_files

### modify input-detail here
### AMP-NODE job
amp_data_type={'interval': '-nt 4000 -ntr 100 -dtype int -dl 1000 1100 1200',
                'interval-te': '-nt 4000 -ntr 100 -dtype int -dl 1000 1100',
                'division': '-nt 4000 -ntr 1000 -dtype div -dl 4 0'}
amp_gaussian_function_param={'log': ' -des gs -pf log10 -pmm 0.05 200 -pn 10 -pmod del', 'pow': ' -des gs -pf powNN -pn 5'}

amp_node_str={ 'noforce': ' -j tr -hl 4 -el 0.001 -fl -0.1  -nc 10',
                'force' : ' -j tr -hl 4 -el 0.001 -fl 0.1 0.04 -nc 10'
                }
amp_node_str_suff=' -hf -g'     # hf: hack force, -g: suppress drawing graph
### AMP-SERVER job
amp_qsub_str={ 'force' : "-pe numa 10 -l mem=3G -v fname=OUTCAR -v pyjob=tr -v hl=4 -v el=0.001 -v tef=force -v fl='0.1 0.04'",
               'noforce': "-pe numa 10 -l mem=3G -v fname=OUTCAR -v pyjob=tr -v hl=4 -v el=0.001 -v tef=force -v fl='-0.1' "
               }
amp_qsub_data_type={'interval':  "-v nt=4000 -v ntr=100 -v dtype=int -v dlist='1000 1100 1200'",
                    'division': "-v nt=4000 -v ntr=1000 -v dtype=div -v dlist='3 0'"
                    }
amp_qsub_gaussian_param={'log': "-v des=gs -v pf=log10 -v pmod=del -v pmm='0.05 200.0' -v pn=10",
                        'pow': "-v des=gs -v pf=powNN -v pn=5"
                        }
amp_qsub_str_suff=" /home/joonho/sandboxg/pypbs/sge_amp.csh"                        

### AMP-GET-QSUB job
amp_getqsub_str = { 'noforce': ' -j tr -hl 4 -el 0.001 -fl -0.1', 'force': ' -j tr -hl 4 -el 0.001 -fl 0.1 0.04 -nc 10'}
amp_getqsub_datatype = amp_data_type
amp_getqsub_symmfuncttype = amp_gaussian_function_param

def show_command(job, job_submit, qname, keyvalues, nodename, sftype, dtype):
    
    if keyvalues:
        kw1 = keyvalues[0]

    if job == 'amp':
        ### modify setting here
        settings = {'force':'force',
                    'dtype': 'interval',
                    'dtype-te': 'interval-te',
                    'gs'   : 'log'
                    }
        print(f"{job}:: in {job_submit}: w. {settings['force']} data:{settings['dtype']} gs:{settings['gs']}-scale ")
        if job_submit == 'node':
            print("Train:\namp_run.py -f OUTCAR", amp_node_str[settings['force']], amp_data_type[settings['dtype']], amp_gaussian_function_param[settings['gs']], amp_node_str_suff)
            print("Test:\namp_run.py -f OUTCAR -j te -nc 10", amp_data_type[settings['dtype-te']], amp_node_str_suff)
            print("\n=== New Method")
            print("amp_wrapper.py -js node [-te] [-c| &] ")
            print("\t -c: for check string but not to submit jobs")
            print("GA")
            print("ga_amp_run.py -js node -nch 8 -hl 10 -nn 10 -np 5 [-c] &")
        elif job_submit == 'qsub':
            print("Structure of script to run multiple amp jobs::\
            \n    ampdir.sh wrapper_subdir [| sh]; which runs\
            \n    amp_wrapper.py -js qsub -qn Npn${n}hl20 -sf NN${n} -hl 20 20 -dl 1000 1300 3500 3600 &; which runs\
            \n    make mlet_tr(te).csh and run qsub\
            \n    qsub -N HL44tr -pe numa 8 -l mem=12G mlet_tr.csh\
            \n\tmlet_tr.csh should have -f force to make force_derivative directory\
            ")
            print("AMP JOBS")
            print("    ampdir.sh job[wrapper wrapper_subdir grep sh etc]")
            print("\nAMP WRAPPER:: train and test with failed-convergence")
            print(f"\tamp_wrapper.py -js qsub [-j [te|tr|trte]] -qn HL44 [-hl 4 4] [-dl int int]  [-s] [-c] &")
            print(f"\tamp_wrapper.py -js qsub &  !!! all options in the script")
            print("\tamp_wrapper.py -js qsub -j te -qn check &")
            print("\tamp_wrapper.py -js qsub -j tr -qn NN9p2 -dl 360 720 &")
            print("\tamp_wrapper.py -js qsub -sf logpn10max200 -qn HL1010 -hl 10 10 -dl 1000 1100 3500 3600 &")
            print("\tqsub -N HL44tr -pe numa 8 -l mem=12G mlet_tr.csh")
            print("    Water molecule:\
                    \n\tamp_wrapper.py -f w1.extxyz -js qsub -qn wat1pn5 -sf NN5 -hl 10 10 -dl 0 300 900 1000 -nt 1000 &\
                    \n\tamp_wrapper.py -f w1.extxyz -qn wat1pn5 -sf NN5 -hl 10 10 -dl 0 700 900 1000 -nt 1000 &\
                    ")

            print("    DB folder")
            print("\tamp_wrapper.py -js qsub -j tr -qn hl4 -dl 0 400 &")
            print("    DB")
            print("\tdiramp.sh fp | sh")
            print("    TRAIN")
            print("\tampdir.sh wrapper")
            print("\tamp_wrapper.py -js qsub -qn HL1010 -hl 10 10 &")
            print("\t    : queue submit job HL44tr-HL44te consecutively")

            print("GA")
            print("ga_amp_run.py -js qsub -nch 12 -hl 7 -nn 15 -np 4")
        elif job_submit == 'getqsub':
            print(f"sge_amp.py -qj {qname} -m 3G", amp_getqsub_str[settings['force']], amp_getqsub_datatype[settings['dtype']], amp_getqsub_symmfuncttype[settings['gs']]) 
        ### this part is for common
        print("\nCLEAN")
        print("clean.py -w amp pbs -a -y")
        print("\nSTATISTICS")
        print("stat_check.py -st [e|f|ef]")
        if 'kw1' in locals():
            d_pre=kw1
        else:
            d_pre='ND'
        print("\nPLOT")
        print("    TRAIN & TEST")
        print("\tampdir.sh grep")
        print("\tampplot_dir.py -p NN -t 'Ndata 100 (Gs-pow)' ")
        print("\tampplot_dir.py -p NN -t 'Training Set: Gs-pow' -y tr -yd . Nd300hl2020 Nd500hl2020")
        print("\tampplot_dir.py -p NN -t 'Test Set    : Gs-pow' -y te -yd . Nd300hl2020 Nd500hl2020")
        print("    TEST")
        print(f"\tampplot_stat_dir.py -p {d_pre} -le mse -t 'Energy Training' -k hl")
        print("\t    in case dir scan is hl: -k hl")
        print(f"\tampplot_stat_dir.py -p {d_pre} -le mse -t 'Force (Energy-only training)' -yl '(eV/A)^2' Correlation")
        print(f"\tampplot_stat_dir.py -p {d_pre} -le maxres -t 'Force (Energy-only training)' -yl eV/A")
        print(f"\tampplot_stat_dir.py -p {d_pre} -le mse -t 'Energy (Force-Training)'")
        print(f"\tampplot_stat_dir.py -p {d_pre} -le mse -t 'Force (Force-Training)' -yl '(eV/A)^2' Correlation")
        print(f"\tampplot_stat_dir.py -p {d_pre} -le maxres -t 'Force (Force-Training)' -yl eV/A ")
        print("\tampplot_stat_dir.py -p {d_pre} -le maxres -t 'Force (Force-Training)' -yl \"F\$_{rmse}\$ (eV/A)\" \"F\$_{maxres}\$ (eV/A)\"")
        print("\tampplot_stat_dir.py -p {d_pre} -le maxres -xl Nparam -t 'Force (Force-Training)' -yl \"F\$_{rmse}\$ (eV/A)\" \"F\$_{maxres}\$ (eV/A)\"")
        print("\tampplot_stat_dir.py -p NN -le maxres -sd Ndata300 -f f -xl Nparam -t 'Ndata 300 (Gs-pow)' -yl \"F\$_{rmse}\$ (eV/A)\" \"F\$_{maxres}\$ (eV/A)\"")
        print("\tampplot_stat_dir.py -p nd100i -le mse -t 'Energy (ND=100)' -xl \"Data index\"")
        print("\tampplot_stat_dir.py -p nd100i -le maxres -t 'Force (ND=100)' -yl \"F\$_{rmse}\$ (eV/A)\" \"F\$_{maxres}\$ (eV/A)\"  -xl 'Data index'")
    elif job == 'qchem':
        if job_submit == 'qsub':
            print("Run Q-Chem")
            print("    N.B.:")
            print("\tv5.1 parallel runs only at skylake@node!!")
            if not nodename:
                print("    qsub_server.py qchem -qj CC6 -i 6-CC-NiFe.in -n 6 -m 6 -no skylake@node14")
            else:
                print(f"    qsub_server.py qchem -qj CC6 -i 6-CC-NiFe.in -n 6 -m 6 -no {nodename}")
            print("\tN.B.: qchem parallel 5.1p doesnot run at Node of (haswell)opt, (sandy)slet -> use skylake@node00")
            print("\t-n for nproc, -m for mem")
            print("NBO analysis")
            print("    qcout_nbo.py nbo -f 5-NiFe.out -a C O O Ni P P Fe P P N -g C O O")

    else:
        print("build more jobs")
    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="show command Amp/Qchem/ etc ")
    parser.add_argument('job', nargs='?', default='amp', choices=['amp','qchem'],  help="one of amp, qchem")
    parser.add_argument('-js','--job_submit', default='qsub', choices=['chi','qsub','getqsub', 'node'],  help="where the job running ")
    parser.add_argument('-qn', '--qname', default='amptest', help="queue name for qsub shown by qstat")
    parser.add_argument('-k', '--keyvalues', nargs='*', help='change a keyword in print')
    parser.add_argument('-no', '--nodename', help='if needed, specify nodename')
    ampg = parser.add_argument_group(title = 'AMP args')
    ampg.add_argument('-dt', '--data_type', default='int', choices=['int','div'], help="data selection type")
    ampg.add_argument('-ft', '--func_type', default='log', choices=['log','pow'], help="gaussian symmetry function parameter type")

    args = parser.parse_args()

    show_command(args.job,args.job_submit,args.qname,args.keyvalues, args.nodename, args.func_type,args.data_type)

if __name__ == "__main__":
    main()
