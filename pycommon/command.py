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

def show_command(job, job_submit, qname, sftype, dtype, keyvalues):
    
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
            print("Shell script to run multiple jobs")
            print(f"dir_cli.sh ampwrap hl 4 5 6 7 8 9 10  | sh")
            print("\nAMP WRAPPER:: train and test with failed-convergence")
            print(f"amp_wrapper.py -js qsub [-j [te|tr|trte]] -qn HL44 [-hl 4 4] [-dl int int]  [-s] [-c] &")
            print(f"amp_wrapper.py -js qsub &  !!! all options in the script")
            print("amp_wrapper.py -js qsub -j te -qn check &")
            print("amp_wrapper.py -js qsub -j tr -qn NN9p2 -dl 360 720 &")
            print("amp_wrapper.py -js qsub -sf logpn10max200 -qn HL1010 -hl 10 10 -dl 1000 1100 3500 3600 &")
            print("\tqsub -N HL44tr -pe numa 8 -l mem=12G mlet_tr.csh")
            print("DB folder")
            print("amp_wrapper.py -js qsub -j tr -qn hl4 -dl 0 400 &")
            print("DB")
            print("diramp.sh fp | sh")
            print("TRAIN")
            print("diramp.sh wrapper")
            print("amp_wrapper.py -js qsub -qn HL1010 -hl 10 10 &")
            print("\t : queue submit job HL44tr-HL44te consecutively")

            print("for GA")
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
        print(f"ampplot_stat_dir.py -p {d_pre} -le mse -t 'Energy Training' -k hl")
        print("    in case dir scan is hl: -k hl")
        print(f"ampplot_stat_dir.py -p {d_pre} -le mse -t 'Force (Energy-only training)' -yl '(eV/A)^2' Correlation")
        print(f"ampplot_stat_dir.py -p {d_pre} -le maxres -t 'Force (Energy-only training)' -yl eV/A")
        print(f"ampplot_stat_dir.py -p {d_pre} -le mse -t 'Energy (Force-Training)'")
        print(f"ampplot_stat_dir.py -p {d_pre} -le mse -t 'Force (Force-Training)' -yl '(eV/A)^2' Correlation")
        print(f"ampplot_stat_dir.py -p {d_pre} -le maxres -t 'Force (Force-Training)' -yl eV/A ")
        print("ampplot_stat_dir.py -p {d_pre} -le maxres -t 'Force (Force-Training)' -yl \"F\$_{rmse}\$ (eV/A)\" \"F\$_{maxres}\$ (eV/A)\"")
        print("ampplot_stat_dir.py -p {d_pre} -le maxres -t 'Force (Force-Training)' -yl \"F\$_{rmse}\$ (eV/A)\" \"F\$_{maxres}\$ (eV/A)\" -k hl -xl 'HL (2layers)'")
        print("ND100, x = 'data index'")
        print("ampplot_stat_dir.py -p nd100i -le mse -t 'Energy (ND=100)' -xl \"Data index\"")
        print("ampplot_stat_dir.py -p nd100i -le maxres -t 'Force (ND=100)' -yl \"F\$_{rmse}\$ (eV/A)\" \"F\$_{maxres}\$ (eV/A)\"  -xl 'Data index'")
    else:
        print("build more jobs")
    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="show command Amp/Qchem/ etc ")
    parser.add_argument('job', nargs='?', default='amp', choices=['amp'],  help="one of amp,  ")
    parser.add_argument('-js','--job_submit', default='qsub', choices=['chi','qsub','getqsub', 'node'],  help="where the job running ")
    parser.add_argument('-qn', '--qname', default='amptest', help="queue name for qsub shown by qstat")
    parser.add_argument('-dt', '--data_type', default='int', choices=['int','div'], help="data selection type")
    parser.add_argument('-ft', '--func_type', default='log', choices=['log','pow'], help="gaussian symmetry function parameter type")
    parser.add_argument('-k', '--keyvalues', nargs='*', help='change a keyword in print')

    args = parser.parse_args()

    show_command(args.job,args.job_submit,args.qname,args.func_type,args.data_type, args.keyvalues)

if __name__ == "__main__":
    main()
