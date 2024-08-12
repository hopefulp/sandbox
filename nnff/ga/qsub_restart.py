import time
import glob
import os

def qsub(num_pop):
    done = glob.glob('*done')
    check = []  
    cnt = 0 
    for i in range(num_pop):
        check.append('pop{:02d}'.format(i))

    for i in range(len(done)):
        done[i] = done[i][:-4]
        check.remove(done[i])


    os.system('qstat > temp.txt')
    with open('temp.txt', 'r') as f:
        flines = f.readlines()
    for j in range(len(flines)):
        if 'w22pot.sh' in flines[j]:
            cnt += 1    
    if cnt == 0:    
        for i in check:
            os.chdir(i) 
            isscore = os.path.isfile('score')
            errfin = glob.glob('w22pot.sh.e*')
            if len(errfin) >= 5 :
                with open('score','w') as f:
                    f.write('99999999\n')
                os.system('cat HL score > pop_fit.txt')
                os.system('cat pop_fit.txt >> ../pop_fit.txt')
                os.remove('HL')
                with open('../{}done'.format(i),'w') as f:
                    f.write('/{}\ncalculation is done'.format(i))
                    os.chdir('../')
                     
            elif len(errfin) < 5:
                if isscore == False:
                    os.system('qsub -cwd w22pot.sh')
                    print('\n---------------------AMP RERUN at '+os.getcwd()+'--------------------\n')
                    os.chdir('../')
    

def csh(num_pop):
    for i in range(num_pop):
        os.chdir('./pop{:02d}'.format(i))
        with open('pid','r') as f:
            fline = f.readline()
            
        run_jobs = os.system('pidof w22pot.sh')
        isscore = os.path.isfile('score')
        if fline not in run_jobs:
            if isscore == False:
                os.system('csh w22pot.sh &')
                os.system("pidof w22pot.sh > pid")
                with open('pid','r') as f:
                    flines = f.readlines()
                flines = ' '.join(flines)
                flines = flines.split(' ')
                flines[-1]
                with open('pid','w') as f:
                    f.write(flines[-1]) 
                print('\n---------------------AMP RERUN at '+os.getcwd()+'--------------------\n')
        os.chdir('../')
