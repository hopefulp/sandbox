#!/home/joonho/anaconda3/bin/python

import argparse
import subprocess 
import re
import os

kw_jid = "Job Id"
kw_fin_jid = "Submit_Host"
kw_uid  = 'Job_Owner'
user = os.environ['USER']

def qstat(options):
    
    cmd = f"qstat -{options}"
    output = subprocess.check_output(cmd, shell=True)
    kws   = [] # the feature names
    vals   = [] # the jobs including all features
    job    = [] # the features
    haskw = 0  # whether I have the extra 'state' kw
    lens   = [] # str length of different field
    states = {} # counting jobs with different state
    #print(output.decode('utf-8').split("\n"))
    print(f"{'Jobid':^10}{'User':^10}{'State'}    {'Job Name':20}")
    Ljobid      = 0
    job_block   = []
    for i, line in enumerate(output.decode('utf-8').split("\n")):
        if Ljobid == 0:
            if re.match(kw_jid, line):
                Ljobid = 1
                job_block.append(line)   # 1st line
                #print(f"{i:>10}: Ljobid = {Ljobid}")
        ### on if Lsave == 'ON'
        ### need condition for 'OFF'
        elif Ljobid == 1:
            job_block.append(line)
            if re.search(kw_uid, line):
                ### set 0 in case not my job to restart
                if not re.search(user, line):
                    Ljobid = 0
                    job_block = []          # initialize again
            ### set 0 in case my job to restart
            if re.search(kw_fin_jid, line):
                Ljobid = 0
                ### Analize Job-Block here & print
                for l in job_block:
                    if "Job Id" in l:
                        jid = re.split('[ :.]+', l.strip())[2]
                    elif "Job_Owner" in l:
                        uname = re.split('[ @=]+',l.strip())[1]
                    elif "Job_Name" in l:
                        jname = re.split('[ =]+', l.strip())[1]
                    elif "job_state" in l:
                        jstate = re.split('[ =]+', l.strip())[1]
                print(f"{jid:>10}{uname:>10}{jstate:^5}    {jname:<25}")

                #print(f"{i:>10}: job_block ended")
    '''    
        if m1:
            # Get the state
            state = m1.group(1)
            job = [state]
            # Count
            if state not in states:
                states[state] = 0
            states[state] += 1
            # Add the job
            vals.append(job)

            # If 'state' kw is not added
            if haskw == 0:
                kws = ['state']
                haskw = 1
            # If it's already added
            elif haskw == 1:
                haskw = 2
        elif m2:
            # The feature name and value
            title = m2.group(1)
            value = m2.group(2)
            # Make sure 'state' is the last feature
            col   = 0 if len (job) == 1 else -1
            # job is a reference, referred to the job in vals
            job.insert (col, value)

            if haskw == 1:
                # Make sure 'state' is the last kw
                col = 0 if len (kws) == 1 else -1
                kws.insert (col, title)

    # calculate the length
    for i, kw in enumerate(kws):
        klen = len (kw)
        vlen = max ([len(x[i]) for x in vals])
        lens.append (max(klen, vlen))

    print("  ".join ([x.upper().ljust(lens[i]) for i,x in enumerate(kws)]))
    print("+".join (['-'*(lens[i]+1) for i,_ in enumerate(kws)]))
    for val in vals:
        print("  ".join ([x.ljust(lens[i]) for i,x in enumerate(val)]))

    print("\nTotal jobs: %s%s" % (len(vals), "".join([", " + k + ": " + str(v) for k,v in list(states.items())])))
    '''
def main():
    parser = argparse.ArgumentParser(description="qstat with option")
    parser.add_argument('-o', '--qst_opts', default='f',  help="the 1st INCAR file ")
    args = parser.parse_args()

    qstat(args.qst_opts) 

if __name__ == "__main__":
    main()
