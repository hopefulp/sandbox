#!/home/noische/Enthought/Canopy_64bit/User/bin/python

import sys
import os
import re
import datetime
import math
import glob
import nutils as nu

# parse LAMMPS log file
logfile = sys.argv[1]
pat1 = re.compile(r"Step[-= ]*(\d*)")
pat2 = re.compile(r"CPU[-= ]*(\d*.\d*)")

result_log = os.popen("tail -n 100 " + logfile).read()

elapsed_step = re.findall(pat1, result_log)
elapsed_step = [int(i) for i in elapsed_step if i]
elapsed_step = elapsed_step[-1]

elapsed_time = re.findall(pat2, result_log)
elapsed_time = [float(i) for i in elapsed_time if i]
elapsed_time = elapsed_time[-1]

restart_step = 0

# parse LAMMPS input file
infile = glob.glob('in.*')
if infile:
    in_log = os.popen("grep -H run in*").read()
    if in_log == "":
        # if there is no run command, follow include link
        include_log = os.popen("grep -h include in*").read()
        _ = include_log.replace('include', '').replace(' ', '').split('\n')
        _ = [i for i in _ if i != ''] 
        __ = ''
        for i in _:
            __ += i + ' '   
        in_log = os.popen('grep -H run ' + __).read()

    pat = re.compile(r"(\S*in\.\S*):run\s*(\d*)")
    requested = sorted(re.findall(pat, in_log), key=lambda x:int(x[1]))
    max_requested = requested[-1]

    if len(infile) == 1:
        with open(infile[0]) as f:
            _ = f.readlines()
            #rst = [i for i in _ if "read_restart" in i][0]
            rst = [i for i in _ if "read_restart" in i]
        try:
            rst = rst[0]
        except:
            pass;
        else:
            pat_restart_step = re.compile(r"\.(\d+)\.")
            restart_step = int(sorted(re.findall(pat_restart_step, rst), key=lambda x:int(x[0]))[0])
            print("- read_restart found in the LAMMPS input script: %d" % restart_step)
            elapsed_step -= restart_step

    try:
        target_step = raw_input("- Enter the requested timestep [%d in %s]: " % (int(max_requested[1]), max_requested[0])) or int(max_requested[1])
    except ValueError:
        nu.die("Wrong value. Quit.")
else:
    try:
        target_step = raw_input("- Enter the requested timestep: ")
    except ValueError:
        nu.die("Wrong value. Quit.")

if target_step == elapsed_step:
    print("ERROR: Job already finished.")
    sys.exit(0)
elif target_step == "":
    sys.exit(0)
else:
    target_step = int(target_step)

now_time = datetime.datetime.now()
job_time = int(target_step * elapsed_time / elapsed_step)
left_time = job_time - elapsed_time
left_time_min = int(math.ceil(left_time/60.0))
left_time_hour = float(left_time/3600.0)
end_time = now_time + datetime.timedelta(seconds=left_time)

if left_time < 0:
    nu.die("Job already finished.")

print("Now on: %d step in %10.2f CPU time" % (elapsed_step, elapsed_time))
print("ETA: %s (about %7.0f seconds / %d minutes / %5.2f hours left)\n" % (end_time, left_time, left_time_min, left_time_hour))
