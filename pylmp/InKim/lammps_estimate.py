#!/home/noische/python

import sys
import os
import re
import datetime

usage = """
LAMMPS_estimateTime.py log_file
"""

# usage
if len(sys.argv) < 2:
    print usage
    sys.exit(0)

# parse LAMMPS log file
logfile = sys.argv[1]
pat1 = re.compile(r"Step[-= ]*(\d*)")
pat2 = re.compile(r"CPU[-= ]*(\d*.\d*)")

result_log = os.popen("tail -n 100 " + logfile).read()

elapsed_step = int(re.findall(pat1, result_log)[-1])
elapsed_time = float(re.findall(pat2, result_log)[-1])

# parse LAMMPS input file
in_log = os.popen("grep -h run in*").read()
pat3 = re.compile(r"run\s*\d*")
requested_step = [int(i.split()[-1]) for i in re.findall(pat3, in_log) if i.split()[-1] != 'run']
requested_step = max(requested_step)

# calculate estimated time
try:
    target_step = raw_input("- Enter the requested timestep in the LAMMPS input file [%d]: " % requested_step) or requested_step
except ValueError:
    sys.exit(0)

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
end_time = now_time + datetime.timedelta(seconds=left_time)

print("Estimated time left: " + str(leftTime) + " seconds (about " + str(int(math.ceil(leftTime/60.0))) + " minutes)")
print("Estimated End time: " + str(endTime))
