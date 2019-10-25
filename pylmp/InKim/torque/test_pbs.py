from torque import *

pbs = PBS()
pbs.parse_qstat()
print(pbs)
#test = pbs.parse_qstat()
#print(test)