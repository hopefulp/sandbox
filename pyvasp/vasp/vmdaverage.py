#!/usr/bin/python2.7
#heejin
import os
import sys
import re
import subprocess
import operator
import math

if len(sys.argv) > 1:
	nstep = int(sys.argv[1])
else:
	nstep = int(500)

pwd = "".join(os.getcwd().split('/')[-1:])

def getsh(command):
	pp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	stdout = pp.stdout
	result = []
	while 1:
		line = stdout.readline()
		if not line:
			break
		result.append(line)
	return result

#get temperature
tmin = float(getsh("grep TEBEG INCAR | awk '{print $3}'")[0].strip())
tmax = float(getsh("grep TEEND INCAR | awk '{print $3}'")[0].strip())
tt = (tmin + tmax)/2
print("-------------------------------------------")
print(" - Target Temperature : %6.1f ~ %6.1f" % (tmin, tmax))

#get output results
if os.path.isfile(pwd+'.log'):
	command = "grep E0 %s.log|tail -n %d" % (pwd, nstep)
elif os.path.isfile(pwd+'.out'):
	command = "grep E0 %s.out|tail -n %d" % (pwd, nstep)
else:
	command = "grep E0 OSZICAR|tail -n %d" % nstep
greposz = getsh(command)

#grepout = getsh("grep 'ETOTAL' OUTCAR| tail -n %d" % nstep)
grepvol = getsh("grep 'volume' OUTCAR| tail -n %d" % nstep)

#calc. average for samples within +-10% target temperature
sumt = float(0)
sume = float(0)
sumf = float(0)
sum0 = float(0)
sump = float(0)
sumk = float(0)
#sumo = float(0)
sumv = float(0)
nimg = float(0)
for i in range(len(greposz)):
	nt = float(greposz[i].strip().split()[2].replace('-.','-0.'))
	ne = float(greposz[i].strip().split()[4].replace('-.','-0.'))
	nf = float(greposz[i].strip().split()[6].replace('-.','-0.'))
	n0 = float(greposz[i].strip().split()[8].replace('-.','-0.'))
	np = float(greposz[i].strip().split()[12].replace('-.','-0.'))
	nk = float(greposz[i].strip().split()[14].replace('-.','-0.'))
	sumt = sumt + nt
	sume = sume + ne
	sumf = sumf + nf
	sum0 = sum0 + n0
	sump = sump + np
	sumk = sumk + nk

#for i in range(len(grepout)):
#	no = float(grepout[i].strip().split()[4])
#	sumo = sumo + no

for i in range(len(grepvol)):
	nv = float(grepvol[i].strip().split()[4])
	sumv = sumv + nv

avgt = sumt / len(greposz)
avge = sume / len(greposz)
avgf = sumf / len(greposz)
avg0 = sum0 / len(greposz)
avgp = sump / len(greposz)
avgk = sumk / len(greposz)
#avgo = sumo / len(grepout)
avgv = sumv / len(grepvol)

#	if ((nt > tt*0.99) and (nt < tt*1.01)):
#		nimg = nimg + 1
#		sumt = sumt + nt
#		sume = sume + ne
#		sumf = sumf + nf
#		sum0 = sum0 + n0

#avgt = sumt / nimg
#avge = sume / nimg
#avgf = sumf / nimg
#avg0 = sum0 / nimg

#calc. standard deviation and standard error
dt = float(0)
de = float(0)
df = float(0)
d0 = float(0)
dp = float(0)
dk = float(0)
#do = float(0)
dv = float(0)
nimg = float(0)
for i in range(len(greposz)):
	nt = float(greposz[i].strip().split()[2].replace('-.','-0.'))
	ne = float(greposz[i].strip().split()[4].replace('-.','-0.'))
	nf = float(greposz[i].strip().split()[6].replace('-.','-0.'))
	n0 = float(greposz[i].strip().split()[8].replace('-.','-0.'))
	np = float(greposz[i].strip().split()[12].replace('-.','-0.'))
	nk = float(greposz[i].strip().split()[14].replace('-.','-0.'))

	dt = (nt-avgt)*(nt-avgt) + dt
	de = (ne-avge)*(ne-avge) + de
	df = (nf-avgf)*(nf-avgf) + df
	d0 = (n0-avg0)*(n0-avg0) + d0
	dp = (np-avgp)*(np-avgp) + dp
	dk = (nk-avgk)*(nk-avgk) + dk

#for i in range(len(grepout)):
#	no = float(grepout[i].strip().split()[4])
#	do = (no-avgo)*(no-avgo) + do

for i in range(len(grepvol)):
	nv = float(grepvol[i].strip().split()[4])
	dv = (nv-avgv)*(nv-avgv) + dv

sdt = math.sqrt(df/len(greposz))
sde = math.sqrt(de/len(greposz))
sdf = math.sqrt(df/len(greposz))
sd0 = math.sqrt(d0/len(greposz))
sdp = math.sqrt(dp/len(greposz))
sdk = math.sqrt(dk/len(greposz))
#sdo = math.sqrt(do/len(grepout))
sdv = math.sqrt(dv/len(grepvol))

set = sdt / math.sqrt(len(greposz))
see = sde / math.sqrt(len(greposz))
sef = sdf / math.sqrt(len(greposz))
se0 = sd0 / math.sqrt(len(greposz))
sep = sdp / math.sqrt(len(greposz))
sek = sdk / math.sqrt(len(greposz))
#seo = sdo / math.sqrt(len(grepout))
sev = sdv / math.sqrt(len(grepvol))

#	if ((nt > tt*0.99) and (nt < tt*1.01)):
#		nimg = nimg + 1
#		dt = (nt-avgt)*(nt-avgt) + dt
#		de = (ne-avge)*(ne-avge) + de
#		df = (nf-avgf)*(nf-avgf) + df
#		d0 = (n0-avg0)*(n0-avg0) + d0

#sdt = math.sqrt(df/nimg)
#sde = math.sqrt(de/nimg)
#sdf = math.sqrt(df/nimg)
#sd0 = math.sqrt(d0/nimg)

#set = sdt / math.sqrt(nimg)
#see = sde / math.sqrt(nimg)
#sef = sdf / math.sqrt(nimg)
#se0 = sd0 / math.sqrt(nimg)

#print
print(" - Last %d images" % (len(greposz)))
print("-------------------------------------------")
#print "Total %d images from the last %d samples" % (nimg, len(greposz))
print("  T  = %9.4f (sd:%8.4f, se:%8.4f)" % (avgt, sdt, set))
print("  E  = %9.4f (sd:%8.4f, se:%8.4f)" % (avge, sde, see))
print("  E0 = %9.4f (sd:%8.4f, se:%8.4f)" % (avg0, sd0, se0))
print("  P  = %9.4f (sd:%8.4f, se:%8.4f)" % (avgp, sdp, sep))
print("  K  = %9.4f (sd:%8.4f, se:%8.4f)" % (avgk, sdk, sek))
#print "  Et = %9.4f (sd:%8.4f, se:%8.4f)" % (avgo, sdo, seo)
print("  V  = %9.4f (sd:%8.4f, se:%8.4f)" % (avgv, sdv, sev))
