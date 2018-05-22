#!/usr/bin/python2.7
#heejin
import os
import sys
import re

#usage description
if len(sys.argv)<5:
	print "Usage: [prefix] [{numbers}] [output file name]"
	print "ex: afmn3-dos-DOS 1 3 12 asdf"
	exit()

#assigns file names
filenames=[]
prefix=sys.argv[1]
outputfilename=sys.argv[-1]

for i in xrange(2, len(sys.argv)-1):
	filename=prefix+sys.argv[i]
	filenames.append(filename)
	
#total data list
totaldata=[]

#stores data into totaldata list
#totaldata[i][j][k]: a datum of ith file, jth line, kth column.
for dosfile in filenames:
	f=open(dosfile)
	data=[]
	while 1:
		line=f.readline()
		if not line:
			break
			f.close()
		r=re.compile("\s?(-?[0-9 ]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s+(-?[0-9]+.[0-9]+)\s?")
		a=r.match(line)
		if a!=None:
			databuffer=[float(a.group(1)), float(a.group(2)), float(a.group(3)), float(a.group(4)), float(a.group(5)), float(a.group(6)), float(a.group(7)), float(a.group(8)), float(a.group(9)), float(a.group(10)), float(a.group(11)), float(a.group(12)), float(a.group(13)), float(a.group(14)), float(a.group(15)), float(a.group(16)), float(a.group(17)), float(a.group(18)), float(a.group(19))]
			data.append(databuffer)
	totaldata.append(data)

#generates a string for output
outstr=""
for j in xrange(0, len(totaldata[0])):
	databuffer=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	for k in xrange(1, 19):
		for i in xrange(0, len(filenames)):
			databuffer[0]=totaldata[i][j][0]
			databuffer[k]+=totaldata[i][j][k]
	outstr+="%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n"%(databuffer[0], databuffer[1], databuffer[2], databuffer[3], databuffer[4], databuffer[5], databuffer[6], databuffer[7], databuffer[8], databuffer[9], databuffer[10], databuffer[11], databuffer[12], databuffer[13], databuffer[14], databuffer[15], databuffer[16], databuffer[17], databuffer[18], databuffer[1]+databuffer[3]+databuffer[5]+databuffer[7]+databuffer[9]+databuffer[11]+databuffer[13]+databuffer[15]+databuffer[17], databuffer[2]+databuffer[4]+databuffer[6]+databuffer[8]+databuffer[10]+databuffer[12]+databuffer[14]+databuffer[16]+databuffer[18])

#generates output file
output=open(outputfilename, 'w')
output.write('E-Ef\ts-up\ts-down\tpy-up\tpy-down\tpz-up\tpz-down\tpx-up\tpx-down\tdxy-up\tdxy-down\tdyz-up\tdyz-down\tdz2-up\tdz2-down\tdxz-up\tdxz-down\tdx2-up\tdx2-down\t'+outputfilename+'\t'+outputfilename+'\n')
output.write(outstr)
output.close()
