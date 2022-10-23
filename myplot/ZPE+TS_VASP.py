#! /usr/bin/python3

### Update : 2019 / 03 / 08
# Min Jong Noh #
# starnmj@kaist.ac.kr

##################################################################################
# This is a script for ZPE & TS calculation to energy profile (HER, ORR, etc.)   #
##################################################################################

###
import os, sys, math   
###

#### THz Data mining in OUTCAR ####

os.system('grep THz OUTCAR > THz.txt')

###################################

#### Get a total number of lines in read-file ####

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i+1

TN_lines = file_len('THz.txt') # Total Number of Lines

f = open('THz.txt')
THz2meV = []
for i in range(TN_lines):
    line = f.readline()
    words = line.split()
    THz2meV.append((int(words[0]), 0.001*float(words[-2])))  # Serial, eV
f.close()

####################################################


#### Zero Point Energy (ZPE) ####
#### += 0.5 * Summation of Frequency Energy [eV] ####

ZPE = 0.0
for i in range(len(THz2meV)):
    E = 0.5 * THz2meV[i][1]
    ZPE += E
##################################

#### Entropy Term ####
#              x             
# T*dS  =  ----------  -  ln(1 - exp(-x))
#          [exp(x)-1]
# 
#  x  = energy / kT
#
######################
kB = 0.0000861733576020577 # eV K-1
T = 298.15 # Temeprature, K
kT = kB * T

TS = 0.0
for i in range(len(THz2meV)):
    x = THz2meV[i][1] / kT
    v1 = x / (math.exp(x) - 1)
    v2 = 1 - math.exp(-x)
    E = v1 - math.log(v2)
    TS += kT*E

#### Print Output ###
print('----------------------------------------------')
print('Zero Point Energy (ZPE) : %15.9f  eV' % ZPE)
print('Entropy Energy (TS)     : %15.9f  eV' % TS)
print('----------------------------------------------')
