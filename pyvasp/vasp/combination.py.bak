#!/usr/bin/python2.7
#heejin
import sys
import os
import re

#usage description
if len(sys.argv)<3:
	print "Usage: [a] [b] of aCb"
	sys.exit()

tot=int(sys.argv[1])
num=int(sys.argv[2])

def xcombinations(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xcombinations(items[:i]+items[i+1:],n-1):
                yield [items[i]]+cc

def xuniqueCombinations(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations(items[i+1:],n-1):
                yield [items[i]]+cc\

def xselections(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for ss in xselections(items, n-1):
                yield [items[i]]+ss

def xpermutations(items):
    return xcombinations(items, len(items))


totset = []
for j in range (tot):
	k = str(j + 1)
	totset.append(k)

for c in xuniqueCombinations(totset,num):
	print ' '.join(c)
