#!/usr/bin/python

import sys
import os
import string
import random
import time
import getopt
import glob
import multiprocessing, Queue
 
def worker_bad(jobs):
	while True:
		try:
			tmp = jobs.get(block=False)
		except Queue.Empty:
			break
 
def worker_good(jobs):
	while True:
		tmp = jobs.get()
		if tmp == None:
			break
 
def main(items, size, num_workers, good_test):
	if good_test:
		func = worker_good
	else:
		func = worker_bad
	jobs = multiprocessing.Queue()
	for i in range(items):
		jobs.put(range(size))
	workers = []
	for i in range(num_workers):
		if good_test:
			jobs.put(None)
		tmp = multiprocessing.Process(target=func, args=(jobs,))
		tmp.start()
		workers.append(tmp)
	for worker in workers:
		worker.join()
	return jobs.empty()
 
if __name__ == '__main__':
	workers = 4
	items = workers * 2
	size = 10000
 
	for good_test in [False, True]:
		passed = 0
		for i in range(100):
			passed += main(items, size, workers, good_test=good_test)
		print '%d%% passed (Good Test: %s)' % (passed, good_test)

