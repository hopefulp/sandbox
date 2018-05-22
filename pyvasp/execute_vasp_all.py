#!/home/jackjack5/epd/bin/python


import argparse
import os




def allvasp():
	host = os.getenv('HOST')
	lists = os.listdir('.')
	for dir in lists:
		if not os.path.isdir(dir):
			continue
		os.chdir(dir)
		print 'now on '+dir
		if os.path.isfile('INCAR') and os.path.isfile('POSCAR') and not os.path.isfile('OUTCAR'):
			os.system('csh /qcfs/jackjack5/vasp/vaspenv.sh')
		os.chdir('..')
	
			
		


def main():
	parser = argparse.ArgumentParser(description='execution of all vaspallallall')
	args = parser.parse_args()
	allvasp()
	return



if __name__=='__main__':
	main()	


