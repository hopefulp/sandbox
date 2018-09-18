#!/home/jackjack5/epd/bin/python

import argparse

# args: OUTCAR energy(sigma->0) -1

def find_energy(outcarfilename,phrase,index):
	f = open(outcarfilename)
	fsp = f.readlines()
	nfsp = len(fsp)
	for i in range(nfsp):
		line = fsp[-i-1]
		if phrase in line:
			parsese = line.split()
			e = float(parsese[index])
			break
	return e


def main():
	parser = argparse.ArgumentParser(description='vasp output energy splitter')
	parser.add_argument('outcar', type=str,  help='vasp outcar file: OUTCAR')
	parser.add_argument('phrase', type=str,  help='phrase for searching line: energy(sigma->0)')
	parser.add_argument('index', type=int,  help='wanted number index: -1')
	args = parser.parse_args()
#	print args.outcar
#	e = find_energy(args.outcar,args.phrase,args.index)
	print args.outcar, e
	return

if __name__=='__main__':
	main()
