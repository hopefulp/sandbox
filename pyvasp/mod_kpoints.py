
import re

def find_kdirection(klines):
    ## as for 2 lines: gives one direction
    return 0

def read_kpoints(kpoints):
    with open(kpoints, 'r') as f:
        lines = f.readlines()
        npoints = int(lines[1].strip())
        klines = lines[4:]
        check_kline = 0
        nkline = 0
        for line in klines:
            if re.search('0.', line):
                if check_kline == 0:
                    check_kline = 1
                elif check_kline == 1:
                    check_kline = 0
                    nkline += 1
                continue
    total_nk = npoints * nkline
    return npoints, nkline
            
