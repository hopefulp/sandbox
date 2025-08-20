
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
            
def mod_kpoints(kpoints, outf='KPOINTS.new', kmulti=3):
    newlines=[]
    with open(kpoints, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        ### read kp's line
        if i == 3:
            newkps = []
            eles = list(map(int, line.strip().split()))
            for k in eles:
                if k != 1:
                    newk = k * kmulti
                else:
                    newk = k
                newkps.append(newk)
            newline = "  ".join(map(str, newkps)) + "\n"
            newlines.append(newline)
        else:
            newlines.append(line)

    with open(outf, 'w') as f:
        f.writelines(newlines)
        