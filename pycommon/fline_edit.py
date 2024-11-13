#!/usr/bin/python

import server_env
import argparse
import re
import os


def inf_print(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
        for line in lines:
            print(line)
        return 0


def modify_file(fname, jobtype, option):

    if jobtype == "python":
        kw = "#!"
        if option == 3:
            pyth = server_env.python
        else:
            pyth = "/usr/bin/python"
        one_line = "#!"+pyth+"\n"
        ind = 0
    elif jobtype == "extxyz":
        kw = "energy"

    with open(fname, "r+") as fp:
        lines = fp.readlines()
        if jobtype == "python":
            if re.match(kw, lines[ind]):
                lines[0] = one_line
                """ rewrite one line in shell
                com = "sed -i '1s[.*[%s[' %s" % (one_line, fname)
                os.system(com)
                """
            else:
                lines.insert(ind, one_line)
                # change file mode
                com = "chmod 755 %s" % fname
                os.system(com)
        elif jobtype == "extxyz":
            i = 0
            while i < len(lines):
                if re.search(kw, lines[i]):
                    line = lines[i].split()
                    for ele in line:
                        if re.match(kw, ele):
                            ene = ele
                            break
                    lines[i]=ele + "\n"
                    print(lines[i])
                else:
                    line = lines[i].split()
                    if 3 < len(line):
                        new_line = line[0]+" "+line[1]+" "+line[2]+" "+line[3]+" "+line[4]+"\n"
                        lines[i] = new_line
                i += 1
        # locate pointer to the 1st byte before write
        fp.seek(0)  
        fp.writelines(lines)

    return

def main():
    parser = argparse.ArgumentParser(description='substitution a line: sed - substitution 1st line, python - add 1st line')
    parser.add_argument('fname', help='input file')
    parser.add_argument('jobtype', choices=['python', 'extxyz'], help='job type')
    parser.add_argument('-o', '--option', default=3, help='options which depend on jobtype')
    args = parser.parse_args()

    modify_file(args.fname, args.jobtype, args.option) 
    return 0

if __name__ == '__main__':
    main()
