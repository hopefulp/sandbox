#!/usr/bin/python

import server_env
import argparse
import re
import os

def modify_file(fname, jobtype, option):

    if jobtype == "python":
        if option == 3:
            pyth = server_env.python
        else:
            pyth = "/usr/bin/python"
        one_line = "#!"+pyth+"\n"
        ind = 0

    with open(fname, "r+") as fp:
        lines = fp.readlines()
        if re.match("#!", lines[ind]):
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
        # locate pointer to the 1st byte before write
        fp.seek(0)  
        fp.writelines(lines)

    return

def main():
    parser = argparse.ArgumentParser(description='substitution a line: sed - substitution 1st line, python - add 1st line')
    parser.add_argument('fname', help='input file')
    parser.add_argument('jobtype', choices=['python'], help='job type')
    parser.add_argument('-o', '--option', default=3, help='options which depend on jobtype')
    args = parser.parse_args()

    modify_file(args.fname, args.jobtype, args.option) 
    return 0

if __name__ == '__main__':
    main()
