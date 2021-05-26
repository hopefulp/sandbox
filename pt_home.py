#!/home/joonho/anaconda3/bin/python
### written by J. Park and J. Lee
import argparse
import os
import re

Nserial = 8     # number of digits

def yes_or_no(question):
    reply = str(input(question+' (y/n): ')).lower().strip()     # raw_input is renamed in v3.6
    print(reply)
    if re.match('y', reply):
        return True
    else:
        return False

def get_serial(dir1):
    dirs = os.listdir(dir1)
    #print(dirs)
    serials=[]
    for d in dirs:
        serial = d.split('-')[-1]
        if len(serial) == Nserial and serial.isdigit():
            serials.append(serial)

    if serials:
        max_serial = max(serials)
        s_num = int(serial)
        s_num += 1
        return  str(s_num).zfill(Nserial)
    else:
        return '1'.zfill(Nserial)

def data_moving(job, dir_home, dir_target, author, dftpackage, dftname, model, kpoints, run):
    if job == None:
        print(f"Usage:: python {__file__} -j [cp|mv|test] -a author -dp [vasp|siesta], -dft [pbe| ... etc ")
        print(f"        {__file__} -h for more info")

    dir1 = dir_home
    dir2 = dir_target + '/' + author + '-' + model + '-' + get_serial(dir_target)
    com  = job + ' ' + dir1 + ' ' + dir2        # cp dirname dirname
    #print(com)
    q = "will you run?: %s" % com
    if run or yes_or_no(q):
        os.system(com)

    return 0 
    

def main():

    parser = argparse.ArgumentParser(description="move data from /home, /home2 to /home3  ")
    parser.add_argument('-j','--job', choices=['cp', 'mv', 'test'], help="mv or cp for data moving ")
    parser.add_argument('-dh', '--dir_home', default='/home', choices=['/home', '/home2'], help="source directory")
    parser.add_argument('-dt', '--dir_target', default='/home3/DB_Material', help='target directory')
    parser.add_argument('-a', '--author', default='jjj', help='author name in 3-letters')
    parser.add_argument('-dp', '--dft_p', choices=['vasp', 'siesta'], help='dft package')
    parser.add_argument('-dft', '--dft_name', default='pbe', help='dft name')
    parser.add_argument('-m', '--model', default='pt', help='name for material system')
    parser.add_argument('-kp', '--kpoints', default='1 1 1', help='kpoints')
    parser.add_argument('-r', '--run', action='store_true', help="execute the command")
    args = parser.parse_args()

    data_moving(args.job, args.dir_home, args.dir_target, args.author, args.dft_p, args.dft_name, args.model, args.kpoints, args.run)

if __name__ == "__main__":
    main()
