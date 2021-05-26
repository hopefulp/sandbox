#!/home/joonho/anaconda3/bin/python
### written by J. Park and J. Lee
import argparse
import os
import re
import json
from collections import OrderedDict
import subprocess

Nserial = 8     # number of digits of serial number

def make_info(dir1, info):
    fname = f"{dir1}/information.json"
    print(json.dumps(info, ensure_ascii=False, indent='\t'))
    with open(fname, 'w', encoding='utf-8') as make_file:
        json.dump(info, make_file,  ensure_ascii=False, indent='\t')
    return 0

def dir_search(dir_target, info):
    dirs = os.listdir(dir_target)
    print(f"Search {dir_target}")
    for d in dirs:
        print(f"search {d}")
    return 0

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
    prefix = 'SN'
    serials=[]
    for d in dirs:
        if len(d) >= 10:
            serial = d[2:10]
            if len(serial) == Nserial and serial.isdigit():
                serials.append(serial)

    if serials:
        max_serial = max(serials)
        s_num = int(serial)
        s_num += 1
        return prefix+str(s_num).zfill(Nserial)
    else:
        return prefix+'1'.zfill(Nserial)

def data_moving(job, dir_home, dir_target, Lrun, info):
    if job == None:
        print(f"Usage:: python {__file__} -j [cp|mv|info|search] -a author -dp [vasp|siesta], -dft [pbe| ... etc ")
        print(f"        {__file__} -h for more info")
        print(f"\tmake information file: {__file__} -j info ")
        print(f"\tsearch directory     : {__file__} -j se")
        print(f"\tmove   directory     : {__file__} -j mv -dh dirname")
        return

    if job == 'rev':

        command = "ls -al "+ dir_home

        home3 = subprocess.check_output([command], shell=True)
        home3 = str(home3)
        home3 = home3.split("->")[1].strip()
        home3 = home3.replace("\\n'", "")

        command = "rm "+ dir_home + "; mv "+ home3 + " "+ dir_home
        os.system(command)

        return 0


    if job == 'info':
        make_info(dir_home, info)
    elif job == 'search':
        dir_search(dir_target, info)
    ### move or copy
    else:
        coms = []
        dir1 = dir_home
        dir1ln = dir1           # + 'l'
        dir2 = dir_target + '/' + get_serial(dir_target) + '-' + info['name'] 
        ### you may change dirname in database with following 
        #dir2 +=
        com  = f"{job} {dir1} {dir2}"        # cp dirname dirname
        coms.append(com)
        if job == 'mv':
            com  = f"ln -s {dir2} {dir1ln}"
            coms.append(com)
        #print(com)
        q = "will you run?: %s" % coms[0]
        if Lrun or yes_or_no(q):
            for com in coms:
                os.system(com)

    return 0 

def author():

    dummy_string = ""
    command = "whoami"
    dummy_string = subprocess.check_output([command], shell=True)
    dummy_string = str(dummy_string)
    dummy_string = dummy_string.replace("b'", "")
    dummy_string = dummy_string.replace("\\n'", "")

    return dummy_string

def main():

    parser = argparse.ArgumentParser(description="move data from /home, /home2 to /home3  ")
    parser.add_argument('-j','--job', choices=['cp', 'mv', 'info', 'search', 'rev'], help="mv or cp for data moving ")
    parser.add_argument('-w', '--write', action='store_true', help='read/write information file')
    parser.add_argument('-dh', '--dir_home', default='/home', help="source directory")
    parser.add_argument('-dt', '--dir_target', default='/home3/DB_Material', help='target directory')
    parser.add_argument('-r', '--run', action='store_true', help="execute the command")
    info = parser.add_argument_group(title='information')
    info.add_argument('-a', '--author', default='', help='author name in 3-letters')
    info.add_argument('-y', '--year', default='20000101', help='data production year')
    info.add_argument('-dp', '--dft_package', choices=['vasp', 'siesta'], help='dft package')
    info.add_argument('-dft', '--dft_name', default='pbe', help='dft name')
    info.add_argument('-m', '--model', help='name for material system')
    info.add_argument('-kp', '--kpoints', help='kpoints')
    info.add_argument('-c', '--comment', help='add additional comment to explain data more')
    args = parser.parse_args()

    ### treat information file here
    file_info = OrderedDict()

    
    if args.author == '': dummy_author = author()    
    else                : dummy_author = args.author
        
    file_info["name"] = dummy_author
    file_info['year'] = args.year
    file_info['package'] = args.dft_package
    file_info['dft'] = args.dft_name
    file_info['model'] = args.model
    file_info['kpoints'] = args.kpoints
    file_info['comment'] = args.comment

    data_moving(args.job, args.dir_home, args.dir_target, args.run, file_info)

if __name__ == "__main__":
    main()
