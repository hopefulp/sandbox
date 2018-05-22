#!/home/jackjack5/epd/bin/python

import argparse
import subprocess as subp
import re
import os
import sys

def cmd(job, m_tag, patterns, exceptions, dir_destination, new_name, run):
    str1 = 'ls '
    p = subp.Popen(str1, shell=True, stdout=subp.PIPE, stderr=subp.PIPE)
    out, err = p.communicate()
    all_files = re.split("\n", out)
    # make list of target, file_target=[]
    file_target=[]
    my_regex_l=[]
    n=0
    if prefixes and not suffixes:
        for prefix in prefixes:
            #my_regex = r"^" + prefix
            my_regex_l.append(prefix)
            #print my_regex
    elif suffixes and not prefixes:
        for suffix in suffixes:
            print suffixes
            my_regex = suffix + "$"
            my_regex_l.append(my_regex)

    for file in all_files:
        if prefixes and not suffixes:
            for prefix in prefixes:
                #prefix = r"^CO2"
                if re.match(prefix, file):
                    file_target.append(file)
        elif suffixes and not prefixes:
            for suffix in suffixes:
                if re.search(suffix, file):
                    file_target.append(file)
        elif prefixes and suffixes:
            prefix = prefixes[0]
            suffix = suffixes[0]
            if re.search(prefix, file) and re.search(suffix, file):
                file_target.append(file)
        elif matches:
            pass
    #print("\n".join(file_target))                       
    if exceptions:
        for fex in exceptions:
            file_target.remove(fex)
    # run job for the file list
    #print("\n".join(file_target))
    print len(file_target), " files are selected"
    if job:
        if job == "ls":
            if file_target:
                print("\n".join(file_target))
            else:
                if prefix:
                    print "No matching files then show all the files"
                print("\n".join(files))
        elif job == "mv" :
            if not os.path.isdir(dir_destination):
                ret = subp.call(['mkdir', dir_destination])
                if ret == 1:
                    sys.exit(1)
            for file in file_target:
                command = job + " " + file + " " + dir_destination
                if run:
                    subp.Popen(command, shell=True, stdout=subp.PIPE)
                else:
                    print command
        elif job == 'rename':
            if not prefixes or len(prefixes) != 1 or not new_name:
                print "cannot rename without 1 prefix and new_name"
                print "try rename -p prefix -rn new_name"
                exit(10)
            else:
                n_pre=len(prefixes[0])
                for file in file_target:
                    f_body = file[n_pre:]
                    f_name  = prefixes[0] + f_body
                    fn_name = new_name + f_body
                    command  = "mv   %-30s %-30s" % (f_name, fn_name)
                    if run:
                        subp.Popen(command, shell=True, stdout=subp.PIPE)
                    else:
                        print command
        elif job == 'rm':
            for file in file_target:
                command = job + " " + file
                if run :
                    subp.Popen(command, shell=True, stdout=subp.PIPE)
                else:
                    print command
    
    if(job == 'rm' or job == 'mv' or job == 'rename') and run == False:
        print "use '-r' for run"
    return                    

#def run_bool(v):
#    if v.lower() in ('yes', '1', 'y', 't', 'true'):
#        return True

def main():
    parser = argparse.ArgumentParser(description='directory management')
    parser.add_argument( 'job', choices=['ls', 'mv', 'rm', 'rename'],  help='shell command')
    group = parser.add_mutually_exclusive_group()
    group.add_argument( '-p', '--prefix', nargs='*', help='prefix of filename')
    group.add_argument( '-s', '--suffix', nargs='*', help='list several suffixes')
    group.add_argument( '-m', '--match', nargs='*', help='find matching string')
    parser.add_argument( '-xcl', '--excluded', nargs='*', help='filename to be excluded')
    parser.add_argument( '-d', '--directory', type=str, default='tmppy', help='target directory to move files')
    parser.add_argument( '-rn', '--rename', help='rename filename')
    parser.add_argument( '-r', '--run', action='store_true', help='run or not-False')
    args = parser.parse_args()

    if args.prefix:
        matching=args.prefix
        m_tag = 'p'
    elif args.suffix:
        matching=args.suffix
        m_tag = 's'
    elif args.matching:
        matching=args.matching
        m_tag = 'm'
    else:
        print "matching should be given"
        return 1

    cmd(args.job, m_tag, matching, args.excluded, args.directory,args.rename,args.run)

if __name__ == "__main__":
    main()
