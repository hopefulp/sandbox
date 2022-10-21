#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import sys
from common import *

def make_newfname(fname, pattern, style, new_word):
    if style == 'ap':
        if re.search('\.',fname):
            fn = re.split('\.',fname)
            newf = fn[0] + new_word + '.' + fn[1]
        else:
            newf = fname + new_word
    elif style == 'rp':
        newf = fname.replace(pattern, new_word)
    return newf


def ch_fname(job, m_tag, pattern, Linverse, style, new_word, L_dir, exceptions,ex_opt, dir_out, Lparents,run):
    pwd = os.getcwd()
    # pattern should have 1 element
    print(pattern)
    ### pattern f: specify filename
    if m_tag != 'f':
        l_file = get_files_patterns(m_tag, pattern, pwd, Ldir=L_dir, Linverse=Linverse, Lparents = Lparents)
    else:
        l_file = pattern
    print("\n".join(l_file))
    #sys.exit(1)
    n=0
    if exceptions:
        ### using whole word or matching word
        if ex_opt == 'm':
            for fex in exceptions:
                ### exclude by matching
                for i in range(len(l_file)-1,-1,-1):
                    if re.match(fex, l_file[i]):
                        l_file.remove(l_file[i])
        # remove by whole name
        else:
            for fex in exceptions:
                l_file.remove(fex)         
    # run job for the file list
    #print("\n".join(l_file))
    print(len(l_file), " files are selected")
    commands=[]
    if job == "ls":
        if l_file:
            print(("\n".join(l_file)))
    elif job in "mvdir" :
        if not os.path.isdir(dir_out):
            comm = "mkdir " + dir_out
            os.system(comm)
        for f in l_file:
            comm = f'{job} {f} {dir_out}")'
            print(comm)
    elif job == 'rename':
        #if not new_name: 
        #    print("add new name with -a or ")
        #    exit(10)
        #else:
        for fname in l_file:
            ### 
            new_name = make_newfname(fname, pattern[0], style, new_word)
            #new_name = fname.replace(patt, new_patt)
            ### rename file with extension
            '''
            if re.search('\.', fname):
                f_name = re.split('\.', fname)
                #print f_name
                new_file=f_name[0]+new_name+'.'+f_name[1]
            else:
                new_file=fname+new_name
            '''
            if os.path.isfile(new_name):
                print(f"{new_name} will be overwritten. Stop!")
                sys.exit(1)
            comm  = "mv  " + fname + " " + new_name
            commands.append(comm)
    else:
        ### rewrite mode
        for f in l_file:
            if job == 'rm':
                com = job + " " + f
            elif job == 'cp':
                for d in dir_out:
                    com = 'cp %s %d' % (f, d)
            elif job == 'chmod':
                com = f'{job} {mode} {f}'
            #elif job == 'mv':
            #    com = f'mv {f} {dir_out}
            commands.append(com)

                
    
    #if(job == 'rm' or job == 'mvdir' or job == 'rename') and run == False:
    #    print "use '-r' for run"
    for com in commands:
        print(com)
    q = "will you run? "
    if yes_or_no(q):
        for com in commands:
            os.system(com)
    return                    


def main():
    parser = argparse.ArgumentParser(description='Command Line Interface to deal with directory')
    parser.add_argument( 'job', choices=['ls', 'mvdir', 'rm', 'rename', 'cp', 'chmod', 'mv'],  help='shell command')
    group = parser.add_mutually_exclusive_group()
    group.add_argument( '-p', '--prefix', nargs='*', help='prefix of filename')
    group.add_argument( '-s', '--suffix', nargs='*', help='list several suffixes')
    group.add_argument( '-m', '--match', nargs='*', help='find matching string')
    group.add_argument( '-i', '--infiles', nargs='*', help='input file list')
    parser.add_argument( '-v', '--inverse', action='store_true', help='after find matching, inverse the selection')
    parser.add_argument('-st', '--style', choices=['ap', 'rp', 'mo'], help='fname changing style: append, replace, mode')
    parser.add_argument( '-rw', '--replace_word', help='string for replacement')
    #change.add_argument('-a', '--append', help="add suffix by -a to the original filename without extension")
    #change.add_argument('-mo', '--mode', default='755', choices=['755','644'], help="input chmod")
    parser.add_argument( '-id', '--include_dir', action='store_true', help='include dirname to filename')
    parser.add_argument( '-ip', '--include_parents', action='store_true', help='include dirname to filename')
    parser.add_argument( '-e', '--excluded', nargs='*', help='filename to be excluded')
    parser.add_argument( '-eo', '--excluded_opt', default='m',help='excluded fname option: matching or fullname')
    parser.add_argument( '-d', '--directory', type=str, default='tmppy', help='target directory to move files')
    parser.add_argument( '-r', '--run', action='store_true', help='run or not-False')
    args = parser.parse_args()

    ### select matching style
    if args.prefix:
        matching=args.prefix
        m_tag = 'p'
    elif args.suffix:
        matching=args.suffix
        m_tag = 's'
    elif args.match:
        matching=args.match
        m_tag = 'm'
    ### specify filename
    elif args.infile:
        matching=args.infiles
        m_tag = 'f'
    else:
        print("matching should be given")
        return 1

    #if args.job == "rename":
    #    if args.append:
    #        new_name = args.append
    #        rn_type = 'a'   # for append except extension
    #cmd(args.job, m_tag, matching, args.excluded, args.directory,args.rename,rn_type,new_name,args.run)
    ch_fname(args.job, m_tag, matching, args.inverse, args.style, args.replace_word, args.include_dir, args.excluded, args.excluded_opt, args.directory,args.include_parents,args.run)

if __name__ == "__main__":
    main()
