#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
import sys
from common import *

def make_newfname(fname, pattern, dic_writing):
    ### ap for append; rp for replace
    print(f"{fname} {pattern} {dic_writing}")
    style = dic_writing['style']
    if dic_writing['word']:
        ch_word = dic_writing['word']

    if  style == 'ap':
        ### rename root or just append suffix
        if dic_writing['opt'] == 'root' and re.search('\.',fname):
            fn = re.split('\.',fname)
            new_fname = fn[0] + dic_writing['word'] + '.' + fn[1]
        else:
            new_fname = fname + dic_writing['word']
    elif style == 'rp':
        if not 'ch_word' in locals():
            ch_word = ""
        new_fname = fname.replace(pattern, ch_word)
    return new_fname


def ch_fname_recursive(dir1, job, m_tag, pattern, Linverse, dic_writing, L_dir, exceptions,ex_opt, dir_out, Lparents,Lrecur, run):
    print(f"{dir1}")
    if not dir1:
        pwd = os.getcwd()
        dir1 = pwd
    print(f"####.... enter {dir1} directory")
    os.chdir(dir1)
    ch_fname(job, m_tag, pattern, Linverse, dic_writing, L_dir, exceptions,ex_opt, dir_out, Lparents,Lrecur, run)

    p_files = os.listdir('.')
    for f in p_files:
        ### without check not link-directory, this makes error
        if os.path.isdir(f) and not os.path.islink(f):
            ch_fname_recursive(f, job, m_tag, pattern, Linverse, dic_writing, L_dir, exceptions,ex_opt, dir_out, Lparents,Lrecur, run)
    os.chdir('..')
    print(f"####....... exit {dir1} directory")
    return 0


def ch_fname(job, m_tag, pattern, Linverse, dic_writing, L_dir, exceptions, ex_opt, dir_out, Lparents,Lrecur, Lrun):
    pwd = os.getcwd()

    # pattern should have 1 element
    print(pattern)
    ### pattern f: specify filename
    if m_tag != 'f':
        l_file = get_files_patterns(m_tag, pattern, pwd, Ldir=L_dir, Linverse=Linverse, Lparents = Lparents)
    else:
        l_file = pattern
    print(f"{l_file} in 1st selection: {len(l_file)} files" )
    #sys.exit(1)
    n=0
    if exceptions:
        ### using whole word or matching word
        print(f"there is exceptions {exceptions} with exopt {ex_opt}")
        if ex_opt == 'm':
            for fex in exceptions:
                ### exclude by matching
                for i in range(len(l_file)-1,-1,-1):
                    if re.match(fex, l_file[i]):
                        l_file.remove(l_file[i])
        ### exception with matching, for whole matching: l_file.remove(fex) 
        else:
            for fex in exceptions:
                l_file = [ x for x in l_file if fex not in x ]
    print(f"{l_file} in 2nd selection: {len(l_file)} files" )
                        
    # run job for the file list
    #print("\n".join(l_file))
    print(len(l_file), " files are selected")
    commands=[]
    if job == "ls":
        if l_file:
            print(("\n".join(l_file)))
    ### mv is move to dir
    elif job == "mv" :
        ### secure dir_out directory exists not to delete files
        if not os.path.isdir(dir_out):
            comm = "mkdir " + dir_out
            os.system(comm)
        for f in l_file:
            comm = f'{job} {f} {dir_out}'
            #print(comm)
            if not os.path.exists(f'{dir_out}/{f}'):
                commands.append(comm)
            else:
                print(f"there exists target file {dir_out}/{f}")

    elif job == 'rename' or job == 'cp':
        if job == 'rename':
            job = 'mv'
        for fname in l_file:
            ### 
            new_name = make_newfname(fname, pattern[0], dic_writing)
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
            comm  = f"{job}  " + fname + " " + new_name
            commands.append(comm)
    else:
        ### rewrite mode
        for f in l_file:
            if job == 'rm':
                com = job + " " + f
            elif job == 'cp':
                for d in dir_out:
                    com = 'cp %s %s' % (f, d)
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
    if Lrun or yes_or_no(q):
        for com in commands:
            os.system(com)
    return 0              


def main():
    parser = argparse.ArgumentParser(description='Command Line Interface to deal with directory')
    parser.add_argument( 'job', choices=['ls', 'mv', 'rm', 'rename', 'cp', 'chmod'],  help='shell command: ?mvdir')
    gmatch = parser.add_mutually_exclusive_group()
    gmatch.add_argument( '-p', '--prefix', nargs='*', help='prefix of filename')
    gmatch.add_argument( '-s', '--suffix', nargs='*', help='list several suffixes')
    gmatch.add_argument( '-m', '--match', nargs='*', help='find matching string')
    gmatch.add_argument( '-mt', '--match_type', default='m', choices=['p', 'm', 's'], help='find matching string in prefix, middle, suffix')
    parser.add_argument( '-ms', '--match_string', nargs='*', help='input matching string')
    parser.add_argument( '-i', '--infiles', nargs='*', help='input file list')
    parser.add_argument( '-v', '--inverse', action='store_true', help='after find matching, inverse the selection')
    grewrite = parser.add_argument_group(title='file name rewriting method')
    grewrite.add_argument('-rs', '--rewrite_style', default='rp', choices=['ap', 'rp', 'mo'], help='fname changing style: append, replace, mode')
    grewrite.add_argument( '-rw', '--rewrite_word', help='string for replacement, if not empty')
    grewrite.add_argument( '-ro', '--rewrite_option', help='which word to be replaced: use period or not')
    #change.add_argument('-a', '--append', help="add suffix by -a to the original filename without extension")
    #change.add_argument('-mo', '--mode', default='755', choices=['755','644'], help="input chmod")
    parser.add_argument( '-id', '--include_dir', action='store_true', help='include dirname to filename')
    parser.add_argument( '-ip', '--include_parents', action='store_true', help='include dirname to filename')
    parser.add_argument( '-e', '--excluded', nargs='*', help='filename to be excluded')
    parser.add_argument( '-eo', '--excluded_opt', help='default==m, excluded fname option: matching or fullname')
    parser.add_argument( '-d', '-nd', '--directory', type=str, default='tmppy', help='target directory to move files')
    parser.add_argument( '-re', '--recursive', action='store_true', help='try subdirectories')
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
    elif args.match_type:
        matching = args.match_string
        m_tag = args.match_type
    elif args.infile:
        matching=args.infiles
        m_tag = 'f'
    else:
        print("matching should be given")
        return 1

    dic_rewrite = {'style': args.rewrite_style, 'word':args.rewrite_word, 'opt': args.rewrite_option}

    pwd = os.getcwd()
    if args.recursive:
        ch_fname_recursive(pwd, args.job, m_tag, matching, args.inverse, dic_rewrite, args.include_dir, args.excluded, args.excluded_opt, args.directory,args.include_parents,args.recursive, args.run)
    else:
        ch_fname(args.job, m_tag, matching, args.inverse, dic_rewrite, args.include_dir, args.excluded, args.excluded_opt, args.directory,args.include_parents,args.recursive, args.run)

if __name__ == "__main__":
    main()
