#!/home/joonho/anaconda3/bin/python
'''
    change into getting only 1 dir
'''    
import argparse
import os
import re
import sys
from common import *
from amp_ini import ampdb
q_list=[]

def dir_clean(d,works,linux_job,prefix, suffix, matches, exclude,excl_fnames,new_dir,Lshowmatch,Lall_rm, Lyes):

    pwd = os.getcwd()
    if not d:
        d = pwd
    print(f"{d} directory in {whereami()}")
    matches=[]
    f_list_all=[]
    for work in works:
        if work == None:
            if prefix:
                f_list=get_files_prefix(prefix, d)
            if suffix:
                f_list=get_files_suffix(suffix, d)
            if matches:
                f_list=get_files_match(matches, d,Lshowmatch)
            if exclude:
                f_list=get_files_exclude(exclude, d)
            if f_list and excl_fnames:
                for f in excl_fnames:
                    if f in f_list:
                        f_list.remove(f)
        elif work == 'qchem':
            pass
        elif work == 'vasp':
            f_list=os.listdir(d)
            if excl_fnames:
                for f in f_list:
                    if os.access(d+'/'+f, os.X_OK):
                        excl_fnames.append(f)
                for efile in excl_fnames:
                    if efile in f_list:
                        f_list.remove(efile)

        elif work == 'pbs':
            #matches=['\.e\d', '\.o\d', '\.pe\d', '\.po\d', 'PI', 'sge']
            matches=['\.e\d', '\.o\d', '\.pe\d', '\.po\d', 'PI']
            f_list = get_files_match(matches, d, Lshow=Lshowmatch)
            f_list_all.extend(f_list)
        elif work == 'amp':
            if linux_job == 'ln':
                for ampdir in ampdb:
                    if os.path.islink(ampdir):
                        f_list_all.append(ampdir)
            else:
                fmatches=['amp','pdf','dat', 'ga', 'GA' ]
                ### if Lall_rm: remove directory also
                f_list = get_files_match(fmatches, d, Lshow=Lshowmatch)
                f_list_all.extend(f_list)
                #dmatches=['ch']
                #f_list = get_files_match(dmatches, d, Lshowmatch, Ldir=Lall_rm)
                #f_list_all.extend(f_list)
                
        elif work == 'lmp':
            matches=['trj', 'log']
            f_list = get_files_match(matches, d, Lshowmatch)

    ### Make directory for 'cp', 'mv'        
    if linux_job == 'mv' or linux_job == 'cp':
        ndir = pwd + '/' + new_dir
        if os.path.isdir(ndir):
            print(f"Dir {new_dir} exists")
        else:
            os.mkdir(new_dir)
            print(f"Dir {new_dir} was made")
    f_list_all.sort()
    ### Make command list
    for f in f_list_all:
        #fname = d+'/'+f
        if Lall_rm:
            comm = "%s -r %s" % (linux_job, f)
        elif linux_job == 'ln':
            comm = "rm %s\n" % f
            comm += "%s -s ../%s %s" % (linux_job, f, f)
        else:
            comm = "%s %s" % (linux_job, f)
        ### extra job for save files
        if linux_job == 'mv' or linux_job == 'cp':
            comm += " %s" % new_dir
        print(comm)
        q_list.append(comm)
        
    #print "all %s files" % len(f_list)
    ### Show command list and run
    ### change f_list to f_list_all
    if f_list_all:
        q = "will you %s %s files? " % (linux_job, len(f_list_all))
        if Lyes or yes_or_no(q):
            i = 0
            for comm in q_list:
                os.system(comm)
                i += 1
            if linux_job == 'rm':
                job_str = 'removed'
            elif linux_job == 'mv':
                job_str = 'moved'
            elif linux_job == 'cp':
                job_str = 'saved'
            elif linux_job == 'ln':
                job_str = 'removed and linked again'
            print(f"{i} files will be {job_str}")         

    return 0


def main():
    parser = argparse.ArgumentParser(description='to clean directory in qchem')
    parser.add_argument('dir1', default=os.getcwd(), help='input one directory')
    parser.add_argument('-w', '--works', nargs='+', choices=['qchem','amp','vasp','pbs','lmp'],help='remove depending on job')
    parser.add_argument('-j', '--job', default='rm', choices=['rm','mv','cp','ln'], help='how to treat files [rm|cp|mv]')
    parser.add_argument('-p', '--prefix', nargs='*', help='remove with prefix')
    parser.add_argument('-s', '--suffix', nargs='*', help='remove with suffix')
    parser.add_argument('-m', '--match', nargs='*', help='remove matching file')
    parser.add_argument('-e', '--exclude', nargs='*', help='remove all files except list') 
    parser.add_argument('-ef', '--excluded_files', nargs='*', help='save this file') 
    parser.add_argument('-jd', '--new_dir', default='tmp', help='directory where files to move')
    parser.add_argument('-ms', '--match_show', action='store_true')
    parser.add_argument('-a', '--all_remove', action='store_true', help='remove all the files')
    parser.add_argument('-y', '--yes', action='store_true', help='execute command')
    args = parser.parse_args()

    if args.works==None and args.prefix==None and args.suffix==None and args.match==None:
        print("input -w|-p|-s|-m")
        print("use %s -h for help" % os.path.basename(__file__))
        sys.exit(0)
    if 'vasp' in args.works and not args.excluded_files:
        args.excluded_files=['POSCAR','POTCAR','KPOINTS','INCAR']
    #if args.work == 'amp' and not args.excluded_files:
    #    args.excluded
    dir_clean(args.dir1,args.works,args.job,args.prefix,args.suffix,args.match,args.exclude,args.excluded_files,args.new_dir,args.match_show,args.all_remove, args.yes)
    return 0

if __name__ == '__main__':
    main()
