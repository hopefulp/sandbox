#!/home/joonho/anaconda3/bin/python
'''
    modified at 02/27/2025
    change one dir and by one option for work [...] -prefix, -suffix, -match
'''    
import argparse
import os
import re
import sys
from common         import *
#from nnff.amp_ini   import ampdb
from mod_vas        import vasf_ini, vasf_out

q_list=[]
def dir_clean(   pwd,      selection, works, subwork,  linux_job, Linc_dir, exclude,excl_fnames,new_dir,Lshowmatch,Lall_rm, Lyes, Loutf):
    print(f"{pwd} directory in {whereami()}")
    matches=[]
    f_list_all=[]
    if selection == 'p':
        f_list=get_files_prefix(works, pwd, Lshow=Lshowmatch, Ldir=Linc_dir)
    elif selection == 's':
        f_list=get_files_suffix(works, pwd, Lshow=Lshowmatch, Ldir=Linc_dir)
    elif selection == 'm':
        f_list=get_files_match(works, pwd,  Lshow=Lshowmatch, Ldir=Linc_dir)
    elif selection == 'w':
        for work in works:
            if work == 'qchem':
                pass
            elif work == 'vasp' or work =='nc':
                fkeep = vasf_ini[:]
                f_list=os.listdir(pwd)
                print(f"total: {f_list}, {len(f_list)} files")
                if Loutf:    # remove only files in vasf_default
                    matches = []
                    for f in f_list:
                        if f in vasf_out:
                            matches.append(f)
                else:
                    if work == 'nc':
                        fkeep=[]
                        exc_suff=['sh', 'py']
                        exc_sufile=get_files_suffix_list(exc_suff, f_list)
                        #print(f"exc_sufile {exc_sufile}")
                        if excl_fnames:
                            excl_fnames.extend(exc_sufile)
                        else:
                            excl_fnames = exc_sufile
                    #print(f"excl files {excl_fnames}")
                    if excl_fnames: 
                        for efile in excl_fnames: # intactable list outside
                            fkeep.append(efile)
                    #print(f"keeps {fkeep}")
                    for f in f_list:
                        if f not in fkeep:
                            matches.append(f)
                    #print(f"matches: {matches}")

                f_list_all = matches

            elif work == 'pbs':
                #matches=['\.e\d', '\.o\d', '\.pe\d', '\.po\d', 'PI', 'sge']
                #matches=['\.e\d', '\.o\d', '\.pe\d', '\.po\d', 'PI', '\.pbs']  #'^\d'
                matches=['\.e\d', '\.o\d', '\.pe\d', '\.po\d', 'PI']  #'^\d'
                f_list = get_files_match(matches, pwd, Lshow=Lshowmatch)
                f_list_all.extend(f_list)
                print(f'{f_list_all}')
            elif work == 'slurm':
                matches=['\.X\d', 'slurm-\d+\.out']
                pres=['\d']
                f_list = get_files_match(matches, pwd, Lshow=Lshowmatch)
                f_list2= get_files_prefix(pres, pwd, Lshow=Lshowmatch)
                f_list_all.extend(f_list)
                f_list_all.extend(f_list2)
            elif work == 'amp':
                if subwork == 'ini':
                    fsuffix = ['ampdb', 'amp-log.txt']
                    f_list = get_files_suffix(fsuffix, pwd, Lshow=Lshowmatch, Ldir=True)
                    f_list_all.extend(f_list)
                else:
                    if linux_job == 'ln':
                        for ampdir in ampdb:
                            if os.path.islink(ampdir):
                                f_list_all.append(ampdir)
                    else:
                        fmatches=['amp','pdf','dat', 'ga', 'GA' ]
                        ### if Lall_rm: remove directory also
                        #dmatches=['ch']
                        #f_list = get_files_match(dmatches, pwd, Lshowmatch, Ldir=Lall_rm)
                        #f_list_all.extend(f_list)
                        f_list = get_files_match(fmatches, pwd, Lshow=Lshowmatch)
                        f_list_all.extend(f_list)

            elif work == 'lmp':
                matches=['trj', 'log']
                f_list = get_files_match(matches, pwd, Lshowmatch)


    if selection == 'p' or selection =='s' or selection == 'm':
        if exclude:
            f_list=get_files_exclude(exclude, d)
        if f_list and excl_fnames:
            for f in excl_fnames:
                if f in f_list:
                    f_list.remove(f)
        if f_list:
            f_list_all.extend(f_list)

    ### Make directory for 'cp', 'mv'        
    if linux_job == 'mv' or linux_job == 'cp':
        ndir = pwd + '/' + new_dir
        if os.path.isdir(ndir):
            print(f"Dir {new_dir} exists")
        else:
            os.mkdir(new_dir)
            print(f"Dir {new_dir} was made")
    
    f_list_all.sort()
    if exclude:
        for f in exclude:
            f_list_all.remove(f)
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
        if Linc_dir:
            comm += " -r"
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
    parser.add_argument('dir1', nargs='?', default=os.getcwd(), help='input one directory')
    gselect = parser.add_mutually_exclusive_group()
    gselect.add_argument('-w', '--works', nargs='*', choices=['qchem','amp','vasp','pbs','slurm','lmp','nc'],help='remove depending on job')
    gselect.add_argument('-p', '--prefix', nargs='*', help='remove with prefix')
    gselect.add_argument('-s', '--suffix', nargs='*', help='remove with suffix')
    gselect.add_argument('-m', '--match', nargs='*', help='remove matching file')
    parser.add_argument('-sw', '--subwork', choices=['ini','ga'], default='ini', help='subwork condition for works')
    parser.add_argument('-j', '--job', default='rm', choices=['rm','mv','cp','ln'], help='how to treat files [rm|cp|mv]')
    parser.add_argument('-e', '--exclude', nargs='*', help='remove all files except list') 
    parser.add_argument('-ef', '--excluded_files', nargs='*', help='save this file') 
    parser.add_argument('-nd', '--new_dir', default='tmp', help='directory where files to move')
    parser.add_argument('-id', '--include_dir', action='store_true', help='include directory')
    parser.add_argument('-ms', '--match_show', action='store_true')
    parser.add_argument('-a', '--all_remove', action='store_true', help='remove all the files')
    parser.add_argument('-y', '--yes', action='store_true', help='execute command')
    parser.add_argument('-o', '--outf', action='store_true', help='remove default files')
    args = parser.parse_args()

    if args.works:
        selection = 'w'
        slist = args.works
    elif args.prefix:
        selection = 'p'
        slist = args.prefix
    elif args.suffix:
        selection = 's'
        slist = args.suffix
    elif args.match:
        selection = 'm'
        slist = args.match
    else:
        print("input on of these options: -w|-p|-s|-m")
        print("use %s -h for help" % os.path.basename(__file__))
        sys.exit(0)
    #if args.work == 'amp' and not args.excluded_files:
    #    args.excluded
    dir_clean(args.dir1, selection, slist, args.subwork, args.job, args.include_dir, args.exclude,args.excluded_files,args.new_dir,args.match_show,args.all_remove, args.yes, args.outf)
    return 0

if __name__ == '__main__':
    main()
