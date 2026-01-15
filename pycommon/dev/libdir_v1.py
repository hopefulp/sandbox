import os
import re
import sys
from common import *
import time

from dataclasses import dataclass

@dataclass
class CleanConfig:
    '''
    used for dir_clean
    '''
    selection: str        # 'p', 's', 'm', 'w'
    slist: list
    works: list
    shelljob: str
    exclude: list
    excl_fnames: list
    new_dir: str
    show_match: bool
    all_rm: bool
    yes: bool



def walk_dirs( root, dir_work, *, follow_symlinks=False, exclude_dirs=None, **work_kwargs ):
    """
    Generic recursive directory walker.
    Calls dir_work(path, **work_kwargs) in every directory.
    """
    if not root:
        root = os.getcwd()

    if exclude_dirs is None:
        exclude_dirs = set()

    root = os.path.abspath(root)

    print(f"#### enter {root}")

    # DO WORK IN THIS DIRECTORY
    dir_work(root, **work_kwargs)

    # RECURSE
    for name in os.listdir(root):
        path = os.path.join(root, name)

        if not os.path.isdir(path):
            continue
        if not follow_symlinks and os.path.islink(path):
            continue
        if name in exclude_dirs:
            continue

        walk_dirs(
            path,
            dir_work,
            follow_symlinks=follow_symlinks,
            exclude_dirs=exclude_dirs,
            **work_kwargs,
        )

    print(f"#### exit  {root}")

### >>> =========================================
### file selection by work should be defined
def select_pbs_files(path, config):
    matches=['\.e\d', '\.o\d', '\.pe\d', '\.po\d', 'PI', 'pbs']  #'^\d'
    f_list = get_files_match(matches, path, Lshow=config.show_match)
    print(f'{f_list}')
    return f_list

def select_slurm_files(path, config):
    '''
    pres: prefixes, slurm files start with digits (job number)
    '''
    matches=['\.X\d', 'slurm-\d+\.out']
    pres=['\d']         
    f_list = get_files_match(matches, path, Lshow=config.show_match)
    f_list2= get_files_prefix(pres, path, Lshow=config.show_match)
    f_list.extend(f_list2)
    return f_list

def select_vasp_files(path, config):
    f_list = []
    return f_list

def select_amp_files(path, config):
    f_list = []
    return f_list

def select_by_work(path, config):
    files = []
    for work in config.works:
        if work == 'vasp':
            files.extend(select_vasp_files(path, config))
        elif work == 'pbs':
            files.extend(select_pbs_files(path, config))
        elif work == 'amp':
            files.extend(select_amp_files(path, config))
        elif work == 'slurm':
            files.extend(select_slurm_files(path, config))
        else:
            print(f"{work} is not defined in {whereami()}() in {__file__}")
            sys.exit(100)
    return files
### <<<

def select_files(path, config):
    if config.selection == 'p':
        return get_files_prefix(config.slist, path)
    elif config.selection == 's':
        return get_files_suffix(config.slist, path)
    elif config.selection == 'm':
        return get_files_match(config.slist, path)
    elif config.selection == 'w':
        return select_by_work(path, config)

def policies():
    '''
    policies: vasp / pbs / slurm / amp / lmp / etc
    '''
    pass

def apply_excludes(files, excl_files):
    for f in excl_files:
        files.remove(f)
    return files

def build_commands(files, config):
    cmds = []
    for f in files:
        if config.shelljob == 'rm':
            cmds.append(f"rm {'-r ' if config.all_rm else ''}{f}")
        elif config.shelljob in ('mv', 'cp'):
            cmds.append(f"{config.shelljob} {f} {config.new_dir}")
        elif config.shelljob == 'ln':
            cmds.append(f"rm {f}; ln -s ../{f} {f}")
    return cmds

def execute_commands(cmds, config):
    if not cmds:
        return
    ### show commands: job number start from 1
    for i, com in enumerate(cmds):
        print(f"{i+1:>3d} {com}")
    if config.yes or yes_or_no(f"Execute {len(cmds)} commands?"):
        for cmd in cmds:
            os.system(cmd)

### DIR JOB: clean from above
def dir_clean(path, config):
    print(f"[CLEAN] {path}")

    files = select_files(path, config)
    if config.excl_fnames:
        files = apply_excludes(files, config.excl_fnames)

    cmds = build_commands(files, config)
    execute_commands(cmds, config)

    return 0
###### DIR JOB: touch ###########################################################
def dir_touch(path):
    print(f"[TOUCH] {path}")

    for name in os.listdir(path):
        fpath = os.path.join(path, name)
        if os.path.isfile(fpath) and not os.path.islink(fpath):
            os.utime(fpath, None)


