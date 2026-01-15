"""
Directory operations library.

Supports:
- recursive directory walking
- clean policies (pbs, slurm, vasp, amp, ...)
- touch job

This can be split into: walker.py, policies.py, clean.py



Adding Policy example in DIR_CLEAN
    def select_cp2k_files(path, config):
        return get_files_match(['*.out', '*.wfn'], path)

    POLICIES['cp2k'] = select_cp2k_files
"""

# ============================================================
# Imports
# ============================================================

import os
from dataclasses import dataclass
from common import *


# ============================================================
# Data structures
# ============================================================
# 

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

# ============================================================
# Generic utilities
# ============================================================

def walk_dirs( root, dir_work, *, follow_symlinks=False, exclude_dirs=None, **work_kwargs ):
    """
    Generic recursive directory walker.
    Calls dir_work(path, **work_kwargs) in every directory.
        dir_clean
        dir_touch
        etc
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

# ============================================================
# Policy helpers (file selection rules)
# ============================================================

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

POLICIES = {
    'pbs': select_pbs_files,
    'slurm': select_slurm_files,
    'vasp': select_vasp_files,
    'amp': select_amp_files,
    # 'lmp': select_lmp_files,
}

def policies(path, config):
    files = []
    for work in config.works:
        if work not in POLICIES:
            raise ValueError(
                f"{work} policy not defined in {__file__}"
            )
        files.extend(POLICIES[work](path, config))
    return files

# ============================================================
# Selection & execution helpers
# ============================================================

def select_files(path, config):
    if config.selection == 'p':
        return get_files_prefix(config.slist, path)
    elif config.selection == 's':
        return get_files_suffix(config.slist, path)
    elif config.selection == 'm':
        return get_files_match(config.slist, path)
    elif config.selection == 'w':
        return policies(path, config)
    else:
        raise ValueError(f"Unknown selection mode: {config.selection}")

def apply_excludes(files, excl_files):
    return [f for f in files if f not in excl_files]


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

# ============================================================
# Directory jobs (public API)
# ============================================================

def dir_clean(path, config):
    print(f"[CLEAN] {path}")

    files = select_files(path, config)
    if config.excl_fnames:
        files = apply_excludes(files, config.excl_fnames)

    cmds = build_commands(files, config)
    execute_commands(cmds, config)

    return 0

def dir_touch(path):
    print(f"[TOUCH] {path}")

    for name in os.listdir(path):
        fpath = os.path.join(path, name)
        if os.path.isfile(fpath) and not os.path.islink(fpath):
            os.utime(fpath, None)


