# ============================================================
# Data structures
# ============================================================
# 
import os
from dataclasses import dataclass
from common import *
from .policy import *

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
# Selection & execution helpers
# ============================================================

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

def dir_clean(path, config):
    print(f"[CLEAN] {path}")

    files = select_files(path, config)
    if config.excl_fnames:
        files = apply_excludes(files, config.excl_fnames)

    cmds = build_commands(files, config)
    execute_commands(cmds, config)

    return 0