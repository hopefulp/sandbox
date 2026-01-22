# ============================================================
# Policy helpers (file selection rules)
# ============================================================

from common import *

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
