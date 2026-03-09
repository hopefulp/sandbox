#!/home/joonho/anaconda3/bin/python

import argparse
import re
import importlib

### not used
def format_text(text, mod):
    if isinstance(text, str) and '{POSCAR}' in text:
        return text.format(POSCAR=getattr(mod, 'POSCAR', '{POSCAR}'))
    return text


def jobs(mod_comm, job_att, subkey=None, poscar=None):

    obj = getattr(mod_comm, job_att, None)
    if obj is None:
        print("No such job")
        return 0

    print(job_att)
    print('---------------------------------------------')

    # =========================
    # SUMMARY MODE
    # =========================
    if subkey is None:

        for name, value in vars(obj).items():

            if name.startswith("_"):
                continue

            # Case 1: direct documentation string
            if isinstance(value, str):
                print(f"{job_att}.{name}")

            # Case 2: nested namespace (e.g. vasp.make)
            elif hasattr(value, "__dict__"):

                for subname, subvalue in vars(value).items():

                    if subname.startswith("_"):
                        continue

                    if isinstance(subvalue, str):
                        print(f"{job_att}.{name}.{subname}")

        print("\nUse -k for detail")
        return 0

    # =========================
    # DETAIL MODE
    # =========================
    value = getattr(obj, subkey, None)
    if value is None:
        print("No such subkey")
        return 0

    print(f"{job_att}.{subkey}")
    print('---------------------------------------------')

    # If leaf string → print doc
    if isinstance(value, str):
        print(value.format(POSCAR=poscar) if poscar else value)
        return 0

    # If namespace → print nested docs
    for name, subvalue in vars(value).items():

        if name.startswith("_"):
            continue

        if isinstance(subvalue, str):
            print(f"\n{name}")
            print(subvalue.format(POSCAR=poscar) if poscar else subvalue)

    return 0



def main():
    parser = argparse.ArgumentParser(description="shows dictionary for all: work, system, package  ")
    #parser.add_argument('-m', '--mod', default='sys', choices=['sys', 'sub'], help='which branch: system|subject')
    parser.add_argument('-s', '--switch', action='store_true', help='choose module comm_sub')
    parser.add_argument('-j', '--job', help='select one attribute')
    parser.add_argument('-p', '--poscar', help='args such as POSCAR name')
    parser.add_argument('-k', '--subkey', help='select one key for subkeys')
    parser.add_argument('-u', '--usage', action='store_true', help='print first keys')
    args = parser.parse_args()
    

    regex = re.compile('__')    # only detect it starts with lower case
    #if args.switch==False and args.mod == 'sys':
    if args.switch==False:
        mod_name = 'comment_subj'
    else:
        mod_name = 'comment_sys'
    my_module = importlib.import_module(mod_name)
    my_module.POSCAR = args.poscar        # module level variable: add new attribute 'POSCAR' to mod my_module
    ### same as
    #my_module.__dict__['POSCAR'] = "anything"
    ### try to pass args 
    if not args.job or args.usage:
        print(my_module.__file__)
        my_module.print_obj( job = args.job, poscar=args.poscar )
        if mod_name == 'comment_subj':
            print(f"\t    -s for other attributes in module 'comment_sys.py' ")
    else:
        jobs(my_module, args.job, args.subkey, args.poscar)

if __name__ == "__main__":
    main()
