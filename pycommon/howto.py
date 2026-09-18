#!/home/joonho/anaconda3/bin/python

import argparse
import re
import importlib
import sys
from pathlib import Path


def job_names_from_source(mod_name):
    """Return top-level MyClass names without importing the info module.

    Keeping help generation import-free means ``showall.py -h`` still works
    when an optional module used by one of the documentation files is absent.
    """
    source = Path(__file__).with_name(f"{mod_name}.py").read_text()
    pattern = re.compile(r"^(\w+)\s*=\s*MyClass\(", re.MULTILINE)
    return list(dict.fromkeys(pattern.findall(source)))


def subkey_names(obj):
    return [
        name for name, value in vars(obj).items()
        if not name.startswith("_") and name != "name"
        and (isinstance(value, str) or hasattr(value, "__dict__"))
    ]

def poscar_dirname(poscar):
    """Return the POSCAR filename without its POSCAR/CONTCAR prefix."""
    filename = Path(poscar).name
    return re.sub(r'^(?:POSCAR|CONTCAR)\.?', '', filename)


def format_text(text, poscar=None):
    if isinstance(text, str) and poscar:
        return text.format(POSCAR=poscar, DIRNAME=poscar_dirname(poscar))
    return text


def first_doc_line(text, poscar=None):
    if not isinstance(text, str):
        return ""
    text = format_text(text, poscar)
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return ""


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

        keys = subkey_names(obj)
        if keys:
            print("Available -k values: " + " ".join(keys))
            print()

        for name, value in vars(obj).items():

            if name.startswith("_"):
                continue
            if name == "name":
                continue

            # Case 1: direct documentation string
            if isinstance(value, str):
                print(f"{job_att}.{name}: {first_doc_line(value, poscar=poscar)}")

            # Case 2: nested namespace (e.g. vasp.make)
            elif hasattr(value, "__dict__"):

                for subname, subvalue in vars(value).items():

                    if subname.startswith("_"):
                        continue
                    if subname == "name":
                        continue

                    if isinstance(subvalue, str):
                        print(f"{job_att}.{name}.{subname}: {first_doc_line(subvalue, poscar=poscar)}")

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
        print(format_text(value, poscar))
        return 0

    # If namespace → print nested docs
    for name, subvalue in vars(value).items():

        if name.startswith("_"):
            continue

        if isinstance(subvalue, str):
            print(f"\n{name}")
            print(format_text(subvalue, poscar))

    return 0



def main():
    # ``-s`` selects system documentation. Detect it before parsing so that
    # argparse help can show the matching top-level ``-j`` values.
    switch_value = None
    for index, arg in enumerate(sys.argv[1:]):
        if arg in ('-s', '--switch'):
            next_index = index + 2
            if next_index < len(sys.argv) and not sys.argv[next_index].startswith('-'):
                switch_value = sys.argv[next_index]
            else:
                switch_value = True
            break
    system_mode = switch_value is True
    help_module = 'comment_sys' if system_mode else 'comment_subj'
    job_names = job_names_from_source(help_module)

    parser = argparse.ArgumentParser(
        description="shows dictionary for all: work, system, package",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Available -j values:\n  " + " ".join(job_names),
    )
    #parser.add_argument('-m', '--mod', default='sys', choices=['sys', 'sub'], help='which branch: system|subject')
    parser.add_argument(
        '-s',
        '--switch',
        nargs='?',
        const=True,
        default=False,
        help='choose module comment_sys; optionally pass POSCAR after -s',
    )
    parser.add_argument('-j', '--job', help='select one attribute from the list below')
    parser.add_argument('-p', '--poscar', help='args such as POSCAR name')
    parser.add_argument('-k', '--subkey', help='select one key for subkeys')
    parser.add_argument('-u', '--usage', action='store_true', help='print first keys')
    args = parser.parse_args()
    

    regex = re.compile('__')    # only detect it starts with lower case
    if isinstance(args.switch, str) and not args.poscar:
        args.poscar = args.switch

    #if args.switch==False and args.mod == 'sys':
    if args.switch is False or isinstance(args.switch, str):
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
        my_module.print_obj(job=args.job)
        if mod_name == 'comment_subj':
            print(f"\t    -s for other attributes in module 'comment_sys.py' ")
    else:
        jobs(my_module, args.job, args.subkey, args.poscar)

if __name__ == "__main__":
    main()
