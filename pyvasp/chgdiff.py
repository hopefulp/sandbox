#!/home/joonho/anaconda3/bin/python

import argparse
import os
import shutil
import subprocess
import sys


def find_script(script):
    found = shutil.which(script)
    if found:
        return found

    base = os.path.dirname(os.path.abspath(__file__))
    for relpath in ('vtstscripts', 'vasp/vtstscripts'):
        candidate = os.path.join(base, relpath, script)
        if os.path.isfile(candidate):
            return candidate

    print(f"ERROR: cannot find {script} in PATH or pyvasp/vtstscripts", file=sys.stderr)
    sys.exit(1)


def chgcar_path(dirname, chgcar='CHGCAR'):
    path = os.path.join(dirname, chgcar)
    if not os.path.isfile(path):
        print(f"ERROR: cannot find {path}", file=sys.stderr)
        sys.exit(1)
    return path


def copy_chgcar(src, dst):
    shutil.copyfile(src, dst)
    print(f"{src} -> {dst}")


def run_stdout(cmd, outfile, workdir=None):
    print(' '.join(cmd) + f" > {outfile}")
    with open(outfile, 'w') as fout:
        subprocess.run(cmd, stdout=fout, check=True, cwd=workdir)


def make_cdd(ab_dir, a_dir, b_dir, workdir='CDD', chgcar='CHGCAR'):
    chgsum = find_script('chgsum.pl')
    chgdiff = find_script('chgdiff.pl')

    os.makedirs(workdir, exist_ok=True)

    chgcar_ab = os.path.join(workdir, 'CHGCAR_AB')
    chgcar_a = os.path.join(workdir, 'CHGCAR_A')
    chgcar_b = os.path.join(workdir, 'CHGCAR_B')
    chgcar_sum = os.path.join(workdir, 'CHGCAR_sum')
    chgcar_diff = os.path.join(workdir, 'CHGCAR_diff')

    copy_chgcar(chgcar_path(ab_dir, chgcar), chgcar_ab)
    copy_chgcar(chgcar_path(a_dir, chgcar), chgcar_a)
    copy_chgcar(chgcar_path(b_dir, chgcar), chgcar_b)

    run_stdout([chgsum, 'CHGCAR_A', 'CHGCAR_B'], chgcar_sum, workdir=workdir)
    run_stdout([chgdiff, 'CHGCAR_AB', 'CHGCAR_sum'], chgcar_diff, workdir=workdir)

    print(f"CDD file: {chgcar_diff}")
    return chgcar_diff


def main():
    parser = argparse.ArgumentParser(
        description='Make charge density difference: CHGCAR_AB - CHGCAR_A - CHGCAR_B, by reading three directories'
    )
    parser.add_argument('ab_dir', help='directory for combined A+B system')
    parser.add_argument('a_dir', help='directory for isolated A system')
    parser.add_argument('b_dir', help='directory for isolated B system')
    parser.add_argument('-w', '--workdir', default='CDD', help='directory for intermediate/output CHGCAR files')
    parser.add_argument('-c', '--chgcar', default='CHGCAR', help='CHGCAR filename in each input directory')
    args = parser.parse_args()

    try:
        make_cdd(args.ab_dir, args.a_dir, args.b_dir, args.workdir, args.chgcar)
    except subprocess.CalledProcessError as err:
        print(f"ERROR: command failed with exit code {err.returncode}: {' '.join(err.cmd)}", file=sys.stderr)
        sys.exit(err.returncode)


if __name__ == '__main__':
    main()
