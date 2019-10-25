import sys
import random
import getopt
import logging

import tqdm
import bgf
import nutils as nu

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(name)s: %(message)s", datefmt="%y-%m-%d %H:%M")
logger = logging.getLogger('add_functional.py')

help = """
Remarks
=======

1) residue name(rName) and id(rNo) is important since this determines the functional group orientation.
    - bottom layer should have: resname GRA and resid 1
    - top layer should have: resname GRA and resid 2

2) BGF file for a functional group is required.
    - alignment: the molecule should be aligned along z-axis.
    - anchoring atom: the anchoring atom should be marked as X on its chain.
    - partial charge: all partial charges should be assigned.
                      countercharge of the functional group should be assigned to the anchoring atom.

This script is only tested for attaching hydroxyl groups on the graphene surface.
Use this script at your own risk.
"""

usage = """Usage: %s -b filename -f filename -r fraction [-o filename]

Add functional groups on the graphene wall.

    -b filename     Graphitic structure in BGF format.
    -f filename     Functional group file in BGF format. Read the following remarks.
    -r fraction     Carbon-to-functional group ratio. (0 < r <= 1)
    -o filename     Filename to save in BGF format.

Please report any bugs to in.kim@kaist.ac.kr
""" % sys.argv[0]


def attach_functional(bgf_file, func_file, out_file="", set_ratio=0.125, silent=False):
    # init
    mybgf = bgf.BgfFile(bgf_file)
    sel_gra = "'C' in atom.ffType and 'GRA' in atom.rName and atom.rNo == 1"  # TODO
    sel_grb = "'C' in atom.ffType and 'GRA' in atom.rName and atom.rNo == 2"  # TODO

    # choose C site
    gra = [atom for atom in mybgf.a if eval(sel_gra)]  # list of bottom C atoms
    grb = [atom for atom in mybgf.a if eval(sel_grb)]  # list of top C atoms
    if not gra or not grb:
        nu.die("No suitable C atoms in the BGF file. Check residue names.")
    candidates = gra + grb
    target_func_n = int(set_ratio * len(candidates))
    logger.info('The script will try to attach around %s functional groups. Calculating..' % target_func_n)

    pool = []
    ratio = 0.0
    while ratio < set_ratio:
        is_good = True
        pick = random.choice(candidates)

        # position check
        if "G" in pick.chain:
            continue

        # neighbors check
        neigh_atoms = [mybgf.getAtom(i) for i in pick.CONECT]
        for i in neigh_atoms:
            if "G" in i.chain:
                is_good = False
                break

        # entrance check
        if pick.y < 28.0 or pick.y > 84.0:
            continue

        if is_good:
            pick.chain = "G"  # position mark
            pool.append(pick)
            ratio = float(len(pool)) / float(len(candidates))
            sys.stdout.write("\rGenerating Positions: %s/%s" % (len(pool), target_func_n))
            sys.stdout.flush()

    print("")
    logger.info("%s positions are successfully chosen." % len(pool))

    # attach functional groups to the site
    for atom in tqdm.tqdm(pool, desc="Attaching", ncols=80):
        # load a group
        func = bgf.BgfFile(func_file)

        # align
        if atom in grb:
            for i in func.a:
                i.x = -i.x
                i.y = -i.y
                i.z = -i.z
                i.rNo = 2

        # translate
        for i in func.a:
            i.x += atom.x
            i.y += atom.y
            i.z += atom.z

        # merge
        mybgf = mybgf.merge(func, False)

        # delete anchoring atom
        head = [a for a in mybgf.a if a.chain == 'X']
        atom.charge = head[0].charge
        head_next = mybgf.getAtom(head[0].CONECT[0])
        mybgf.delAtom(mybgf.a2i[head[0].aNo])

        # connect
        mybgf.connect(mybgf.a2i[atom.aNo], mybgf.a2i[head_next.aNo])
        mybgf.renumber()

    # save
    if out_file:
        fname = out_file
    else:
        fname = bgf_file.split('.bgf')[0] + ".mod.bgf"

    logger.info("Saving file to %s" % fname)
    mybgf.saveBGF(fname)

    logger.info("Done.")


if __name__ == '__main__':
    bgf_file = ""; func_file = ""; out_file = ""; ratio = 0.0

    print("Requested options: %s" % sys.argv)
    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    options, args = getopt.getopt(sys.argv[1:], 'hb:f:r:o:', ['help', 'bgf=', 'func=', 'ratio=', 'out='])
    for option, value in options:
        if option in ('-h', '--help'):
            print(usage)
            print(help)
            sys.exit(0)
        elif option in ('-b', '--bgf'):
            bgf_file = str(value)
        elif option in ('-f', '--func'):
            func_file = str(value)
        elif option in ('-r', '--ratio'):
            ratio = float(value)
        elif option in ('-o', '--out'):
            out_file = str(value)
        else:
            nu.die("Option %s is not recognizable." % option)

    attach_functional(bgf_file, func_file, out_file=out_file, set_ratio=ratio, silent=False)
