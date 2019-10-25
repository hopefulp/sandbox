#!/home/noische/python

import sys
import os
import math
import random
import getopt
import copy
import logging
import random
import time
import shutil

import bgf
import bgftools
import numpy as np
import nutils as nu
import networkx as nx
import LAMMPS_trj2bgf


def calculate_wiener_index(body):
    """
    Calculates Wiener index of the generated polymer.
    :BgfAtom body
    :returns int(index): Wiener index 

    Method from from PEI_calcWienerIndex_parallel.py
    """

    myBGF = copy.deepcopy(body)

    # remove hydrogens 
    pei = bgftools.getBackbone(myBGF, 0)

    # convert connection into dictionary
    d_graph = bgftools.getConnectionDict(pei)

    # convert dictionary into Graph
    G = nx.Graph(d_graph)

    # calculate the all shortest length path
    d_dist = nx.all_pairs_shortest_path_length(G)

    # get Wiener Index
    index = 0;
    for key in d_dist.iterkeys():
        for key2 in d_dist[key].iterkeys():
            index += d_dist[key][key2];
    index = index / 2.0

    return int(index)


def refresh_branch(body_bgf):

    branches = []

    def countHydrogen(bgf_atom):
        n_Hatom = 0
        for ano in bgf_atom.CONECT:
            atom = body_bgf.getAtom(ano)
            if "H" in atom.ffType or "H" in atom.aName:
                n_Hatom += 1

        return n_Hatom

    for atom in body_bgf.a:
        if "B" in atom.chain:
            # remark: no aNo will be appended for primary -- no branching points    

            if "SEC" in atom.rName:
                n = countHydrogen(atom)
                # add one ano if two are hydrogens
                # else this branch is fully occupied already
                if n == 2:
                    branches.append(atom.aNo)
                    
            elif "TER" in atom.rName:
                # count hydrogen atoms
                n = countHydrogen(atom)

                # if two are hydrogens, add two anos
                # else if one is hydrogen and the other is head atom, add one ano
                # else if all is connected with head atoms, then this branch is fully occupied already.
                if n == 2:
                    branches.append(atom.aNo)
                    branches.append(atom.aNo)
                elif n == 1:
                    branches.append(atom.aNo)

    return branches


def check_monomer(monomer, ff_file):
    """
Check whether monomer has proper number of branching points

    :str monomer:
    :str ff_file:
    :bool : object
    """
    try:
        f = bgf.BgfFile(monomer)
    except IOError:
        return False

    branch = 0
    head = 0

    for atom in f.a:
        atom.rNo = 0	# this residue number will be the key to distinguish body and the new block 
        if atom.chain == "B":
            branch += 1
        if atom.chain == "H":
            head += 1
    if branch != 1:
        raise ValueError('More than two branch atoms inside ' + str(monomer))
        return False
    if head != 1:
        raise ValueError('More than two head atoms inside ' + str(monomer))
        return False

    mass = bgftools.getMass(f, ff_file)
    f.saveBGF(monomer)
    return mass


def get_head_ano_of_monomer(monomer):
    """
    """
    ano = 0;
    f = bgf.BgfFile(monomer)
    for atom in f.a:
        if atom.chain == "B":
            ano = atom.aNo

    return ano


def main(primary, secondary, tertiary, ratio, mw, suffix, ff, n_requested_gen, isForce, silent, isPreserveCharge=False):
    """
    generate polymer by slow addition method

    :str primary: primary monomer BGF location
    :str secondary: secondary monomer BGF location
    :str tertiary: tertiary monomer BGF location
    :list ratio: primary, secondary, tertiary amines ratio
    :float mw: desired molecular weight
    :str suffix: used to make a log filename and generated structure filename
    :str ff: force field file
    :n_requested_gen: requested number of polymer generation
    :bool silent:
    :bool isPreserveCharge:
    input: - primary, secondary, tertiary BGF file marked with branching point B and head H
           - the amine ratio
           - the target molecular weight

    remarks: - connect B and H
             - hydrogens attached on H or B will be automatically removed with respect to their role.
             - iteration will go on until we spend all designated monomers.

    writes an output: an n-oligomer BGF file

    """

    ### 0. initialize
    # log: for different format for levels, see this: http://stackoverflow.com/questions/28635679/python-logging-different-formatters-for-the-same-log-file
    log = logging.getLogger(__name__)
    log.propagate = 0
    logging.basicConfig(level=logging.DEBUG)

    formatter_file = logging.Formatter("[%(asctime)s] [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")
    formatter_debug = logging.Formatter("[%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")

    if not silent:
        streamhandler = logging.StreamHandler(sys.stdout)
        streamhandler.setFormatter(formatter_debug)
        streamhandler.setLevel(logging.INFO)

        log.addHandler(streamhandler)

    filehandler = logging.FileHandler(suffix + ".log")
    filehandler.setFormatter(formatter_file)
    filehandler.setLevel(logging.DEBUG)

    log.addHandler(filehandler)

    log.info("***** Welcome to Hyperbranched Polymer Generator by Slow Addition Method *****")
    log.info(str(sys.argv))


    # load primary, secondary, tertiary amines
    primary_mass = check_monomer(primary, ff)
    secondary_mass = check_monomer(secondary, ff)
    tertiary_mass = check_monomer(tertiary, ff)
    if not primary_mass or not secondary_mass or not tertiary_mass:
        log.error("Error found on reading monomers: Exiting.")
        sys.exit(0)


    # parameters for LAMMPS
    #mpi_command = os.environ["MPI"]
    #lammps_parallel_command = mpi_command + " " + lammps_command   # parallel not supported yet
    lammps_command = os.environ["EXEC"]
    curr_dir = os.path.abspath(".")
    temp_dir = curr_dir + ("/scratch/")
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
    shutil.copy(primary, temp_dir)
    shutil.copy(secondary, temp_dir)
    shutil.copy(tertiary, temp_dir)
    os.chdir(temp_dir)


    # determine numbers of required amine blocks from requested Mw
    ratio = np.array(ratio)
    log.info("Requested ratio: " + str(ratio))

    masses = np.array([primary_mass, secondary_mass, tertiary_mass])
    log.info("Calculated monomer mass (T/L/D): " + str(masses))

    blocks = np.ceil(mw * ratio / masses)
    log.info("Calculated required monomer blocks (T/L/D): " + str(blocks))

    estimated_mw = np.sum(masses*blocks)
    log.info("Estimated hyperbranched polymer molecular weight: " + str(estimated_mw))

    monomers = [primary, secondary, tertiary]
    blocks = list(blocks); blocks = [int(i) for i in blocks]
    log.debug(str(blocks))

    log.info("The script will generate " + str(n_requested_gen) + " polymers.")


    n_generation = 0;
    n_trial = 0;
    branch_pt_profile = [];
    while n_generation < n_requested_gen:

        n_trial += 1
        # make a pool
        # 0 = primary, 1 = secondary, 2 = tertiary
        pool = []
        for index, i in enumerate(blocks):
            log.debug("index, i: " + str(index) + " " + str(i))
            for j in range(i):
                pool.append(index)


        # the first block should not be primary(terminal)
        log.info("<<< Trial " + str(n_trial) + " >>>")
        log.info("Shuffling monomer blocks..")
        random.shuffle(pool)
        while pool[-1] == 0:
            log.debug("pool: " + str(pool))
            log.info("\tResuffle..")
            random.shuffle(pool)

        log.debug("pool: " + str(pool))


        # pick the initial block
        n_addition = 1
        init_block = pool.pop()	# random pick: the last item of the shuffled list will be chosen for the first pick.
        body = bgf.BgfFile(monomers[init_block])	# the first monomer
        for atom in body.a:
            atom.rNo = n_addition    # mark

        branches = refresh_branch(body)
        branch_pt_profile.append([n_addition, len(branches)])

        log.info("*** iteration " + str(n_addition) + " ***")
        log.debug("init block: " + str(monomers[init_block]))

        is_generation_success = False;
        tester = 0;
        # addition iteration
        #while tester == 0:
        while len(pool) > 0:
            # tester ######
            tester = 1

            # counter
            n_addition += 1
            log.info("*** iteration " + str(n_addition) + " ***")

            ### find(update) branching points of the body
            branches = refresh_branch(body)
            branch_pt_profile.append([n_addition, len(branches)])
            log.debug("aNo for Possible branches: " + str(branches))
            if len(branches) == 0:
                log.warn("No possible branching points. Failed to generate a polymer. Retrying..")
                os.system("rm in* dat* *trj *pbs *log")
                is_generation_success = False
                break

            ### if trials too many then quit the script.
            if not isForce and n_trial >= 3*n_requested_gen:
                log.error("The generator does not seem to generate proper structures. Quit.")
                # branch point print
                output = "*** number of branching points ***"
                for i in branch_pt_profile:
                    output += i[0] + '\t' + i[1] + '\n'

                log.debug(output)
                sys.exit(0)

            if n_trial > 1000 and n_generation == 0:
                log.error("No single proper structure though many attemption. Quit.")
                # branch point print
                output = "*** number of branching points ***"
                for i in branch_pt_profile:
                    output += str(i[0]) + '\t' + str(i[1]) + '\n'

                log.debug(output)
                sys.exit(0)



            ### randomly pick a branch point
            random.shuffle(branches)
            branch_pt_body = branches.pop()
            log.debug("Selected branch point aNo: " + str(branch_pt_body))


            ### randomly pick a block (pop)
            random.shuffle(pool)
            monomer_type = pool.pop()
            monomer = bgf.BgfFile(monomers[monomer_type])
            for atom in monomer.a:
                atom.rNo = n_addition   # mark

            log.debug("Selected block: " + str(monomer_type))	


            ### pick an atom from the branch atom of the body to make a bond
            branch_atom = body.getAtom(branch_pt_body)
            log.debug(str(branch_atom))
            bonding_candidate_body = []
            for ano in branch_atom.CONECT:
                atom = body.getAtom(ano)
                if "H_" in atom.ffType or "H" in atom.aName:
                    bonding_candidate_body.append(ano)
            #log.debug(str(bonding_candidate_body))
            bonding_candidate_body = random.choice(bonding_candidate_body)
            bonding_candidate_body_atom = body.getAtom(bonding_candidate_body)
            log.debug(str(bonding_candidate_body))
            log.debug(str(bonding_candidate_body_atom))


            ### translate a block
            (x, y, z) = (bonding_candidate_body_atom.x, bonding_candidate_body_atom.y, bonding_candidate_body_atom.z)
            log.debug("script will translate the block by -" + str((x, y, z)))
            for atom in monomer.a:
                atom.x += x
                atom.y += y
                atom.z += z


            ### merge
            body = body.merge(monomer) 
            body.renumber()
            #log.debug(str(body))


            ### pick an atom from the head atom of a monomer
            head_atom_ano = 0
            for atom in body.a:
                if "H" in atom.chain and atom.rNo == n_addition:
                    head_atom_ano = atom.aNo
            head_atom = body.getAtom(head_atom_ano)
            #log.debug(str(head_atom))
            bonding_candidate_monomer = []
            for ano in head_atom.CONECT:
                atom = body.getAtom(ano)
                if "H" in atom.ffType and "H" in atom.aName:
                    bonding_candidate_monomer.append(ano)
            #log.debug(str(bonding_candidate_monomer))
            bonding_candidate_monomer = random.choice(bonding_candidate_monomer)
            bonding_candidate_monomer_atom = body.getAtom(bonding_candidate_monomer)
            #log.debug(str(bonding_candidate_monomer))
            #log.debug(str(bonding_candidate_monomer_atom))


            ### connect head(C) to tail(N)
            # branch_atom -- bonding_candidate_body ---- bonding_candidate_monomer -- head_atom
            # connect
            #log.debug("Connecting " + str(branch_atom.aNo) + " and " + str(head_atom.aNo))
            branch_atom.connect(head_atom)
            head_atom.connect(branch_atom)
            #log.debug(str(body))


            # charge
            branch_atom.charge += bonding_candidate_body_atom.charge
            head_atom.charge += bonding_candidate_monomer_atom.charge


            # remove
            delatoms = [body.a2i[bonding_candidate_body], body.a2i[bonding_candidate_monomer]]
            body.delAtoms(delatoms)
            body.renumber()
            #log.debug(str(body))


            # translate
            _ = bgftools.getCom(body, ff)
            for atom in body.a:
                atom.x -= _[0]
                atom.y -= _[1]
                atom.z -= _[2]


            # refresh
            branches = refresh_branch(body)
            branch_pt_profile.append([n_addition, len(branches)])
            #log.debug(str(branches))


            # save
            temp_suffix = "_polymer_" + str(n_addition)
            temp_file = temp_suffix + ".bgf"   # ex) _polymer_1.bgf
            body.CRYSTX = [50.0, 50.0, 50.0, 90.0, 90.0, 90.0]
            body.PERIOD = "111"
            body.AXES = "ZYX"
            body.SGNAME = "P 1                  1    1"
            body.CELLS = [-1, 1, -1, 1, -1, 1]
            body.saveBGF(temp_file)


            # Minimization on LAMMPS
            createLammpsInput = "~tpascal/scripts/createLammpsInput.pl" + " -b " + temp_file + " -f " + ff + " -s " + temp_suffix + " -o 'no shake' -t min " + " > /dev/null"
            os.system(createLammpsInput)
            in_file = "in." + temp_suffix
            data_file = "data." + temp_suffix

            # LAMMPS input patch
            #os.system("sed -i 's/' " + in_file)
            os.system("sed -i 's/dielectric      1/dielectric      72/' " + in_file)
            os.system("sed -i 's/kspace_style    pppm 0.0001/kspace_style    none/' " + in_file)
            os.system("sed -i 's/boundary        p p p/boundary        s s s/' " + in_file)
            os.system("sed -i 's/lj\/charmm\/coul\/long\/opt 7.5 8.50000/lj\/cut\/coul\/debye 0.142 10/' " + in_file)

            # LAMMPS data patch
            os.system("sed -i 's/0.000000  50.000000/-50.000000  50.000000/' " + data_file)
            os.system("sed -i 's/0 # X/0 0 # X/' " + data_file)
            os.system("sed -i 's/Impropers//' " + data_file)

            t1 = time.time()
            runLammps = lammps_command + " -in in." + temp_suffix + " -log " + temp_suffix + ".log " + "-screen none"
            log.debug("Running " + runLammps)
            os.system(runLammps)
            t2 = time.time()
            log.debug("Elapsed time for minimization: " + str(t2 - t1) + " sec")


            # update coordinates
            trj_file = temp_suffix + ".min.lammpstrj"
            LAMMPS_trj2bgf.getLAMMPSTrajectory(temp_file, trj_file, temp_file, -1, False, True)
            body = bgf.BgfFile(temp_file)

            # end of a block addition
            is_generation_success = True


        # generation success
        if is_generation_success:
            n_generation += 1

            # amine group check
            _ = bgftools.getAmineGroupInfo(body)
            if _ == blocks:
                body.REMARK.append("Polymer generated with number of requested blocks.")
            else:
                log.warn("Requested monomer block does not match with the generated polymer.")
                break
            
            # charge check
            charge = 0.0
            for atom in body.a:
                charge += atom.charge
            if abs(charge) >= 0.00001:
                log.warn("Charge is not neutral: " + str(charge))

            # save the structure
            filename = suffix + "_" + str(n_generation) + ".bgf"
            body.DESCRP = "Polymer generated using " + os.path.basename(sys.argv[0]) + " by " + os.environ["USER"] + " on " + time.asctime(time.gmtime())
            body.REMARK.insert(0, "Dendritic(tertiary) monomer with " + os.path.abspath(tertiary))
            body.REMARK.insert(0, "Linear(secondary) monomer with " + os.path.abspath(secondary))
            body.REMARK.insert(0, "Terminal(primary) monomer with " + os.path.abspath(primary))
            body.saveBGF(filename)
            shutil.copy(filename, curr_dir)
            log.info("\tSaving the structure to the filename " + filename)

            # wiener index
            idx = calculate_wiener_index(body)
            log.info("Wiener index: " + str(idx))
            body.REMARK.append("Requested Mw: " + str(mw))
            body.REMARK.append("Actual mass:" +  str(bgftools.getMass(body, ff)))
            body.REMARK.append("Wiener index: " + str(idx))


    log.info(str(n_generation) + " structure generated.")
    log.info("Number of trials: " + str(n_trial))


    # branch point print
    output = "*** number of branching points ***"
    for i in branch_pt_profile:
        output += i[0] + '\t' + i[1] + '\n'
    output += '\n'


    log.debug(output)

    # done


if __name__ == '__main__':
    option = ""
    args = ""

    usage = """
Usage: PEI_generate.py -p primary -s secondary -t tertiary -r "ratio" -m Mw -o logfile -f ff_file -a
    -a Ignores failure and continuously try to generate the structure
    """

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(0)

    # Defaults
    silent = False
    out_file = ""
    T = ""
    L = ""
    D = ""
    r = []
    M = 0.0
    ff = ""
    n = 1
    isForce = False

    options, args = getopt.getopt(sys.argv[1:], 'hp:s:t:r:m:o:f:n:a',
                                  ['help', 'primary=', 'secondary=', 'tertiary=', 'ratio=', 'mw=', 'output=', 'ff=', 'n=', 'force='])
    for option, value in options:
        if option in ('-h', '--help'):
            print usage
            sys.exit(0)
        elif option in ('-p', '--primary'):
            T = str(value)
        elif option in ('-s', '--secondary'):
            L = str(value)
        elif option in ('-t', '--tertiary'):
            D = str(value)
        elif option in ('-r', '--ratio'):
            r = [float(i) for i in str.split(value)]
        elif option in ('-m', '--mw'):
            M = float(value)
        elif option in ('-o', '--output'):
            out_file = value
        elif option in ('-f', '--ff'):
            ff = str(value)
        elif option in ('-n', '--n'):
            n = int(value)
        elif option in ('-a', '--force'):
            isForce = True
        elif option in (''):
            print usage
            sys.exit(0)

    # check ratio
    if np.array(r).sum() == 100:
        r = [float(i/100) for i in r]
    elif np.array(r).sum() != 1:
        nu.die("Ratio sum should be 1 (fraction) or 100 (percentage)")

    # defaults
    if ff == "":
        ff = "/home/noische/ff/DREIDING2.21.ff"
    if n == 0:
        n = 1

    # main run
    main(T, L, D, r, M, out_file, ff, n, isForce, silent)
