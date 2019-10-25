"""A class for monomer
20170213: class is modified for python 3
"""
import os
import random
import math
import subprocess

import networkx as nx
import tqdm
import bgf
import bgftools
import nutils as nu
import LAMMPS_trj2bgf


class BuildError(Exception):
    def __str__(self):
        return "Cannot build a monomer."


class RandomPolymer(nx.Graph):
    def __init__(self, max_branch=2):
        super(RandomPolymer, self).__init__()
        self.monomers = {}
        self.num_monomer = 0
        self.max_num_branch = max_branch  # number of allowed branchs per a node
        self.target_mw = 0  # target molecular weight
        self.num_target_node = 0  # number of nodes required to build a random polymer
        self.Terminal = ""
        self.Linear = ""
        self.Dendron = ""
        self.bgfmodel = bgf.BgfFile()
        self.ff = ""

    def check_branching_condition(self, monomer_id, branch_type):
        # number of branches check
        num_current_branch = len(self.neighbors(monomer_id))
        if num_current_branch > self.max_num_branch:
            nu.warn("Cannot add a monomer: branch already full.")
            return False

        # type check
        for i in self[monomer_id]:
            if monomer_id < i and self[monomer_id][i]['branch'] == branch_type:
                nu.warn("Cannot add a monomer: branch already exists.")
                return False

        return True

    def add_monomer(self, monomer_id, branch_type):
        """
        branch_type:  1: left,  2: right
        """
        if not self.check_branching_condition(monomer_id, branch_type):
            nu.warn("Failed to add monomer.")
            return False

        self.num_monomer += 1
        new_monomer_id = self.num_monomer

        self.add_node(self.num_monomer)  # add a new monomer
        self.add_edge(monomer_id, new_monomer_id, branch=branch_type)

        return new_monomer_id

    def add_focal_point(self):
        """
        focal point id == 0
        """
        self.add_node(0)

        return 0

    def del_monomer(self, monomer_id):
        self.remove_node(monomer_id)

    def update_branch(self):
        # REMARK: [i, j] -- monomer i, branch j
        # possible branching points over monomers
        possible_branch = [[i, j] for i in self.nodes() for j in range(1, self.max_num_branch + 1)]

        # existing branching points
        occupied_branch = []
        for i in self.nodes():
            for j in self.neighbors(i):
                if i < j:
                    occupied_branch.append([i, self[i][j]['branch']])

        # substract and return
        return [i for i in possible_branch if i not in occupied_branch]

    def build_random(self):
        if self.num_target_node == 0:
            nu.warn("Number of target monomers not specified.")
            raise BuildError

        # while the polymer reach to the target nodes
        while self.num_monomer < self.num_target_node:

            # update_branch
            avail_branch = self.update_branch()
            if len(avail_branch) == 0:
                nu.warn("No more available branch sites.")
                break

            # pick one site
            pick = random.choice(avail_branch)

            # attach
            self.add_monomer(pick[0], pick[1])  # pick = [monomer_id, monomer_type]

            # return: do nothing

    def calculate_distance(self):
        """
        Calculates distance from a focal point to every node
        Return a dictionary
        """
        distance = {}
        for i in self.nodes():
            distance[i] = nx.shortest_path_length(self, source=0, target=i)

        return distance

    def calculate_wiener_index(self):
        """
        Calculates Wiener index of the graph
        """
        d_dist = nx.all_pairs_shortest_path_length(self)

        index = 0.0
        for key in d_dist:
            for key2 in d_dist[key]:
                index += d_dist[key][key2]
        index /= 2.0

        return index

    def set_nodes_degree_of_branching(self):
        """
        Set each node its branch attribute (T/L/D)
        This is used to color the nodes
        """

        for i in self.nodes():
            n_branch = 0
            for j in self[i]:
                if i < j:
                    n_branch += 1

            if n_branch == 0:
                self[i]['n_branch'] = 0
            elif n_branch == 1:
                self[i]['n_branch'] = 1
            elif n_branch == 2:
                self[i]['n_branch'] = 2

    def remove_nodes_degree_of_branching(self):
        """
        After drawing networkx nodes, remove n_branch to function nx.draw() function.
        (An error occurs if other keys are defined in nodes while drawing)
        """
        for i in self.nodes():
            if 'n_branch' in self[i]:
                del self[i]['n_branch']

    def calculate_db(self):
        """
        Calculates degree of branching (DB) = 2D / (T+L+D)
        """
        t = l = d = 0.0

        for i in self.nodes():
            n_branch = 0
            for j in self[i]:
                if i < j:
                    n_branch += 1

            if n_branch == 0:
                t += 1
            elif n_branch == 1:
                l += 1
            elif n_branch == 2:
                d += 1

        return float((2 * d) / (t + l + d))

    def compile(self, filename, fastmode=True):
        """
        generate random hyperbranched polymer structure according to the graph.
        returns a BgfFile object.
        """
        # REMARK: make sure self.Terminal/Linear/Dendron is properly assigned.
        if self.Terminal == "" or self.Linear == "" or self.Dendron == "":
            nu.die("BGF file for monomers are not properly set.")
            return False

        # os.command environments
        curr_dir = os.path.abspath('.')
        temp_dir = curr_dir + "/scratch/"
        if not os.path.isdir(temp_dir):
            os.makedirs(temp_dir)
        os.chdir(curr_dir)

        for i in tqdm.tqdm(self.nodes(), ncols=120, desc="compile"):
            # introduce a monomer
            # count a number of branches connected to a node
            n_branch = 0
            parent_monomer = 0
            for j in self[i]:
                if i < j:
                    n_branch += 1
                elif i > j:
                    parent_monomer = j
            if n_branch == 0:
                new_monomer = bgf.BgfFile(self.Terminal)
            elif n_branch == 1:
                new_monomer = bgf.BgfFile(self.Linear)
            elif n_branch == 2:
                new_monomer = bgf.BgfFile(self.Dendron)

            for atom in new_monomer.a:
                atom.rNo = i  # set the residue number to the node id
                # at this moment, need to move CM of the new monomer too

            # connect the monomer to the body
            if i == 0:
                self.bgfmodel = self.bgfmodel.merge(new_monomer)
                self.bgfmodel.renumber()
            else:
                # head: in monomer (C) // tail: in the body (N)
                head_atom_ano = 0
                tail_atom_ano = 0
                for atom in self.bgfmodel.a:
                    if ("T" in atom.chain or "B" in atom.chain) and atom.rNo == parent_monomer:
                        tail_atom_ano = atom.aNo
                for atom in new_monomer.a:
                    if "H" in atom.chain and atom.rNo == i:
                        head_atom_ano = atom.aNo

                head_atom = new_monomer.getAtom(head_atom_ano)
                tail_atom = self.bgfmodel.getAtom(tail_atom_ano)

                # choose a random H atom to detach from the tail (N)
                bonding_candidate_body = []
                for ano in tail_atom.CONECT:
                    atom = self.bgfmodel.getAtom(ano)
                    if "H" in atom.ffType and "H" in atom.aName:
                        bonding_candidate_body.append(ano)
                bonding_candidate_body = random.choice(bonding_candidate_body)
                bonding_candidate_body_atom = self.bgfmodel.getAtom(bonding_candidate_body)

                # translate
                x = bonding_candidate_body_atom.x
                y = bonding_candidate_body_atom.y
                z = bonding_candidate_body_atom.z
                for atom in new_monomer.a:
                    atom.x += x
                    atom.y += y
                    atom.z += z

                # merge
                self.bgfmodel = self.bgfmodel.merge(new_monomer)
                self.bgfmodel.renumber()

                # choose a random H atom to detach from the head (C)
                bonding_candidate_monomer = []
                for ano in head_atom.CONECT:
                    atom = self.bgfmodel.getAtom(ano)
                    if "H" in atom.ffType and "H" in atom.aName:  # any H atoms are exposed to detachment
                        bonding_candidate_monomer.append(ano)
                bonding_candidate_monomer = random.choice(bonding_candidate_monomer)
                bonding_candidate_monomer_atom = self.bgfmodel.getAtom(bonding_candidate_monomer)

                # connect
                tail_atom.connect(head_atom)
                head_atom.connect(tail_atom)
                tail_atom.charge += bonding_candidate_body_atom.charge
                head_atom.charge += bonding_candidate_monomer_atom.charge
                delatoms = [self.bgfmodel.a2i[bonding_candidate_body], self.bgfmodel.a2i[bonding_candidate_monomer]]
                self.bgfmodel.delAtoms(delatoms)
                self.bgfmodel.renumber()

                # translate
                _ = bgftools.getCom(self.bgfmodel, self.ff)
                for atom in self.bgfmodel.a:
                    atom.x -= _[0]
                    atom.y -= _[1]
                    atom.z -= _[2]

                # save
                os.chdir(temp_dir)
                temp_suffix = "_polymer_" + str(i)
                temp_file = temp_dir + temp_suffix + ".bgf"  # ex) _polymer_1.bgf
                self.bgfmodel.CRYSTX = [50.0, 50.0, 50.0, 90.0, 90.0, 90.0]
                self.bgfmodel.PERIOD = "111"
                self.bgfmodel.AXES = "ZYX"
                self.bgfmodel.SGNAME = "P 1                  1    1"
                self.bgfmodel.CELLS = [-1, 1, -1, 1, -1, 1]
                self.bgfmodel.saveBGF(temp_file)

                # minimization is not performed for every attachment process
                log_i = math.ceil(math.log(i))
                if not log_i:
                    continue
                if i % log_i != 0:
                    continue

                # LAMMPS environments
                n_atoms = len(self.bgfmodel.a)
                if n_atoms < 1000:
                    mpi_command = "mpirun -n 1 "
                elif n_atoms < 2000:
                    mpi_command = "mpirun -n 2 "
                else:
                    mpi_command = "mpirun -n 4 "

                hostname = os.environ['HOSTNAME']
                lammps_command = ''
                tpascal_script_path = ''
                if "kdft" in hostname or "psi" in hostname or "rho" in hostname:
                    lammps_command = os.environ['EXEC'] + " -screen none "
                    tpascal_script_path = "/home/tpascal/scripts/"
                elif "out" in hostname:
                    lammps_command = "lammps -screen none "
                    tpascal_script_path = "/home/noische/tod_scripts/"
                elif "in" in hostname:
                    lammps_command = "lammps -screen none "
                    tpascal_script_path = "/Users/noische/codes/tod_scripts/"

                # Minimization on LAMMPS
                with open(os.devnull, 'wb') as devnull:
                    create_lmp_input = tpascal_script_path + "createLammpsInput.pl " + " -b " + temp_file + " -f " \
                                        + self.ff + " -s " + temp_suffix + " -o 'no shake' -t min "
                    subprocess.call(create_lmp_input, stdout=devnull, shell=True)
                    in_file = "in." + temp_suffix
                    data_file = "data." + temp_suffix

                    # LAMMPS input patch
                    subprocess.call("sed -i 's/dielectric      1/dielectric      72/' " + in_file, stdout=devnull, shell=True)
                    subprocess.call("sed -i 's/kspace_style    pppm 0.0001/kspace_style    none/' " + in_file, stdout=devnull, shell=True)
                    # subprocess.call("sed -i 's/boundary        p p p/boundary        s s s/' " + in_file, stdout=devnull, shell=True)
                    subprocess.call(
                        "sed -i 's/lj\/charmm\/coul\/long\/opt 7.5 8.50000/lj\/cut\/coul\/debye 0.142 10/' " + in_file, stdout=devnull, shell=True)

                    # LAMMPS data patch
                    subprocess.call("sed -i 's/0.000000  50.000000/-50.000000  50.000000/' " + data_file, stdout=devnull, shell=True)
                    subprocess.call("sed -i 's/0 # X/0 0 # X/' " + data_file, stdout=devnull, shell=True)
                    subprocess.call("sed -i 's/Impropers//' " + data_file, stdout=devnull, shell=True)

                    run_lmp = mpi_command + lammps_command + " -in in." + temp_suffix + " -log " + temp_suffix + ".log " + "-screen none"
                    subprocess.call(run_lmp, stdout=devnull, shell=True)

                    # update coordinates
                    trj_file = temp_dir + temp_suffix + ".min.lammpstrj"
                    LAMMPS_trj2bgf.getLAMMPSTrajectory(temp_file, trj_file, temp_file, -1, False, silent=True)
                    self.bgfmodel = bgf.BgfFile(temp_file)
                # end of minimization loop

        # check
        charge = 0.0
        for atom in self.bgfmodel.a:
            charge += atom.charge
        if abs(charge) >= 0.00001:
            nu.warn("Charge is not neutral: " + str(charge))
            self.bgfmodel.REMARK.append("Partial charge not neutral: " + "{0:8.5f}".format(charge))

        # compile success!
        os.chdir(curr_dir)
        self.bgfmodel.REMARK.append(str(filename) + " WI " + str(self.calculate_wiener_index()))
        self.bgfmodel.REMARK.append(str(filename) + " DB " + str(self.calculate_db()))
        self.bgfmodel.saveBGF(filename)
