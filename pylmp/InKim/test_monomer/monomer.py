import sys
import os
import time
import random
import itertools
import random

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import bgf
import bgftools
import nutils as nu
import LAMMPS_trj2bgf

class Monomer:
    def __init__(self, id):
        self.id = int(id)
        self.link = {};


    def __str__(self):
        return str(self.id) + " connected to: " + str([x.id for x in self.link])


    def add_link(self, target_monomer):
        self.link[target_monomer] = 1


class RandomPolymer(nx.Graph):
    def __init__(self, max_branch=2):
        super(RandomPolymer, self).__init__()
        self.monomers = {}
        self.num_monomer = 0
        self.max_num_branch = max_branch        # number of allowed branchs per a node
        self.target_mw = 0          # target molecular weight 
        self.num_target_node = 0    # number of nodes required to build a random polymer
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
        #links = self[monomer_id]
        #if links.has_key(branch_type):
        #    nu.warn("Cannot add a monomer: branch already occupied.")
        #    return False
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

        self.add_node(self.num_monomer)     # add a new monomer
        self.add_edge(monomer_id, new_monomer_id, branch = branch_type)

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

            ## update_branch
            avail_branch = self.update_branch()
            #print(avail_branch)
            if len(avail_branch) == 0:
                nu.warn("No more available branch sites.")
                break

            ## pick one site
            random.shuffle(avail_branch)
            pick = random.choice(avail_branch)
            #print(pick)
            
            ## attach
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
        index = 0.0;
        for key in d_dist.iterkeys():
            for key2 in d_dist[key].iterkeys():
                index += d_dist[key][key2]
        index /= 2.0

        return index


    def compile(self, filename="compile.bgf"):
        """
        generate random hyperbranched polymer structure according to the graph.
        returns a BgfFile object.
        """
        # LAMMPS preparation
        lammps_command = os.environ["EXEC"]
        curr_dir = os.path.abspath(".")
        temp_dir = curr_dir + ("/scratch/")
        if not os.path.isdir(temp_dir):
            os.makedirs(temp_dir)
        os.chdir(temp_dir)

        # REMARK: make sure self.Terminal/Linear/Dendron is properly assigned.
        if self.Terminal == "" or self.Linear == "" or self.Dendron == "":
            nu.die("BGF file for monomers are not properly set.")
            return False

        # for i in nodes:
        for i in self.nodes():
            ## introduce a monomer
            ### count a number of branches connected to a node
            n_branch = 0;
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
                atom.rNo = i    # set the residue number to the node id

            ## connect the monomer to the body
            if i == 0:
                self.bgfmodel = self.bgfmodel.merge(new_monomer)
                self.bgfmodel.renumber()
                continue    ### for headnode, no connection is required.
            else:
                ### head: in monomer (C) // tail: in the body (N)
                for atom in self.bgfmodel.a:
                    if ("B" in atom.chain or "T" in atom.chain) and atom.rNo == parent_monomer:
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
                (x, y, z) = (bonding_candidate_body_atom.x, bonding_candidate_body_atom.y, bonding_candidate_body_atom.z)
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
                    if "H" in atom.ffType and "H" in atom.aName:    # any H atoms are exposed to detachment
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
                temp_suffix = "_polymer_" + str(i)
                temp_file = temp_suffix + ".bgf"   # ex) _polymer_1.bgf
                self.bgfmodel.CRYSTX = [50.0, 50.0, 50.0, 90.0, 90.0, 90.0]
                self.bgfmodel.PERIOD = "111"
                self.bgfmodel.AXES = "ZYX"
                self.bgfmodel.SGNAME = "P 1                  1    1"
                self.bgfmodel.CELLS = [-1, 1, -1, 1, -1, 1]
                self.bgfmodel.saveBGF(temp_file)

                # Minimization on LAMMPS
                createLammpsInput = "~tpascal/scripts/createLammpsInput.pl" + " -b " + temp_file + " -f " + self.ff + " -s " + temp_suffix + " -o 'no shake' -t min " + " > /dev/null"
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
                #print("Running " + runLammps)
                os.system(runLammps)
                t2 = time.time()
                #print("Elapsed time for minimization: " + str(t2 - t1) + " sec")

                # update coordinates
                trj_file = temp_suffix + ".min.lammpstrj"
                LAMMPS_trj2bgf.getLAMMPSTrajectory(temp_file, trj_file, temp_file, -1, False, True)
                self.bgfmodel = bgf.BgfFile(temp_file)

        # compile success!
        self.bgfmodel.saveBGF(curr_dir + '/' + filename)
        # check?





