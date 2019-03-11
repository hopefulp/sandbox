# -*- coding: utf-8 -*-
'''
This is a python class for the QEq method
according to Rappe et al, JPC 1991, 95, 3358-3363
Note: alist has to be an uppercase string array of width 2
Added initial support for ACKS2 in V0.5
(c) 2010 Christian Spickermann, R. Schmid, D. Pache

V 0.2 removing callback and add a hook for the QEq call in the
      overall energy routine called by all MD, optimization routines (RS)
      
v 0.3 major rewrite: instead of subclassing for the different codes (with explicit
      reference to these codes (like importing _pydlpoly within this module) we expect a certain API
      for the charge equlibration. this means upon instantiation of the QEq an instance
      of the MM code needs to be passed to the init. with this we have access to certain
      methods that need to be supplied by the MM engine.
      
      for the ease of comparison with the dl_poyl part we convert all from kJ/mol
      to kcal/mol (it would be better the other way round: use SI units in dl_poly,
          but this is much more difficult to do so i go the lazy way :-) RS
      
      These methods are:
          get_natoms       -> return the number of atoms in the system
          get_elements     -> return a list of element keys to init the params
          get_atomtypes    -> alternative to elements (in case we need different params for the same element)
          get_charges      -> get the charges as a numpy array
          set_charges      -> set the charges from a numpy array
          calc_coul_energy -> compute just the Coulombic part of the energy of the MM system
                              (no other contribs/forces etc.) and return the force on the
                              charges (potential) in a numpy array
          set_onsite_charge_energy -> a way to hand back the energy of the onsite terms in QEq to the Coul energy of the MM engine
          set_recalc       -> tell qeq that the geometry has changed and a recalcualtion is needed
      optional helpers/flags related to parallel runs:
          nodes    -> number of nodes (if 1 then not parallel)
          idnode   -> id of the node
          
      optional methods:
          get_xyz          -> get atomic coordinates (should not be necessary if the potential is computed in the MM engine)
                              but might be needed for some additional calculations
          calc_Jmat         -> get the matrix of all interactions, which was originally called Hessian. It is exactly the matrix of
                              second derivatives in the charge space but in order not confuse it with the xyz Hessian we call it
                              J-matrix (see definitions document on the meaning of J)
          
      as a convention all the further data to control QEq will be read from the control file AFTER the "finish" keyword
      which makes pydlpoly stop parsing its input. This file is always present local and is always called CONTROL.
      
V 0.4 We will have to rename this thing because it will be more then just QEq in the future. Currently just updates to match with
      IO handeled by :class:`pdlpio.pdlpio` and molecules in :class:`pdlpmol.pdlpmol`.
      Also a per molecule constraint is added
     
     
V 0.5 initial ACKS2 version added (keep the name) by a named parameter in the __init__ (R. Schmid)

V 0.6 Full ACKS2 implementation via matrix multiplication (see q_solver). DP
'''

import sys
import numpy
import string
import timer
#import inter
#import nlist

from erf import erf

import scipy.optimize as sciopt

pi4e0 = 1389.354566140940121 # returns Coulomb energy in kJ/mol if distance is in Angstrom and "chemical" charges are employed
factor = 96.48534            # from ev to kJ/mol
j2cal = 1.0/4.184            # from kJ/mol to kcal/mol
#dlp2kJ = 1.0E-2              # dlpoly energy unit is 10 J/mol
#hart2kJ = 2625.4996356       # from a.u. to kJ/mol
#ang2bohr = 1.0/0.52917721    # from Angstrom to a.u.

# do all in kcal/mol for comparison with dl_poly energies
pi4e0 *= j2cal
factor *= j2cal

numpy.set_printoptions(threshold=5000)

class qeq:

    def __init__(self, pd, input_file=None, exclude=None, acks2=False):
        """
        Initalizes qeq object
        
        :Parameters:
            
            - pd         : reference to parent pydlpoly instance
            - input_file : name of the input file, if `None` (default) a file named <runname>.qeq is expected
            - exclude    : list of molecules excluded from polarizability
            - acks2      : add s-basis ACKS2 variant (pure python)
        """
        self.pd = pd
        self.natoms = self.pd.get_natoms()
        self.tstep = self.pd.get_tstep()
        self.exclude = exclude
        # acks2
        self.acks2 = acks2
        # set up parallel info
        self.idnode = self.pd.idnode
        self.nodes  = self.pd.nodes
        if input_file:
            self.input_file = input_file
        else:
            self.input_file = self.pd.start_dir + "/" + self.pd.name + ".qeq"
        self.mark_line = 80*"&"
        self.print_header()
        self.par_with_sigma = False
        self.molconst = False
        # set up core charge usage
        self.use_core_charge=self.pd.inquire_core_charge_use()
        if self.use_core_charge:
            self.core_charges = self.pd.get_core_charges()
            self.tot_core_charge = numpy.sum(self.core_charges)
        self.parse_input()
        self.allocate_arrays()
        self.assign_atomparams()
        self.init_q()
        self.init_invH = False
        self.report_step = 20
        self.timer = timer.timer("QEq Timer")
        #
        self.last_E = 0.0
        self.E = 0.0
        self.step = 0
        # defaults
        self.qmass = 0.0006
        # set the qeq_iter flag to True in pydlpoly to compute qpot 
        self.pd.dlp_setup.qeq_iter = True
        ###Charge constraint
        self.limit_q = True
        self.Jij3    = 400.0
        self.q_thresh= 2.0
        self.E_ii3 = 0.0
        #extended lagrangian
        # the flag ext_lagr defines whether we are in the extended lagrangian mode or not
        self.ext_lagr = False
        #NEWNEW
        if self.acks2:
            acks2_params = inter.read_Xij_params("water.acks2")
            print acks2_params
            self.Xij = inter.Xij(self.pd, 8.0, 1.0, acks2_params, sys.stdout)
            self.Xij.nlist.exclude_mols("water")
            self.Xij.nlist.set_inject_excluded(True)
        return


    def pprint(self,s):
        if self.idnode ==  0 : print(s)
        return

    def print_header(self):
        self.pprint(self.mark_line+"\n")
        self.pprint("QEq Charge Equilibration initialized"+"\n")
        if self.acks2:
            self.pprint("   USING INITIAL ACKS2!")
        self.pprint(self.mark_line)
        return
        
    
    def parse_input(self):
        infile = open(self.input_file,'r')
        stop = False
        while not stop:
            line = infile.readline()
            if len(line) == 0:
                raise IOError, "QEq block in input file is missing"
            line = string.split(line)
            if len(line) > 0:
                if line[0] == "QEq": stop = True
        # reached QEq params in CONTROL
        self.pprint("parsing QEq params in the input file")
        # first set the default values
        self.Qtot = 0.0
        self.prntlvl = 1
        self.unit = 1.0
        self.conv = 0.01
        self.params = {}
        self.strategy = "lbfgs"
        self.maxiter = 5000
        self.sd_step = 0.0015
        # check if mandatory params are read in
        params_read = False
        # now parse until QEq_finish
        stop = False
        while not stop:
            line = infile.readline()
            if len(line)== 0:
                self.pprint("QEq_fnish marker is missing in file")
                stop=True
            line = string.split(line)
            if len(line) > 0:
                if line[0] == "params":
                    if line[1] == "element":
                        self.assign_by_elements = True
                    elif line[1] == "atomtypes":
                        self.assign_by_elements = False
                    else:
                        raise IOError, "Unknown parameter option in CONTROL file for QEq"
                    if len(line) == 3:
                        if line[2] == "withsigma":
                            self.par_with_sigma = True
                    if len(line) > 3:
                        # ok, we have a parfile read => read this instead
                        try:
                            parfile = open(line[3],'r')
                        except:
                            raise IOError, "Could not read the QEq parfile %s" % line[3]
                        content = parfile.readlines()
                        parfile.close()
                    else:
                        # read params directly from CONTROL
                        content = []
                        parstop = False
                        parline = infile.readline()
                        while not parstop:
                            if string.split(parline)[0] == "end":
                                parstop = True
                            else:
                                content.append(parline)
                            parline = infile.readline()
                    for line in content:
                        if self.par_with_sigma:
                            atom, En, Hard, sigma = line.split()
                            atom = string.upper(atom)
                            if len(atom)<2: atom = string.ljust(atom,2)
                            self.params[atom] = (float(En)*factor, float(Hard)*factor, float(sigma))  # all values are stored in kcal/mol
                            self.pprint("   & parameter for %4s in kcal/mol EN: %12.6f   Jii: %12.6f  sigma [A] %10.5f" % \
                                 (atom,self.params[atom][0],self.params[atom][1],self.params[atom][2]) )
                        else:
                            atom, En, Hard = line.split()
                            atom = string.upper(atom)
                            if len(atom)<2: atom = string.ljust(atom,2)
                            self.params[atom] = (float(En)*factor, float(Hard)*factor)  # all values are stored in kcal/mol
                            self.pprint("   & parameter for %4s in kcal/mol EN: %12.6f   Jii: %12.6f" % \
                                 (atom,self.params[atom][0],self.params[atom][1]) )
                    params_read = True
                elif line[0] == "Qtot":
                    self.Qtot = string.atof(line[1])
                    self.pprint("   & Qtot = %12.6f " % self.Qtot)
                elif line[0] == "molconst":
                    if len(line) > 1:
                        if line[1] == "on":
                            self.molconst = True
                            self.pprint("   & Constraining Q per molecules")
                elif line[0] == "unit":
                    self.unit = string.atof(line[1])
                    self.pprint("   & unit = %12.6f " % self.unit)
                elif line[0] == "conv":
                    self.conv = string.atof(line[1])
                    self.pprint("   & conv = %12.6f " % self.conv)
                elif line[0] == "strategy":
                    self.strategy = line[1]
                    self.pprint("   & strategy = %s " % self.strategy)
                elif line[0] == "printlevel":
                    self.prntlvl = string.atof(line[1])
                    self.pprint("   & printlevel = %12.6f " % self.prntlvl)
                elif line[0] == "maxiter":
                    self.maxiter = string.atoi(line[1])
                    self.pprint("   & maxiter = %12.6f " % self.sd_step)
                elif line[0] == "sd_step":
                    self.sd_step = string.atof(line[1])
                    self.pprint("   & steepest descent stepsize = %12.6f " % self.sd_step)
                #Extended lagrangian
                elif line[0] == "qmass":
                    self.params_extlagr = {}
                    all_qmass = []
                    qmass_stop = False
                    qmass_line = infile.readline()
                    while not qmass_stop:
                        if string.split(qmass_line)[0] == "endqmass":
                            qmass_stop = True
                        else:
                            all_qmass.append(qmass_line)
                            qmass_line = infile.readline()
                    for line in all_qmass:
                        qmass_atom, q_mass = line.split()
                        qmass_atom = string.upper(qmass_atom)
                        q_mass = float(q_mass)
                        if len(qmass_atom)<2: qmass_atom = string.ljust(qmass_atom,2)
                        self.params_extlagr[qmass_atom] = q_mass
                        self.pprint("   & qmass for %5s : %10.5f" % (qmass_atom, q_mass))
                # ACKS2 params quick and dirty ....
                elif line[0] == "acks2":
                    self.params_acks2 = {}
                    all_acks2 = []
                    acks2_stop = False
                    acks2_line = infile.readline()
                    while not acks2_stop:
                        if string.split(acks2_line)[0] == "endacks2":
                            acks2_stop = True
                        else:
                            all_acks2.append(acks2_line)
                            acks2_line = infile.readline()
                    for line in all_acks2:
                        sline = line.split()
                        atoms = map(string.strip, sline[:2])
                        atoms = map(string.upper, atoms)
                        pars  = map(float, sline[2:])
                        self.params_acks2[string.join(atoms, ":")] = pars
                        atoms.reverse()
                        self.params_acks2[string.join(atoms, ":")] = pars
                        self.pprint("   & ACSK Xij params %s - %s : %s" % (atoms[1], atoms[0], pars))
                elif line[0] == "QEq_finish":
                    stop = True
                else:
                    self.pprint ("Unknown keyword in QEq params: %s" % line[0])
        # all read in
        infile.close()
        # check mandatory params
        if not params_read:
            raise IOError, "No QEq atom parameters in your CONTROL file. Aborting!"
        # setup Qtot in case of core charge usage
        if self.use_core_charge:
            self.Qtot -= self.tot_core_charge
            self.pprint("   & USING CORE CHARGES!! total core charge is   %12.6f" % self.tot_core_charge)
            self.pprint("   & in QEq the valence charge is constrained to %12.6f" % self.Qtot)
        self.pprint("\nFINISHED READING QEQ PARAMS")
        self.pprint(self.mark_line+"\n")
        return

    def allocate_arrays(self):
        self.En   = numpy.zeros([self.natoms], dtype="float64")
        self.Jii  = numpy.zeros([self.natoms], dtype="float64")
        if self.molconst:
            self.mol_pot = numpy.zeros([self.pd.molecules.nmolecules], dtype="float64")
            self.mol_Qtot = numpy.zeros([self.pd.molecules.nmolecules], dtype="float64")
            self.mol_ones = numpy.zeros([self.natoms, self.pd.molecules.nmolecules], dtype="float64")
            self.mol_natoms = numpy.zeros([self.pd.molecules.nmolecules], dtype="int32")
            for i, m in enumerate(self.pd.molecules.mlist):
                self.mol_natoms[i] = m.natoms
                self.mol_Qtot[i] = m.Q
                for j in m.matoms:
                    self.mol_ones[j, i] = 1.0
        if self.acks2:
            self.u = numpy.zeros([self.natoms], dtype="float64")
            self.Xii  = numpy.zeros([self.natoms], dtype="float64")
        return

    def assign_atomparams(self):
        # assigns atomic parmeters form parmeter dictonary to the corresponding arrays
        if self.assign_by_elements:
            alist = self.pd.get_elements()
        else:
            alist = self.pd.get_atomtypes()
        if self.par_with_sigma:
            sigmas = numpy.zeros([self.natoms], dtype="float64")
        for i in xrange(self.natoms):
            ali = alist[i]
            if len(ali)<2: ali = string.ljust(ali,2)
            try:
                p = self.params[ali.upper()]
            except KeyError:
                self.pprint("No QEq paramter for atom %s" % ali.upper())
                raise IOError
            self.En[i]  = p[0]
            self.Jii[i] = p[1]
            if self.par_with_sigma: sigmas[i] = p[2]
        if self.par_with_sigma:
            self.pd.set_sigmas(sigmas)
        ### ACKS2 ###
        # we just test if all params are there
        if self.acks2:
            for a1 in alist:
                for a2 in alist:
                    #print string.join([a1.strip().upper(), a2.strip().upper()], ":")
                    assert string.join([a1.strip().upper(), a2.strip().upper()], ":") in self.params_acks2.keys(), "ACKS2 parameter for %s %s is missing" % (a1, a2)
        return

    def init_q(self, Zero=False):
        if Zero:
            self.q = numpy.zeros([self.natoms], dtype="float64")
            self.q[:] = self.Qtot/self.natoms
        else:
            self.q = self.pd.get_charges()
        self.set_recalc()
        return

    def set_recalc(self):
        self.recalc_qpot = True
        self.recalc_Jij  = True
        return

    def update_q(self):
        # update charges in the mm engine
        self.pd.set_charges(self.q)
        self.recalc_qpot = True
        self.recalc_Jij  = True
        return

    def set_Eii(self):
        # update the coulomb energy in the mm engine
        self.pd.set_onsite_charge_energy(self.E_ii/self.unit)
        return

    def zero_Eii(self):
        # set the one center term in the MM engine to zero in order to get the rigth E_ij
        self.pd.set_onsite_charge_energy(0.0)
        return

    def calc_energy_and_potential(self, Jij=False, fill_Jii=True, core_pot=False, constrain=False, debug=None):
        """
        Compute energies and potentials for QEq_finish

        :Parameters:

            - Jij (boolean) : set to `True` to call the calc_Jmat method of the mme
            - fill_Jij (boolean) : set to `True` to fill the diagonal elements of Jij
            - core_pot (boolean) : use core potentials
            - constrain (boolean) : constrain charge to Qtot (if "molconst on" is specified in the input a per molecule contraint is used)

        """
        #flag for the calculation of degrees of freedom for get_dof()
        if constrain:
            self.constraint = True
        # first compute the twocenter part with the MM engine
        if Jij:
            if self.recalc_Jij:
                self.timer.switch_to("CORE calc_Jmat")
                if core_pot:
                    self.E_ij, self.qpot_ij, self.Jij, self.core_qpot = self.pd.calc_Jmat(get_core_pot=True)
                    #print "DEBUG DEBUG DEBUG"
                    #print "E_ij from fortran code", self.E_ij
                    #print "Eval   ", self.pd.get_val_energy()
                    #print "Ecore  ", self.pd.get_core_energy()
                    #print "Ecore_val  ", self.pd.get_coreval_energy()
                    #print "Ecoreval from core_pot",  numpy.sum(self.core_qpot*self.q)
                    #print "Eij_val_only from qpot", numpy.sum(self.qpot_ij*self.q)
                else:
                    self.E_ij, self.qpot_ij, self.Jij = self.pd.calc_Jmat()
                self.Jij *= self.unit
                self.E_ij    *= self.unit
                self.qpot_ij *= self.unit
                self.recalc_Jij = False
                self.recalc_qpot = False
        else:
            if self.recalc_qpot:
                self.timer.switch_to("CORE calc_Jij")
                if core_pot:
                    self.E_ij, self.qpot_ij, self.core_qpot = self.pd.calc_coul_energy(get_core_pot=True)
                    #print "DEBUG DEBUG DEBUG"
                    #print "E_ij from fortran code", self.E_ij
                    #print "Eval   ", self.pd.get_val_energy()
                    #print "Ecore  ", self.pd.get_core_energy()
                    #print "Ecore_val  ", self.pd.get_coreval_energy()
                    #print "Ecoreval from core_pot",  numpy.sum(self.core_qpot*self.q)
                    #print "Eij_val_only from qpot", numpy.sum(self.qpot_ij*self.q)
                else:
                    self.E_ij_for, self.qpot_ij = self.pd.calc_coul_energy()
                    #print "DEBUG DEBUG DEBUG"
                    #print "E_ij from fortran code", self.E_ij_for
                    #print "Eij_val_only from qpot", 0.5*numpy.sum(self.qpot_ij*self.q)
                    # FIX this ... fortran gives us wrong E_ij ... on the second call. must be an initalization error.
                    #self.E_ij = 0.5*numpy.sum(self.qpot_ij*self.q)
                    self.E_ij = self.E_ij_for
                self.E_ij    *= self.unit
                self.qpot_ij *= self.unit
                self.recalc_qpot = False
        self.timer.switch_to("CORE calc qpot")
        # now add the one-center terms, note that following the original implementation the Jii is really 2*Jii
        # as we can add it directly (times the charge) to the potential
        # => this convention means we need to take it 0.5 for the energy
        qpot_ii = self.Jii*self.q
        self.qpot = self.En + qpot_ii + self.qpot_ij
        # add core_pot if core charges are used. Note: the flag core_pot triggers a recalc but if
        #     core charges are used we add it in any case .. so it needs to be generated in the first call
        if self.use_core_charge:
            self.qpot += self.core_qpot
        self.E_ii = numpy.sum(self.q*(self.En + 0.5*qpot_ii))
        if self.limit_q:
            q_large = numpy.clip(numpy.fabs(self.q)-self.q_thresh, 0.0, 1.0e6)
            q_large2 = q_large*q_large
            qpot_ii3 = numpy.copysign(3.0*self.Jij3*q_large2, self.q)
            self.qpot += qpot_ii3
            self.E_ii3 = numpy.sum(self.Jij3*q_large*q_large2)
            self.E_ii += self.E_ii3
        self.E    = self.E_ii+self.E_ij
        if debug:
            self.pprint("DEBUG This is QEq calc_energy_force! called from %s \n Eii = %12.6f Eij = %12.6f E = %12.6f" % (debug, self.E_ii, self.E_ij, self.E))
        if Jij:
            if fill_Jii:
                for i in xrange(self.natoms):
                    self.Jij[i,i] = self.Jii[i]
        # charge conservation constraint
        # NOTE: this is currently terribly non-parrallel
        # could be done in two loops first summing up and then correcting
        # with a allreduce sum inbetween.
        if constrain:
            self.timer.switch_to("CORE constrain")
            if self.molconst:
                #for m in self.pd.molecules.mlist:
                #    molpot = self.qpot[m.matoms]
                #    tot_mpot = molpot.sum()
                #    molpot -= tot_mpot/m.natoms
                #    self.qpot[m.matoms] = molpot[:]
                # NOTE this is a bit of a hack ... the fortran constrain is private to QEq, but I did not
                #      want to have another module for it, so i put it into the coulomb module of pydlpoly ... so this is not a "generic mme"
                #      interface any more
                self.pd.dlp_coul.constrain_qpot_permol(self.idnode, self.nodes, self.qpot, self.natoms,
                                                    self.mol_pot, self.pd.molecules.nmolecules, self.mol_natoms)
            else:
                # constrain for the total system
                tot_pot = numpy.sum(self.qpot)
                self.qpot -= tot_pot/self.natoms

        #NEWNEW
        if self.acks2:
            self.calc_acks2_X()
            self.calc_acks2_u()
            self.E += self.calc_acks2_energy(self.u)
            self.qpot += self.u

        if self.exclude:
            for i in self.exclude:
                self.qpot[self.pd.molecules.mlist[i].matoms]=0
        #self.pprint("###Forces from calc_energy_and_potential###")
#        self.pprint("Identical with qpot sh, only here to show initial lbfgs forces")
        #self.pprint(self.qpot_ij)
        #self.pprint(self.qpot)
#        self.print_charges()
        return self.E

    def calc_Eii(self):
        self.E_ii = numpy.sum(self.q*(self.En + 0.5*(self.Jii*self.q)))
        if self.acks2:
            self.E_ii += numpy.sum(self.u*(-self.q+0.5*(self.Xij.diag*self.u)))
#        if self.limit_q:
#            self.E_ii += self.E_ii3
        return self.E_ii


    def chk_conv(self):
        rms = numpy.sqrt(numpy.vdot(self.qpot,self.qpot)/self.natoms)
        if rms < self.conv:
            # if self.prntlvl > 0: self.pprint("QEq converged .. rms is %12.6f" % rms)
            return rms, True
        else:
            return rms, False

    def qopt_callback(self, new_q):
        """
        callback function for iterative charge optimizers
        it is called with the new charge and returns current energy and cosntrained qpot
        """
        self.q[:] = new_q[:]
        self.update_q()
        e = self.calc_energy_and_potential(constrain=True, debug=False)
        return e, self.qpot

    def analyze_J_eigenvalues(self, negative_only=True):
        evals,vect = numpy.linalg.eigh(self.Jij)
        for i,e in enumerate(evals):
            v = vect[:,i]
            if e < 0.0:
                self.pprint("%10.3f  : %s" % (e, numpy.array_str(v, precision=2, suppress_small=True)))
            else:
                if not negative_only:
                    self.pprint("%10.3f  : %s" % (e, numpy.array2string(v, precision=2, suppress_small=True)))
        return

    def q_solver(self, test_convergence=True, check_eigenvalues=True):
        if self.molconst:
            nconst = self.pd.molecules.nmolecules
        else:
            nconst = 1
        # solve equations .. in this case no intial solution is needed
        self.zero_Eii()
        self.calc_energy_and_potential(Jij=True,core_pot=self.use_core_charge)
        if check_eigenvalues: self.analyze_J_eigenvalues()
        self.H = numpy.zeros([self.natoms+nconst,self.natoms+nconst],dtype="float64")
        if self.acks2: self.H = numpy.zeros([(self.natoms+nconst)*2,(self.natoms+nconst)*2],dtype="float64")
        self.x = numpy.zeros([self.natoms+nconst],dtype="float64")
        if self.acks2: self.x = numpy.zeros([(self.natoms+nconst)*2],dtype="float64")
        if self.use_core_charge:
            self.x[0:self.natoms] = -self.En[:]-self.core_qpot[:]
        else:
            self.x[0:self.natoms] = -self.En[:]
        # set up H matrix
        # NOTE: maybe it is better to keep this matrix and just rewrite the Jij
        self.H[:self.natoms,:self.natoms] = self.Jij[:,:]
        if self.molconst:
            self.x[self.natoms:self.natoms+nconst]  = self.mol_Qtot[:]
            self.H[self.natoms:self.natoms+nconst,:self.natoms] = self.mol_ones.transpose()[:,:]
            self.H[:self.natoms,self.natoms:self.natoms+nconst] = self.mol_ones[:,:]
        else:
            self.x[self.natoms] = self.Qtot
            self.H[self.natoms,:self.natoms] = 1.0
            self.H[:self.natoms,self.natoms] = 1.0
        if self.acks2:
            self.calc_acks2_X(non_sparse=False)
            self.H[:self.natoms,self.natoms+nconst:2*self.natoms+nconst] = -numpy.eye(self.natoms)
            self.H[self.natoms+nconst:2*self.natoms+nconst,:self.natoms] = -numpy.eye(self.natoms)
            self.H[self.natoms+nconst:2*self.natoms+nconst,self.natoms+nconst:2*self.natoms+nconst] = self.X[:]
            if self.molconst:
                self.H[2*self.natoms+nconst:2*(self.natoms+nconst),self.natoms+nconst:2*self.natoms+nconst] = self.mol_ones.transpose()[:,:]
                self.H[self.natoms+nconst:2*self.natoms+nconst,2*self.natoms+nconst:2*(self.natoms+nconst)] = self.mol_ones[:,:]
            else:
                self.H[2*self.natoms+nconst,self.natoms+nconst:2*self.natoms+nconst] = 1.0
                self.H[self.natoms+nconst:2*self.natoms+nconst,2*self.natoms+nconst] = 1.0
        q = numpy.linalg.solve(self.H,self.x)
        self.q = q[0:self.natoms]
        self.lamb = q[self.natoms:self.natoms+nconst]
        if self.acks2:
            self.u = q[self.natoms+nconst:2*self.natoms+nconst]
            self.ny = q[2*self.natoms+nconst:2*(self.natoms+nconst)]
        self.update_q()
        if test_convergence:
            E = self.calc_energy_and_potential(constrain=True,core_pot=self.use_core_charge)
            self.set_Eii()
            rms, converged = self.chk_conv()
            if self.prntlvl > 0:
                if self.molconst:
                    self.pprint("  q_solver : E = %12.6f, rms = %12.6f, sum of q: %12.6f" % \
                            (E, rms, numpy.sum(self.q)))
                else:
                    self.pprint("  q_solver : E = %12.6f, rms = %12.6f, lambda = %12.6f, sum of q: %12.6f" % \
                            (E, rms, self.lamb[0], numpy.sum(self.q)))
        else:
            # compute and set final onecenter terms (if no converegnce test is made)
            self.calc_Eii()
            self.set_Eii()
            converged = True
        #if self.acks2:
            #self.calc_and_print_acks2_pot()
            #print "###### self.u#####"
            #print self.u

        #delete delete
        self.manual_output_selection()
        return converged


    def steepest_descent(self):
        self.timer.switch_to("sd_solver pre loop")
        step = 0
        self.zero_Eii()
        E = self.calc_energy_and_potential(constrain=True)
        rms, converged = self.chk_conv()
        while (not converged and (step < self.maxiter)):
            self.timer.switch_to("sd_solver propagate q")
            step += 1
            self.q -= self.sd_step*self.qpot
            self.update_q()
            self.timer.switch_to("sd_solver calc_E")
            E = self.calc_energy_and_potential(constrain=True)
            self.timer.switch_to("sd_solver check conv")
            rms, converged = self.chk_conv()
            if (self.prntlvl > 1) and (step%self.report_step == 0) :
                self.pprint("  step %d : E = %12.6f, rms = %12.6f, sum of q: %12.6f" % \
                            (step, E, rms, numpy.sum(self.q)))
                self.pprint("            E_ii: %12.6f E_ij: %12.6f  E: %12.6f" %
                    (self.E_ii, self.E_ij, self.E))
        self.timer.switch_to("sd_solver post loop")
        if self.prntlvl > 0:
            self.pprint("  sd converged after %d steps : E = %12.6f, rms = %12.6f, sum of q: %12.6f" % \
                            (step, E, rms, numpy.sum(self.q)))
        self.set_Eii()
        return converged

    def lbfgs_scipy(self):
        self.timer.switch_to("lbfgs solver")
        self.zero_Eii()
        q_init  = self.q.copy()
        q_final, E, optstat = sciopt.fmin_l_bfgs_b(self.qopt_callback, q_init, pgtol=self.conv, disp=0, factr=1e5, iprint=0)
        if optstat["warnflag"] > 2:
            self.pprint("SOEMTHIGN WENT WRONG!!")
            self.pprint(optstat["task"])
        self.timer.switch_to("lbfgs check_conv")
        # not sure if this is necessary ... probaly the values are already pushed to the MME
        self.q[:] = q_final
        self.update_q()
        E = self.calc_energy_and_potential(constrain=True, debug=None)
        rms, converged = self.chk_conv()
        #qpot_final = optstat["grad"]
        #rms = numpy.sqrt(numpy.vdot(qpot_final,qpot_final)/self.natoms)
        if self.prntlvl > 0:
            self.pprint("  lbfgs converged after %d energy calcs: E = %12.6f, rms = %12.6f, sum of q: %12.6f" % \
                            (optstat["funcalls"], E, rms, numpy.sum(self.q)))
            if self.E_ii3 > 1.0e-6:
                self.pprint("cubic penalty term: %12.6f" % self.E_ii3)
        self.set_Eii()
        return converged

    def preopt_q(self):
        # if we are in the ext_lagr mode skip this
        if self.ext_lagr:
            return 0.0
        self.timer.switch_on("preopt_q")
        self.step += 1
        if self.strategy == "q_solver":
            self.q_solver()
        elif self.strategy == "sd":
            self.steepest_descent()
        elif self.strategy == "lbfgs":
            self.lbfgs_scipy()
        else:
            self.pprint("Unknown strategy specified in your control file: %s" % self.strategy)
        self.timer.switch_to("after solve")
        if self.step > 1:
            if self.E < self.last_E - 100.0:
                # this is a suspicious step
                self.pprint("SUSPICIOUS energy change .. checking eigenvalues of J matrix")
                self.calc_energy_and_potential(Jij=True)
                self.analyze_J_eigenvalues()
        self.last_E = self.E
        self.timer.switch_off()
        return self.E

    def print_charges(self):
        if self.molconst:
            for i,m in enumerate(self.pd.molecules.mlist):
                mq = self.q[m.matoms]
                mQ =mq.sum()
                self.pprint("mol %5d: Q: %10.5f %s" % (i, mQ, numpy.array_str(mq, precision=2, suppress_small=True)))
        else:
            self.pprint(numpy.array_str(self.q, precision=2, suppress_small=True))
        return

############### Extended Lagrangian #################################################
# methods for an extra_system added to pydlpoly

    def set_ext_lagr(self):
        """
        switch to ext_lagr mode
        """
        self.qvel = numpy.zeros([self.natoms], dtype="float64")
        self.qmass = numpy.zeros([self.natoms], dtype="float64")
        self.ekin = 0.0
        self.pd.add_extra_system(self)
        self.preopt_q()
        self.ext_lagr = True
        if self.assign_by_elements:
            alist = self.pd.get_elements()
        else:
            alist = self.pd.get_atomtypes()
        for i in xrange(self.natoms):
            try:
                self.qmass[i] = self.params_extlagr[alist[i]]
            except KeyError:
                self.pprint("No qmass for atom %s in the .qeq File" % alist[i])
                raise IOError
        return

    def unset_ext_lagr(self):
        self.ext_lagr=False
        return


    def get_name(self):
        return "QEq"

    def start_up(self):
        return

    def get_ener_names(self):
        return ["QEq"]

    def get_epot(self):
        """ return total potential energy of QEq which is NOT YET part of the coulomb energy """
        # self.E si the total Qeq energy ... self.E_ij is the "regular" coulomb energy which is computed anyway
        # and self.E_ii are the on-site contributions to be added
        #return self.E_ii, []
        return 0.0, []

    def get_ekin(self):
        """ return kinetic energy of dynamic system """
        # self.ekin = numpy.sum((self.ekin_fh + self.ekin_sh)/2.0)
        return self.ekin

    def get_dof(self):
        if self.constraint:
            self.dof = self.natoms - self.pd.molecules.nmolecules
        else:
            self.dof = self.natoms - 1
        self.constraint = False
        return self.dof

    def vv_fh(self):
        """ do first half of VV step"""
        self.timer.switch_on("vv_fh")
#       NOTE: qpot is the derivative of the QEq-energy with respect to the charge (possibly with the constraint included)
#              in the sciyp.fmin_l_bfgs_b that is what we need (the gradient)
#              in case of steepest descent (see corresp method) or in extended lagrangian we need the force
#              which is the negative of the gradient"""
        self.qvel -= 0.5*(self.tstep/self.qmass)*self.qpot
#       now compute the new charges ("postitions")
        self.q += self.tstep*self.qvel
        self.update_q()
        self.timer.switch_off()
        return

    def calc_energy(self):
        """ compute current energy and forces """
        self.timer.switch_on("calc_energy_extlagr")
        #self.zero_Eii()
        self.calc_energy_and_potential(constrain=True, debug=None)
        self.timer.switch_off()
        return

    def vv_sh(self):
        """ do second half of VV step """
        self.timer.switch_on("vv_sh")
        self.qvel -= 0.5*(self.tstep/self.qmass)*self.qpot
        self.ekin = 0.5*numpy.sum(self.qvel*self.qvel*self.qmass)
        self.timer.switch_off()

        #self.pprint("Qpot Second Half Step")
        #self.pprint(self.qpot)
        #self.print_charges()
        return

############### additional stuff for ACKS2 ################################

    def compute_acks2_X(self, rcut= 6.0):
        rcut2 = rcut*rcut
        assert self.acks2, "Initalize with acks2=True"
        bcond = self.pd.get_bcond()
        assert bcond<3, "Not implemented for tricilinc systems"
        if bcond>0:
            cell=self.pd.get_cell().diagonal()
            icell = 1.0/cell
        self.X = numpy.zeros([self.natoms, self.natoms], dtype="float64")
        if self.assign_by_elements:
            alist = self.pd.get_elements()
        else:
            alist = self.pd.get_atomtypes()
        alist = map(string.strip, alist)
        alist = map(string.upper, alist)
        # simple non-parallel loop ... improve later with a simple python level parallelization
        # this is a full N2 distance calc ... better use neigbor lists here in the future
        xyz = self.pd.get_xyz()
        for i in xrange(self.natoms-1):
            # get distances to all higher index atoms
            ixyz = xyz[i]
            dxyz = ixyz-xyz[i+1:]
            # fix PBC (cubic or orthorombic only)
            if bcond>0:
                dxyz -= cell*numpy.around(icell*dxyz)
            d2 = (dxyz*dxyz).sum(axis=1)
            incutoff = numpy.less(d2, rcut2)
            for j in xrange(self.natoms-i-1):
                if incutoff[j]:
                    d = numpy.sqrt(d2[j])
                    jp = j+i+1
                    par = self.params_acks2[alist[i]+":"+alist[jp]]
                    Xij = par[0]*numpy.exp(-par[1]*d)
                    # Xij goes on the off-diagonals negative and on the diagonal it adds up negative
                    self.X[i,jp] =  Xij
                    self.X[jp,i] =  Xij
                    self.X[i,i]  -= Xij
                    self.X[jp,jp]-= Xij
        #replaces intramolecular parts with matrix given in init
        for i in range(self.natoms):
            if i%self.molatoms == 0:
                for j in range(self.molatoms):
                    self.X[i+j][i:i+self.molatoms] = self.intrax[j][:]
        for i in range(self.natoms):
                self.Xii[i] = self.X[i][i]
        return self.X


    #NEWNEW
    def calc_acks2_energy(self,u,sign=1.0):
        u_all = numpy.sum(u)
        u -= u_all/self.natoms
        self.calc_acks2_upot(u)
        self.E_acks2 = sum(0.5*u*self.upot-0.5*self.q*u)
        return sign*self.E_acks2

    def calc_acks2_upot(self,u,sign=1.0):
        self.upot = numpy.zeros([self.natoms], dtype="float64")
        self.upot[::] = -self.q
        self.upot += self.Xij.mult_vec(u)
        return sign*self.upot
    
    def calc_acks2_u(self):
        
        results = sciopt.fmin_l_bfgs_b(self.calc_acks2_energy, self.u, fprime = self.calc_acks2_upot, args = (-1.0,), maxiter=10000000)

        self.u = results[0]
        
        print "###########################"
        print results
        print self.u
        print sum(self.u)
        print "###########################"

        return self.u
        
    def calc_acks2_X(self,non_sparse=True):
        if non_sparse:
            self.X= self.Xij.calc()
        elif not non_sparse:
            self.Xij.calc()
            self.X = self.Xij.get_nonsparse()
        return self.X

    def sd_acks2(self):
        step = 0
        self.conv = 0.001
        self.maxiter = 10000000
        E = self.calc_acks2_energy(self.u,sign=-1.0)
        rms, converged = self.chk_conv_acks2()
        while (not converged and (step < self.maxiter)):
            step += 1
            self.u += self.sd_step*self.upot
            u_all = numpy.sum(self.u)
            self.u -= u_all/self.natoms
            E = self.calc_acks2_energy(self.u,sign=-1.0)
            rms, converged = self.chk_conv_acks2()
            if (self.prntlvl > 1) and (step%self.report_step == 0) :
                self.pprint("  step %d : E = %12.6f, rms = %12.6f, sum of u: %12.6f" % \
                            (step, E, rms, numpy.sum(self.u)))
        if self.prntlvl > 0:
            self.pprint("  sd converged after %d steps : E = %12.6f, rms = %12.6f, sum of u: %12.6f" % \
                            (step, E, rms, numpy.sum(self.u)))
        print self.u
        print "######sum######"
        print sum(self.u)
        return converged

    def chk_conv_acks2(self):
        rms = numpy.sqrt(numpy.vdot(self.upot,self.upot)/self.natoms)
        if rms < self.conv:
            # if self.prntlvl > 0: self.pprint("QEq converged .. rms is %12.6f" % rms)
            return rms, True
        else:
            return rms, False


############### additional stuff for parametrization ################################

    def add_charge_site(self, Qxyz, Q):
        # this function modifies the current parameter set (En's only)
        # to mimick the presence of a single point charge at point xyz with charge Q
        # this charge can not be changed
        # the original parmeter set is kept and with the remove_charge_site function
        # it can be recovered
        self.En_keep = self.En.copy()
        xyz = self.pd.get_xyz()
        sigma = self.pd.get_sigmas()
        # compute distances
        d = xyz-Qxyz
        d = numpy.sqrt(numpy.sum(d*d,axis=1))
        for i in xrange(self.natoms):
            self.En[i] += pi4e0 * Q/d[i]*erf(d[i]/sigma[i])
        return

    def remove_charge_site(self):
        self.En = self.En_keep
        return

    def solve_charge_site(self, Qxyz, Q):
        # solve qeq equations assuming the Jij matrix has alread been computed
        # under the presence of a charge site
        # first compute the interaction with the extra charge
        xyz = self.pd.get_xyz()
        sigma = self.pd.get_sigmas()
        d = xyz-Qxyz
        d = numpy.sqrt(numpy.sum(d*d,axis=1))
        En = self.En + (pi4e0*Q/d*erf(d/sigma))
        self.H = numpy.ones([self.natoms+1,self.natoms+1],dtype="float64")
        if self.acks2: self.H = numpy.zeros([(self.natoms+1)*2,(self.natoms+1)*2],dtype="float64")
        self.x = numpy.zeros([self.natoms+1],dtype="float64")
        if self.acks2: self.x = numpy.zeros([(self.natoms+1)*2],dtype="float64")
        self.x[0:self.natoms] = -En
        if self.use_core_charge:
            self.x[0:self.natoms] -= self.core_qpot
        self.x[self.natoms] = self.Qtot
        self.H[0:self.natoms,0:self.natoms] = self.Jij[:,:]
        self.H[self.natoms,self.natoms] = 0.0
        self.H[self.natoms,:self.natoms] = 1.0
        self.H[:self.natoms,self.natoms] = 1.0
        if self.acks2:
            self.compute_acks2_X()
            self.H[:self.natoms,self.natoms+1:2*self.natoms+1] = -numpy.eye(self.natoms)
            self.H[self.natoms+1:2*self.natoms+1,:self.natoms] = -numpy.eye(self.natoms)
            self.H[self.natoms+1:2*self.natoms+1,self.natoms+1:2*self.natoms+1] = self.X[:]
            self.H[2*self.natoms+1,self.natoms+1:2*self.natoms+1] = 1.0
            self.H[self.natoms+1:2*self.natoms+1,2*self.natoms+1] = 1.0
        sol = numpy.linalg.solve(self.H,self.x)

        q = sol[0:self.natoms]
        Eij = 0.5*numpy.sum(self.Jij*numpy.outer(q,q))
        Eii = numpy.sum(En*q)
        if self.use_core_charge:
            # in this case Eij is Eij_val only -> add the other contribs
            Eij_coreval = numpy.sum(self.core_qpot*q)
            Eij_core = self.pd.get_core_energy()
            Eij += Eij_coreval +Eij_core
            # in this case the interaction of the extra charge with the valence distribs is in Eii
            #        however, the extra charge <-> core interaction (does not depend on qval) needs to be
            #        added
            Eii += pi4e0*Q*numpy.sum(self.core_charges/d)
        if self.acks2:
            u = sol[self.natoms+1:2*self.natoms+1]
            Eacks2 = numpy.sum(u*(-q+0.5*self.X*u))
            Eii += Eacks2
        return Eii+Eij


    def manual_output_selection(self):
        negative_only=True
        evals,vect = numpy.linalg.eigh(self.Jij)
        for i,e in enumerate(evals):
            print "Eigenvalue %d" % (i)
            print e
            v = vect[:,i]
            if e < 0.0:
                self.pprint("%10.3f  : %s" % (e, numpy.array_str(v, precision=2, suppress_small=True)))
            else:
                if not negative_only:
                    self.pprint("%10.3f  : %s" % (e, numpy.array2string(v, precision=2, suppress_small=True)))

        self.analyze_J_eigenvalues()
        print "###self.En###"
        print  self.En/23.06
        print "###self.Jii###"
        print self.Jii/23.06

        return
