# -*- coding: utf-8 -*-
import numpy
import string
import subprocess
import shutil
import os
import copy

import simple_dialog

PDLPDIR = os.environ["PDLPDIR"]

HA2KCAL = 627.5095
BOHR2A  = 0.5291772

PBC_array = []    
for x in xrange(-1,2):
    for y in xrange(-1,2):
        for z in xrange(-1,2):
            PBC_array.append([x,y,z])
PBC_array = numpy.array(PBC_array, "d")


class turbo_interface:

    def __init__(self, turbodir, defprot, sed=None, restartdir=None, EE=False, alink_assignment = False):
        self.EE = EE
        self.turbodir = turbodir
        self.restartdir=restartdir
        if os.path.exists(defprot):
            self.defprot = os.getcwd()+"/"+defprot
        else:
            print 'No define protocol found'
            raise IOError
#            self.defprot  = PDLPDIR+"/QMMM/"+defprot
        self.sed = sed
        self.do_cosmo = False
        self.skip_QM_energy = False
        self.alink_assignment = alink_assignment
        return
        
    def init(self, pd, QMmol, links):
        self.set_pydlpoly_instance(pd)
        self.set_QMmol(QMmol)
        if self.alink_assignment:
            self.get_links(links)
        else:
            self.set_links(links)
        self.cycle = 0
        self.turbodir = self.pd.rundir+"/"+self.turbodir
        return
        	
    def use_cosmo(self, cosmo_prot):
        self.do_cosmo = True
        self.cosmo_prot = os.getcwd()+"/"+cosmo_prot
        return
	
    def set_pydlpoly_instance(self, pd):
        self.pd = pd
        self.pd.pprint("QMMM Calculation")
        return
        
    def set_QMmol(self, QMmol):
        self.QMmol = QMmol
        self.QMatoms = self.pd.mol.mols[QMmol]
        self.nqmatoms = len(self.QMatoms)
        # set the QMelems
        self.QMelems = []
        for a in self.QMatoms:
            self.QMelems.append(self.pd.mol.elems[a])
        self.MMatoms = range(self.pd.get_natoms())
        return
    
    def set_links(self, links):
        self.nlinks = len(links)
        self.link_scales          = []
        self.link_MMatoms         = []
        self.link_QMatoms         = []
        self.link_QMatoms_i       = [] 
        self.link_elems           = []
        for l in links:
            qm, mm, scale, elem = l
            qm -= 1
            mm -= 1
            self.link_scales.append(scale)
            self.link_elems.append(elem)
            if not self.QMatoms.count(qm):
                raise ValueError, "The QM atom %d defined in the link %s is not in the QM set" % (qm+1, str(l))
            if self.QMatoms.count(mm):
                raise ValeuError, "The MM atom %d defiend in the link %s is a QM atom" % (mm+1, str(l))
            self.link_MMatoms.append(mm)
            self.link_QMatoms.append(qm)
            self.link_QMatoms_i.append(self.QMatoms.index(qm))
        # append link elems to the QMelems list
        for i in xrange(self.nlinks):
            self.QMelems.append(self.link_elems[i])
        return

    def get_links(self, links):
        print 'Performing automatic link assignment'
        scale = links[0]
        elem  = links[1]
        self.nlinks = 0
        self.link_scales          = []
        self.link_MMatoms         = []
        self.link_QMatoms         = []
        self.link_QMatoms_i       = [] 
        self.link_elems           = []
        for i in self.QMatoms:
            for j in self.pd.mol.cnct[i]:
                if j not in self.QMatoms:
                    self.nlinks += 1
                    self.link_scales.append(scale)
                    self.link_elems.append(elem)
                    self.link_MMatoms.append(j)
                    self.link_QMatoms.append(i)
                    self.link_QMatoms_i.append(self.QMatoms.index(i))
        for i in xrange(self.nlinks):
            self.QMelems.append(self.link_elems[i])
        return

    def setup(self):
        if os.path.isdir(self.turbodir):
            raise IOError, "TURBODIR exists!!"
        self.pd.pprint("QMMM: starting setup!")
        os.mkdir(self.turbodir)
        self.get_QMxyz()
        self.write_coord()
        self.MMatoms = range(self.pd.get_natoms())
        for i in xrange(self.pd.get_natoms()):
            if i in self.QMatoms:
                self.MMatoms.remove(i)
        if (self.EE==True):
            self.find_frags()
        os.chdir(self.turbodir)
        if self.restartdir:
            self.pd.pprint("QMMM: doing a restart and taking all files from %s" % self.restartdir)
            #os.system("cp %s/* ." % self.restartdir)
            os.system("cp %s/control ." % self.restartdir)
            os.system("cp %s/basis ." % self.restartdir)
            os.system("cp %s/auxbasis ." % self.restartdir)
            os.system("cp %s/alpha ." % self.restartdir)
            os.system("cp %s/beta ." % self.restartdir)
            os.system("cp %s/mos ." % self.restartdir)
            os.system("cp %s/master ." % self.restartdir)
            #os.system("cp %s/energy ." % self.restartdir)
            #os.system("cp %s/gradient ." % self.restartdir)        
        else:
            self.pd.pprint("QMMM: start talking with define ...")
            #setup = simple_dialog.simple_dialog(self.defprot, "define.log")
            #setup.talk_to("define")
            os.system('define < %s > define.log' % self.defprot)
	    if self.do_cosmo:
		    cosmo_setup = simple_dialog.simple_dialog(self.cosmo_prot, "cosmo.log") 
		    cosmo_setup.talk_to("cosmoprep")
            shutil.copy("control", "control.from_define")
            self.pd.pprint("QMMM: done with define. wrap up")
            if self.sed:
                #f = open("control.sed","w")
                #retcode = subprocess.call(["sed", self.sed, "control"], stdout=f)
                #f.close()
                os.system("sed '%s' control > control.sed" % self.sed)
                shutil.copy("control.sed", "control")
        os.chdir(self.pd.start_dir)        
        self.pd.pprint("QMMM: setup done")
        return

            
    def __call__(self, force):
        e = 0.0
        if self.skip_QM_energy :
            self.get_QMxyz()
            self.write_coord()
            #self.write_charges()  ### Testzwecke
            return e
        self.cycle += 1
        self.get_QMxyz()
        self.write_coord()
        if (self.EE==True):
            self.write_charges() 
        pddir = os.getcwd()
        os.chdir(self.turbodir)
        if os.path.exists("control"):
            #self.cycle += 1
            ferr = open("turbo_blurb.dump","w")
            f = open("turbo_%d.out" % self.cycle, "w")
            retcode = subprocess.call(["ridft"], stdout=f, stderr=ferr)
            if retcode != 0: raise IOError, "QMMM Error: RIDFT failed!!!"
            if force:
                retcode = subprocess.call(["rdgrad"], stdout=f, stderr=ferr)
                if retcode != 0: raise IOError, "QMMM Error: RDGRAD failed!!!"
            f.close()
            ferr.close()
            if force:
                self.read_gradient()
                self.set_MM_force()
            else:
                self.read_energy()
            e = self.QMenergy
        os.chdir(self.pd.start_dir)        
        return e
                
    def get_QMxyz(self):
        """ extract xyz coords and pick the proper image for all atoms (take first atom as refernce) """
        xyz = []
        all_xyz = self.pd.get_xyz()
        # first collect all true QM atoms
        for a in self.QMatoms:
            xyz.append(all_xyz[a].tolist())
        # now add the link MM atoms (first we need to get them into the local image before we scale the position
        for i in xrange(self.nlinks):
            xyz.append(all_xyz[self.link_MMatoms[i]].tolist())
        xyz = numpy.array(xyz)
        # compute distances of all atoms to first
        if self.pd.mol.boundarycond == 3:
            raise valueError, "QMMM for triclinic cells currently not implemented"
        cell = self.pd.get_cell().diagonal()
        d = xyz - xyz[0]
        too_small = numpy.less(d, -0.5*cell)
        xyz = numpy.where(too_small, xyz+cell, xyz)
        too_large = numpy.greater_equal(d, 0.5*cell)
        xyz = numpy.where(too_large, xyz-cell, xyz)
        # now all atoms are in a common frame - now we need to scale the nlink last atoms
        for l in xrange(self.nlinks):
            link_xyz = xyz[self.nqmatoms+l]
            qm_xyz   = xyz[self.link_QMatoms_i[l]]
            link_xyz[:] = qm_xyz+(link_xyz-qm_xyz)/self.link_scales[l]
        self.QMxyz = xyz
        return xyz
        
        
    def write_coord(self):
        #if os.path.isfile("coord"): ### um Änderungen der Ladungsverteilung zu berechenen
         #   os.system("cp coord coord_%d" % self.cycle)
        f = open(self.turbodir+"/coord", "w")
        f.write("$coord\n")
        for i,e in enumerate(self.QMelems):
            c = self.QMxyz[i]/BOHR2A
            f.write("%15.10f %15.10f %10.10f %s\n" % (c[0], c[1], c[2], e))
        f.write("$end\n")
        f.close()
        return
        
    def read_gradient(self):
        f = open("gradient", "r")
        line = string.split(f.readline())
        if line[0] != "$grad":
            raise IOError, "Not a TURBOMOLE gradient file!"
        line = string.split(f.readline())
        cycle = string.atoi(line[2])
        while cycle < self.cycle:
            for i in xrange(2*(self.nqmatoms+self.nlinks)):
                f.readline()
            line = string.split(f.readline())
            cycle = string.atoi(line[2])
        # now we are at the right spot .. get the energy first
        self.QMenergy = string.atof(line[6])*HA2KCAL
        # skip the coords
        for i in xrange(self.nqmatoms+self.nlinks):
            f.readline()
        gradient = []
        for i in xrange(self.nqmatoms+self.nlinks):
            rawline = f.readline().replace("D", "E")
            line = string.split(rawline)
            gradient.append(map(string.atof, line[:3]))
        gradient = numpy.array(gradient)
        # convert to angstrom and kcal/mols
        gradient *= HA2KCAL/BOHR2A
        # redistribute the link atoms gradient to the qm and mm atoms
        #    the MM atoms gets 1/scale times the gradient whereas the QM atom gets 1-1/scale
        for l in xrange(self.nlinks):
            s = 1.0/self.link_scales[l]
            i = self.nqmatoms+l
            qmgrad = (1-s)*gradient[i]
            gradient[i] *= s
            gradient[self.link_QMatoms_i[l]] += qmgrad
        # remove the net force (does not remove the roational net force!)
        gradient -= gradient.sum(axis=0)/(self.nqmatoms+self.nlinks)
        self.QMgradient = gradient
        return
        
    def read_energy(self):
        f = open("energy", "r")
        line = string.split(f.readline())
        if line[0] != "$energy":
            raise IOError, "Not a TURBOMOLE energy file!"
        line = string.split(f.readline())
        cycle = string.atoi(line[0])
        while cycle < self.cycle:
            line = string.split(f.readline())
            cycle = string.atoi(line[0])
        # now we are at the right spot .. get the energy
        self.QMenergy = string.atof(line[1])*HA2KCAL
        return

    def set_MM_force(self):
        force = numpy.zeros([self.pd.natoms,3], "d")
        for i in xrange(self.nqmatoms):
            force[self.QMatoms[i],:] = -self.QMgradient[i,:]
        for l in xrange(self.nlinks):
            force[self.link_MMatoms[l],:] = -self.QMgradient[self.nqmatoms+l,:]
        self.pd.add_force(force)
        return
        
    ######################################################################
    #             2nd derivatives added by JPD (turbo part)
    ######################################################################    
            
    def aoforce(self): ### Aenderungen in Zeile 285 notwendig in Bezug auf $noproj
        self.pd.pprint("preparing the control file")
        os.system("sed -e '/^$ricore/d' -e '/^$rij/d' -e '/^$jbas/d' -e 's/$marij/$noproj/g'  control > control.aoforce")
        shutil.copy("control.aoforce", "control")
        ferr = open("force_blurp.dump","w")
        f = open("force.out","w")
        self.pd.pprint("aoforce is running")
        retcode = subprocess.call(["aoforce"], stdout=f, stderr=ferr)
        if retcode != 0: raise IOError, "QMMM Error: aoforce failed!!!"
        f.close()
        ferr.close()
        self.pd.pprint("aoforce done")
        return
        
    def read_hessian(self): ### Turbomole input für QMMM, schon symmetrisiert, nicht massengewichtet ###
        self.aoforce()
        f = open("control", "r") 
        n = self.nqmatoms + self.nlinks
        found = False
        while not found:
          line = string.split(f.readline())
          if line[0] == "$nprhessian": found = True 
        hessian = numpy.zeros((3*n,3*n), dtype = "float64")
        line = string.split(f.readline()) 
        r = int(round(((3*n/5.0 - 3*n/5)*5)))
        if r == 0:
            a = 0
        else:
            a = 1
        for x in xrange(3*n): 
            for y in xrange((3*n/5) + a): ### a reinbasteln
                l = len(line)
                if y+1 <= 3*n/5:
                    values = numpy.array(map(string.atof, line[l-5:]))
                    hessian[x,(y*5):(y*5)+5] = values
                else:
                    values = numpy.array(map(string.atof, line[l-r:]))
                    hessian[x,(y*5):(y*5)+r] = values ###????
                line = string.split(f.readline())
        f.close
        hessian = (hessian * 2240.877106) ### scale factor from a.u. to kcal/(A**2 * mol)  ###
        return hessian    

        
    def build_transhessian(self):
        # go to turbo rundirectory 
        #pddir = os.getcwd()
        os.chdir(self.turbodir)
        qm_hessian = self.read_hessian()
        os.chdir(self.pd.start_dir)        
        #os.chdir(pddir)
        # done ...
        self.pd.pprint("translation of QM hessian started")
        ndeg = 3 * self.pd.get_natoms()
        trans_hess = numpy.zeros((ndeg, ndeg), dtype="float64") 
        ### Übergabe der reinen QM Felder ###
        for i in xrange (self.nqmatoms):
          ii = self.QMatoms[i]
          for j in xrange (self.nqmatoms):
            jj = self.QMatoms[j]
            trans_hess[3*jj:(3*jj)+3,3*ii:(3*ii)+3] = qm_hessian[3*j:(3*j)+3,3*i:(3*i)+3]
        ### Übergabe und Umrechnung der QM-link Felder ###
        for k in xrange (self.nlinks):
          s = 1.0/self.link_scales[k]
          kk = self.link_QMatoms[k]
          aa = self.link_MMatoms[k]
          for l in xrange (self.nqmatoms):
            ll = self.QMatoms[l]
            ### QMatoms ### # korrekt
            trans_hess[3*ll:(3*ll)+3,3*kk:(3*kk)+3]+=(1-s)*qm_hessian[3*l:(3*l)+3,(3*self.nqmatoms+(3*k)):(3*self.nqmatoms+(3*k))+3]
            trans_hess[3*kk:(3*kk)+3,3*ll:(3*ll)+3]+=(1-s)*qm_hessian[(3*self.nqmatoms+(3*k)):(3*self.nqmatoms+(3*k))+3,3*l:(3*l)+3]
            ### MMatoms ### #korrekt
            trans_hess[3*ll:(3*ll)+3,3*aa:(3*aa)+3]+=s*qm_hessian[3*l:(3*l)+3,(3*self.nqmatoms+(3*k)):(3*self.nqmatoms+(3*k))+3] 
            trans_hess[3*aa:(3*aa)+3,3*ll:(3*ll)+3]+=s*qm_hessian[(3*self.nqmatoms+(3*k)):(3*self.nqmatoms+(3*k))+3,3*l:(3*l)+3]
      ### Übergabe und Umrechnung der link-link Felder
        for m in xrange (self.nlinks):
          s = 1.0/self.link_scales[m]
          mm = self.link_QMatoms[m]
          oo = self.link_MMatoms[m]
          for n in xrange (self.nlinks):
            nn = self.link_QMatoms[n]
            pp = self.link_MMatoms[n]
            ### QMatoms ###
            trans_hess[3*nn:(3*nn)+3,3*mm:(3*mm)+3] += ((1-s)**2)*qm_hessian[3*(self.nqmatoms + n):3*(self.nqmatoms + n)+3,3*(self.nqmatoms + m):3*(self.nqmatoms + m)+3] 
           ### MMatoms ###
            trans_hess[3*pp:(3*pp)+3,3*oo:(3*oo)+3] += (s**2)*qm_hessian[3*(self.nqmatoms + n):3*(self.nqmatoms + n)+3,3*(self.nqmatoms + m):3*(self.nqmatoms + m)+3] 
          ### QM-MMatoms ###
            trans_hess[3*pp:(3*pp)+3,3*mm:(3*mm)+3] += s*(1-s)* qm_hessian[3*(self.nqmatoms + n):3*(self.nqmatoms + n)+3,3*(self.nqmatoms + m):3*(self.nqmatoms + m)+3]
            trans_hess[3*nn:(3*nn)+3,3*oo:(3*oo)+3] += s*(1-s)* qm_hessian[3*(self.nqmatoms + n):3*(self.nqmatoms + n)+3,3*(self.nqmatoms + m):3*(self.nqmatoms + m)+3]
        self.pd.pprint("translation of QM hessian done")
        return trans_hess   
          
          
    ######################################################################
    #             electrostatic embedding added by JPD
    ######################################################################
    
    def get_shortest(self,atomnumber):
        cell = self.pd.get_cell().diagonal()
        xyz = self.pd.get_xyz()
        xyz_mm_images = xyz[atomnumber,:] + (PBC_array*cell) ### 27 Abbilder eines MMatoms
        d = xyz_mm_images[:,numpy.newaxis,:] - self.QMxyz
        r = numpy.sqrt(numpy.sum(d*d,axis = 2))
        shortest = numpy.unravel_index(numpy.argmin(r),r.shape)
        shortest_distance = numpy.amin(r)
        shortest_xyz = xyz_mm_images[shortest[0],:]
        return shortest_distance, shortest_xyz
       
    def get_chargeatoms(self,cutoff): ###entfernen nur für hättig
        cell = self.pd.get_cell().diagonal()
        #cutoff = 12.0
        realcutoff = cutoff
        if cutoff >= cell[0]/2:
            raise ValueError, "cutoff > half times the lattice constant"
        xyz = self.pd.get_xyz()
        chargelist = []
        charge_xyz = [] 
        for i in self.MMatoms:
            shortest_distance, shortest_xyz = self.get_shortest(i)
            if shortest_distance <= cutoff and (i not in chargelist) and (i not in self.link_MMatoms):
                chargelist.append(i)
                charge_xyz.append(shortest_xyz.tolist())               
        ### fragmentatome mit herrein nehmen ### an welcher Stelle der loop?
                for j in range(len(self.fragments)):
                    if i in self.fragments[j]:
                        for k in range(len(self.fragments[j])):
                            if self.fragments[j][k] not in chargelist and (self.fragments[j][k] in self.MMatoms):
                                shortest_distance, shortest_xyz = self.get_shortest(self.fragments[j][k])
                                if shortest_distance >= cell[0]/2:
                                    raise ValueError, "cutoff > half times the lattice constant"
                                else:
                                    chargelist.append(self.fragments[j][k])
                                    charge_xyz.append(shortest_xyz.tolist())
                                    if shortest_distance > realcutoff:
                                        realcutoff = shortest_distance
        charge_xyz = numpy.array(charge_xyz, dtype="float64")
        return chargelist, charge_xyz, realcutoff
      
    def write_charges(self): 
        #if os.path.isfile("charges"): ### um Änderungen der Ladungsverteilung zu berechenen
         #   os.system("cp charges charges_%d" % self.cycle)
        for cutoff in range(8,78):
            l,xyz,realcutoff = self.get_chargeatoms(cutoff)
            xyz /= BOHR2A
            charges = self.pd.get_charges()
            sigmas = self.pd.get_sigmas()
            sigmas /= BOHR2A
            a = 1/(sigmas**2)
            f = open(self.turbodir+"/charges_%d" % cutoff, "w")
            f.write("$point_charges gaussians list\n")
            charge = 0.0
            for i in xrange(len(l)):
                f.write("%15.10f %15.10f %15.10f %15.10f %15.10f\n" % (xyz[i,0], xyz[i,1], xyz[i,2], charges[l[i]], a[l[i]]))
                charge += charges[l[i]]
            f.write("$end\n")
            f.write ("sum")
            f.write("%15.10f\n" % (charge))
            f.write ("cutoff")
            f.write("%15.10f" % (realcutoff))        
            f.close()
        return   
        
    def read_frags(self):
        f = open("../neutral_groups", "r")
        self.frags = {} ### dictionary
        for line in f:
            linecontent = string.split(line)
            self.frags[linecontent[0]] = linecontent[1:]
        print self.frags
        return
    
    def bond_walker(self, atom, l, typelist, remain_list, atypes, conn):
        bonds = conn[atom]
        # check atoms for still being in the remain_list
        # print "bonds for atom %d :  %s" % (atom, str(bonds))
        for a in copy.copy(bonds):
            # print "atom %d is %d times in remaining_list" % (a, remain_list.count(a))
            if not remain_list.count(a):
                # ok this atom is no longer in the remaining lists ... we have it
                # already and can remove it
                bonds.remove(a)
        # print "bonds after removal  %s" % (str(bonds))
        # now check if it is part of the fragment
        for a in bonds:
            at = atypes[a]
            if typelist.count(at):
                # print ("atom %d with type %d is part of the fragment" % (a, at))
                # print "current fragmentlist:"
                # print l
                if l.count(a):
                    pass
                    # print "atom %d is already in the fragment list!!" % a
                    #print "found a ring closure"
                else:
                    # this atom is part of the fragment
                    # remove from the remain_list and add to l
                    # then walk on
                    remain_list.remove(a)
                    l.append(a)
                    self.bond_walker(a, l, typelist, remain_list, atypes, conn)
        return 
    
    
    def find_frags(self):
        self.read_frags()
        self.atypes = self.pd.get_atomtypes()
        nnets=1
        self.frag_names = []
        self.frag_atoms = int(self.pd.get_natoms())*[0,]
        self.fragments  = []
        self.mols       = []
        self.frag_count = {}
        c_fnumber = 0
        for fn in self.frags.keys(): self.frag_count[fn] = 0
        remaining_atoms = range(self.pd.natoms)
        remaining_atoms.reverse()
        # now do this as long as we have remaining atoms
        while (len(remaining_atoms) > 0):
            c_atom = remaining_atoms.pop()
            c_atype = self.atypes[c_atom]
            # now we need to find the fragment type to which this atom belongs
            c_fname = None
            c_mol   = 0
            for fn in self.frags.keys():
                if self.frags[fn].count(c_atype):
                    c_fname = fn
                    if (nnets>1):
                        # find out which net/mol this is
                        found = 0
                        for i in xrange(nnets):
                            if (c_atom < mol_index[i]) and not found:
                                found = 1
                                c_mol = i
            if not c_fname:
                print "Could not find atomtype %d for atom %d in frags" % (c_atom, c_atype)
                return
            c_flist = self.frags[c_fname]
            self.frag_names.append(c_fname)
            self.frag_count[c_fname] += 1
            self.frag_atoms[c_atom] = c_fnumber
            # now follow bonds and check if the atom is also part of the current frag
            other_atoms = []
            self.bond_walker(c_atom, other_atoms, c_flist, remaining_atoms, self.atypes, self.pd.mol.cnct)
            c_fragment = [c_atom] + other_atoms
            self.fragments.append(c_fragment)
            self.mols.append(c_mol)
            for a in other_atoms: self.frag_atoms[a] = c_fnumber
            c_fnumber += 1
        return
            
            
    
    def do_scan(self): ### Methode, die nur für das Hättig Vertiefungspraktikum wichtig war!
        a = 14# atomnumber des Cs des CO Moleküls
        b = 15# atomnumber des Os des CO Moleküls
        c = 1# atomnumber des Cus des CO Moleküls
        r = self.QMxyz[a] - self.QMxyz[c]
        d = numpy.sqrt(numpy.sum(r*r))
        rn = (0.05/d)*r
        xyz_orig = self.QMxyz
        xyz = xyz_orig.copy()
        for i in range(41,51):
            xyz = xyz_orig.copy()
            xyz[a] = xyz[a] + (i)*rn
            xyz[b] = xyz[b] + (i)*rn
            nr = xyz[a] - xyz[c]
            nd = numpy.sqrt(numpy.sum(nr*nr))
            f = open(self.turbodir+"/coord_+%d" % i, "w")
            f.write("$coord\n")
            for k,e in enumerate(self.QMelems):
                d = xyz[k]/BOHR2A
                f.write("%15.10f %15.10f %10.10f %s\n" % (d[0], d[1], d[2], e))
            f.write("$end\n")
            f.write ("distance")
            f.write("%15.10f" % (nd))        
            f.close()
#        for j in range(20):
#            xyz = xyz_orig.copy()
#            xyz[a] = xyz[a] - (j+1)*rn
#            xyz[b] = xyz[b] - (j+1)*rn
#            nr = xyz[a] - xyz[c]
#            nd = numpy.sqrt(numpy.sum(nr*nr))
#            z = j+1
#            f = open(self.turbodir+"/coord_-%d" % z, "w")
#            f.write("$coord\n")
#            for k,e in enumerate(self.QMelems):
#                d = xyz[k]/BOHR2A
#                f.write("%15.10f %15.10f %10.10f %s\n" % (d[0], d[1], d[2], e))
#            f.write("$end\n")
#            f.write ("distance")
#            f.write("%15.10f" % (nd))        
#            f.close()
        return
        
        
