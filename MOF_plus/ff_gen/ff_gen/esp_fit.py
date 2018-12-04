#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 12:00:57 2017

@author: johannes
"""

from horton import load_h5
import numpy
import molsys
import molsys.stow as stow
import click
import copy
import cPickle as pickle
import cma

class espfit(object):
    
    def __init__(self, mol, cost):
        self._m    = mol
        self._cost = cost
        self._A    = copy.copy(cost._A)
        self._B    = copy.copy(cost._B)
#        self._V    = vol
        self.equivs_act = False
        self.fixed_act  = False
        self.equivs = {}
        self.fixed  = {}
        self.var_map = range(self._m.natoms)
        return
    
    
    @staticmethod
    def calculate_dipole(mol, charges):
        mol.set_real_mass()
        com = mol.get_com()
        d = (com-mol.xyz)*charges[:,numpy.newaxis]
        d = numpy.sum(d, axis = 0)
        m = numpy.linalg.norm(d)
        return m/0.2082, d/m

    def apply_fixed(self):
        assert self.equivs_act == False
        fixes = self.fixed
        A = self._A
        B = self._B
        for i,c in fixes.items(): B-= c*A[:,i]
        keep = []
        maplist = range(len(A))
        count = 0
        for i in range(len(A)):
            if fixes.keys().count(i) == 0:      
                keep.append(i)
                maplist[i]=count
                count += 1
            else:
                maplist[i]=None 
        self.fixed_act = True
        self._A = A[:,keep]
        self._B = B
        self.var_map = maplist
        return
    
    def apply_equivs(self):
        A = self._A
        B = self._B
        equivs = self.equivs
        maplist = self.var_map
        for i,j in equivs.items(): A[:,maplist[j]]+=A[:,maplist[i]]
        keep = []
        for i in range(len(A)):
            if maplist[i] != None and i not in equivs.keys():
                keep.append(maplist[i])
        self._A = A[:,keep]
        self._B = B
        self.equiv_act = True
        return 
    
    def apply_lagrangian(self, qtot = 0.0):
        assert self.equivs_act == self.fixed_act == False
        A = self._A
        B = self._B
        B = numpy.append(B,qtot)
        A = numpy.append(A, numpy.ones([1,numpy.shape(A)[1]]), axis = 0)
        A = numpy.append(A, -numpy.ones([numpy.shape(A)[0],1]), axis = 1)
        A[-1,-1] = 0.0
        self._A = A
        self._B = B
        return
    
    def reconstruct(self, variable, qtot = 0.0):
        equivs = self.equivs
        fixed = self.fixed
        natoms = self._m.natoms
        charges = numpy.zeros(natoms)
        variable = variable.tolist()
        variable.reverse()
        lvar = []
        for i in xrange(natoms):
            if i in fixed.keys():
                charges[i]=fixed[i]
            else:
                if i not in equivs.keys(): 
                    charges[i] = variable.pop()
                    lvar.append(i)
        ### apply equivs
        for i,j in equivs.items():
            lvar.append(i)
            charges[i]=charges[j]
        ### apply total carge constraint
        if qtot != None:        
            aid = numpy.zeros(natoms)
            qdiff = numpy.sum(charges)-qtot
            for i in lvar: aid[i] = qdiff/len(lvar)
            charges -= aid
        return charges
    
    @staticmethod
    def format_charges(mol, cost, vol, **kwargs):
        buffer = "                "
        buffer += (len(kwargs)) * "%12s" % tuple(kwargs.keys()) +"\n"
        elems = mol.elems
        atypes = mol.atypes
        for i in range(len(kwargs.values()[0])):
            c = []
            for k,j in kwargs.items(): c.append(j[i])
            buffer += ("%3i %3s %10s" + len(kwargs)*"%12.6f" +"\n") % tuple([i+1, elems[i], atypes[i]]+c)
        buffer += "worst: %12.6f\n" % cost.worst()
        for k,v in kwargs.items():
            val, worst, dip = espfit.get_stats(mol, cost, v)
            buffer +="-----------%s------------\n" % k
            buffer += "cost:  %12.6f\n" % val
            buffer += "rmsd:  %12.6f\n" % (val/vol)**0.5
            buffer += "rrmsd: %12.6f\n" % (val/cost.worst())**0.5
            buffer += "dipm:  %12.6f\n" % dip[0]
            buffer += "dipv:  %12.6f %12.6f %12.6f\n" % (dip[1][0],dip[1][1],dip[1][2])
            buffer += "qtot:  %12.9f\n" % numpy.sum(v)
        return buffer
    
    @staticmethod
    def get_stats(mol, cost, charges):
        return cost.value(charges), cost.worst(), espfit.calculate_dipole(mol, charges)
    
    def todb(self, parnames, verbose = True):
        result = {}
        for i,p in enumerate(parnames):
            ### check if parname was fitted
            if i not in self.fixed.keys():
                if i not in self.equivs.keys():
                    result[p] = self.charges[i]
        if verbose:
            for k,v in result.items(): print "%-35s %12.6f" % (k,v)
        return result
      
    def anafit(self, qtot = None, reconstruct = True):
        res = numpy.linalg.lstsq(self._A, self._B)[0]
        self.charges = self.reconstruct(res,qtot)
        if reconstruct:
            return self.charges
        else:
            return res
        
    def numfit(self, guess, qtot = 0.0, reconstruct = True):
        es = cma.CMAEvolutionStrategy(guess,0.01, {"maxiter":500})
        es.optimize(self.concost, verb_disp = 0)
        res = es.result[0]
        self.charges = self.reconstruct(res,0.0)
        if reconstruct:
            return self.charges
        else:
            return res
        
    def concost(self, params, qtot = 0.0):
        charges = self.reconstruct(params,qtot = qtot)
        return self._cost.value(charges)

### for fitting missing charges with new assignment scheme
@click.command()
@click.argument("mfpx", type = click.Path(exists=True))
@click.argument("fcost", type = click.Path(exists=True))
@click.argument("const", type = click.Path(exists=True))
def get_charges(mfpx, fcost, const):
    m      = molsys.mol.fromFile(mfpx)
    cost   = load_h5(fcost)['cost']
    vol    = load_h5(fcost)['used_volume']
    fitdat = pickle.load(open(const, "rb"))
    esp1 = espfit(m,cost)
    esp1.apply_lagrangian()
    esp1.anafit()
    esp2 = espfit(m,cost)
    esp2.equivs = fitdat["equivs"]
    esp2.fixed  = fitdat["fixed"]
    esp2.apply_fixed()
    esp2.apply_equivs()
    guess = esp2.anafit(qtot = 0.0, reconstruct = False)
    anacharges = copy.copy(esp2.charges)
    esp2.numfit(guess)
    print espfit.format_charges(m, cost,vol,free = esp1.charges, anafit = anacharges, numfit = esp2.charges)
    return esp2.todb(fitdat["parnames"])
    

if __name__ == "__main__":
    get_charges()
    


#get_charges("7.mfpx","cost.h5","espfit_azone.pickle")
#get_charges("7.mfpx",pname = "espfit.pickle")
#m      = molsys.mol.fromFile("7.mfpx")
#cost   = load_h5("cost.h5")['cost']
#vol    = load_h5("cost.h5")['used_volume']
#fitdat = pickle.load(open("espfit_azone.pickle", "rb"))
#esp2 = espfit(m,cost)
#esp2.equivs = fitdat["equivs"]
#esp2.fixed  = fitdat["fixed"]
#esp2.apply_fixed()
#esp2.apply_equivs()
#sol = esp2(qtot = None, reconstruct = False)
#print sol
#print espfit.format_charges(m, cost,vol,const = esp2.charges)

#es = cma.CMAEvolutionStrategy(sol,0.01, {"maxiter":500})

#es.optimize(esp2.concost, verb_disp = 1)

#sol2 = es.result()[0]

#print espfit.format_charges(m, cost,vol, cma = esp2.reconstruct(sol2, qtot = 0.0))
#esp2.todb(fitdat["parnames"])

#esp1 = espfit(m,cost)
#esp1()

#es = cma.CMAEvolutionStrategy(esp1.charges,0.01, {"maxiter":500})

#es.optimize(esp1.concost, verb_disp = 1)

#sol2 = es.result()[0]
#print espfit.format_charges(m, cost,vol,const = esp1.charges)
#get_charges("7.mfpx")

#cost = load_h5("cost.h5")['cost']
#vol  = load_h5("cost.h5")['used_volume']

#fixed = dict(zip(range(22,32), 10*[0.12]))
#fixed.update(dict(zip([0,2,3,4,5,13,14,15,16,17],10*[-0.12])))

#equivs = {19:18,
#          21:20,
#          8:7,
#          12:1,
#          11:10,
#          6:1,
#          9:1
#          }
#m = molsys.mol.fromFile("final.xyz")

#esp1 = espfit(m,cost)
#esp1()

#esp2 = espfit(m,cost)
#esp2.equivs = equivs
#esp2.fixed  = fixed
#esp2.apply_fixed()
#esp2.apply_equivs()
#esp2()

#esp3 = espfit(m,cost)
#esp3.equivs = equivs
#esp3.fixed  = fixed
#esp3.apply_fixed()
#esp3.apply_equivs()
#esp3(qtot = 0.0)

#kollman = load_h5("7.hdf5")["hessians"]["primary"]["charges"]


#print espfit.format_charges(m, cost,vol,f = esp1.charges, c = esp2.charges,c0 = esp3.charges, k = kollman)

