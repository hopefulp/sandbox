#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import numpy
import string
import assign_FF
#import bb
import tools as tools
from itertools import combinations
import pydlpoly
import unit_cell
import molsys.util.elems as elements
import copy

class topoFF:

    def __init__(self, net):
        self.mol = assign_FF.mol(is_master = True)
        self.mol.verbose = False
        self.net = net
        self.conn = copy.deepcopy(net.net.conn)
        self.mol.pass_molsys(net.net)
        self.bbs = self.net.bbs
        self.mol.find_internals(do_smallring_check = False, do_lin_check = False,
                sqp = False, ltorsion = False)
        return

    def get_bonds(self, k = 0.1, ccdist=1.5):
        for a in self.mol.bonddata:
            bbs = []
            for i in string.split(a,':'):
                bbs.append(self.bbs[i])
            #bbs = string.split(a,':')
            d = 0.0 
            for i,s in enumerate(bbs):
                if s.specific_conn != None:
                    if i == 0:
                        type = s.specific_conn.index(bbs[1].label)
                    else:
                        type = s.specific_conn.index(bbs[0].label)
                    ind  = s.connectors_type.index(type)
                    #d += elements.call(s.conn_elems[ind], 'cov_radii') + s.conn_dist[ind]
                    d += elements.cov_radii[s.conn_elems[ind]] + s.conn_dist[ind]
                else:
                    d += elements.cov_radii[s.conn_elems[0]]+\
                            numpy.mean(s.conn_dist)
                    #d += elements.call(s.conn_elems[0], 'cov_radii') +\
                    #        numpy.mean(s.conn_dist)
#            self.mol.bonddata[a] = [k, numpy.mean(self.bbs[bbs[0]].conn_dist)+
#                    numpy.mean(self.bbs[bbs[1]].conn_dist)+ccdist]
            self.mol.bonddata[a] = [k,d]
        return

    def get_combos(self,l1,l2):
        combos = []
        for i in l1:
            for j in l2:
                combos.append([i,j])
        return combos


    def get_angles(self,k = 1.0, noangle = False):
        if noangle:
            for a in self.mol.angledata:
                self.mol.angledata[a]=[0.0,0.0]
        for a in self.mol.angledata:
            types = string.split(a,':')
            #bb = self.bbs[string.split(a,':')[1]]
            bb = self.bbs[types[1]]
            ctypes = [types[0],types[2]]
            geom = tools.geometry(numpy.append(bb.connector_xyz,[[0,0,0]],axis =0))
            angles = []
            if bb.specific_conn != None:
                connecs = []
                for i in ctypes:
                    type = bb.specific_conn.index(i)
                    #connecs += numpy.where(numpy.array(bb.connectors_type)==type)[0].tolist()
                    connecs.append(numpy.where(numpy.array(bb.connectors_type)==type)[0].tolist())
                if ctypes[0]==ctypes[1]:
                    combos = list(combinations(connecs[0],2))
                else:
                    combos = self.get_combos(connecs[0],connecs[1])
                    ### build combos direct
            else:
                connecs = range(len(bb.connectors))
                combos = list(combinations(connecs,2))
            nconnecs = len(bb.connectors)
            #for i in combinations(connecs,2):
            for i in combos:
                atoms = [i[0],nconnecs,i[1]]
                angle = geom.get_angle(atoms)
                angles.append(angle)
                #angles.append(geom.get_angle([i[0],nconnecs,i[1]]))
            if numpy.std(angles) < 5.0:
                self.mol.angledata[a] = [k, numpy.mean(angles)]
            else:
                self.mol.angledata[a] = map(float,bb.angleterm[1:])

    def get_torsions(self):
        for a in self.mol.dihedraldata:
            self.mol.dihedraldata[a]=[0.0,0.0,0.0,0.0]

    def get_oops(self):
        for a in self.mol.oopdata:
            self.mol.oopdata[a]=[0.0]

    def __call__(self):
        key = tools.keycreator()
        key.pass_mol(self.mol)
        key.get_atoms()
        #key.get_vdws()
        self.get_bonds()
        self.get_angles(noangle = True)
        self.get_torsions()
        self.get_oops()
        key.write_key(self.net.name+'.key', [key.buffer_nonbonded(),
            key.buffer('bond', 'bonddata'),
            key.buffer('angle', 'angledata'),
            key.buffer('torsion','dihedraldata'),
            key.buffer('opbend','oopdata')
            ])
        return

    def pd_init(self):
        self.net.net.write(self.net.name+'_init.txyz', ftype='txyz')
        self.pd = pydlpoly.pydlpoly(self.net.name)
        self.pd.control["cut"] = 0.5
        self.pd.control["no link"] = " "
        self.pd.setup(xyz = self.net.name+"_init.txyz", key = self.net.name+".key", 
                local = True, bcond = 3)
        return

    def prescale(self):
        self.pd.LATMIN_sd(0.1,0.1,3000, latonly = True)
        self.pass_coordinates()
        return

    def optimize(self):
        self.pd.LATMIN_sd(0.1,0.1,3000, latonly = False)
        self.pass_coordinates()
        return

    def pass_coordinates(self):
        self.net.net.xyz = self.pd.get_xyz()
        self.net.net.cell = self.pd.get_cell()
        self.net.net.cellparams = unit_cell.abc_from_vectors(self.pd.get_cell())
        return














