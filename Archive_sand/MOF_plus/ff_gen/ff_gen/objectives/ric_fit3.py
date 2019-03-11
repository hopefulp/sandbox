import numpy as np
import os
import copy
#import hessutils
import pdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from molsys.util.timing import timer, Timer
import logging
import string
from mpi4py import MPI
from ff_gen import hessutils
from ff_gen.refclass import refclass2 as refclass
from ff_gen.objectives.base import base


logger = logging.getLogger("FFgen.ric_fit")

rad2deg = 180.0/np.pi

mapper = {"bnd": "str",
        "ang": "ibe",
        "dih": "tor",
        "oop": "obe",
        }

def scatter(optdats, refdats, labels, xlabel, ylabel, fname, 
        writedat = True, legend = False, confidence=None):
    """
    Method to create scatter plots with matiplotlib comparing
    the ff performance to the corresponding reference data. 
    In addition plain ascii files are created with the corresponding
    data.

    :Parameters:
        - optdat(numpy.ndarray): ff data
        - refdat(numpy.ndarray): ref data
        - xlabel(str): label on x axis
        - ylabel(str): label on y axis
        - fname(str): filename of scatterplot, has to end with .png
    """
    assert type(xlabel) == type(ylabel) == type(fname) == str
    assert len(refdats) == len(optdats)
#    assert np.shape(refdat) == np.shape(optdat)
    assert fname.split(".")[-1] == "png"
    plt.clf()
    for r,o,l in zip(refdats,optdats,labels):
        plt.plot(r, o, linestyle = "none", marker = "o", label = l)
    xmin = np.amin([np.amin(refdats[0]),np.amin(optdats[0])])*0.9
    xmax = np.amin([np.amax(refdats[0]),np.amax(optdats[0])])*1.1
    plt.plot([xmin,xmax],[xmin,xmax])
    if confidence is not None:
        c = confidence
        plt.plot([xmin,xmax],[xmin+c,xmax+c])
        plt.plot([xmin,xmax],[xmin-c,xmax-c])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if legend: plt.legend()
    plt.savefig(fname)
    if writedat:
        with open(fname.split(".")[0]+".dat", "w") as f:
            for i in range(len(refdats[0])):
                f.write("%5d %12.6f %12.6f\n" % (i, refdats[0][i], optdats[0][i]))
    return


class ric_fit(base):
    """
    class to compute redundant internal coordinates (ric)
    by using the ric addon of the molclass and to compute deviations from a corresponding
    reference system.

    :Parameters:
        - tag (str): tag of the reference info in the corresponding\
            hdf5 file, it is assumed that there is Group called "hessians" in the root\
            subgroup of the hdf5 file
        - disp (float, optional): displacement used for the finite difference\
            hessian, defaults to 0.002A
        - weights (dict, optional): dictionary holding the weights for the\
            individual ric contributions, defaults all besides obe, lbe and eck to 1.0
        - only_diag (bool, optional): Flag to specify if only the diagonal elements of the hessian\
            should be included to the objective function, defaults to True
        - catchneg (bool, optional): Flag to specify if models with negtive frequencies should be
            penalized, defaults to True
        - doublefd (bool, optional): Flag to specify if double or single sided FD hessian should
            be computed, defaults to True
"""

    
    def __init__(self, pd, reffile, tag, start_dir, disp = 0.001, doublefd = True,tresh = 0.001,
            wstr = 1.0, wibe = 1.0, wobe=0.0, wlbe = 0.0, wtor = 1.0, weck = 0.0, whes = 1.0, wstress = 0.0,
            only_diag = True, catchneg = False, maxneg = 0, latmin = False, fullric = True, 
            log=True, vweights = None, lindict = {}, mpi_comm = None, out = None):
        assert type(disp)    == float
        assert type(tresh)    == float
        # init objective base class
        super(ric_fit,self).__init__(pd,reffile,tag,start_dir,mpi_comm,out)
        self.name = 'RicFit'
        # set log and and timer settings
        self.log = log
        self.timer = Timer()
        # set optimizer settings
        self.tresh = tresh
        if self.backend == "lammps":
            self.opt = self.pd.MIN_cg
            self.latopt = self.pd.LATMIN_boxrel
            if doublefd == True:
                self.hess  = hessutils.doublehessian(self.pd, lmp=True)
            else:
                self.hess  = hessutils.singlehessian(self.pd, lmp=True)
        else:
            self.opt = self.pd.MIN_lbfgs
            self.latopt = self.pd.LATMIN_sd
            if doublefd == True:
                self.hess  = hessutils.doublehessian(self.pd, lmp=False)
            else:
                self.hess  = hessutils.singlehessian(self.pd, lmp=False)
        # set finite difference hessian settings
        self.disp = disp
        self.pd.mol.set_real_mass()
        self.massmatrix = hessutils.hessutils.get_massmatrix(self.pd.get_natoms(),
                self.pd.mol.get_mass())
        self.catchneg = catchneg
        self.maxneg = maxneg
        # init ric addon
        self.pd.mol.addon('ric')
        self.pd.mol.ric.setup_rics(full = fullric, lindict=lindict)
        # set unit prefactors for the different contributions
        self.fact = {"str": 1.0, "ibe": 1.0, "obe":1.0, "lbe":1.0, "tor": 1.0, 
                "eck":1.0, "hes": 1.0/143.88, "stress":1.0/143.88}
        # set the weights for the individual contributions to the overall objective function
        self.weights = {}
        self.weights['str']=wstr
        self.weights['ibe']=wibe
        self.weights['obe']=wobe
        self.weights['lbe']=wlbe
        self.weights['tor']=wtor
        self.weights['eck']=weck
        self.weights['hes']=whes
        self.weights['stress']=wstress
        self.oweight = np.zeros(len(self.pd.mol.ric.active_rics)+2)
        for i,k in enumerate(self.pd.mol.ric.active_rics+["hes", "stress"]):
            self.oweight[i] = self.weights[k]
        # handle the inclusion of stress in the objective function for periodic systems
        if wstress > 0.0:
            self.include_stress = True
        else:
            self.include_stress = False
        self.latmin = latmin
        if self.include_stress: assert self.latmin == False
        if self.latmin: assert self.include_stress == False
        self.cellforce = np.zeros([3,3],dtype = float)
        # generate the reference information
        self.generate_reference(self.pd.get_cell())
        # handle setup of the individual weights per ric
        if vweights is None: 
            self.generate_vweights()
        else:
            self.read_vweights(os.path.join(self.start_dir,vweights))
        self.set_weights(norm = True, only_diag = only_diag)
        self.set_multiplicities()
        return

    def generate_reference(self, cell):
        """
        Method to setup the internal data  structures concerning the reference information.

        :Parameters:
            - cell(numpy.ndarray): cell tensor
        """
        self.refxyz = copy.deepcopy(self.ref(info = 'coord', branch = 'hessians', tag = self.tag))
        self.refhes = copy.deepcopy(self.ref(info = 'hessian', branch = 'hessians', tag = self.tag))
        # check if structure is periodic
        try:
            self.refcell= copy.deepcopy(self.ref(info = 'cell', branch = 'hessians', tag = self.tag))
            self.periodic = True
            if 'periodicity' in self.ref.f['hessians'][self.tag]['cell'].attrs:
                self.periodicity = self.ref.f['hessian'][self.tag]['cell'].attrs['periodicity']
            else:
                self.periodicity = 3
#            import pdb; pdb.set_trace()
        except:
            logger.info('No reference cell given --> not a periodic system')
            self.periodic = False
            self.periodicity = 0
            assert self.include_stress == False, "Fitting the stress is not possible for a nonperiodic system"
            self.refcell = cell
        self.construct_and_project(self.refcell, self.refxyz, self.refhes)
        self.ref_val = {}
        for tr,lk in self.pd.mol.ric.val_dict.items():
            self.ref_val[tr] = copy.deepcopy(lk())
        # if structure is periodic and lammps is used as backend refcell and refhess has to be
        # overwritten due to the rotation of the box in the triclinic case
        if self.backend == 'lammps' and self.periodic:
            self.refcell = self.pd.ff2lmp.rotate_cell(self.refcell)
            invcell = np.linalg.inv(self.refcell)
            frac = np.dot(self.refxyz, invcell)
            self.refxyz = np.dot(frac, self.refcell)
        return

    @timer("construct and project")
    def construct_and_project(self, cell, xyz, hes = None):
        """
        Method to construct the Wilson B matrix for a given geometry and
        to invert it, if hessian is given also the hessian is projected

        :Parameters:
            - cell(numpy.ndarray): cell tensor
            - xyz(numpy.ndarray): cartesian coordinates
            - hes(numpy.ndarray, optional): hessian of the system, defaults to None
        """
        with self.timer("construct"): self.pd.mol.ric.construct_b_matrix(cell, xyz)
        with self.timer("invert"): 
            invb = self.pd.mol.ric.invert_b_matrix()
            #if rankb < self.pd.mol.natoms: raise ValueError("Rank too small!")
        with self.timer("project"):
            if type(hes) != type(None): self.pd.mol.ric.project_hessian(hes)
        return

    def set_multiplicities(self):
        """
        Method to set the multiplicities for the torsion values.
        """
        self.pd.mol.ff.ric.compute_rics()
        self.multiplicities = np.ones(self.pd.mol.ric.num_torsion)
        for i,r in enumerate(self.pd.mol.ff.ric_type["dih"]):
            j, j_glob = self.pd.mol.ric.map_ric("tor",r)
            if j == None: continue
            m = r.value[1]
            if m == None:
                self.pprint("Multiplicity for dihedral %s is not defined" % r)
                self.multiplicities[j] = 1
            else:
                self.multiplicities[j] = r.value[1]
        return

    def get_symmetry(self):
        ff = self.pd.mol.ff
        self.syms = {
                'str':{},
                'ibe':{},
                'obe':{},
                'tor':{},
                'lin':{},
                }
        ics = ["bnd", "ang", "dih", "oop"]
        #ics = ["bnd"]
        for ic in ics:
            for i,pi in enumerate(ff.parind[ic]):
                # get corresponding ric
                ric = ff.ric_type[ic][i]
                # we ignore linear bends
                if ic == "ang" and ric.lin == True: continue
                # map this ric
                j, j_glob = self.pd.mol.ric.map_ric(mapper[ic],ric)
                # get potential name
                name = pi[0]
                # add j, jglob to the corresponding point in syms dict
                if name not in self.syms[mapper[ic]]:
                    self.syms[mapper[ic]][name]= {
                            'j':[j],
                            'j_glob':[j_glob],
                            'opt':{'val_avg':None,
                                'val_std':None,
                                'hes_avg':None,
                                'hes_std':None,},
                            'ref':{'val_avg':None,
                                'val_std':None,
                                'hes_avg':None,
                                'hes_std':None,}
                            }
                else:
                    self.syms[mapper[ic]][name]['j'].append(j)
                    self.syms[mapper[ic]][name]['j_glob'].append(j_glob)
        return

    def compute_averages(self, datatype = 'ref'):
        assert datatype in ['ref', 'opt']
        if datatype == 'ref': 
            data = self.ref_val
        else:
            data = {}
            for tr,lk in self.pd.mol.ric.val_dict.items():
                data[tr] = copy.deepcopy(lk())
        sym = self.syms
        self.sym_couple = {}
        for ic in sym.keys():
            for name in sym[ic].keys():
                j = sym[ic][name]['j']
                j_glob = sym[ic][name]['j_glob']
                if ic =='tor':
                    per = 2*np.pi/np.mean(self.multiplicities[j])
                    sym[ic][name][datatype]['val_avg'] = np.mean(abs(data[ic][j])%per)
                    sym[ic][name][datatype]['val_std'] = np.std(abs(data[ic][j])%per)
                else:
                    sym[ic][name][datatype]['val_avg'] = np.mean(data[ic][j])
                    sym[ic][name][datatype]['val_std'] = np.std(data[ic][j])
                sym[ic][name][datatype]['hes_avg'] = np.mean(data['hes'].diagonal()[j_glob]) 
                sym[ic][name][datatype]['hes_std'] = np.std(data['hes'].diagonal()[j_glob]) 
        return

    def print_averages(self, ic, dt = 'hessian', prefix = 'ref'):
        assert dt in ['hessian','geom']
        assert ic in ['str','ibe', 'tor', 'obe']
        assert prefix in ['ref','opt']
        if self.mpi_rank != 0: return
        fact = {"str": 1.0, "ibe": rad2deg, "obe":rad2deg, "lbe":rad2deg, "tor": rad2deg, 
                "eck":1.0, "hes": 1.0/143.88}
        buffer = ''
        sym = self.syms
        buffer += '### %s ###\n' % ic
        dtype = [('name', 'S200'),('avg',float),('sig',float)]
        data = []
        namefmt = 0
        for name in sorted(sym[ic].keys()):
            if len(name) > namefmt: namefmt=len(name)
            j = sym[ic][name]['j']
            j_glob = sym[ic][name]['j_glob']
            hes_avg = sym[ic][name][prefix]['hes_avg']  
            hes_std = sym[ic][name][prefix]['hes_std'] 
            val_avg = sym[ic][name][prefix]['val_avg'] * fact[ic]
            val_std = sym[ic][name][prefix]['val_std'] * fact[ic]
            if dt == 'hessian':
                data.append((name, hes_avg,hes_std))
            else:
                data.append((name, val_avg,val_std))
        nfmt = '%-'+str(namefmt)+'s'
        data = np.array(data,dtype=dtype)
        np.savetxt('%s-%s-%s.csv' % (prefix, ic, dt),data,delimiter=';', fmt=[nfmt,'%12.6f','%12.6f' ])
        return data

    def write_vweights(self):
        """
        Method to write individual weights per ric to a file called
        weights.dat
        """
        with open('weights.dat', 'w') as f:
            varpotnums = self.pd.mol.ff.par.variables.varpotnums(self.pd.mol.ff)
            ics = ["bnd", "ang", "dih", "oop"]
            for ic in ics:
                f.write("%s\n" % ic)
                for v, w in self.vweights.items():
                    if v[0] == ic:
                        f.write("%-80s %3i %12.6f %12.6f\n" % (v[1],  varpotnums[v], w[0], w[1]))
        return

    def read_vweights(self, fname):
        """
        Method to read the indiviual weights per ric from a file

        :Parameters:
            - fname(str): Name of the file
        """
        ics = ["bnd", "ang", "dih", "oop"]
        self.vweights = {}
        with open(fname, 'r') as f:
            for line in f.readlines():
                sline = line.split()
                if len(sline) == 1 and sline[0] in ics:
                    ic = sline[0]
                else:
                    self.vweights[(ic,sline[0])] = [float(sline[2]), float(sline[3])]


    def generate_vweights(self):
        """
        Method to generate the individual weights per varied ric type,
        every term get total weight of zero.
        """
        self.vweights = {}
        varpots    = self.pd.mol.ff.par.variables.varpots
        for i in varpots:
            self.vweights[i] = [1.0,1.0]
        return

    def set_weights(self, norm = True, only_diag = True, equal = True):
        """
        Method to set the weights.

        :Parameters:
            - norm(bool, optional): if true use normalized weights, if false
                use for every ric a weight of 1.0, optional to True
            - only_diag(bool, optional): only diag flag, defaults to True 
        """
        # allocate individual weight arrays
        self.wgt   = {}
        self.wgt["str"] = np.zeros(self.pd.mol.ric.num_stretch)
        self.wgt["ibe"] = np.zeros(self.pd.mol.ric.num_in_bend)
        self.wgt["obe"] = np.zeros(self.pd.mol.ric.num_out_bend)
        self.wgt["lbe"] = np.zeros(self.pd.mol.ric.num_lin_bend)
        self.wgt["tor"] = np.zeros(self.pd.mol.ric.num_torsion)
        self.wgt["eck"] = np.zeros(self.pd.mol.ric.num_eckart)
        self.wgt["hes"] = np.zeros([self.pd.mol.ric.num_ric, self.pd.mol.ric.num_ric])
        # compute average values in the hessian per ic type
        avg_bnd = np.average(self.ref_val['hes'].diagonal()[self.pd.mol.ric.first_str:self.pd.mol.ric.first_ibe])
        avg_ang = np.average(self.ref_val['hes'].diagonal()[self.pd.mol.ric.first_ibe:self.pd.mol.ric.first_obe])
        avg_oop = np.average(self.ref_val['hes'].diagonal()[self.pd.mol.ric.first_obe:self.pd.mol.ric.first_tor])
        avg_dih = np.average(self.ref_val['hes'].diagonal()[self.pd.mol.ric.first_tor:self.pd.mol.ric.first_lbe])
        if equal:
            hess_facts = {'bnd':1.0,
                    'ang':1.0,
                    'dih':1.0,
                    'oop':1.0}
        else:
            hess_facts = {'bnd':1.0,
                    'ang': avg_bnd/avg_ang,
                    'dih': avg_bnd/avg_dih,
                    'oop': avg_bnd/avg_oop}
        varpotnums = self.pd.mol.ff.par.variables.varpotnums(self.pd.mol.ff)
        ics = ["bnd", "ang", "dih", "oop"]
        for ic in ics:
            for i,r in enumerate(self.pd.mol.ff.ric_type[ic]):
                j, j_glob = self.pd.mol.ric.map_ric(mapper[ic],r)
                ### check for linear bend case
                if ic == "ang" and r.lin == True: 
                    j, j_glob = self.pd.mol.ric.map_ric("lbe",r)
                if j == None: continue 
                pi = self.pd.mol.ff.parind[ic][i]
                gweight = 0.0
                hweight = 0.0
                for p in pi:
                    if (ic, p) in self.vweights:
                        gweight = self.vweights[(ic,p)][0]
                        hweight = self.vweights[(ic,p)][1]
                        if norm: 
                            gweight/=varpotnums[(ic, p)]
                            hweight/=varpotnums[(ic, p)]
                        break
#                for p in pi:
#                    if (ic, p) in varpots:
#                        weight = 1.0
#                        if norm: weight/=varpotnums[(ic, p)]
#                        break
                if ic == "oop":
                    self.wgt[mapper[ic]][j:j+3]  = gweight/3.0
                    self.wgt["hes"][j_glob, j_glob]     = hweight/3.0*hess_facts[ic]
                    self.wgt["hes"][j_glob+1, j_glob+1] = hweight/3.0*hess_facts[ic]
                    self.wgt["hes"][j_glob+2, j_glob+2] = hweight/3.0*hess_facts[ic]
                    #import pdb; pdb.set_trace()
                elif ic == "ang" and r.lin == True:
                    self.wgt["lbe"][j]   = gweight/2.0
                    self.wgt["lbe"][j+1] = gweight/2.0
                    self.wgt["hes"][j_glob,j_glob]     = hweight/2.0*hess_facts[ic]
                    self.wgt["hes"][j_glob+1,j_glob+1] = hweight/2.0*hess_facts[ic]
                else:
                    self.wgt[mapper[ic]][j] = gweight
                    self.wgt["hes"][j_glob,j_glob] = hweight*hess_facts[ic]
        ### only for strbnd
        if only_diag == False:
            for i,r in enumerate(self.pd.mol.ff.ric_type["ang"]):
                j, j_glob = self.pd.mol.ric.map_ric("ibe",r)
                if j == None: continue 
                pi = self.pd.mol.ff.parind["ang"][i]
                for p in pi:
                    if ((("ang", p) in self.vweights) and (self.pd.mol.ff.par["ang"][p][0] == "strbnd")):
                        weight = self.vweights[('ang',p)]/6.0
                        if norm: weight/=varpotnums[("ang", p)]
                        s1, s1_glob = self.pd.mol.ric.map_ric("str", r[:2])
                        s2, s2_glob = self.pd.mol.ric.map_ric("str", r[1:])
                        self.wgt["hes"][j_glob,  s1_glob] += weight
                        self.wgt["hes"][s1_glob, j_glob]  += weight
                        self.wgt["hes"][j_glob,  s2_glob] += weight
                        self.wgt["hes"][s2_glob, j_glob]  += weight
                        self.wgt["hes"][s1_glob, s2_glob] += weight
                        self.wgt["hes"][s2_glob, s1_glob] += weight
        return

    @timer("call")
    def __call__(self):
        """
        Method to compute the fitness function as msd between ric geometry and hessian
        between model and reference.
        """
        ### write variables to the calculators internal data structure
        with self.timer("write params"): self.writer(self.pd)
        ### optimze and calc hessian
        if self.periodic: self.pd.set_cell(self.refcell)
        self.pd.set_xyz(self.refxyz)
        with self.timer("optimize structure"):
            if self.latmin: self.latopt(0.1,0.01)
            self.opt(self.tresh)
            if self.include_stress: self.cellforce = self.pd.get_cellforce()
            #self.opt(self.tresh, verbose = False)
        with self.timer('compute hessian'):
            self.hessian = self.hess(disp = self.disp)
        ### constuct and project
        self.construct_and_project(self.pd.get_cell(), self.pd.get_xyz(), self.hessian)
        msd, msds = self.msd
        if self.catchneg:
            val, vec = np.linalg.eigh(self.hessian*self.massmatrix)
            if np.sum(np.less(val,-9.0e18)) > self.maxneg: 
                self.pprint("WARNING")
                msd = np.inf 
        self.cycle += 1
        if self.log:
            self.pprint(len(msds)*"%12.6f " % tuple(msds))
        return msd

    @property
    def msd(self):
        msd_keys = self.pd.mol.ric.active_rics + ["hes", "stress"]
        msds     = np.zeros(len(msd_keys))
        for i,k in enumerate(msd_keys[:-1]):
            ### for torsions the multiplicity has to be included
            if k == "tor":
                delt  = self.ref_val[k] - self.pd.mol.ric.val_dict[k]()
                per   = 2*np.pi/self.multiplicities 
                wdelt = (self.pd.mol.ric.val_dict[k]()-(self.ref_val[k] - per*np.around(delt/per)))*self.wgt[k]*self.fact[k]
#                wdelt = (abs(abs(self.pd.mol.ric.val_dict[k]())-abs(self.ref_val[k]))%(2*np.pi/self.multiplicities)) \
#                        *self.wgt[k]*self.fact[k]
            else:
                wdelt = (self.pd.mol.ric.val_dict[k]()-self.ref_val[k])*self.wgt[k]*self.fact[k]
            if self.wgt[k].sum() > 0.0:
                msds[i] = (wdelt**2).sum()/self.wgt[k].sum()
        if self.include_stress: 
            msds[-1] = 0.01 * np.sqrt(np.sum(self.cellforce**2)/9.)**4
        self.pprint(np.sqrt(np.sum(self.cellforce**2/9.)))
        msd = np.sum(self.oweight*msds)
        return msd, msds

    def finish(self):
        """
        Method to call after the optimization is finished. It writes out the scatter
        plot and the comparisons between reference and model data.
        """
        ### call after the last __call__ with final params
        ### move to rundir of the calculator
        if self.mpi_rank != 0: return
        retval = os.getcwd()
        os.chdir(self.pd.rundir)
        # write optimized params
        self.pd.mol.ff.write('opt')
        # write vweights
        self.write_vweights()
        #self.pd.mol.ff.write('opt')
        #hessian = self.hess(disp = self.disp)
        if ((self.backend=='lammps') and (self.periodicity < 3)):
            xyz = self.pd.mol.map2image(self.pd.get_xyz())
        else:
            xyz = self.pd.get_xyz()
        ### write molden input files for freq and opt
        hu = hessutils.hessutils(xyz, self.hessian, self.pd.get_elements())
        hu.write_molden_input("opt.freq")
        wn_opt = copy.deepcopy(hu.wn)
        hu = hessutils.hessutils(self.refxyz, self.refhes, self.pd.get_elements())
        hu.write_molden_input("ref.freq")
        wn_ref = copy.deepcopy(hu.wn)
        ### write opt.mfpx
        if self.periodic: self.pd.mol.set_cell(self.pd.get_cell())
        self.pd.mol.set_xyz(xyz)
        self.pd.mol.write("opt.mfpx")
        ### print results in rics
        self.print_all()
        ### print information on averages
        self.get_symmetry()
        self.compute_averages('ref')
        self.compute_averages('opt')
        for ic in ['str', 'ibe', 'tor', 'obe']:
            for t in ['hessian','geom']:
                ref = self.print_averages(ic,t,'ref')
                opt = self.print_averages(ic,t,'opt')
                scatter([opt['avg']], [ref['avg']], ["None"],'ref', 'opt', 'avg_%s-%s.png' % (ic,t))
        ### make compare freqs plot
        scatter([wn_opt], [wn_ref],"None","ref freq [cm$^{-1}$]", "opt freq [cm$^{-1}$]",
                "vib.png")
#        self.analyze_couplings(['n3_c2zn1', 'zn4_n4', 'n3_c2zn1'])
        ### write timings
        self.timer.write_logger(logger.info)
        os.chdir(retval)
        return

        
    def print_geometry(self,plt_scatter = True):
        """
        Method to generate buffer with the difference between the ric and models geometries

        :Parameters:
            -plt_scatter(bool, optional): if true the corresponding scatter plots are created,
            defaults to True
        """
        msd_keys = self.pd.mol.ric.active_rics
        rics     = self.pd.mol.ric.all_rics
        count = 0
        buffer = ""
        buffer += "!!!!!!!!!!!GEOM!!!!!!!!!!!!\n"
        fact = {"str": 1.0, "ibe": rad2deg, "obe":rad2deg, "lbe":rad2deg, "tor": rad2deg, 
                "eck":1.0, "hes": 1.0/143.88}
        scatterfiles = {"str": "all_str-geom.png", "ibe": "all_ibe-geom.png", 
                "tor": "all_tor-geom.png", "obe": "all_obe-geom.png"}
        scatteraxes  = {"str": ("ref bnd length", "opt bnd length"), 
                "ibe": ("ref ang value", "opt ang value"),
                "tor": ("ref dih value", "opt dih value"),
                "obe": ("ref oop value", "opt oop value")}
        for k in msd_keys:
            if k == "tor":
                delt   = self.ref_val[k] - self.pd.mol.ric.val_dict[k]()
                per    = 2*np.pi/self.multiplicities 
                mapped = self.ref_val[k] - per*np.around(delt/per)
                diff   = (self.pd.mol.ric.val_dict[k]()-mapped)*fact[k]
                if plt_scatter and k in scatterfiles.keys():
                    scatter([self.pd.mol.ric.val_dict[k]()*fact[k]], [mapped*fact[k]],["None"],
                            scatteraxes[k][0], scatteraxes[k][1], scatterfiles[k])
            else:
                diff  = (self.ref_val[k] - self.pd.mol.ric.val_dict[k]())*fact[k]
                if plt_scatter and k in scatterfiles.keys():
                    scatter([self.pd.mol.ric.val_dict[k]()*fact[k]], [self.ref_val[k]*fact[k]],["None"],
                            scatteraxes[k][0], scatteraxes[k][1], scatterfiles[k])

                #ratio = abs(1-(self.ref_val[k]/self.pd.mol.ric.val_dict[k]()))*100
            for i,r in enumerate(rics[k]):
                rstring = "%4d %3d %5s" % (count,i, k)+len(r)*" %4d" % tuple(r)
                buffer += rstring + "%12.4f %12.4f %12.4f %12.4f\n" % \
                        (self.ref_val[k][i]*self.fact[k], self.pd.mol.ric.val_dict[k]()[i]*self.fact[k],
                                diff[i],  self.wgt[k][i])
                count += 1
        return buffer

    def print_diagonals(self, plt_scatter=True):
        """
        Method to create buffer with the difference in the diagonal elements of the hessian in rics.

        :Returns:
            -buffer(str): string holding the information
        """
        buffer   = ""
        buffer += "!!!!!!!!!!!DIAG!!!!!!!!!!!!\n"
        msd_keys = self.pd.mol.ric.active_rics
        rics     = self.pd.mol.ric.all_rics
        count = 0
        k = "hes"
        diff  = (self.ref_val[k] - self.pd.mol.ric.val_dict[k]())*self.fact[k]
        ratio = abs(1-(self.ref_val[k]/self.pd.mol.ric.val_dict[k]()))*100
        for k in msd_keys:
            for i,r in enumerate(rics[k]):
                rstring = "%4d %3d %5s" % (count,i,k)+len(r)*" %4d" % tuple(r)
                buffer += rstring + "%12.4f %12.4f %12.4f %12.4f %12.4f\n" % \
                        (self.ref_val["hes"][count,count]*self.fact["hes"], self.pd.mol.ric.val_dict["hes"]()[count,count]*self.fact["hes"],
                                diff[count,count], ratio[count,count], self.wgt["hes"][count,count])
                count += 1
        k = 'hes'
        if plt_scatter:
#            import pdb;pdb.set_trace()
            rang = np.array(range(self.pd.mol.ric.first_str,self.pd.mol.ric.first_ibe))
            scatter([self.pd.mol.ric.val_dict[k]().diagonal()[rang]],
                    [self.ref_val[k].diagonal()[rang]], ["None"], 'ref', 'opt', 'all_str-hessian.png')
            rang = np.array(range(self.pd.mol.ric.first_ibe,self.pd.mol.ric.first_obe))
            scatter([self.pd.mol.ric.val_dict[k]().diagonal()[rang]],
                    [self.ref_val[k].diagonal()[rang]], ["None"], 'ref', 'opt', 'all_ibe-hessian.png')
            rang = np.array(range(self.pd.mol.ric.first_obe,self.pd.mol.ric.first_tor))
            scatter([self.pd.mol.ric.val_dict[k]().diagonal()[rang]],
                    [self.ref_val[k].diagonal()[rang]], ["None"], 'ref', 'opt', 'all_obe-hessian.png')
            rang = np.array(range(self.pd.mol.ric.first_tor,self.pd.mol.ric.first_lbe))
            scatter([self.pd.mol.ric.val_dict[k]().diagonal()[rang]],
                    [self.ref_val[k].diagonal()[rang]], ["None"], 'ref', 'opt', 'all_tor-hessian.png')
        return buffer

    def print_full(self):
        """
        Method to create buffer with the difference of the full hessian in rics.

        :Returns:
            -buffer(str): string holding the information
        """
        buffer = ""
        buffer += "!!!!!!!!!!!FULL!!!!!!!!!!!!\n"
        k = "hes"
        nrics = self.pd.mol.ric.num_ric
        diff  = (self.ref_val[k] - self.pd.mol.ric.val_dict[k]())*self.fact[k]
        ratio = abs(1-(self.ref_val[k]/self.pd.mol.ric.val_dict[k]()))*100
        for i in range(nrics):
            for j in range(nrics):
                buffer += "%3d %3d %12.6f %12.6f %12.6f %12.6f\n" % \
                        (i,j,self.ref_val[k][i,j]*self.fact["hes"], self.pd.mol.ric.val_dict[k]()[i,j]*self.fact[k],
                                diff[i,j], ratio[i,j])
        return buffer
        
    def print_msd(self):
        """
        Method to create buffer with all msd values (bnd, ibe, obe, lbe, tor, hess)

        :Returns:
            -buffer(str): string holding the information
        """
        buffer = ""
        buffer += "!!!!!!!!!!!MDSs!!!!!!!!!!!!\n"
        msd_keys = self.pd.mol.ric.active_rics + ["hes", "stress"]
        msd, msds = self.msd
        for i,k in enumerate(msd_keys):
            buffer+= "%3s:  %12.8f\n" % (k, msds[i])
        buffer += "msd:  %12.8f\n" % msd 
        return buffer

    def print_all(self, file = "result.dat"):
        """
        Method to write information from the print_msd, print_geometry, print_hessian
        and print_diagobals information to a file.

        :Parameters:
            - file(str,optional): name of output file, defaults to result.dat
        """
        if self.mpi_rank != 0: return
        with open(file, "w") as f:
            f.write(self.print_msd())
            f.write(self.print_geometry())
            f.write(self.print_diagonals())
            f.write(self.print_full())

    def analyze_couplings(self, atypes):
        with open('couplings.dat', 'w') as f:
            for i,r in enumerate(self.pd.mol.ff.ric_type["ang"]):
                if map(lambda a: self.pd.mol.atypes[a], r) == atypes: 
                    j, j_glob = self.pd.mol.ric.map_ric('ibe',r)
                    aval = self.ref_val['ibe'][j]
                    js1, js_glob1 = self.pd.mol.ric.map_ric('str',r[:2])
                    js2, js_glob2 = self.pd.mol.ric.map_ric('str',r[1:])
                    sval1 = self.ref_val['str'][js1]
                    sval2 = self.ref_val['str'][js2]
                    f.write("%12.6f %12.6f %12.6f\n" % (aval*rad2deg, sval1, sval2))
