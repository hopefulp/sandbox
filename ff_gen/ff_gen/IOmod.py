import numpy
import string
import unit_cell
import tools
import copy
import elements as elemclass
import itertools

try:
    import pandas
except:
    print("pandas package could not be imported, zmat conversion will not work")
    print("try installing pandas package with: pip install pandas")

try:
    import chemcoord
except:
    print("chemcoord package could not be imported, zmat conversion will not work")
    print("try installing chemcoord package with: pip install chemcoord")

class io:
    """
    class which is used for reading geometric information without setting up an pd instance.
    In the near future this should be replaced by an overall mol or net class which serves as backend
    for reading, storing and generation of gemoetric data and is used in weaver, pydlpoly and ffgen
    """

    def __init__(self, master = True):
        """
        initializes an io instance
        """
        self.is_master = master
        self.cell = None
        return

    ####### Read Routines #######

    def read_xyz(self, fname):   
        """
        reads geometric information from an xyz file

        :Parameters:
            -fname      (str) : name of the xyz file

        :Returns:
            -elements   (list): list of element symbols
            - xyz       (ndarray): N,3 numpy array with the geometry
        """

        with open(fname, 'r') as f:
            natoms = int(f.readline().split()[0])
            f.readline() #comment
            xyz = numpy.zeros((natoms, 3))
            elements = []
            for i in range(natoms):
                line = f.readline().split()
                elements.append((line[0]).upper())
                xyz[i,:] = line[1:4]
            #if len(f.readlines()) > 0:
            #    raise Exception('File '+fname+' longer than natoms+2: xyz movie?')
        #f = open(fname, 'r')V
        #natoms = string.atoi(string.split(f.readline())[0])
        #f.readline()
        #xyz = numpy.zeros((natoms, 3))
        #elements = []
        #for i in range(natoms):
        #    line = string.split(f.readline())
        #    elements.append(string.upper(line[0]))
        #    xyz[i,:] = line[1:4]
        for i in range(len(elements)):
            if len(elements[i]) == 1:
                elements[i] += ' '
        elements = map(string.lower, elements)
        elements = map(string.capitalize, elements)
        self.set_elements(elements)
        self.set_natoms(natoms)
        self.set_masses()
        self.set_vdwr()
        self.set_xyz(xyz)
        return elements, xyz


    def read_tinker_xyz(self, name):
        elems = []
        xyz = []
        atypes = []
        cnct = []
        f = open(name, "r")
        lbuffer = string.split(f.readline())
        natoms = string.atoi(lbuffer[0])
        if len(lbuffer) > 1 and lbuffer[1] != 'molden':
            boundarycond = 3
            if lbuffer[1] == "#":
                # read full cellvectors
                celllist = map(string.atof,lbuffer[2:11])
                cell = numpy.array(celllist)
                cell.shape = (3,3)
                cellparams = unit_cell.abc_from_vectors(cell)
            else:
                cellparams = map(string.atof, lbuffer[1:7])
                cell = unit_cell.vectors_from_abc(cellparams)
            if ((cellparams[3]==90.0) and (cellparams[4]==90.0) and (cellparams[5]==90.0)):
                boundarycond=2
                if ((cellparams[0]==cellparams[1])and(cellparams[1]==cellparams[2])and\
                    (cellparams[0]==cellparams[2])):
                        boundarycond=1
        for i in xrange(natoms):
            lbuffer = string.split(f.readline())
            xyz.append(map(string.atof, lbuffer[2:5]))
            elems.append(string.lower(lbuffer[1]))
            t = lbuffer[5]
            atypes.append(t)
            cnct.append((numpy.array(map(string.atoi, lbuffer[6:]))-1).tolist())
        # done: wrap up
        xyz = numpy.array(xyz)
        for i in range(len(elems)):
            if len(elems[i]) == 1:
                elems[i] += ' '
        self.set_elements(elems)
        self.set_natoms(natoms)
        self.set_masses()
        self.set_vdwr()
        self.set_atomtypes(atypes)
        self.set_xyz(xyz)
        self.set_cnct(cnct)
        if 'cell' not in locals():
            cell = None
            boundarycond = None
        self.set_cell(cell)
        self.set_boundarycond(None)
        return elems, atypes, xyz, cnct, cell, boundarycond

    def add_atom(self, elem, atype, coords, connectors):
        self.elements.append(elem)
        xyz = list(self.get_xyz())
        xyz.append(coords)
        self.set_xyz(numpy.array(xyz))
        self.atypes.append(atype)
        self.natoms += 1
        self.set_masses()
        self.set_vdwr()
        for i in connectors:
            self.cnct[i].append(self.natoms-1)
        self.cnct.append(connectors)
        return

    def add_virtual(self, connectors):
        xyz = self.get_com(connectors)
        self.add_atom('xx', '100', xyz, connectors)
        return

    def get_com(self, idx):
        com = 0.0
        mass = 0.0
        for i in idx:
            com += self.masses[i]*self.xyz[i,:]
            mass += self.masses[i]
        com /= mass
        return com

    def set_masses(self):
        self.masses = []
        for i in range(self.natoms): 
            self.masses.append(float(elemclass.call(self.elements[i], 'mass')))
        self.masses = numpy.array(self.masses)
        return

    def get_masses(self):
        return self.masses
    
    def set_vdwr(self):
        self.vdwr = []
        for i in range(self.natoms): 
            self.vdwr.append(float(elemclass.call(self.elements[i], 'vdw_radii')))
        self.vdwr = numpy.array(self.vdwr)
        return

    def get_vdwr(self):
        return self.vdwr

    def set_natoms(self, natoms):
        self.natoms = natoms

    def get_natoms(self):
        return self.natoms

    def get_elements(self):
        return self.elements

    def set_elements(self, elements):
        self.elements = elements

    def get_atomtypes(self):
        return self.atypes

    def set_atomtypes(self, atypes = None):
        if atypes == None:
            self.atypes = numpy.shape(self.xyz)[0]*['0']
        else:
            self.atypes = atypes

    def get_xyz(self):
        return self.xyz

    def set_xyz(self, xyz):
        self.xyz = xyz

    def get_cnct(self):
        return self.cnct

    def set_cnct(self, cnct):
        self.cnct = cnct

    def get_cell(self):
        return self.cell

    def set_cell(self, cell):
        self.cell = cell

    def get_boundarycond(self):
        return self.boundarycond

    def set_boundarycond(self, boundarycond):
        self.boundarycond = boundarycond


    ####### Write Routines #######

    def write_xyz(self, fname, elems = None, xyz=None):
        """ writing out a simple xyz file for DEBUG purposes
        """
        if elems == None: elems = self.elements
        if xyz  == None: xyz = self.xyz
        natoms = numpy.shape(xyz)[0]
        if self.is_master:
            fxyz = open(fname,"w")
            fxyz.write("%d\n" % natoms)
            fxyz.write("PyDLP xyz file\n")
            for i in xrange(natoms):
                fxyz.write("%s %12.6f %12.6f %12.6f\n" % (elems[i], xyz[i,0], xyz[i,1], xyz[i,2]))
            fxyz.close()
        return
    
    def write_tinker_xyz(self, fname, elems=None, atypes=None, xyz=None, 
            cnct = None, cell=None,fullcell=False, mode="w"):
        """ if fullcell is not set we write the cell parameters (a,b,c,alpha,beta,gamma) and rotate the system properly
             if fullcell is true we write the complete cell vectors (9 values) to the first line
             
             :Parameters:
                 - moldenr : replaces atom type strings with integers to ensure readability by molden when set to True
             """
        if elems == None: elems = self.elements
        if atypes== None: atypes= self.atypes
        if xyz  == None: xyz = self.xyz
        if cnct == None: cnct= self.cnct
        #if cell == None: cell= self.cell
        if numpy.equal(cell, None): cell = self.cell
        natoms = numpy.shape(xyz)[0]
        if type(cell) != type(None):
            cellparams = unit_cell.abc_from_vectors(cell)
        if type(cell) != type(None):
            if fullcell:
                pass
            else:
                rotcell = unit_cell.vectors_from_abc(cellparams)
                # compute fractional coords in old cell
                inv_cell = numpy.linalg.inv(cell)
                abc = numpy.dot(xyz, inv_cell)
                # compute cartesian coordinates in new cellvectors from fractional coords
                xyz = numpy.dot(abc,rotcell)
        if self.is_master:
            f = open(fname, mode)
            if type(cell) != type(None):
                if fullcell:
                    f.write(("%5d # "+9*"%10.4f "+"\n") % tuple([natoms]+cell.ravel().tolist()))                
                else:
                    f.write("%5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n" % tuple([natoms]+cellparams))
            else:
                f.write("%5d \n" % natoms)
            for i in xrange(natoms):
                line = ("%3d %-3s" + 3*"%12.6f" + " %5s") % \
                    tuple([i+1]+[elems[i]]+ xyz[i].tolist() + [atypes[i]])
                conn = (numpy.array(cnct[i])+1).tolist()
                if len(conn) != 0:
                    line += (len(conn)*"%7d") % tuple(conn)
                f.write("%s \n" % line)
            f.close()
        return

    def generate_cnct(self, xyz = None, elements=None, tresh = 0.1, cell = None, remove_duplicates = False):
        if xyz == None:
            xyz = self.xyz
        if elements  == None:
            elements = self.elements
        if type(cell) == type(None): 
            cell = self.cell
        if type(cell) != type(None):
            converted_cell = unit_cell.abc_from_vectors(cell)
            cell_abc = converted_cell[:3]
            cell_angles = converted_cell[3:]
            if cell_angles[0] != 90.0 or cell_angles[1] != 90.0 or cell_angles[2] != 90.0:
                inv_cell = numpy.linalg.inv(cell)
        natoms = numpy.shape(xyz)[0]
        cnct = []
        for i in xrange(natoms):
            a = xyz - xyz[i]
            if type(cell) != type(None):
                if cell_angles[0] == 90.0 and cell_angles[1] == 90.0 and cell_angles[2] == 90.0:
                    a -= cell_abc * numpy.around(a/cell_abc)
                else:
                    frac = numpy.dot(a, inv_cell)
                    frac -= numpy.around(frac)
                    a = numpy.dot(frac, cell)
            dist = ((a**2).sum(axis=1))**0.5 # distances from i to all other atoms
            cnct_local = []
            for j in xrange(natoms):
                if i != j and dist[j] <= self.get_covdistance([elements[i],elements[j]])+tresh:
                    cnct_local.append(j)
            cnct.append(cnct_local)
        self.set_cnct(cnct)
        return cnct

    def get_covdistance(self, elems):
        return elemclass.call(elems[0],'cov_radii')+elemclass.call(elems[1], 'cov_radii')
 
    def remove_duplicates(self, thresh = 0.1, cell = None):
        """
        removes duplicates in orthorombic systems
        Parameters:
          - thresh : distance in A, below which the atoms will be regarded as duplicates.
        """
        xyz = self.xyz
        elements = self.elements
        natoms = xyz.shape[0]
        if type(cell) == type(None): 
            cell = self.cell
        if type(cell) != type(None):
            converted_cell = unit_cell.abc_from_vectors(cell)
            cell_abc = converted_cell[:3]
            cell_angles = converted_cell[3:]
        if cell_angles[0] != 90.0 or cell_angles[1] != 90.0 or cell_angles[2] != 90.0:
            inv_cell = numpy.linalg.inv(cell)
        # find duplicates
        duplicates = []
        for i in xrange(natoms):
            a = xyz - xyz[i]
            if type(cell) != type(None):
                if cell_angles[0] == 90.0 and cell_angles[1] == 90.0 and cell_angles[2] == 90.0:
                    a -= cell_abc * numpy.around(a/cell_abc)
                else:
                    frac = numpy.dot(a, inv_cell)
                    frac -= numpy.around(frac)
                    a = numpy.dot(frac, cell)
            dist = ((a**2).sum(axis=1))**0.5
            for j in xrange(i, natoms):
                if i != j and dist[j] < thresh:
                    duplicates.append(j)
        if duplicates:
            print("Number of duplicates: %d" % len(duplicates))
            new_xyz = numpy.delete(xyz, duplicates, 0)
            new_elements = numpy.delete(elements, duplicates)
            print("Number of atoms before removal: %d, after removal: %d" % (xyz.shape[0], new_xyz.shape[0]))
            self.xyz = new_xyz
            self.elements = new_elements
        return
        
    def to_Cartesian(self):
        """
        transforms the molecule into a chemcoord.Cartesian object
        """
        natoms = self.xyz.shape[0]
        elems = copy.deepcopy(self.elements)
        for i, j in enumerate(elems): 
            elems[i] = j.strip().capitalize()
        xyz = pandas.DataFrame(self.xyz, columns=["x","y","z"], dtype='float64')
        elems = pandas.DataFrame(elems, columns=["atom"], dtype='str')
        output = chemcoord.Cartesian(pandas.concat([elems, xyz], axis=1))
        output.index = range(1, natoms+1)
        return output
        
    def from_Cartesian(self, cartesian):
        """
        loads molecule data from a chemcoord.Cartesian object
        """
        self.xyz = cartesian[:, ['x', 'y', 'z']].as_matrix()
        elements = cartesian[:, 'atom'].as_matrix()
        for i in range(len(elements)):
            if len(elements[i]) == 1:
                elements[i] += ' '
        self.natoms = self.xyz.shape[0]
        self.elements = elements
        self.set_masses()
        self.set_vdwr()
        return
    
    def rotate_dihedral(self, idx, deg):
        """ 
        Rotates a dihedral angle
        Parameters:
          - idx : List of atom indices of the atoms spanning the dihedral
          - deg : target angle in degrees
        """
        if self.xyz.shape[0] < 4:
            raise IOError('The amount of atoms in the molecule is smaller than 4!')
        if len(idx) != 4:
            raise IOError('The amount of indices is not 4!')
        xyz = self.to_Cartesian()
        idx_array = [[idx[0],      0,      0,      0], \
                     [idx[1], idx[0],      0,      0], \
                     [idx[2], idx[1], idx[0],      0], \
                     [idx[3], idx[2], idx[1], idx[0]]]
        idx_array = numpy.array(idx_array)
        buildlist = xyz._get_buildlist(fixed_buildlist = idx_array)
        zmat = xyz.to_zmat(buildlist)
        zmat[idx[3], 'dihedral'] = deg
        xyz = zmat.to_xyz()
        self.from_Cartesian(xyz)
        return
    





