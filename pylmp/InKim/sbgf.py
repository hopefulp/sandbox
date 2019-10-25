#!/home/noische/program/python27/bin/python
"""
BGF.py
Original: Aug 12 2004 -(c)- caglar tanrikulu
Modified: Jan 01 2011 In Kim

Module containing BGF-file specific tools including:
(*) BgfAtom class:
        Stores atom information
(*) BgfFile class:
        Stores all components of a bgf file
        The structural information is stored as a sequence of 'BgfAtom's
"""

import os, sys, re, string, math, copy
import cutils as cu

version = '110808'

#
# - BGF file constants:
#


# default values for the optional fields:
defaultBgfFields = {'fixed':0, 'field1':0, 'fsmA':0.0, 'fsmT':0 }


#
# - General BGF file utilities:
#

def getAtomNumberFromLine(bgfline):
    """
getAtomNumberFromLine(bgfline):
    returns the atom number form a BGF ATOM/HETATM/CONECT/ORDER line
    """
    return int(bgfline[7:12])


# ----------------------------
#
# - Define the BgfAtom object:
#

class BgfAtom:
    """
    BgfAtom object stores atom information contained in the
    ATOM/HETATM, CONECT and ORDER lines
    """

    def __init__(self, *lines):
        """
    __init__(self, *lines):
        initialize atom with empty fields, or data from atom lines
        as provided

        attributes of the BgfAtom class are:
        self.aTag   = 0       # atom tag, 0->'ATOM  ', 1-> 'HETATM'
        self.aNo    = 0       # atom number
        self.aName  = ''      # atom name
        self.rName  = ''      # residue name
        self.chain  = ''      # chain name
        self.rNo    = 0       # residue number
        self.x      = 0.0     # x coordinate
        self.y      = 0.0     # y coordinate
        self.z      = 0.0     # z coordinate
        self.ffType = '     ' # forcefield type
        self.bonds  = 0       # number of connected bonds
        self.lpair  = 0       # number of lone pairs
        self.charge = 0.0     # atom charge
    
        the following fields have -ve default values in order to 
        differentiate between their initialized and unitintialized states
        self.fixed  = -1      # movable rec., 0-> free, 1-> fixed
        self.field1 = -1      # unknown field 1
    
        self.fsmA   = -1.0    # fsm atomic weight
        self.fsmT   = -1      # fsm solvation type
        
        self.CONECT =[]       # list of connected atom numbers
        self.ORDER  =[]       # bond order for connected atoms """
        
        # atom attributes
        # initialized to default values:
        self.aTag   = 0       # atom tag, 0->'ATOM  ', 1-> 'HETATM'
        self.aNo    = 0       # atom number
        self.aName  = ''      # atom name
        self.rName  = ''      # residue name
        self.chain  = ''      # chain name
        self.rNo    = 0       # residue number
        self.x      = 0.0     # x coordinate
        self.y      = 0.0     # y coordinate
        self.z      = 0.0     # z coordinate
        self.ffType = '     ' # forcefield type
        self.bonds  = 0       # number of connected bonds
        self.lpair  = 0       # number of lone pairs
        self.charge = 0.0     # atom charge

        # the following fields have -ve default values in order to 
        #  differentiate between their initialized and unitintialized states
        self.fixed  = -1      # movable rec., 0-> free, 1-> fixed
        self.field1 = -1      # unknown field 1

        self.fsmA   = -1.0    # fsm solvation per area (?)
        self.fsmT   = -1      # fsm solvation type (?)
        
        self.CONECT = []       # list of connected atom numbers
        self.ORDER  = []       # bond order for connected atoms

        # the following fields are added due to velocity
        self.vx     = 0.0
        self.vy     = 0.0
        self.vz     = 0.0
        self.fx     = 0.0
        self.fy     = 0.0
        self.fz     = 0.0

        # if BGF lines are provided, read them:
        if lines:
            for line in lines:  self.readBgfLine(line)


    def ATOMline(self):
        """
    ATOMline(self):
        returns a string: the BGF ATOM line """
        
        if self.aTag:
            output = 'HETATM'
        else:
            output = 'ATOM  '

        item = (self.aNo, self.aName, self.rName, self.chain, self.rNo, float(self.x), float(self.y), float(self.z), self.ffType, self.bonds, self.lpair, float(self.charge))
        output += ' {0:>5} {1:<5} {2:3} {3:<1} {4:>5}{5:>10.5f}{6:>10.5f}{7:>10.5f} {8:<5}{9:3}{10:2} {11:>8.5f}'.format(*item)
        if self.fixed >= 0:
            output += '%2d%4d' % (self.fixed, self.field1)

            if self.fsmA > 0:
                output += '  %8.4f %3d' % (self.fsmA, self.fsmT)

        output += '\n'
        
        # return
        return output


    def CONECTline(self):
        """
    CONECTline(self):
        returns a string: the BGF CONECT line """
        output = ''
        
        if self.CONECT:
            output  = 'CONECT %5d' + ('%6d' * len(self.CONECT)) + '\n'
            output %= tuple([self.aNo] + self.CONECT)

        # return
        return output


    def ORDERline(self):
        """
    ORDERline(self):
        returns a string: the BGF ORDER line """
        output = ''

        if self.ORDER:
            output  = 'ORDER  %5d' + ('%6d' * len(self.ORDER)) + '\n'
            output %= tuple([self.aNo] + self.ORDER)
        
        # return
        return output


    def __str__(self):
        """
    __str__(self):
        prints atom lines """
        return self.ATOMline() + self.CONECTline() + self.ORDERline()

        
    def readBgfLine(self,line):
        """
    readBgfLine(self,line):
        enters info from ATOM/HETATM, CONECT and ORDER lines into a BgfAtom object """

        line = line.strip()
	line = line.replace("-", " -")
        
        if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':

            if line[0:6] == 'HETATM':  self.aTag = 1
            l_line = line.split()
            self.aNo    =   int(l_line[1])
            self.aName  =     l_line[2]
            self.rName  =     l_line[3]
            self.chain  =     l_line[4]
            self.rNo    =   int(l_line[5])
            self.x      = float(l_line[6])
            self.y      = float(l_line[7])
            self.z      = float(l_line[8])
            self.ffType =       l_line[9]
            self.bonds  =   int(l_line[10])
            self.lpair  =   int(l_line[11])
            self.charge = float(l_line[12])

        elif line[0:6] == 'CONECT' or line[0:6] == 'ORDER ':
            if not self.aNo:
                cu.die('readBgfLine:  Cannot read CONECT/ORDER data prior to ATOM data!')

            elif int(line[7:12]) != self.aNo:
                cu.die('readBgfLine:  The line read does not refer to the atom, "%d":\n%s' % (self.aNo, line))

            else:       # using re.split is OK in the following lines since aNo cannot be 6 digits
                if line[0:6] == 'CONECT':
                    conectInfoString = line[12:].strip()
                    if conectInfoString:
                        self.CONECT = cu.list2intList( re.split(r'\s+', conectInfoString) )
                
                else:   #elif line[0:6] == 'ORDER ':
                    if not self.CONECT:
                        cu.die('readBgfLine:  Cannot read ORDER data prior to CONECT data!')                        
                    else:
                        self.ORDER  = cu.list2intList( re.split(r'\s+', line[12:].strip()) )

                        if len(self.CONECT) != len(self.ORDER):
                            cu.die('readBgfLine:  CONECT and ORDER line elements do not match:\n%s%s' % \
                                (self.CONECTline(), line) )
            
        else:
            cu.die('readBgfLine:  Line provided does not carry any atom information:\n%s' % (line,))


    """
    def connect(self,other,order=1):
        """
    connect(self, other, order=1):
        forms a bond between the atoms 'self' and 'other' with the specified 'order' """
        if other.aNo not in self.CONECT:
            self.CONECT.append(other.aNo)
            if self.ORDER:
                self.ORDER.append(order)
            elif order != 1:
                self.ORDER = [1]*len(self.CONECT)
                self.ORDER[-1] = order

        if self.aNo not in other.CONECT:
            other.CONECT.append(self.aNo)
            if other.ORDER:
                other.ORDER.append(order)
            elif order != 1:
                other.ORDER = [1]*len(other.CONECT)
                other.ORDER[-1] = order


    def disconnect(self,other,strict=False):
        """
    disconnect(self, other, strict=False):
        removes the bond between the atoms 'self' and 'other'
        if strict is True, dies when no bond is found between the atoms"""
        # remove other from self's list
        for n in range(len(self.CONECT)):
            if self.CONECT[n] == other.aNo:
                del self.CONECT[n]
                if self.ORDER:
                    del self.ORDER[n]
                break
        else:
            if strict:
                cu.die("atoms %d is not connected to %d"%(self.aNo,other.aNo))
        
        # remove self from other's list
        for n in range(len(other.CONECT)):
            if other.CONECT[n] == self.aNo:
                del other.CONECT[n]
                if other.ORDER:
                    del other.ORDER[n]
                break
        else:
            if strict:
                cu.die("atoms %d is not connected to %d"%(other.aNo,self.aNo))
    # done
    """

            
    # ... checking atom properties:

    """
    def is_hydrogen(self):
        """
    is_hydrogen(self):
        returns true if atom is a hydrogen atom """
        if self.ffType[0:2] == 'H_':
            return True
        else:
            return False

    def is_carbon(self):
        """
    is_carbon(self):
        returns true if atom is a carbon atom """
        if self.ffType[0:2] == 'C_':
            return True
        else:
            return False

    def is_free(self):
        """
    is_free(self):
        returns true if atom is free"""
        if not self.fixed:   return True
        else:                return False
    
    def is_fixed(self):
        """
    is_fixed(self):
        returns true if atom is fixed"""
        if self.fixed:   return True
        else:            return False
        
    def move(self, dx=0, dy=0, dz=0):
        """
    move(atom, dx=0):
        move the atom position by adding dx, dy, and dz to the original x, y, and z coordinates respectively.
        """
        self.x += float(dx)
        self.y += float(dy)
        self.z += float(dz)
    """

    # end of BgfAtom class #

        
#
# - Functions/Operations using the atom object:
#

# ... changing atom attributes:

"""
def setFree(atom):
    """
setFree(atom):
    sets the movable rec. of the atom to 0
    """
    atom.fix()
    

def setFixed(atom):
    """
setFixed(atom):
    sets the movable rec. of the atom to 1
    """
    atom.free()


# ... checking atom properties:


def is_free(atom):
    """
is_free(atom):
    returns true if atom is free
    """
    return atom.is_free()
    

def is_fixed(atom):
    """
is_fixed(atom):
    returns true if atom is fixed
    """
    return atom.is_fixed()
"""


# ... changing atom formatting/BGF version:
    
"""
def guessBGFformat(atom):
    """
guessBGFformat(atom):
    returns the first integer of the BIOGRF version, based on the line length
    """
    if atom.fixed < 0:
        return 2
    elif atom.fsmA < 0:
        return 3
    else:
        return 4
"""


# ... evaluating properties of two atoms:

"""
def is_bonded(atom1, atom2):
    """
is_bonded(atom1, atom2):
    returns true if atoms are bonded
    """
    return is_12(atom1,atom2)


def is_12(atom1, atom2):
    """
is_12(atom1, atom2):
    returns true if atoms are bonded
    """
    for i in atom1.CONECT:
        if i == atom2.aNo:
            return 1
    else:
        return 0
"""


def sqrDistance(atom1, atom2):
    """
sqrDistance(atom1, atom2):
    returns the distance squared (r^2) between two atoms
    """
    #return ( (atom2.x - atom1.x) ** 2 +\
    #         (atom2.y - atom1.y) ** 2 +\
    #         (atom2.z - atom1.z) ** 2  )
    # define vector a:
    a = ( (atom2.x-atom1.x), (atom2.y-atom1.y), (atom2.z-atom1.z) ) #1->2
    return a_dot_b(a,a)


def distance(atom1, atom2):
    """
distance(atom1, atom2):
    returns the distance r (in A) between two atoms
    """
    return math.sqrt( sqrDistance(atom1, atom2) )


def angle(atom1, atom2, atom3, radians=False):
    """
angle(atom1, atom2, atom3, radians=False):
    returns the angle theta (in degrees unless radians=True) between three 
    atoms around central atom, atom2
    """
    # define vectors a and b
    a = ( (atom1.x-atom2.x), (atom1.y-atom2.y), (atom1.z-atom2.z) ) #2->1
    b = ( (atom3.x-atom2.x), (atom3.y-atom2.y), (atom3.z-atom2.z) ) #2->3
    
    #                   a . b
    #  theta = arccos ( ------ )
    #                   |a||b|
    theta = math.acos( a_dot_b(a,b) / (math.sqrt(a_dot_b(a,a))*math.sqrt(a_dot_b(b,b))) )
    if radians:
        return theta
    else:
        return rad2deg(theta)     


def dihedral(atom1, atom2, atom3, atom4, radians=False):
    """
dihedral(atom1, atom2, atom3, atom4, radians=False):
    returns the torsion angle theta (in degrees unless radians=True) around 
    the bond between atom2 and atom3, assuming that the atoms are bonded 
    as 1-2-3-4
    """
    # define vectors a, b, and c between the atoms
    a = ( (atom2.x-atom1.x), (atom2.y-atom1.y), (atom2.z-atom1.z) ) #1->2
    b = ( (atom3.x-atom2.x), (atom3.y-atom2.y), (atom3.z-atom2.z) ) #2->3
    c = ( (atom4.x-atom3.x), (atom4.y-atom3.y), (atom4.z-atom3.z) ) #3->4

    # normal vectors to the planes defined by vectors a&b and b&c:
    bc = a_cross_b(a,b)
    ab = a_cross_b(b,c)
    
    # calculate magnitude of the angle between these two normal vectors:
    #                   ab . bc
    #  theta = arccos ( -------- )
    #                   |ab||bc|
    theta = math.acos( a_dot_b(ab,bc) / (math.sqrt(a_dot_b(ab,ab))*math.sqrt(a_dot_b(bc,bc))) )
    
    # to get the sign, we look at the direction of the vector abXbc =(ab x bc) and compare
    # it with vector b.  if the two vectors are in the same direction, the sign is negative
    abXbc = a_cross_b(ab,bc)
    if a_dot_b(b, abXbc) > 0:
        theta = -theta
    if radians:
        return theta     
    else:
        return rad2deg(theta)     


#
# mathematical functions used in the measurements between atoms
#

def a_dot_b(a,b): 
    """ 
a_dot_b(a,b): 
    returns the dot product of vector a=(ax,ay,az) with vector b=(bx,by,bz) 
    """  
    return float( (a[0]*b[0]) + (a[1]*b[1]) + (a[2]*b[2]) ) 


def a_cross_b(a,b): 
    """ 
a_cross_b(a,b): 
    calculates c, the cross product of vector a=(ax,ay,az) with b=(bx,by,bz) 
    returns vector c = a x b = (cz,cy,cz) as a tuple
    """  
    return ( (a[1]*b[2])-(a[2]*b[1]), \
             (a[2]*b[0])-(a[0]*b[2]), \
             (a[0]*b[1])-(a[1]*b[0]) )


def rad2deg(radians):
    """
rad2deg(radians):
    converts radians to degrees
    """
    return 180*(radians/math.pi)



# ---------------------------
#
# - Define the BgfFile class:
#

class BgfFile:
    """
    Class to contain collections of BgfAtoms.
    Each object contains a list of atoms, which are indexed consecutively
    starting from zero, as well as indexing information that relates chains
    and residues to these atoms.
    """

    def __init__(self, file=0, safe=0):
        """
    __init__(self, file=0, safe=0):
        creates a new, empty bgf-file object and, if a file is specified,
        populates this object with the contents of the file.
        if 'safe' is set, information after the charges is not read. 
        self.BIOGRF = ''    # Biograf version 
        self.DESCRP = 'bgf' # file descriptor
        self.REMARK = []    # REMARK list
        self.FF     = 'FORCEFIELD DREIDING'
        self.FORMAT = {'ATOM'  :'FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,i2,i4,f10.5)',
                       'CONECT':'FORMAT CONECT (a6,14i6)',
                       'ORDER' :'FORMAT ORDER (a6,i6,13f6.3)'}
        self.PERIOD = ''
        self.CRYSTX = []
        self.OTHER  = []    # a list of unresolved/unknown lines

        # file info:
        self.source  = ''   # file source
        self.nlines  = 0    # number of lines in file
        self.format  = 0    # bgf format: 2/3/4

        # atoms:
        self.a = []         # index of first atom in file is 0
        self.natoms  = 0    # number of atoms

        # bookkeeping:
        self.a2i = {}       # translates actual atomNo (a) to internal atomNo (i) 
                            # the following information is obtained from self.a2i: 
        self.max_aNo = 0    #  largest atom number in the file
                            # smallest atom number in the file
        self.min_aNo = sys.maxint

        # compilation populates the data structures below:
          this allows the information to be reached in the order

              chains --> resIDs ---,-> internal atom numbers --,-> atom numbers
               (c)        (r)     /       (i)                 /      (a)
              resIDs ------------'                           /
              internal atom numbers ------------------------'

          however it is not possible to reach atom numbers directly.
    
        self.compiled=False # tells if the res/chain range info is compiled or not
        self.r2i = {}       # translates resID to atom ranges
        self.c2i = {}       # translates chain to atom ranges

                            # for .chains and .residues below: 
                            # (I'm assuming that chains may be broken in the BGF file,
                            # but that residues won't be broken) 

        self.chains = []    # keeps a list of chains
        self.residues = {}  # unlike .chains, .residues is a dictionary, where each 
                            #  chain name points to an array of resID's in that chain
        self.resIDs  = []   # keeps a list of residues
        self.nres    = 0    # number of residues
        self.nchains = 0    # number of chains """
        
        # bgf lines with default values:
        self.BIOGRF = ''    # Biograf version 
        self.DESCRP = 'bgf' # file descriptor
        self.REMARK = []    # REMARK list
        self.FF     = ''
        self.FORMAT = {'ATOM'  :'FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,i2,i4,f10.5)', \
                       'CONECT':'FORMAT CONECT (a6,14i6)', \
                       'ORDER' :'FORMAT ORDER (a6,i6,13f6.3)'}
        self.PERIOD = ''
        self.AXES = ''
        self.SGNAME = []
        self.CRYSTX = []
        self.CELL = []
        self.OTHER  = []    # a list of unresolved/unknown lines

        # file info:
        self.source  = ''   # file source
        self.nlines  = 0    # number of lines in file
        self.format  = 0    # bgf format: 2/3/4
        
        # atoms:
        self.a = []         # index of first atom in file is 0
        self.natoms  = 0    # number of atoms

        # bookkeeping:
        self.a2i = {}       # translates actual atomNo (a) to internal atomNo (i) 
                            # the following information is obtained from self.a2i: 
        self.max_aNo = 0    #  largest atom number in the file
                            # smallest atom number in the file
        self.min_aNo = sys.maxint
        
        # compilation populates the data structures below:
        self.compiled=False # tells whether the res/chain range information is compiled or not
        self.r2i = {}       # translates resID to atom ranges
        self.c2i = {}       # translates chain to atom ranges

                            # for .chains and .residues below: 
                            # ( I'm assuming that chains may be broken in the BGF file,
                            # but that residues won't be broken) 

        self.chains = []    # keeps a list of chains
        self.residues = {}  # unlike .chains, .residues is a dictionary, where each chain
                            #  name points to an array of resID's in that chain
        self.resIDs  = []   # keeps a list of residues
        self.nres    = 0    # number of residues
        self.nchains = 0    # number of chains
        self.cmatrix = []   # connection matrix

        # if a filename is specified, read the file in:
        if file:  self.readBGF(file,safe)


    def __len__(self):
        """
    __len__(self):
        returns the number of atoms in the BgfFile object"""
        return self.natoms
    

    def __str__(self):
        """
    __str__(self):
        prints self as a .bgf file """
        lines = self.writeBgfToLines()
        return string.join(lines,'')


    """
    def copy(self):
        """
    copy(self):
        returns a (deep)copy of the BgfFile object """
        return copy.deepcopy(self)
        

    def compile(self, forceCompile=True):
        """
    compile(self, forceCompile=True):
        compiles the chain and residue information on a BgfFile object
          (warns if recompiling without forceCompile set to true)
        by doing this, we have all the information on:
        - which residues are in the file       : self.r2i.keys()
        -   and where they start and end       : self.r2i
        - which chains are in the file         : self.chains
        -   and where they start and end       : self.c2i
        -   and which residues they contain    : self.resIDs
        -   and residues by chain              : self.residues
        - how many residues and chains we have : self.nres and self.nchains"""
    """

        # done
        

    def readBgfFromLines(self,lines,safe=0):
        """
    readBgfFromLines(self,lines,safe=0):
        reads a list of BGF lines and populates basic fields
        if 'safe' is set, data after the charges are not read """

        aCount = 0          # counts atoms
                
        for line in lines:
            tag = re.split(r'\s+', line)[0]
            parse = line.split()

            if tag == 'ATOM' or tag == 'HETATM':

                # REMARK: maybe safe mode won't work since line is parsed
                if safe:
                    line = line[:80] # do not read anything after 80 chars
                
                atom = BgfAtom()
                atom.readBgfLine(line.strip())
                self.format = max(self.format, guessBGFformat(atom))
                self.a.append(atom)
                
                self.a2i[atom.aNo] = aCount
                self.max_aNo = max(atom.aNo, self.max_aNo)
                self.min_aNo = min(atom.aNo, self.min_aNo)

                aCount += 1

            elif tag == 'PERIOD':
                if len(parse) == 2:
                    self.PERIOD = parse[1]
                else:
                    self.PERIOD = ""

            elif tag == 'CRYSTX':
                if len(parse) > 2:
                    self.CRYSTX = [ float(parse[1]), float(parse[2]), float(parse[3]), float(parse[4]), float(parse[5]), float(parse[6]) ]
                else:
                    self.CRYSTX = []

            elif tag == 'FORCEFIELD':
                if len(parse) == 2:
                    self.FF = str.strip(parse[1])
                else:
                    self.FF = ''

            elif tag == 'AXES':
                if len(parse) == 2:
                    self.AXES = str.strip(parse[1])

            elif tag == 'SGNAME':
                if len(line) > 8:
                    self.SGNAME = line[7:]

            elif tag == 'CELLS':
                if len(parse) > 2:
                    self.CELLS = parse[1:]

            elif tag == 'CONECT' or tag == 'ORDER':
                aNo = getAtomNumberFromLine(line)
                self.a[ self.a2i[aNo] ].readBgfLine(line)

            elif tag == 'REMARK':
                self.REMARK.append( line.strip()[7:] )

            elif tag == 'FORMAT':
                if re.search(r'^FORMAT ATOM', line):
                    self.FORMAT['ATOM'] = line.strip()                    
                elif re.search(r'^FORMAT CONECT', line):
                    self.FORMAT['CONECT'] = line.strip()                    
                elif re.search(r'^FORMAT ORDER', line):
                    self.FORMAT['ORDER'] = line.strip()                    
                else:
                    cu.die("readBgfFromLines: Unkown FORMAT line:\n%s" % line)

            elif tag == 'DESCRP':
                self.DESCRP = line[7:].strip()

            elif tag == 'BIOGRF':
                self.BIOGRF = line[7:].strip()

            # appended to prevent the errors from XTLGRP keyword. @ 110808
            elif tag == 'XTLGRF':
                self.BIOGRF = line[7:].strip()

            elif tag == 'FORCEFIELD':
                self.FF = line.strip()

            elif tag == 'REM' or \
                 tag == 'PHI' or \
                 tag == 'PSI':       # these are fields in the rotamer library
                self.OTHER.append(line)

            elif tag == 'END':
                break

            else:
                cu.warn("readBgfFromLines: Encountered unknown field:\n%s" % line)
                line = line.strip()
                if line:
                    self.OTHER.append(line)
        else:
            if self.source[-4:] != '.lib':
                cu.warn("readBgfFromLines:  Did not see the 'END' tag at the end of file, %s!" % self.source)

        self.natoms  = aCount


    def readBGF(self,file,safe=0):
        """
    readBGF(self,file,safe=0):
        reads a .bgf file into self
        if 'safe' is set, data after the charges are not read """

        # read bgf file:
        lines = cu.readLinesFromFile(file)
        
        # modify bgf info fields:
        self.source = file
        self.DESCRP = file[0:8]    
        self.nlines = len(lines)
        
        # read structure information:
        self.readBgfFromLines(lines,safe)


    def writeBgfToLines(self):
        """
    writeBgfToLines(self):
        writes self into a list of lines and returns it """
        lines = []

        lines.append("BIOGRF %s\n" % self.BIOGRF)    # Biograf version
        lines.append("DESCRP %s\n"  % self.DESCRP)    # File descriptor

        for entry in self.REMARK:                     # Remarks
            lines.append("REMARK %s\n" % entry)

        if self.FF != "":
            lines.append("FORCEFIELD %s\n" % self.FF)

        if self.PERIOD != "":
            lines.append("PERIOD %s\n" % self.PERIOD)
            lines.append("AXES   %s\n" % self.AXES)
            lines.append("SGNAME %s" % self.SGNAME)
            lines.append("CRYSTX {0:<11.5f}{1:<11.5f}{2:<11.5f}{3:<11.5f}{4:<11.5f}{5:<11.5f}\n".format(*self.CRYSTX))
            lines.append("CELLS  {0:<5}{1:<5}{2:<5}{3:<5}{4:<5}{5:<5}\n".format(*self.CELLS))
        """
            elif tag == 'PERIOD':
                self.PERIOD = parse[1]

            elif tag == 'CRYSTX':
                self.CRYSTX = [ float(parse[1]), float(parse[2]), float(parse[3]), float(parse[4]), float(parse[5]), float(parse[6]) ]

            elif tag == 'AXES' or tag == 'CELLS' or tag == 'SGNAME':
                line = line.strip()
                if line:
                    self.OTHER.append(line)
        """

        for entry in self.OTHER:                      # Unknown entries that are not 
            lines.append("%s\n" % entry)              #  empty lines

        lines.append("%s\n" % self.FORMAT['ATOM'])    # FORMAT ATOM line

        for atom in self.a:                           # ATOM lines
            lines.append(atom.ATOMline())

        lines.append("%s\n" % self.FORMAT['CONECT'])  # FORMAT CONECT line
        lines.append("%s\n" % self.FORMAT['ORDER'])   # FORMAT ORDER line        

        for atom in self.a:                           # CONECT and ORDER lines
            lines.append(atom.CONECTline())
            lines.append(atom.ORDERline())

        lines.append('END\n')                         # END

        # return
        return lines
    

    def saveBGF(self,file):
        """
    saveBGF(self,file):
        writes self into a file """
        # write bgf file into text:
        lines = self.writeBgfToLines()    
        
        # save structure information:
        cu.writeLinesToFile(lines,file)



    def addAtom(self,atom,index=None):
        """
    addAtom(self,atom,index=None):
        adds a new atom to self; renumbers atom if necessary
        if index is specified, the new atom is inserted into the atom list at the
        specified index, otherwise atom is appended at the end of the list """
        atom = atom.copy()   # create a copy (not reference)

        if self.a2i.has_key(atom.aNo):   # aNo exists, need to renumber this atom:
            atom.aNo = self.max_aNo + 1

        if index != None:
            self.a.insert(index,atom)    # insert new atom
            self.a2i[atom.aNo] = index
        else:
            self.a.append(atom)          # add new atom
            self.a2i[atom.aNo] = self.natoms

        self.max_aNo = max(atom.aNo, self.max_aNo)  #
        self.min_aNo = min(atom.aNo, self.min_aNo)  # ... fill in the rest ...
                                                    #
        self.natoms += 1                            #

        if self.compiled:
            # need to update the compiled records
            # at this point!!!!!!!!!!!!!!!!!!!!!!!!!!!
            self.compiled = False
        # done

        
    def delAtom(self,index):
        """
    delAtom(self,index):
        deletes one atom from self, and updates all x2i information """
        if index < 0 or index >= self.natoms:  # if atom doesn't exist: 
            cu.die('BgfFile.delAtom: Atom at index "%d" does not exist' % index)

        del self.a2i[self.a[index].aNo]     # delete the a2i record
        delaNo = self.a[index].aNo
        del self.a[index]                   # delete the atom

        self.natoms -= 1                    # update no.of atoms

        for i in range(index, self.natoms): # decrement/update all the a2i 
            atom = self.a[i]                #  records that come after this
            self.a2i[atom.aNo] = i           #  index number

        """
        for i in range(0, len(self.a)):
            atom = self.a[i]
            ano = index + 1
            if delaNo in atom.CONECT:
                atom.CONECT.remove(delaNo)
        """
        if self.compiled:
            # need to update the compiled records
            # at this point!!!!!!!!!!!!!!!!!!!!!!!!!!!
            self.compiled = False

        self.max_aNo = max( self.a2i.keys() ) # ... update constants ...
        self.min_aNo = min( self.a2i.keys() ) #
        # done


    def delRange(self, iRange):
        """
    delRange(self, iRange):
        deletes atoms indicated by a range of indices from self """
        if iRange[1] > self.natoms:
            cu.die('BgfFile.delAtom:  Range (%d,%d) extends outside the current atom list!' % iRange)

        for index in range(iRange):
            del self.a2i[self.a[index].aNo]     # delete the a2i record
            del self.a[index]                   # delete the atom            
            self.natoms -= 1                    # update no.of atoms

        for i in range(iRange[0], self.natoms): # decrement/update all the a2i 
            atom = self.a[i]                    #  records that come after this
            self.a2i[atom.No] = i               #  index number

        if self.compiled:
            # need to update the compiled records
            # at this point!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            self.compiled = False

        self.max_aNo = max( self.a2i.keys() ) # ... update constants ...
        self.min_aNo = min( self.a2i.keys() ) #
        # done


    def addAtoms(self,atoms,index=None):
        """
    addAtoms(self,atoms,index=None):
        adds a list of atoms to self; renumbers new atoms if necessary
        if index is specified, the new atoms are inserted into the atom list at the
        specified index, otherwise atoms are appended at the end of the list """
        if index != None:
            atoms.reverse()
        for atom in atoms:  self.addAtom(atom,index)


    def delAtoms(self,indices,silent=True):
        """
    delAtoms(self,indices):
        deletes atoms indicated by a list of indices from self """
        indices.sort()
        indices.reverse()
        n_atom = len(indices)
        n_progress = 0;
        for index in indices:
            self.delAtom(index)
            n_progress += 1
            if n_atom > 1000:
                if not silent: sys.stdout.write("\rDeleting " + str(n_progress) + " atoms out of " + str(n_atom))
        #print("")

        
    def merge(self,other,dereference=False):
        """
    merge(self,other,dereference=0):
        merges two bgf files by adding <other> at the end of <self>.  Creates a new
        BgfFile object to house the merged data
        if <dereference> is set to True, dereferences all information before moving
        them to the newly created object
        returns the new BgfFile object """
        if dereference:
            result = self.copy()
            other = other.copy()
        else:
            if not self.natoms:                  # if self is empty
                result = self                    #   -> no need to protect it
                result.BIOGRF = other.BIOGRF     #   -> pick up other's information
                result.DESCRP = other.DESCRP
                result.min_aNo= 0
            else:                                # if it is a real file
                result = BgfFile().merge(self)   #   -> we need a new object for the result
                
        
        # figure out whether to renumber other's atom numbers before merging:
        if result.max_aNo < other.min_aNo:
            increment = 0
        else:
            increment = result.max_aNo - other.min_aNo + 1

        # renumber atoms of other, adding them the increment
        #  and add these atoms to sum:
        for atom in other.a:
            if increment:
                atom.aNo += increment
                for i in range(len(atom.CONECT)):   atom.CONECT[i] += increment

            result.a.append(atom)                # not using the addAtom() function here
            result.a2i[atom.aNo] = result.natoms #  so that we don't check for renumbering
            result.natoms += 1                   #  each time we read a line from the file

        # adjust the remaining variables:
        result.max_aNo = max( result.a2i.keys() )
        result.min_aNo = min( result.a2i.keys() )
        
        # return:
        return result

            
    def __add__(self,other):
        """
    __add__(self,other):
        merges two bgf files into one BgfFile object
        does dereference information when moving it to the newly created object
        <other> can be a BgfFile or a BgfAtom object
        returns the new BgfFile object """
        if isinstance(other,BgfFile):      # if second object is a BgfFile
            return self.merge(other,True)  #   merge it with the first
        elif isinstance(other,BgfAtom):    # if second object is a BgfAtom
            result = self.copy()           #   create new object,
            result.addAtom(other.copy())   #   add the atom to it,
            return result                  #   and return result


    def __radd__(self,other):
        """
    __radd__(self,other):
        reverse __add__
        works only if the second object is a BgfFile """
        if isinstance(other,BgfFile):    # execute if second object is a BgfFile
            return BgfFile.__add__(other,self)
        elif isinstance(other,BgfAtom):  # die if second object is a BgfAtom
            cu.die("Cannot add a BgfFile to a BgfAtom object")
            

    def renumber(self,start=1):
        """
    renumber(self,start=1):
        renumbers all atom numbers (.aNo) in the BgfFile starting from 'start' """
        old2new_aNo = {}
        a2i = {}
        
        # renumber ATOM lines:
        for i in range(len(self.a)):
            cur_aNo = i + start

            old2new_aNo[ self.a[i].aNo ] = cur_aNo
            self.a[i].aNo = cur_aNo

            a2i[cur_aNo] = i

        # save new a2i array
        self.a2i = a2i

        # renumber CONECT lines:
        for i in range(self.natoms):
            newConectList = []
            myConectList = self.a[i].CONECT

            for j in range(len( myConectList )):
                connectedAtom = myConectList[j]
                if old2new_aNo.has_key(connectedAtom):
                    newConectList.append(old2new_aNo[connectedAtom])
            self.a[i].CONECT = newConectList
 
            if (len(newConectList) != len(myConectList)) and self.a[i].ORDER:
                newOrderList = []
                myOrderList  = self.a[i].ORDER
                for j in range(len( myConectList )):
                    connectedAtom = myConectList[j]
                    if old2new_aNo.has_key(connectedAtom):
                        newOrderList.append(myOrderList[j])
                self.a[i].ORDER  = newOrderList
        # done

    """
    def connect(self, index1, index2, order=1):
        """
    connect(self, index1, index2, order=1):
        forms a bond of the specified order, between the atoms identified by
        index1 and index2 """
        self.a[index1].connect( self.a[index2], order)

    def connectAtoms(self, aNo1, aNo2, order=1):
        """
    connectAtoms(self, aNo1, aNo2, order=1):
        forms a bond of the specified order, between the atoms identified by
        atom number 1 and atom number 2 """
        i1 = self.getAtomIndex(aNo1)
        i2 = self.getAtomIndex(aNo2)
        self.a[i1].connect( self.a[i2], order)

    def disconnect(self, index1, index2, strict=False):
        """
    disconnect(self, index1, index2, strict=False):
        removes the bond/connectivity between the atoms identified by index1 and index2
        if strict is True, dies with an error when the bond is not found in self"""
        self.a[index1].disconnect( self.a[index2], strict )
    """


    def removeDanglingBonds(self):
        """
    removeDanglingBonds(self):
        removes any CONECT/ORDER entries records that refer to atoms that are no longer
        described in the current BgfFile object """
        
        for i in range(self.natoms):
            iatom = self.a[i]

            itemsToBeDeleted = []
            itemNo = -1
            for ja in self.a[i].CONECT:
                itemNo += 1
                if not self.a2i.has_key(ja):
                    itemsToBeDeleted.append(itemNo)

            if itemsToBeDeleted:
                itemsToBeDeleted.reverse()
                for itemIndex in itemsToBeDeleted:
                    del self.a[i].CONECT[itemIndex]
                    if self.a[i].ORDER:
                        del self.a[i].ORDER[itemIndex]
        # done
        
    """
    def setFormat(self, format=3):
        """
    setFormat(self, format=3):
        sets the BgfFile atom line format to 'format' (2, 3, or 4) for each atom"""
        for i in range(self.natoms):
            self.a[i].setFormat(format)

    def setVersion(self,version):
        """
    setVersion(self,version):
        sets the BIOGRF line to 'version' and formats all atom lines accordingly"""
        version = str(version).strip()
        self.BIOGRF = version
        self.FORMAT = BgfFile().FORMAT        
 
        if version[0] == '4':
            self.setFormat(4)
        elif version[0] == '3':
            self.setFormat(3)
        elif version[0] == '2':
            self.setFormat(2)
        else:
            pass
        # done

    
    def charge(self):
        """
    charge(self):
        returns the sum of the charges of all atoms in the BgfFile object"""
        totalCharge = 0.0
        for i in range(self.natoms):
            totalCharge += self.a[i].charge
        
        return totalCharge


    def getAtomByIndex(self,intAtomNo=0):
        """
    getAtomByIndex(self,intAtomNo=0):
        returns the atom object indicated by the internal atom no (intAtomNo)"""
        return self.a[intAtomNo]


    def getAtomByNumber(self, atomNo):
        """
    getAtomByIndex(self, atomNo):
        returns the atom object indicated by the internal atom no (intAtomNo)"""
        return self.a[self.a2i[atomNo]]


    def getAtomIndex(self, atomNo):
        """
    getAtomIndex(self, atomNo):
        returns the BgfAtom specified by the atom number"""
        index = self.a2i.get(atomNo, (-1))

        # return:
        if index < 0:
            cu.die("getAtom: No such atom number, %s, in BgfFile object!" % atomNo)
        else:
            return index
    """
     
    # functions to return atoms, residues, chains without the need to deal with atom indices

    def getAtom(self, atomNo):
        """
    getAtom(self, atomNo):
        returns the BgfAtom specified by the atom number (internal index)"""
        index = self.a2i.get(atomNo, (-1))

        # return:
        if index < 0:
            cu.die("getAtom: No such atom number, %s, in BgfFile object!" % atomNo)
        else:
            return self.a[index]
        

    """
    def getResidue(self, resID):
        """
    getResidue(self, resID):
        returns the residue specified by the residue ID as a BgfFile object"""
        if not self.compiled:
            self.compile()

        start,stop = self.r2i.get(resID, (-1,-1))

        if start < 0:
            cu.die("getResidue: No such resID, %s, in BgfFile object!" % resID)

        myResidue = BgfFile()
        for index in range(start,stop):
            myResidue.addAtom( self.a[index] )

        myResidue.BIOGRF = self.BIOGRF
        
        # return:
        return myResidue
    
    
    def getChain(self, chainID):
        """
    getChain(self, chainID):
        returns the chain specified by the residue ID as a BgfFile object"""
        if not self.compiled:
            self.compile()

        atomRangeList = self.c2i.get(chainID,[])
        if not atomRangeList:
            cu.die("getChain: No such chainID, %s, in BgfFile object!" % chainID)
                
        myChain = BgfFile()
        for (start,stop) in atomRangeList:
            for index in range(start,stop):
                myChain.addAtom( self.a[index] )

        myChain.BIOGRF = self.BIOGRF

        # return:
        return myChain
    

    def getSidechain(self, resID=''):
        """
    getSidechain(self, resID=''):
        returns the sidechain atoms in the BgfFile as a BgfFile object.
        if resID is defined, only the sidechain atoms in that residue are returned."""
        if not self.compiled:
            self.compile()

        if resID:
            mybgf = self.getResidue(resID)
        else:
            mybgf = self

        sidechainAtoms = BgfFile()
        for index in range(mybgf.natoms):
            if mybgf.a[index].is_sidechain():
                sidechainAtoms.addAtom( mybgf.a[index] )
                
        sidechainAtoms.BIOGRF = self.BIOGRF
        
        # return:
        return sidechainAtoms

    
    def getIndicesByAtomName(self,atomName,indexList=[],expectedHits=None):
        """
    getIndicesByAtomName(self, atomName, indexList=[], expectedHits=None):
        searches for atoms with the given atomName among atoms specified by the list
        of indices, indexList, and returns the matching indices in an array.
        if no indexList is specified the search is carried out among all atoms present.
        if expectedHits is set to an integer, the function will die with an error if
        that many atoms are not identified."""
        atomName = atomName.strip()
        if not indexList:
            indexList = range(self.natoms)

        matchingIndices = []
        for i in indexList:
            if self.a[i].aName.strip() == atomName:
                matchingIndices.append(i)
        
        if expectedHits != None and expectedHits != len(matchingIndices):
            cu.die("getIndicesByAtomName: Only %d %s was expected, but found %d:  %s"%(expectedHits,atomName,len(matchingIndices),cu.flatWrite([ self.a[i].aNo for i in matchingIndices ])))
            
        # return
        return matchingIndices


    def getIndicesByAtomNameWithinRes(self,atomName,resID,expectedHits=None):
        """
    getIndicesByAtomNameWithinRes(self, atomName, resID, expectedHits=None):
        searches for atoms with the given atomName among atoms of the specified residue
        and returns the matching indices in an array.
        if expectedHits is set to an integer, the function will die with an error if
        more than that many atoms are found in the residue."""
        if not self.compiled:
            self.compile()
        start,stop = self.r2i[resID]

        return self.getIndicesByAtomName(atomName,range(start,stop),expectedHits)


    def fix(self):
        """
    fix(self):
        sets all atoms in the BgfFile object fixed """
        for i in range(self.natoms):
            self.a[i].fix()


    def free(self):
        """
    free(self):
        sets all atoms in the BgfFile object movable """
        for i in range(self.natoms):
            self.a[i].free()

    """

    # ... measurements:

    def distance(self, index1, index2):
        """
    distance(self,index1,index2):
        returns the distance r (in A) between the two atoms specified """ 
        return distance(self.a[index1],self.a[index2])

    def angle(self, index1, index2, index3, radians=False):
        """
    angle(self, index1, index2, index3, radians=False):
        returns the angle theta (in degrees unless radians=True) between three 
        atoms around the central atom with the index, index2 """
        return angle(self.a[index1],self.a[index2],self.a[index3],radians)

    def dihedral(self, index1, index2, index3, index4, radians=False):
        """
    dihedral(self, index1, index2, index3, index4, radians=False):
        returns the torsion angle theta (in degrees unless radians=1) around 
        the bond between the atoms specified by index2 and index3, assuming
        that the atoms are bonded as 1-2-3-4 """
        return dihedral(self.a[index1],self.a[index2],self.a[index3],self.a[index4],radians)


    ### end of BgfFile class ###


#
# - functions using BgfFiles:
#

def saveBGF(bgffile, filename):
    """
saveBGF(bgffile, filename):
    same as BgfFile.saveBGF
    """
    bgffile.saveBGF(filename)


def readBGF(filename):
    """
readBGF(filename):
    returns a BgfFile object read from 'filename'
    """
    return BgfFile(filename)


def merge(myBGF1,myBGF2,dereference=False):
    """
merge(myBGF1, myBGF2, dereference=False):
    merges the two BgfFile objects entered.
    if dereference=True is set dereferences both BgfFile objects
    prior to merging them.
    """
    return myBGF1.merge(myBGF2,dereference)


def moveBGF(myBGF, dx, dy, dz):
    """
moveBGF(bgffile, dx, dy, dz):
    move the BgfFile objects.
    """

    for atom in myBGF.a:
        atom.x += dx
        atom.y += dy
        atom.z += dz


"""
def guessFileFormat(myBGF):
    """
guessFileFormat(myBGF):
    based on the values read its atoms, evaulates the format
    of the BgfFile myBGF:  returns 2, 3, or 4
    *** BgfFile(filename).format makes this obsolete ***    
    """
    for i in range(myBGF.natoms):
        if myBGF.a[i].field1 >= 0:
            break
    else:
        return 2

    for i in range(myBGF.natoms):
        if myBGF.a[i].fsmT >= 0:
            break
    else:
        return 3
    return 4
    

def totalCharge(myBGF):
    """
totalCharge(myBGF):
    returns the sum of the charges of all atoms in the BgfFile object
    """
    return myBGF.charge()


def rmsd(myBGF1,myBGF2,safe=1,noH=0,scOnly=0):
    """
rmsd(myBGF1,myBGF2,safe=1,noH=0,scOnly=0):
    returns the RMSD (in A) between the two BgfFile objects, myBGF1 and myBGF2,
    that share the same structure.
    if 'safe' is set to 0, atom-by-atom similarity comparison will not be carried
    out between myBGF1 and myBGF2.
    if 'noH' is set to 1, hydrogens will not be included in the RMSD calculation
    if 'scOnly' is set to 1, sidechain-only RMSD will be returned
    """
    natoms = myBGF1.natoms
    if natoms != myBGF2.natoms:  # check natoms
        cu.die("BGF: rmsd: BgfFile objects don't contain the same number of atoms: %d =/= %d"%(natoms,myBGF2.natoms))

    atomCount = 0
    sumOfSqrDist = 0.0
    for i in range(natoms):
        atom1 = myBGF1.a[i]
        atom2 = myBGF2.a[i]
        if safe and (atom1 != atom2):  # check atoms
            cu.die("BGF: rmsd: Atom number %d refers to different atoms in the BgfFile objects provided."%i)
        if noH and atom1.is_hydrogen():
            continue
        if scOnly and atom1.is_sidechain():
            continue

        atomCount += 1
        sumOfSqrDist += sqrDistance(atom1, atom2)

    myRMSD = math.sqrt(sumOfSqrDist / atomCount)

    return myRMSD


def getBGFSize(myBGF, margin=0):
    """
getBGFSize(myBGF, margin=0):
    returns the list [xlo, xhi, ylo, yhi, zlo, zhi] of the BGF,
    which is a size of the molecule.
    If the margin (unit in Angstrom) is set, then the margin is added to each value.
    """
    bgfsize = [myBGF.a[0].x, myBGF.a[0].x, myBGF.a[0].y, myBGF.a[0].y, myBGF.a[0].z, myBGF.a[0].z];

    for atom in myBGF.a:
        if atom.x < bgfsize[0]:
            bgfsize[0] = atom.x
        elif atom.x > bgfsize[1]:
            bgfsize[1] = atom.x

        if atom.y < bgfsize[2]:
            bgfsize[2] = atom.y
        elif atom.y > bgfsize[3]:
            bgfsize[3] = atom.y

        if atom.z < bgfsize[4]:
            bgfsize[4] = atom.z
        elif atom.z > bgfsize[5]:
            bgfsize[5] = atom.z

    bgfsize[0] -= margin
    bgfsize[1] += margin
    bgfsize[2] -= margin
    bgfsize[3] += margin
    bgfsize[4] -= margin
    bgfsize[5] += margin

    return bgfsize
"""


#-------------------------------------
#
# - Print information
#
if __name__ == '__main__':

    # get directory:
    directory = dir()

    # set imported stuff we don't want to see:
    imported = ['os', 'sys', 'math', 'string', 'copy', 're', 'cu']

    # print __doc__ for the module:
    print("\n")
    print("-"*60)
    if 'version' not in directory:  version = '??????'
    print("%-45s%15s" % (os.path.basename(sys.argv[0]), 'ver: '+version))

    print("-"*60)
    print(__doc__)


    # import types:
    import types

    # create hash-table:
    hashtable = {}
    for item in directory:
        actual_item = eval(item)
        if item in imported:
            # don't show imported stuff:
            pass
        elif type(actual_item) is types.ModuleType:
            # don't discuss other modules:
            pass
        elif type(actual_item) is types.FunctionType:
            # show __doc__s for functions:
            hashtable[item] = actual_item.__doc__
        elif type(actual_item) is types.ClassType:
            # show __doc__s for classes:
            title = item+' class: '
            hashtable[item] = title +  ( '-' * (60-len(title)) )
            hashtable[item] += actual_item.__doc__

            # show __doc__s for class elements:
            for classItem in dir(actual_item):
                actual_class_item = eval(item+'.'+classItem)
                if type(actual_class_item) is types.ModuleType:
                    # don't discuss other modules:
                    pass
                elif type(actual_class_item) is types.UnboundMethodType \
                         or type(actual_class_item) is types.MethodType:
                    # show __doc__s for functions:
                    hashtable[item] += actual_class_item.__doc__ + '\n'

                elif classItem in ['__doc__','__module__']:
                    pass
                else:
                    # for other stuff show the value:
                    hashtable[item] += '\n'+classItem+' = '+str(actual_class_item)+'\n'

            hashtable[item] +=  ( '-'*60 )+'\n\n'
            
        elif item[0] != '_':
            # for other stuff show the value:
            hashtable[item] = '\n'+item+' = '+str(actual_item)+'\n'

    # print info out
    keys = hashtable.keys()
    keys.sort()
    
    print("Contents:")
    print("-"*60 )
    for item in keys:
        print(hashtable[item])

    print("\n")
    print("-"*60 )
    print("contact: noische@kaist.ac.kr")
    print("-"*60 )
    print("\n")
    
    # done!

