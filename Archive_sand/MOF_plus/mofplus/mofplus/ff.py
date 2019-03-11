#!/usr/bin/env python
# -*- coding: utf-8 -*-

import string
import logging
import sys
from decorator import faulthandler, download
import user
from molsys.util.aftypes import aftype, aftype_sort, afdict

logger = logging.getLogger("mofplus")

bodymapping = {1:"onebody", 2:"twobody",3:"threebody",4:"fourbody"}

allowed_ptypes = {1: ["charge", "vdw", "equil"],
        2: ["bnd", "chargemod", "vdwpr"],
        3: ["ang"],
        4: ["dih", "oop"]
        }

allowed_potentials = {"charge": [["point",1], ["gaussian",2], ["slater",2]],
        "equil": [["equil", 1]],
        "vdw": [["LJ",2], ["buck",2], ["buck6d",2]],
        "bnd": [["harm",2], ["mm3",2], ["quartic",5], ["morse",3], ["equiv", 2]],
        "chargemod": [["point",1], ["gaussian",2], ["slater",2]],
        "vdwpr": [["LJ",2], ["buck",2], ["damped_buck",2]],
        "ang": [["harm",2],["mm3",2], ["quartic",5], ["fourier",5],  ["strbnd", 6]],
        "dih": [["harm",2], ["cos3",3], ["cos4",4]],
        "oop": [["harm",2]]}



class FF_api(user.user_api):   
    """
    Via the ff_api class the API routines of MOFplus concerning the retrieval of MOF-FF parameters can be used.
    The credentials can be set either as environment variables MFPUSER and MFPPW or can be given interactively or
    can be stated in ~/.mofplusrc.
    
    The FF_api class inherits from the user_api class.

    :Attrributes:
        - mfp (obj)     : ServerProxy XMLRPCLIB object holding the connection to MOF+
        - username (str): Username on MOF+
        - pw (str)      : Password corresponding to the username

    :Args:
        - local        (bool, optional): Use to connect directly to a MySQL server, defaults to False
        - localhost    (bool, optional): Use to connect to an MFP server running on localhost, defaults to False
        - banner       (bool, optional): If True, the MFP API banner is printed to SDTOUT, defaults to False
    """

    
    def format_atypes(self, atypes, ptype, potential):
        """
        Helper function to extract fragments out of atypes and to
        order atypes and fragments in dependence of the ptype. 
        
        :Parameters:
            - atypes (list): list of atom types in the form "atype@fragtype"
            - ptype (str): ric type
            - potential (str): potential type
        """
        assert type(ptype) == str
        ### split into tuples of aftypes, then sort and then split into 
        ### frags and atypes
        aftypes = []
        for i in string.split(atypes,":"): aftypes.append(aftype(*i.split("@")))
        aftypes = aftype_sort(aftypes, ptype)
        atypes  = [i.atype for i in aftypes]
        fragments = [i.fragtype for i in aftypes]
        if ptype not in allowed_ptypes[len(atypes)]:
            raise ValueError("ptype %s not allowed for %s term" % (ptype, bodymapping[len(atypes)]))
        if potential not in [i[0] for i in allowed_potentials[ptype]]:
            raise ValueError("potential %s not allowed for ptype %s" % (potential, ptype))
        return atypes, fragments

    def get_params_from_ref(self, FF, ref):
        """
        Method to look up all FF parameters that are available for a reference system
        
        :Parameters:
            - FF (str): Name of the FF the parameters belong to
            - ref (str): Name of the reference system the parameters belong to
        """
        paramsets = self.mfp.get_params_from_ref(FF,ref)
        paramdict = {"onebody":{"charge":afdict(),"vdw":afdict(),"equil":afdict()},
                "twobody":{"bnd":afdict(),"chargemod":afdict(), "vdwpr":afdict()},
                "threebody":{"ang":afdict()},
                "fourbody": {"dih":afdict(),"oop":afdict()}}
        # RS (explanation to be improved by JPD)
        # paramset is a nested list of lists provided by MOF+
        # it is resorted here in to a number of nested directories for an easier retrieval of data
        # i loops over the lists from paramset
        # each entry is
        #      i[0] : atomtype (len i[0] determines n-body via gloabl bodymapping)
        #      i[1] : fragment
        #      i[2] : type (e.g. charge, vdw, equiv)   TODO: change euilv -> equiv
        #      i[3] : ptype
        #      i[4] : paramstring
        for i in paramsets:
            #typestr =""
            #for a,f in zip(i[0],i[1]):
            #    typestr+="%s@%s:" % (a,f)
            ## cut off last ":"
            #typestr = typestr[:-1]
            typelist = [aftype(a,f) for a,f in zip(i[0],i[1])]
            typedir = paramdict[bodymapping[len(i[0])]][i[2]]
            tt = tuple(typelist)
            if tt in typedir:
                # another term .. append
                typedir.appenditem(tt, (i[3],i[4]))
            else:
                typedir[tt] = [(i[3],i[4])]
        return paramdict
 
    @faulthandler
    def get_params(self,FF, atypes, ptype, potential,fitsystem):
        """
        Method to look up parameter sets in the DB
        
        :Parameters:
            - FF (str): Name of the FF the parameters belong to
            - atypes (list): list of atypes belonging to the term
            - ptype (str): type of requested term
            - potential (str): type of requested potential
            - fitsystem (str): name of the FFfit/reference system the
              parameterset is obtained from
        """
        assert type(FF) == type(ptype) == type(atypes) == type(potential) == str
        atypes, fragments = self.format_atypes(atypes,ptype, potential)
        params = self.mfp.get_params(FF, atypes, fragments, ptype, potential, fitsystem)
        return params
           

    def list_FFrefs(self,FF):
        """
        Method to list names and meta properties of all available reference systems in the DB
        
        :Parameters:
            - FF (str): Name of the FF the reference systems belong to, give "*" in order to 
                get all available references independent from the FF
        """
        res = self.mfp.list_FFrefs(FF)
        dic = {}
        for i in res:
            dic[i[0]] = i[1:]
        return dic


    @download("FFref")
    def get_FFref_graph(self,name, mol = False):
        """
        Downloads the reference system in mfpx file format
        
        :Parameters:
            -name (str): name of the reference system
            -mol    (bool,optional): if true a mol object is returned, if false
                            fragment is written to a file, defaults to False
        """
        assert type(name) == str
        lines = self.mfp.get_FFref_graph(name)
        return lines

    @download("FFref", binary = True)
    def get_FFref(self,name):
        """
        Method to retrieve an reference file in hdf5 file format from the DB
        
        :Parameters:
            - name (str): name of the entry in the DB
        """
        assert type(name) == str
        bstr = self.mfp.get_FFref(name).data
        return bstr

    @download("FFfrag")
    def get_FFfrag(self,name, mol = False):
        """
        Downloads a FFfrag in mfpx file format
        
        :Parameters:
            -name (str): name of the fragment
            -mol    (bool,optional): if true a mol object is returned, if false
                            fragment is written to a file, defaults to False
        """
        assert type(name) == str
        lines = self.mfp.get_FFfrag(name)
        return lines

    def list_FFfrags(self):
        """
        Method to list names and meta properties of all available FFfrags in the DB
        """
        return self.mfp.list_FFfrags()

    def list_special_atypes(self):
        """
        Method to get a dictionary of aftypes with special properties
        """
        res = self.mfp.list_special_atypes()
        dic = {"linear": []}
        for l in res:
            af = aftype(l[0], l[1])
            dic[l[2]].append(af)
        return dic

    def create_fit(self, FF, ref, azone=None, atfix = None, comment = ""):
        """
        Method to create a FFfit entry in the database which is necessary
        for storing parameters for a predefined FF
        :Parameters:
            - FF (str): name of the FF
            - ref (str): name of the reference system
            - azone (list): list of integers describing the active zone of the fit
            - atfix (dict): dictionary containing special atypes information
            - comment (string): comment
        """
        return self.mfp.create_fit(FF, ref, azone, atfix, comment)

    def set_params(self, FF, atypes, ptype, potential, fitsystem,params):
        """
        Method to upload parameter sets in the DB
        :Parameters:
            - FF (str): Name of the FF the parameters belong to
            - atypes (str): list of atypes belonging to the term
            - ptype (str): type of requested term
            - potential (str): type of requested potential
            - params (list): parameterset
            - fitsystem (str): name of the FFfit/reference system the
              parameterset is obtained from
        """
        #assert type(FF) == type(ptype) == type(potential) == type(atypes) == str
        assert type(params) == list
        atypes, fragments = self.format_atypes(atypes,ptype, potential)
        rl = {i[0]:i[1] for i in allowed_potentials[ptype]}[potential]
        if len(params) != rl:
            raise ValueError("Required lenght for %s %s is %i" %(ptype,potential,rl))
        ret = self.mfp.set_params(FF, atypes, fragments, ptype, potential, fitsystem,params)
        return ret
    
    def set_params_interactive(self, FF, atypes, ptype, potential, fitsystem, params):
        """
        Method to upload parameter sets in the DB interactively
        :Parameters:
            - FF (str): Name of the FF the parameters belong to
            - atypes (str): list of atypes belonging to the term
            - ptype (str): type of requested term
            - potential (str): type of requested potential
            - params (list): parameterset
            - fitsystem (str): name of the FFfit/reference system the
              parameterset is obtained from
        """
        stop = False
        while not stop:
            print "--------upload-------"
            print "FF      : %s" % FF
            print "atypes  : " +len(atypes)*"%s " % tuple(atypes)
            print "type    : %s" % ptype
            print "pot     : %s" % potential
            print "ref     : %s" % fitsystem
            print "params  : ",params
            print "--------options---------"
            print "[s]: skip"
            print "[y]: write to db"
            print "[a]: modify atypes"
            print "[t]: modify type"
            print "[p]: modify pot"
            print "[r]: modify ref"
            x = raw_input("Your choice:  ")
            if x == "s":
                stop = True
                print "Entry will be skipped"
            elif x == "y":
                ret = self.set_params(FF, string.join(atypes,":"), ptype, potential, fitsystem, params)
                print ret
                if type(ret) != int:
                    "Error occurred during upload, try again!"
                else:
                    print "Entry is written to db"
                    stop = True
            elif x == "a":
                inp = raw_input("Give modified atypes:  ")
                atypes = string.split(inp)
            elif x == "t":
                ptype = raw_input("Give modified type:  ")
            elif x == "p":
                potential = raw_input("Give modified pot:  ")
            elif x == "r":
                fitsystem = raw_input("Give modified ref:  ")

if __name__ == '__main__':
    import admin
    if len(sys.argv) > 1:   
        if sys.argv[1] == "user":
            api = user.user_api(banner=True, local =False, localhost = False)
        elif sys.argv[1] == "admin":
            api = admin.admin_api(banner=True, localhost = False)
        elif sys.argv == "ff":
            api = FF_api(banner=True, localhost = False, local = False)
    else:
        api = FF_api(banner=True, localhost = False, local = False)



