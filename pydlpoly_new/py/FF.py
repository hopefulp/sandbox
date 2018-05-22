# -*- coding: utf-8 -*-
"""
   FF.py
   
   implements a class to represent "a force field"
   primary target: read a tinker prm file and be able to write it as various 
   input files for MM codes (primary target tinker and dl_poly2)
   
   Note: this is not a complete implementation of tinker mm3
         all things not used in our forcefields (like pi-bonding) are ignored

   Remark: all settings are converted to lower case upon reading
"""

import numpy as num
import string
import copy

#########  some helpers

def numberify_list(l, formatcode):
    if len(formatcode) == 0:
        return l
    elif formatcode[0] == "*":
        if formatcode[1] == "i":
            try:
                l = map(string.atoi,l)
            except ValueError:
                print "Error converting %s to integer!!" % str(l)
                raise ValueError
        elif formatcode[1] == "f":
            l = map(string.atoi,l)
        elif formatcode[1] == "s":
            pass
        else:
            raise ValueError, "unknown formatcode %s" % formatcode
    else:
        for i in xrange(len(formatcode)):
            if formatcode[i] == "i":
                try:
                    l[i] = string.atoi(l[i])
                except ValueError:
                    print "Error converting %s to integer!!" % str(l[i])
                    raise ValueError
            if formatcode[i] == "f":
                try:
                    l[i] = string.atof(l[i])
                except ValueError:
                    print "Error converting %s to float!!" % str(l[i])
                    raise ValueError
                except IndexError:
                    print "Formatcode: %s" % formatcode
                    print "Data      : %s" % str(l)
                    raise IndexError('Maybe is the last opbend parameter missing? (not sure)')
            if formatcode[i] == "o":
                # optional float .. need to check if l is long enough ... o must be at the end!!
                if (len(l) > i) :
                    l[i] = string.atof(l[i])
                else:
                    l.append(0.0)
    return l

def with_index(seq):
    for i in xrange(len(seq)):
        yield i, seq[i]

def replace_all(seq, obj, replacement):
    for i, elem in with_index(seq):
        if elem == obj:
            seq[i] = replacement

#########


class FF:
    
    def __init__(self, name, verbose):
        self.name = name
        self.verbose = verbose
        self.pterms = {
#            "atom"        : (1, "ssifi") , \
            "atom"        : (1, "so") , \
            "vdw"         : (1, "ff")    ,\
            "vdwpr"       : (2, "ff")    ,\
            "hbond"       : (2, "ff")    ,\
            "bond"        : (2,"ffo")   ,\
            "bondq"       : (2,"ffff")   ,\
            "bond5"       : (2,"ffo")   ,\
            "bond4"       : (2,"ff")   ,\
#            "bond3"       : (2,"ff")   ,\
#            "electneg"    : (3,"f")    ,\
            "angle"       : (3, "ffoo") ,\
            "angleq"       : (3, "ffff") ,\
            "angle5"      : (3, "ffoo") ,\
            "angle4"      : (3, "ffoo") ,\
#            "angle3"      : (3, "ffoo") ,\
            "anglef"      : (3, "ffioo")  ,\
            "anglef-2"    : (3, "ff")   ,\
            "strbnd"      : (1, "fff") ,\
            "angang"      : (4, "f") ,\
#            "opbend"      : (2, "f")   ,\
            "opbend"      : (4, "fo")   ,\
            "opbendq"      : (4, "fff")   ,\
#            "torsion"     : (4, "ffiffiffioo") ,\
            "torsion"     : (4, "fffo") ,\
#            "torsion5"    : (4, "ffiffiffioo") ,\
            "torsion5"    : (4, "fffo") ,\
#            "torsion4"    : (4, "ffiffiffioo") ,\
            "torsion4"    : (4, "fffo") ,\
#            "strtors"     : (2, "fff") ,\
            "charge"      : (1, "f") ,\
            "chargemod"   : (2, "f") ,\
            "chargeadd"   : (1, "f") ,\
#            "dipole"      : (2, "ff") ,\
#            "dipole5"     : (2, "ff") ,\
#            "dipole4"     : (2, "ff") ,\
#            "dipole3"     : (2, "ff") ,\
#            "piatom"      : (1, "fff") ,\
#            "pibond"      : (2, "ff") ,\
#            "pibond5"     : (2, "ff") ,\
#            "pibond4"     : (2, "ff") ,\
            "molname"     : (1, "*s") ,\
            "virtbond"    : (2, "") ,\
            "restrain-distance"  : (2, "ff"),\
            "restrain-angle"     : (3, "ff"),\
            "restrain-torsion"   : (4, "ff"),\
            }
#        self.virtual_atoms = []
        self.string_settings = ["forcefield", "vdwtype", "radiusrule", "radiustype", "radiussize", \
                "epsilonrule", "strbndtype", "spacegroup", "ewald", "parameters", "opbendtype", \
                "opbendpot", "bondtype","chargetype","dispersion", "rigid", "freeze","sqp", "sqpp",\
                "sqp-next", "opbend-next", "dihedral-bridge", "raw_buck"]
        # active settings must be string settings. if a certain switch is set up, a function is called
        #    usually in order to change the format of a potential type
        self.active_settings = {\
            "strbndtype": [(["mmff"], self.set_strbend_mmff)] ,\
            "opbendtype" : [(["mmff"],self.set_opbend_mmff)],\
            "chargetype": [ (["gaussian"],   self.set_charge_gaussian),\
                            (["gauss_core"], self.set_charge_gausscore),\
                            (["gauss_sp"],   self.set_charge_gauss_sp),\
                          ],\
            "dispersion": [(["d3"], self.set_vdw_use_D3)],\
            "raw_buck"  : [(["on"], self.set_use_raw_buckingham)],\
            }
                                  
        # generate a dictionary of dictionaries for the parameters of the potential terms
        self.params = {}
        self.variables = {}
        self.convariables = {}
        self.ptnames = self.pterms.keys()
        for pt in self.ptnames:
            self.params[pt] = {}
        # generate a dictionary for all the other "switches"
        self.settings = {}
        #self.settings['version'] = None
        # set defaults
        self.settings["dispersion"] = [""]
        self.settings["raw_buck"]   = ["off"]
        # generate a dictionary for the atomtype equivalences
        self.equivalences = {}
        return
    
    def read_prm(self, key_filename):
        f = open(key_filename, "r")
        if self.verbose: print("reading %s tinker key file" % key_filename)
        # read this until file ends
        line = f.readline()
        stop = False
        prm_file = None
        while not stop:
            sline = string.split(line)
            if ((len(sline)>0) and (sline[0][0] != "#")):
                keyword = string.lower(sline[0])
                # check if keyword is "parameters" .. then we have to read the file before going on
                # this is some kind of include mechanism
                if (keyword == "parameters") and (sline[1] != "none") :
                    prm_file = sline[1] + ".prm"
                    if self.verbose: print ("including prm file %s" % prm_file)
                    f_key = f
                    # switch to prm file and just keep going
                    f = open(prm_file, "r")
                else:
                    if sline.count("!!"):
                        # there is a note .. cut it off
                        note = sline.index("!!")
                        sline = sline[:note]
                    # check if this is a potential term keyword
                    if self.ptnames.count(keyword):
                        # yes, it is a term
                        natoms, formatcode = self.pterms[keyword]
                        atomkey = string.join(sline[1:natoms+1], ":")
                        self.params[keyword][atomkey] = numberify_list(sline[natoms+1:], formatcode)
                    elif ((keyword == 'var') and (self.ptnames.count(string.lower(sline[1])))):
                        keyword = string.lower(sline[1])
                        natoms, formatcode = self.pterms[keyword]
                        atomkey = string.join(sline[2:natoms+2], ":")
                        paramindex = string.atoi(sline[natoms+2])-1
                        if keyword not in self.variables:
                            self.variables[keyword]={}
                        if atomkey not in self.variables[keyword]:
                            self.variables[keyword][atomkey]=[{},[]]
                        if paramindex not in self.variables[keyword][atomkey][0]:
                            self.variables[keyword][atomkey][0][paramindex]='xxx'
                        else:
                            print 'ERROR: The variable parameter %s %s %s is define twice --> check your keyfile' % (keyword, atomkey, str(paramindex))  
                            raise ValueError
                        if ((formatcode[paramindex] == 'f') or (formatcode[paramindex] == 'o')):
                            self.variables[keyword][atomkey][0][paramindex]=map(string.atof, sline[natoms+3:])
                        else:
                            print 'ERROR: The parameter %s %s %s you have chosen to optimize is not a float Parameter --> not possible to optmize' % (keyword, atomkey, str(paramindex))
                            raise ValueError
                    elif ((keyword == 'var') and (string.lower(sline[1]) == 'con') and (self.ptnames.count(string.lower(sline[2])))):
                        keyword = string.lower(sline[2])
                        natoms, formatcode = self.pterms[keyword]
                        atomkey1 = string.join(sline[3:natoms+3], ":")
                        paramindex1 = string.atoi(sline[natoms+3])-1
                        atomkey2 = string.join(sline[natoms+4:natoms+natoms+4], ":")
                        paramindex2 = string.atoi(sline[2*natoms+4])-1
                        if keyword not in  self.convariables:
                            self.convariables[keyword] = []
                        self.convariables[keyword].append([atomkey1, atomkey2, paramindex1, paramindex2])
                    elif keyword == "equivalent":
                        self.equivalences[sline[1]] = sline[2]
#                    elif keyword == "virtual":
#                        atomlist = []
#                        line = sline[1:]
#                        line.reverse()
#                        while len(line) > 0:
#                            a = line.pop()
#                            if a == "range":
#                                first = string.atoi(line.pop())
#                                last  = string.atoi(line.pop())
#                                atomlist += range(first-1, last)
#                            else:
#                                atomlist.append(string.atoi(line.pop())-1)
#                        self.virtual_atoms.append(atomlist)
                    else:
                        # ok this is just a regular setting
                        if self.string_settings.count(keyword):
                            if len(sline) >= 2:
                                self.settings[keyword] = map(string.lower, sline[1:])
                            else:
                                self.settings[keyword] = "ON"
                            if self.active_settings.keys().count(keyword):
                                active = self.active_settings[keyword]
                                for act in active:
                                    if self.settings[keyword] == act[0]:
                                        # yes, the switch is set. execute the function
                                        act[1]()
                        else:
                            self.settings[keyword] = string.atof(sline[1])
            line = f.readline()
            if len(line) == 0:
                if prm_file:
                    # we are reading parameters file .. switch back to original key and continue
                    f.close()
                    f = f_key
                    prm_file = None
                else:
                    f.close()
                    stop = True
#        if self.variables != {}:
#            self.check_variables()
#            self.sort_variablekeys()
            #self.get_variables()
        if self.settings['version'] != 2.0:
            print 'Error: The actual version of pydlpoly can only handle input files version 2!\n \
                Use the script convert_key to vonvert to the new file format. For helf type in your shell:\
                convert_key -h'
#        else:
#          self.settings['strbndtype'] = 'mmff'
#          self.settings['opbendtype'] = 'mmff'
#          self.set_strbend_mmff()
#          self.set_opbend_mmff()
#        print self.params
        return

                
                
    def get_params(self, term, types, permute=False, reverse=True, verbose=True):
        params = None
        equi = None
        if permute:
            tl = []
            ts = string.split(types,":")
            perm = self.permute3(ts[1:])
            for p in perm:
                tl.append(string.join([ts[0]]+p, ":"))
        elif reverse:
            tl = [types, self.reverse(types)]
        else:
            tl = [types]
        term_params = self.params[term]
        for t in tl:
            try:
                params = term_params[t]
                return params, equi
            except KeyError:
                pass
        if not params and term not in ('chargemod','chargeadd'):
            tl = string.split(types, ":")
            for t in tl:
                if self.equivalences.has_key(t):
#                    if verbose and self.verbose: print (" ----> found equivalence for %s : replacing with %s " % (t, self.equivalences[t]))
                    etl = copy.copy(tl)
                    replace_all(etl, t, self.equivalences[t])
                    # recursively call it again
                    etypes = string.join(etl, ":")
                    params, bla = self.get_params(term,etypes, permute=permute, reverse=reverse, verbose=verbose)
                    if params:
                        return params, etypes 
        return params, equi

    def get_params_and_pottype(self, term_list, types, permute = False):
        for term in term_list:
            params, equi = self.get_params(term, types, permute)
            if params:
                # we found parameters and are done (the first case is used!)
                return (params, term, equi)
        # if we have reached this point then no parameters were found in all terms
        return (None, None, None)
    
    def set_strbend_mmff(self):
        if self.verbose: print ("setting format of strbend cross term to mmff")
        self.pterms["strbnd"] = (3, "fff")
        return
            
    def set_opbend_mmff(self):
        if self.verbose: print ("setting format of opbend to mmff: 4 atomic types")
        self.pterms["opbend"] = (4, "fo")
        return
        
    def set_charge_gaussian(self):
        if self.verbose: print ("setting format of charge, chargemod and chargeadd to gaussian: include sigma")
        self.pterms["charge"]    = (1, "ff")
        self.pterms["chargemod"] = (2, "ff")
        self.pterms["chargeadd"] = (2, "ff")
        return

    def set_charge_gausscore(self):
        if self.verbose: print ("setting format of charge, chargemod and chargeadd to gausscore: include sigma and core charge")
        self.pterms["charge"]    = (1, "fff")
        self.pterms["chargemod"] = (2, "fff")
        self.pterms["chargeadd"] = (2, "fff")
        return 

    def set_charge_gauss_sp(self):
        if self.verbose: print ("setting format of charge, chargemod and chargeadd to gaus_sp: include sigma for s and p")
        self.pterms["charge"]    = (1, "fff")
        self.pterms["chargemod"] = (2, "fff")
        self.pterms["chargeadd"] = (2, "fff")
        return 

    def set_vdw_use_D3(self):
        if self.verbose: print ("setting format of vdwpr to use D3: 1/B, A, C6, C8, Rab")
        self.pterms["vdwpr"]     = (2, "fffff")
        return

    def set_use_raw_buckingham(self):
        if self.verbose: print ("setting vdwpr to read raw buckingham values A, B, C")
        self.pterms["vdwpr"]     = (2, "ffff")
        return
        
    def reverse(self, t):
        l = string.split(t,":")
        l.reverse()
        return string.join(l,":")

    def permute3(self, tl):
        return [[tl[0],tl[2],tl[1]],[tl[1],tl[0],tl[2]],[tl[2],tl[0],tl[1]], \
                [tl[1],tl[2],tl[0]],[tl[2],tl[1],tl[0]],[tl[0],tl[1],tl[2]]]
        

            
if __name__ == "__main__":
    ff = FF("mof-ff")
    ff.read_prm("mof5.key_ref")
    #print ff.settings
            
            
