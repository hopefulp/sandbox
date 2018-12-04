#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xmlrpclib
from xmlrpclib import ServerProxy
import molsys.stow as stow
import sys
import logging
import os
import molsys
from assign_FF import sort_bond, sort_angle, sort_dihedral, sort_oop
import getpass
logger = logging.getLogger("mofplus")
logger.setLevel(logging.DEBUG)
shandler = logging.StreamHandler()
shandler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m-%d %H:%M')
shandler.setFormatter(formatter)
logger.addHandler(shandler)


def download(dtype, binary = False):
    """
    mfp download decorator
    """
    def download_decorator(func):
        def inner(*args, **kwargs):
            try:
                lines = func(*args, **kwargs)
                if "mol" in kwargs.keys():
                    if kwargs["mol"] == True:
                        if dtype == "topology":
                            m = molsys.topo()
                        else:
                            m = molsys.mol()
                        m.fromString(lines)
                        return m
                if binary == False:
                    f=open(str(args[1])+'.mfpx', 'w')
                    f.write(lines)
                    f.close()
                else:
                    with open("%s.hdf5" % str(args[1]), "wb") as handle:
                        handle.write(lines)
                logger.info('%s %s downloaded from mofplus' % (dtype,args[1]))
            except xmlrpclib.Fault:
                logger.error('Requested %s %s not available on mofplus' % (dtype, args[1]))
        return inner
    return download_decorator

def faulthandler(func):
    def inner(*args,**kwargs):
        ret = func(*args, **kwargs)
        if type(ret) == dict:
            logger.error(ret['faultString'])
            return ret["faultString"]
        return ret
    return inner

class user_api(object):
    """
    Via the user_api class the API routines of MOFplus which are accessable for normal users can be used.
    The credentials can be set either as environment variables MFPUSER and MFPPW or can be given interactively.

    :Attrributes[MaÈ[MaÇ[MaÇ:
        - mfp (obj)     : ServerProxy XMLRPCLIB object holding the connection to MOF+
        - username (str): Username on MOF+
        - pw (str)      : Pa[MaÇ[MaÇssword corresponding to the username

    :Args:
        - experimental (bool, option[MaÆal): Use to connect to experimental DB, defaults to False
        - banner       (bool, optional): If True, the MFP API banner is printed to SDTOUT, defaults to False
    """

    def __init__(self,experimental=False, banner = False):
        if banner: self.print_banner()
        try:
            self.username = os.environ['MFPUSER']
            self.pw       = os.environ['MFPPW']
        except KeyError:
            logger.error("Credentials not found!")
            self.username, self.pw = self.credentials_from_cmd()
        if experimental:
	    self.mfp = ServerProxy('https://%s:%s@www.mofplus.org/MFP_JPD/API/user/xmlrpc' % (self.username, self.pw))
        else:
            #self.mfp = ServerProxy('https://%s:%s@www.mofplus.org/API/call/xmlrpc' % (username, pw))
            self.mfp = ServerProxy('https://%s:%s@www.mofplus.org/API/user/xmlrpc' % (self.username, self.pw))
        self.check_connection()
        return

    def credentials_from_cmd(self):
        """
        Method to get the credentials from the command line
        """
        username = raw_input("Email:")
        pw       = getpass.getpass()
        return username, pw

    def check_connection(self):
        """
        Method to check if the connection to MFP is alive
        """
        try:
            self.mfp.add(2,2)
            logger.info("Connection to user API established")
        except xmlrpclib.ProtocolError:
            logger.error("Not possible to connect to MOF+. Check your credentials")
            exit()
        return

    def print_banner(self):
        """
        Prints the MFP banner
        """
        print  ":##::::'##::'#######::'########:::::::::::::::'###::::'########::'####:\n\
:###::'###:'##.... ##: ##.....::::'##::::::::'## ##::: ##.... ##:. ##::\n\
:####'####: ##:::: ##: ##::::::::: ##:::::::'##:. ##:: ##:::: ##:: ##::\n\
:## ### ##: ##:::: ##: ######:::'######::::'##:::. ##: ########::: ##::\n\
:##. #: ##: ##:::: ##: ##...::::.. ##.::::: #########: ##.....:::: ##::\n\
:##:.:: ##: ##:::: ##: ##::::::::: ##:::::: ##.... ##: ##::::::::: ##::\n\
:##:::: ##:. #######:: ##:::::::::..::::::: ##:::: ##: ##::::::::'####:\n\
:..:::::..:::.......:::..:::::::::::::::::::..:::::..::..:::::::::....:"

   
    @download('topology')
    def get_net(self,netname, mol = False):
        """Downloads a topology in mfpx file format
        :Parameters:
            -netname (str): name of the net
            -mol    (bool,optional): if true a mol object is returned, if false
                            topology is written to a file, defaults to False
        """
        lines = self.mfp.get_net(netname)
        return lines

    def get_list_of_nets(self):
        """Returns a list of all topologies in the db"""
        return self.mfp.get_list_of_nets()
    

    def get_list_of_bbs(self):
        """Returns a list of all BBS in the db"""
        return self.mfp.get_list_of_bbs()

    @download('building block')
    def get_bb(self,bbname, mol = False):
        """Downloads a bb in mfpx file format
        :Parameters:
            -bbname (str): name of the bb
            -mol    (bool,optional): if true a mol object is returned, if false
                            bb is written to a file, defaults to False
        """
        lines = self.mfp.get_bb(bbname)
        return lines
    
    @download('MOF')
    def get_mof_structure_by_id(self,strucid, mol = False):
        """Downloads a MOF structure in mfpx file format
        :Parameters:
            -strucid (str): id of the MOF structure in the DB
            -mol    (bool,optional): if true a mol object is returned, if false
                            bb is written to a file, defaults to False
        """
        lines,name = self.mfp.get_mof_structure_by_id(strucid)
        return lines

    def get_cs(self,name):
        """
        Returns the coordinations sequences of a topology as a list of lists.
        :Parameters:
            -name (str): Name of the topology
            
        """
        return self.mfp.get_cs(name)
    
    def get_vs(self,name):
        """
        Returns the vertex symbol of a topology as a list of strings
        :Parameters:
            -name (str): Name of the topology
        """
        return self.mfp.get_vs(name)

    def search_cs(self, cs, vs, cfilter = True):
        """
        Searches nets with a given coordination sequences and given vertex symbols and returns
        the corresponding netnames as a list of strings.
        :Parameters:
            -cs (list): List of the coordination sequences
            -vs (list): List of the vertex symbols
            -cfilter (bool): If true no catenated nets are returned, defaults to True
        """
        assert type(cs) == list
        assert type(vs) == list
        nets = self.mfp.search_cs(cs, vs)
        rl = []
        if cfilter:
            for i,n in enumerate(nets):
                if n.find('-c') != -1: rl.append(n)
            for i in rl: nets.remove(i)
        return nets

class admin_api(user_api):

    """
    Via the admin_api class the API routines of MOFplus which are accessable for normal users and for admin users
    can be used. Class is inherited from the user_api class.
    
    The credentials can be set either as environment variables MFPUSER and MFPPW or can be given interactively.

    :Attrributes:
        - mfp (obj)     : ServerProxy XMLRPCLIB object holding the connection to MOF+
        - username (str): Username on MOF+
        - pw (str)      : Password corresponding to the username

    :Args:
        - experimental (bool, optional): Use to connect to experimental DB, defaults to False
        - banner       (bool, optional): If True, the MFP API banner is printed to SDTOUT, defaults to False
    """

    def __init__(self, experimental = False, banner = False):
        user_api.__init__(self,experimental, banner)
        if experimental:
	    self.mfp = ServerProxy('https://%s:%s@www.mofplus.org/MFP_JPD/API/admin/xmlrpc' % (self.username, self.pw))
        else:
            self.mfp = ServerProxy('https://%s:%s@www.mofplus.org/API/admin/xmlrpc' % (self.username, self.pw))
        self.check_adminconnection()

    def check_adminconnection(self):
        """
        Method to check if the connection to MFP is alive
        """
        try:
            self.mfp.add2(2,2)
            logger.info("Connection to admin API established")
            print """
            We trust you have received the usual lecture from the MOF+ system administrator.
            It usually boils down to these two things:
                #1) Think before you type.
                #2) With great power comes great responsibility.
            """
        except xmlrpclib.ProtocolError:
            logger.error("Not possible to connect to MOF+ admin API. Check your credentials")
            exit()
        return

    def delete_net(self, name):
        """
        Deletes a net from the db
        :Parameters:
            -name (str): name of the net
        """
        assert type(name) == str
        self.mfp.delete_net(name)
   
    def add_bb_penalties(self,data):
        """
        Method to adds penalties to building blocks
        """
        retstring = self.mfp.add_bb_penalties(data)
        print retstring
        return

    def upload_weaver_run(self, fwid, fname):
        """
        Method to upload the results of a weaver run to the db
        :Parameters:
            -fwid: firework id of the job
            -fname: filename of the structure file
        """
        data = {}
        data['fwid'] = str(fwid)
        f = open(fname, 'r')
        data['fmfpx'] = f.read()
        a = self.mfp.upload_weaver_run(data)
        return

    def upload_mof_structure_by_id(self, fname, strucid):
        """
        Method to upload a structure file to the DB
        :Parameters:
            - fname (str): path to the mfpx file
            - strucid (int): id of the structure in the db
        """
        data = {}
        f = open(fname, 'r')
        data['id'] = strucid
        data['fmfpx'] = f.read()
        self.mfp.upload_mof_structure_by_id(data)
        return

    def upload_topo_file_by_name(self, fname, name):
        """
        Method to upload a topo file to the DB
        :Parameters:
            - fname (str): path to the mfpx file
            - name (str): name of the topology
        """
        data = {}
        f = open(fname, 'r')
        data['name'] = name
        data['fmfpx'] = f.read()
        self.mfp.upload_topo_file_by_name(data)
        return

    ### method in principle obsolete
    def upload_pa_run(self,data):
        ret = self.mfp.upload_pa_run(data)
        return

    def upload_bbfile_by_name(self, fname, name):
        """
        Method to upload a bb file to the DB
        :Parameters:
            - fname (str): path to the mfpx file
            - name (str): name of the bb
        """
        data = {}
        f = open(fname, 'r')
        data['name'] = name
        data['fmfpx'] = f.read()
        self.mfp.upload_bbfile_by_name(data)
        return

    def insert_bb(self,name, fname, chemtype, frag = False):
        """
        Method to create a new entry in the bb table.
        :Parameters:
            - name (str): name of the bb
            - fname (str): path to the mfpx file
            - chemtype (str): string describing the character of the bb
            - frag (bool, optional): Option to set a BB as fragment, defaults to False
        """
        data = {}
        data['name'] = name
        data['fmfpx'] = open(fname, 'r').read()
        data['type'] = chemtype
        self.mfp.insert_bb(data)
        return

    def set_cs(self, name, cs):
        """
        Method to set the cs of a topology.
        :Parameters:
            - name (str): name of the topology
            - cs (list): list of lists with the cs
        """
        data = {}
        data['name'] = name
        data['cs'] = cs
        self.mfp.set_cs(data)
        return
    
    def set_vs(self, name, vs):
        """
        Method to set the vs of a topology.
        :Parameters:
            - name (str): name of the topology
            - vs (list): list with the vs
        """
        data = {}
        data['name'] = name
        data['vs'] = vs
        self.mfp.set_vs(data)
        return
    
    def connect_nets(self, pnet, cnet, pattern):
        """
        Method to create relationchips between nets in the DB
        :Parameters:
            - pnet (str): name of the parent net
            - cnet (str): name of the child net
            - pattern (str): derivation type
        """
        assert type(pnet) == str
        assert type(cnet) == str
        assert type(pattern) == str
        assert cnet != pnet
        self.mfp.connect_nets(pnet,cnet,pattern)
        return

    def add_skal_property(self, strucid, ptype, prop):
        """
        Method to add a skalar property to a structure
        :Parameters:
            - strucid (int): id of the structure in the DB
            - ptype (str): name of the property
            - prop (float): property value
        """
        assert type(strucid) == int
        assert type(ptype) == str
        self.mfp.add_skal_property(strucid, ptype, prop)
        return

    def add_xy_property(self,strucid,ptype,data):
        """
        Method to add a dataset as property to the DB
        :Parameters:
            - strucid (int): id of the structure in the DB
            - ptype (str): name of the property
            - data (dict): dataset as dictionary 
        """
        assert type(strucid) == int
        assert type(ptype) == str
        self.mfp.add_xy_property(strucid, ptype,data)
        return

    def fa_finish(self,faid):
        """
        Method to register a fireanalyzer run as finished
        :Parameters:
            - faid (int): id of fireanalyzer run
        """
        assert type(faid) == int
        self.mfp.fa_finish(faid)


class FF_api(admin_api):

    def format_atypes(self, atypes, ptype):
        """
        Helper function to extract fragments out of atypes and to 
        order atypes and fragements in dependence of the ptype.
        """
        assert type(ptype) == str
        if ptype == "bond":
            atypes = sort_bond(atypes)
        elif ptype == "angle":
            atypes = sort_angle(atypes)
        elif ptype == "oop":
            atypes = sort_oop(atypes)
        elif ptype == "dihedral":
            atypes = sort_dihedral(atypes)
        latypes = atypes.split(":")
        atypes = []
        fragments = []
        for at in latypes: 
            atypes.append(at.split("@")[0])
            fragments.append(at.split("@")[1])
        assert len(atypes) == len(fragments)
        return atypes, fragments

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
        atypes, fragments = self.format_atypes(atypes,ptype)
        params = self.mfp.get_params(FF, atypes, fragments, ptype, potential, fitsystem)
        return params

    @faulthandler
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
        assert type(FF) == type(ptype) == type(potential) == type(atypes) == str
        assert type(params) == list
        atypes, fragments = self.format_atypes(atypes,ptype)
        ret = self.mfp.set_params(FF, atypes, fragments, ptype, potential, fitsystem,params)
        return ret

    def list_FFrefs(self):
        """
        Method to list names and meta properties of all available reference systems in the DB
        """
        return self.mfp.list_FFrefs()

    def set_FFref(self, name, hdf5path, mfpxpath, comment=""):
        """
        Method to create a new entry in the FFref table and to upload a file with
        reference information in the hdf5 file format.
        :Parameters:
            - name (str): name of the entry in the DB
            - path (str): path to the hdf5 reference file
        """
        assert type(name) == type(hdf5path) == type(mfpxpath) == type(comment) == str
        with open(hdf5path, "rb") as handle:
            binary = xmlrpclib.Binary(handle.read())
        with open(mfpxpath, "r") as handle:
            mfpx = handle.read()
        self.mfp.set_FFref(name, binary, mfpx, comment)
        return

    def set_FFref_graph(self,name, mfpxpath):
        with open(mfpxpath, "r") as handle:
            mfpx = handle.read()
        self.mfp.set_FFref_graph(name,mfpx)
        return

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

    def set_FFfrag(self,name,path,comment=""):
        """
        Method to create a new entry in the FFfrags table.
        :Parameters:
            - name (str): name of the entry in the db
            - path (str): path to the mfpx file of the fragment
            - comment (str): comment
        """
        assert type(name) == type(path) == type(comment) == str
        with open(path, "r") as handle:
            lines = handle.read()
            m = molsys.mol()
            m.fromString(lines, ftype = "mfpx")
            prio = m.natoms-m.elems.count("x")
        self.mfp.set_FFfrag(name, lines, prio, comment)
        return

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

    def get_parameter_history(self, id):
        assert type(id) == int
        return self.mfp.get_parameter_history(id)

    def get_FFfit(self, id):
        return 

    def set_FFfit(self,id):
        return

if __name__ == '__main__':
    option = [
            ['', 't', 'topology', "Name of topology which is downloaded from mofplus"],
            ['', 'b', 'buildingblock', "Name of building block which is downloaded from mofplus"],
            ]
    shellval = stow.main(stow.sys.argv[1:], option)
    if shellval[0] != '' or shellval[1] != '':
        api = FF_api(banner=False, experimental = False)
        if shellval[0] != '': api.get_net(shellval[0], mol = False)
        if shellval[1] != '': api.get_bb(shellval[1], mol = False)
    else:
        api = FF_api(banner=True, experimental = False)

