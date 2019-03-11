#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xmlrpclib
from xmlrpclib import ServerProxy
import logging
import string
import molsys
from decorator import faulthandler, download, nolocal
import ff
logger = logging.getLogger("mofplus")


class admin_api(ff.FF_api):

    """
    Via the admin_api class the API routines of MOFplus which are accessable for normal users and for admin users
    can be used. Class is inherited from the ff_api class.
    
    The credentials can be set either as environment variables MFPUSER and MFPPW or can be given interactively or
    can be stated in ~/.mofplusrc.

    :Attrributes:
        - mfp (obj)     : ServerProxy XMLRPCLIB object holding the connection to MOF+
        - username (str): Username on MOF+
        - pw (str)      : Password corresponding to the username

    :Args:
        - localhost    (bool, optional): Use to connect to an MFP server running on localhost, defaults to False
        - banner       (bool, optional): If True, the MFP API banner is printed to SDTOUT, defaults to False
    """

    def __init__(self, banner = False, localhost = False):
        ff.FF_api.__init__(self, banner=banner, localhost = localhost)
        if localhost:
            self.mfp = ServerProxy('http://%s:%s@localhost/MOFplus_final2/API/admin/xmlrpc' % (self.username, self.pw))
        else:
            self.mfp = ServerProxy('https://%s:%s@www.mofplus.org/API/admin/xmlrpc' % (self.username, self.pw), allow_none = True)
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

    def upload_weaver_run(self, fwid, scid, fname, energy):
        """
        Method to upload the results of a weaver run to the db
        :Parameters:
            -fwid: firework id of the job
            -fname: filename of the structure file
        """
        data = {}
        data['fwid'] = str(fwid)
        data["scid"] = scid
        data["energy"] = energy
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

    def set_structure(self,scid, name, path, ff, properties = {}):
        """
        Sets the structures for a given supercell db entry id.
        :Parameters:
            - scid(int): id of the supercell entry in the db
            - name(str): name of the structure
            - path(str): path to the mfpx file of the structure
            - ff(str): FF/LOT on which the structure was computed
            - properties(dict, optional): dictionary of properties
            of the structures, defaults to {}

        """
        with open(path , "r") as f:
            fstring = f.read()
        self.mfp.set_structure(scid, name, fstring, ff, properties)
        return

    def RTAfinish(self,jobid):
        """
        Method to mark an RTA job in the mofplus db as finished
        
        :Parameters:
            -jobid(int): id of the job in the db
        """
        self.mfp.RTAfinish(int(jobid))
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
    
    def set_FFfrag(self,name,path,comment=""):
        """
        Method to create a new entry in the FFfrags table.
        :Parameters:
            - name (str): name of the entry in the db
            - path (str): path to the mfpx file of the fragment
            - comment (str): comment
        """
        assert type(name) == type(path) == type(comment) == str
        with open(path, "r") as f:
            lines = f.read()
            m = molsys.mol.fromString(lines)
        m = molsys.mol.fromFile(path)
        prio = m.natoms-m.elems.count("x")
        self.mfp.set_FFfrag(name, lines, prio, comment)
        return
    
    def set_special_atype(self, at, ft, stype = "linear"):
        """
        Method to assign an attribute to an aftype
        :Parameters:
            - at (str): atype
            - ft (str): fragtype
            - stype (str,optional): attribute, defaults to linear
        """
        assert type(at) == type(ft) == type(stype) == str
        self.mfp.set_special_atype(at,ft,stype)
        return

    def set_orients(self, scid, path):
        """
        Method to upload an orients file
        :Parameters:
            - scid: id of the supercell entry
            - path: path to the orientsfile
        """
        with open(path, "r") as f:
            lines = f.read()
        self.mfp.set_orients(scid, lines)
        return

    def set_scaledtopo(self, scid, path):
        """
        Method to upload a scaled topo file
        :Parameters:
            - scid: id of the supercell entry
            - path: path ot the scaled topo file
        """
        with open(path, "r") as f:
            lines = f.read()
        self.mfp.set_scaledtopo(scid,lines)
        return

