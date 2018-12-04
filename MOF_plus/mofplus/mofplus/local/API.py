# -*- coding: utf-8 -*-
import os
from db import db

locpath = os.environ["MFPLOC"]

def search_cs(cs, vs = None):
    import netsearch as ns
    r = ns.search_cs_and_vs(db, cs = cs, vs = vs)
    res = []
    for i in r: res.append(i.name)
    return res

def get_cs(name):
    import netdbhelper as ndb
    return ndb.get_cs(db,name)
#    r = db(db.nets.name== name).select()[0]
#    return [r.cs1, r.cs2, r.cs3, r.cs4, r.cs5, r.cs6, r.cs7, r.cs8, r.cs9, r.cs10]

def get_vs(name):
    import netdbhelper as ndb
    return ndb.get_vs(db,name)

def get_net(name):
    row = db(db.nets.name== name).select()[0]
    f = open(locpath +'/nets/files/'+row.topofile, 'r')
    string = f.read()
    f.close()
    return string

def get_list_of_nets():
    avail = db(db.nets).select(db.nets.name)
    return [i.name for i in avail]

def get_list_of_bbs():
    avail = db(db.bbs).select(db.bbs.name)
    return [i.name for i in avail]

def get_bb(name):
    row = db(db.bbs.name==name).select()[0]
    f = open(locpath +'bbs/files/'+row.xyzfile, 'r')
    string = f.read()
    f.close()
    return string

def get_mof_structure_by_id(strucid):
    row = db(db.structures.id == strucid).select()[0]
    f = open(locpath+'mofs/structures/'+row.filename, 'r')
    string = f.read()
    f.close()
    return string, row.name

def get_params(FF, atypes, fragments, ptype, potential, FFfit):
    import FFhelper
    return FFhelper.get_params(db,FF, atypes, fragments, ptype, potential, FFfit)

def get_params_from_ref(FF, ref):
    import FFhelper
    return FFhelper.get_params_from_ref(db, FF, ref)

def get_FFref(name):
    import FFhelper
    return FFhelper.get_FFref(db,name)

def get_FFref_graph(name):
    import FFhelper
    return FFhelper.get_FFref_graph(db,name)

def list_FFrefs(FF):
    import FFhelper
    return FFhelper.list_FFrefs(db,FF)

def get_FFfrag(name):
    import FFhelper
    return FFhelper.get_FFfrag(db,name)

def list_FFfrags():
    import FFhelper
    return FFhelper.list_FFfrags(db)
