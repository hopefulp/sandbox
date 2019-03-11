#!/usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import xmlrpclib
import os

locpath = os.environ["MFPLOC"]

# methods to retrieve and fill tables with the MOF-FF parameters

def get_params(db,FF, atypes, fragments, ptype, potential, FFfit):
    plen = len(atypes)
    if plen == 1:
        return get_onebody(db, FF, atypes, fragments, ptype, potential, FFfit)
    elif plen == 2:
        return get_twobody(db, FF, atypes, fragments, ptype, potential, FFfit)
    elif plen == 3:
        return get_threebody(db, FF, atypes, fragments, ptype, potential, FFfit)
    elif plen == 4:
        return get_fourbody(db, FF, atypes, fragments, ptype, potential, FFfit)
    else:
        return []
    
def get_params_from_ref(db, FF, ref):
    ffid = db(db.FF.name == FF).select()[0].id
    refid = db(db.FFrefs.name == ref).select()[0].id
    fitid = db((db.FFfits.refID1 == refid) & (db.FFfits.FFID == ffid)).select()[0].id
    params = []
    onebody = db(db.onebody.fitID==fitid).select()
    for r in onebody:
        atypes = []
        fragtypes = []
        atypes.append(db.atypes[r.atype1].name)
        fragtypes.append(db.FFfrags[r.frag1].name)
        params.append([atypes,fragtypes, r.type, r.pot, r.params])
    twobody = db(db.twobody.fitID==fitid).select()
    for r in twobody:
        atypes = []
        fragtypes = []
        atypes.append(db.atypes[r.atype1].name)
        fragtypes.append(db.FFfrags[r.frag1].name)
        atypes.append(db.atypes[r.atype2].name)
        fragtypes.append(db.FFfrags[r.frag2].name)
        params.append([atypes,fragtypes, r.type, r.pot, r.params])
    threebody = db(db.threebody.fitID==fitid).select()
    for r in threebody:
        atypes = []
        fragtypes = []
        atypes.append(db.atypes[r.atype1].name)
        fragtypes.append(db.FFfrags[r.frag1].name)
        atypes.append(db.atypes[r.atype2].name)
        fragtypes.append(db.FFfrags[r.frag2].name)
        atypes.append(db.atypes[r.atype3].name)
        fragtypes.append(db.FFfrags[r.frag3].name)
        params.append([atypes,fragtypes, r.type, r.pot, r.params])
    fourbody = db(db.fourbody.fitID==fitid).select()
    for r in fourbody:
        atypes = []
        fragtypes = []
        atypes.append(db.atypes[r.atype1].name)
        fragtypes.append(db.FFfrags[r.frag1].name)
        atypes.append(db.atypes[r.atype2].name)
        fragtypes.append(db.FFfrags[r.frag2].name)
        atypes.append(db.atypes[r.atype3].name)
        fragtypes.append(db.FFfrags[r.frag3].name)
        atypes.append(db.atypes[r.atype4].name)
        fragtypes.append(db.FFfrags[r.frag4].name)
        params.append([atypes,fragtypes, r.type, r.pot, r.params])
    return params

def get_header_ids(db,FF,atypes,fragments,fitsystem):
    ### first select FF
    try:
        ffid = db(db.FF.name == FF).select()[0].id
    except IndexError:
        return xmlrpclib.Fault("DBError", "FF %s not in DB" % FF)
    try:
        refid = db(db.FFrefs.name == fitsystem).select()[0].id
    except IndexError:
        return xmlrpclib.Fault("DBError", "referencesystem %s not in DB" % fitsystem)
    try:
        fitid = db((db.FFfits.refID1 == refid) & (db.FFfits.FFID == ffid)).select()[0].id
    except IndexError:
        return xmlrpclib.Fault("DBError", "FFfit for FF %s and referencesystem %s not in DB" % (FF,fitsystem))
    ### select atypes and fragments
    atids = []
    fragids = []
    for at, frag in zip(atypes,fragments):
        atrows = db(db.atypes.name == at).select()
        if len(atrows) == 0:
            return xmlrpclib.Fault("DBError", "atype %s not in DB" % at)
        else:
            if db((db.atypes2FFrefs.atypeID == atrows[0].id) & (db.atypes2FFrefs.refID == refid)).isempty():
                return xmlrpclib.Fault("DBError", "atype %s not in reference system %s" % (at,fitsystem)) 
            atids.append(atrows[0].id)
        fragrows = db(db.FFfrags.name == frag).select()
        if len(fragrows) == 0:
            return xmlrpclib.Fault("DBError", "fragment %s not in DB" % frag)
        else:
            if db((db.FFfrags2FFrefs.fragID == fragrows[0].id) & (db.FFfrags2FFrefs.refID == refid)).isempty():
                return xmlrpclib.Fault("DBError", "fragment %s not in reference system %s" % (frag,fitsystem))
            fragids.append(fragrows[0].id)
    return atids, fragids, fitid

def get_onebody(db, FF, atypes, fragments, ptype, potential, FFfit):
    header_ids = get_header_ids(db, FF, atypes, fragments, FFfit)
    if type(header_ids) == tuple:
        atids, fragids, fitid = header_ids
    else:
        return header_ids
    table = db.onebody
    prow = db((table.atype1 == atids[0]) & (table.frag1 == fragids[0])
            & (table.type == ptype) & (table.pot == potential) & (table.fitID==fitid)).select()
    if len(prow) == 0:
        return xmlrpclib.Fault("DBError", "Requested parameter set not available")
    return prow[0].params, prow[0].id

def get_twobody(db,FF, atypes, fragments, ptype, potential, FFfit):
    header_ids = get_header_ids(db, FF, atypes, fragments, FFfit)
    if type(header_ids) == tuple:
        atids, fragids, fitid = header_ids
    else:
        return header_ids
    table = db.twobody
    prow = db((table.atype1 == atids[0]) & (table.frag1 == fragids[0])
            & (table.atype2 == atids[1]) & (table.frag2 == fragids[1])
            & (table.type == ptype) & (table.pot == potential) & (table.fitID==fitid)).select()
    if len(prow) == 0:
        return xmlrpclib.Fault("DBError", "Requested parameter set not available")
    return prow[0].params, prow[0].id

def get_threebody(db,FF, atypes, fragments, ptype, potential, FFfit):
    header_ids = get_header_ids(db, FF, atypes, fragments, FFfit)
    if type(header_ids) == tuple:
        atids, fragids, fitid = header_ids
    else:
        return header_ids
    table = db.threebody
    prow = db((table.atype1 == atids[0]) & (table.frag1 == fragids[0])
            & (table.atype2 == atids[1]) & (table.frag2 == fragids[1])
            & (table.atype3 == atids[2]) & (table.frag3 == fragids[2])
            & (table.type == ptype) & (table.pot == potential) & (table.fitID==fitid)).select()
    if len(prow) == 0:
        return xmlrpclib.Fault("DBError", "Requested parameter set not available")
    return prow[0].params, prow[0].id

def get_fourbody(db,FF, atypes, fragments, ptype, potential, FFfit):
    header_ids = get_header_ids(db, FF, atypes, fragments, FFfit)
    if type(header_ids) == tuple:
        atids, fragids, fitid = header_ids
    else:
        return header_ids
    table = db.fourbody
    prow = db((table.atype1 == atids[0]) & (table.frag1 == fragids[0])
            & (table.atype2 == atids[1]) & (table.frag2 == fragids[1])
            & (table.atype3 == atids[2]) & (table.frag3 == fragids[2])
            & (table.atype4 == atids[3]) & (table.frag4 == fragids[3])
            & (table.type == ptype) & (table.pot == potential) & (table.fitID==fitid)).select()
    if len(prow) == 0:
        return xmlrpclib.Fault("DBError", "Requested parameter set not available")
    return prow[0].params, prow[0].id

def set_params(db, FF, atypes, fragments, ptype, potential, FFfit, params):
    plen = len(atypes)
    if plen == 1:
        return set_onebody(db, FF, atypes, fragments, ptype, potential,FFfit,params)
    elif plen == 2:
        return set_twobody(db, FF, atypes, fragments, ptype, potential,FFfit,params)
    elif plen == 3:
        return set_threebody(db, FF, atypes, fragments, ptype, potential,FFfit,params)
    elif plen == 4:
        return set_fourbody(db, FF, atypes, fragments, ptype, potential,FFfit,params)

    
# user und timestamp fehlen
def set_onebody(db, FF, atypes, fragments, ptype, potential,FFfit, params):
    header_ids = get_header_ids(db, FF, atypes, fragments, FFfit)
    if type(header_ids) == tuple:
        atids, fragids, fitid = header_ids
    else:
        return header_ids
    check = get_onebody(db, FF, atypes, fragments, ptype, potential, FFfit)
    if type(check) == tuple: return xmlrpclib.Fault("DBError", "Parameter set already in DB")
    row = db.onebody.insert(creationtime = datetime.datetime.now(),
            creator = current.session.auth.user.id,
            atype1 = atids[0],
            frag1 = fragids[0],
            type = ptype,
            pot = potential,
            fitID = fitid,
            params = params)
    return row.id

def set_twobody(db, FF, atypes, fragments, ptype, potential,FFfit, params):
    header_ids = get_header_ids(db, FF, atypes, fragments,FFfit)
    if type(header_ids) == tuple:
        atids, fragids, fitid = header_ids
    else:
        return header_ids
    check = get_twobody(db, FF, atypes, fragments, ptype, potential, FFfit)
    if type(check) == tuple: return xmlrpclib.Fault("DBError", "Parameter set already in DB")
    row = db.twobody.insert(creationtime = datetime.datetime.now(),
            creator = current.session.auth.user.id,
            atype1 = atids[0],
            frag1 = fragids[0],
            atype2 = atids[1],
            frag2 = fragids[1],
            type = ptype,
            pot = potential,
            fitID = fitid,
            params = params)
    return row.id

def set_threebody(db, FF, atypes, fragments, ptype, potential,FFfit, params):
    header_ids = get_header_ids(db, FF, atypes, fragments,FFfit)
    if type(header_ids) == tuple:
        atids, fragids, fitid = header_ids
    else:
        return header_ids
    check = get_threebody(db, FF, atypes, fragments, ptype, potential, FFfit)
    if type(check) == tuple: return xmlrpclib.Fault("DBError", "Parameter set already in DB")
    row = db.threebody.insert(creationtime = datetime.datetime.now(),
            creator = current.session.auth.user.id,
            atype1 = atids[0],
            frag1 = fragids[0],
            atype2 = atids[1],
            frag2 = fragids[1],
            atype3 = atids[2],
            frag3 = fragids[2],
            type = ptype,
            fitID = fitid,
            pot = potential,
            params = params)
    return row.id

def set_fourbody(db, FF, atypes, fragments, ptype, potential,FFfit, params):
    header_ids = get_header_ids(db, FF, atypes, fragments,FFfit)
    if type(header_ids) == tuple:
        atids, fragids, fitid = header_ids
    else:
        return header_ids
    check = get_fourbody(db, FF, atypes, fragments, ptype, potential, FFfit)
    if type(check) == tuple: return xmlrpclib.Fault("DBError", "Parameter set already in DB")
    row = db.fourbody.insert(creationtime = datetime.datetime.now(),
            creator = current.session.auth.user.id,
            atype1 = atids[0],
            frag1 = fragids[0],
            atype2 = atids[1],
            frag2 = fragids[1],
            atype3 = atids[2],
            frag3 = fragids[2],
            atype4 = atids[3],
            frag4 = fragids[3],
            type = ptype,
            pot = potential,
            fitID = fitid,
            params = params)
    return row.id

### methods to save and retrieve hdf5 file with reference informations

def get_FFref(db,name):
    row = db(db.FFrefs.name == name).select()[0]
    path = locpath + '/FFs/refs/'+row.reffile
    with open(path, "rb") as handle:
         return xmlrpclib.Binary(handle.read())
        
def get_FFref_graph(db,name):
    row = db(db.FFrefs.name == name).select()[0]
    path = locpath + '/FFs/refs/'+row.reffile
    with open(path, "r") as handle:
        return handle.read()
    
def set_FFref_graph(db,name,mfpxfile):
    from cStringIO import StringIO
    import atypedbhelper as adb
    row = db(db.FFrefs.name == name).select()[0]
    f2 = StringIO(mfpxfile)
    row.update_record(graph = db.FFrefs.graph.store(f2, "%s.mfpx" % name))
    adb.insert_atypes_of_FFref(db, row.id)
    return
        
def set_FFref(db,name,reffile,mfpxfile,comment):
    from cStringIO import StringIO
    import atypedbhelper as adb
    f = StringIO(reffile.data)
    f2 = StringIO(mfpxfile)
    refid = db.FFrefs.insert(name = name,
                    creator = current.session.auth.user.id,
                    uploadtime = datetime.datetime.now(),
                    comment = comment,
                    graph = db.FFrefs.graph.store(f2, "%s.mfpx" % name),
                    reffile = db.FFrefs.reffile.store(f, "%s.hdf5" % name))
    adb.insert_atypes_of_FFref(db, refid)
    return

def list_FFrefs(db,FF):
    try:
        ffid = db(db.FF.name == FF).select()[0].id
    except IndexError:
        return xmlrpclib.Fault("DBError", "FF %s not in DB" % FF)
    fitrows = db(db.FFfits.FFID == ffid).select()
    names = []
    for fr in fitrows:
        ftypes = [db.FFfrags[rf.fragID].name for rf in db(db.FFfrags2FFrefs.refID ==fr.refID1).select()]
        names.append((db.FFrefs[fr.refID1].name,db.FFrefs[fr.refID1].priority, ftypes))
    #rows = db(db.FFrefs.id > 0).select()
    #names = []
    #for r in rows: 
    #    ftypes = [db.FFfrags[rf.fragID].name for rf in db(db.FFfrags2FFrefs.refID ==r.id).select()]
    #    names.append((r.name,r.priority, ftypes ))
    return names

def get_FFfrag(db, name):
    row = db(db.FFfrags.name == name).select()[0]
    path = locpath + '/FFs/frags/'+row.file
    with open(path, "r") as handle:
         return handle.read()
        
def set_FFfrag(db,name, fragfile, prio, comment):
    from cStringIO import StringIO
    import atypedbhelper as adb
    f = StringIO(fragfile)
    fragid = db.FFfrags.insert(name = name,
                     creator = current.session.auth.user.id,
                     creationtime = datetime.datetime.now(),
                     comment = comment,
                     priority = prio,
                     file = db.FFfrags.file.store(f, "%s.mfpx" % name))
    adb.insert_atypes_of_fragment(db,fragid)
    return

def frags2atypes(db,row):
    import molsys
    m = molsys.mol()
    m.read()

def list_FFfrags(db):
    import json
    rows = db(db.FFfrags.id > 0).select()
    names = []
    for r in rows: 
        #[db(db.fireanalyzer.structureID==i.id).select(db.fireanalyzer.ALL) for i in structures]
        atypes = [db.atypes[rf.atypeID].name for rf in db(db.atypes2FFfrags.fragID ==r.id).select()]
        names.append((r.name, r.priority, atypes))
    return names

def create_empty_FFfit(db, refrow):
    ffid = db(db.FF.name == "MOF-FF").select()[0].id
    db.FFfits.insert(creator = 1,
                    creationtime = datetime.datetime.now(),
                    comment = "old fit, no input and output available",
                    refID1 = refrow.id,
                    FFID = ffid)

    
#db.define_table('FFfits',
#    Field('name', 'string', unique = True),
#    Field('creator', 'reference auth_user'),
#    Field('creationtime', 'datetime'),
#    Field('settings', 'json'),
#    Field('input', 'upload', uploadfolder=request.folder+'/static/FFs/inp'),
#    Field('output', 'upload', uploadfolder=request.folder+'/static/FFs/out'),
#    Field('refID1', 'reference FFrefs'),
#    Field('comment', 'text', default=""))
    #'./applications/MFP_JK/temp/netfiles'
    #path = 
    #with open(path, "wb") as handle:
    #    handle.write(reffile.data)

#db.define_table('FFrefs',
#    Field('name', 'string', unique = True),
#    Field('creator','reference auth_user'),
#    Field('uploadtime', 'datetime'),
#    Field('reffile','upload', uploadfolder=request.folder+'/static/FFs/refs'),
#    Field('comment', 'text', default=""),
#    format='%(name)s')
