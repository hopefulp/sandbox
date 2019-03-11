#!/usr/bin/env python
# -*- coding: utf-8 -*-

def get_net_info(db,netname):
    query=(((db.nets.name == netname)&(db.vertex.netID==db.nets.id)))
    vs = db(query).select(db.vertex.ALL)
    conn = sorted([int(i.coordination_number) for i in vs])
    connshort = sorted(list(set(conn)))
    q = (db.nets.name == netname) & (db.vertex.netID==db.nets.id) & ((db.edge.fromID==db.vertex.id) | (db.edge.toID==db.vertex.id))
    es = db(q).select(db.edge.ALL,distinct=True)
    return vs, conn, connshort, es

def delete_net(db, netname):
    #select all related entries
    net = db(db.nets.name==netname).select(db.nets.ALL)[0]
    mof = db(db.mofs.net==net.id).select()
    if len(mof)>0: return 'Not possible to delete! MOF is connected to this net!'
    vertices = db(db.vertex.netID == net.id).select(db.vertex.ALL)
    fedges = [db(db.edge.fromID==i.id).select(db.edge.ALL) for i in vertices]
    tedges = [db(db.edge.toID==i.id).select(db.edge.ALL) for i in vertices]
    fconn = [db(db.conn.fromID==i.id).select(db.conn.ALL) for i in vertices]
    tconn = [db(db.conn.toID==i.id).select(db.conn.ALL) for i in vertices]
    shapes = [db(db.shapes.vertexID==i.id).select(db.shapes.ALL) for i in vertices]
    rels = db((db.net_relations.cID == net.id) | (db.net_relations.pID == net.id)).select(db.net_relations.ALL)
    # then delete
    for i in fedges+tedges+fconn+tconn+shapes:
        for j in i:
            j.delete_record()
    for i in list(vertices) + list(rels):
        i.delete_record()
    net.delete_record()
    db.commit()
    return

def get_cs(db, netname):
    netid = db(db.nets.name == netname).select(db.nets.id)[0]
    verts = db(db.vertex.netID == netid).select()
    cs = []
    for i,v in enumerate(verts):
        assert v.idx == i
        cs.append([v.cs1,v.cs2,v.cs3,v.cs4,v.cs5,v.cs6,v.cs7,v.cs8,v.cs9,v.cs10])
    return cs

def get_vs(db, netname):
    netid = db(db.nets.name == netname).select(db.nets.id)[0]
    verts = db(db.vertex.netID == netid).select()
    vs = []
    for i,v in enumerate(verts):
        assert v.idx == i
        vs.append(v.vs)
    return vs

def create_relation(db, pname, cname, pattern):
    assert type(pname) == str
    assert type(cname) == str
    assert type(pattern) == str
    pid = db(db.nets.name == pname).select()[0].id
    cid = db(db.nets.name == cname).select()[0].id
    db.net_relations.insert(pID = pid, cID = cid, pattern = pattern)
    db.commit()
    return

def get_childs(db, pname):
    assert type(pname) == str
    pid = db(db.nets.name == pname).select()[0].id
    cids = db(db.net_relations.pID==pid).select(db.net_relations.cID)
    cnames = [db(db.nets.id == i.cID).select()[0].name for i in cids]
    return cnames

def get_parent(db, cname):
    assert type(cname) == str
    cid = db(db.nets.name == cname).select()[0].id
    try:
        rel = db(db.net_relations.cID == cid).select()[0]
        pid = rel.pID
        pattern = rel.pattern
        pname = db(db.nets.id == pid).select()[0].name
        return pname, pattern
    except IndexError:
        return None, None
    
def get_siblings(db, cname, pattern):
    assert type(cname) == str
    assert type(pattern) == str
    cid = db(db.nets.name == cname).select()[0].id
    rel = db(db.net_relations.cID == cid).select()
    if len(rel) < 1: return []
    pid = rel[0].pID
    pattern = rel[0].pattern
    sids = db((db.net_relations.pID == pid) & (db.net_relations.pattern == pattern)).select()
    snames = []
    for i in sids:
        sname = db(db.nets.id == i.cID ).select()[0].name
        if sname != cname: snames.append(sname)
    return snames




    
#def add_vs(db,n,vs):
#    verts = db(db.vertex.netID ==n.id).select()
#    assert len(verts) == len(vs)
#    for i,v in enumerate(verts):
#        v.vs = vs[i]
#        v.update_record()
#    db.commit()
