#!/usr/bin/env python
# -*- coding: utf-8 -*-

def search_cs_and_vs(db,cs,vs=None):
    ncs = len(cs)
    nets = db(db.nets.p == ncs).select(db.nets.id) 
    for i, c in enumerate(cs):
        q = (db.vertex.cs10==c[9]) &\
        (db.vertex.cs9==c[8]) &\
        (db.vertex.cs8==c[7]) &\
        (db.vertex.cs7==c[6]) &\
        (db.vertex.cs6==c[5]) &\
        (db.vertex.cs5==c[4]) &\
        (db.vertex.cs4==c[3]) &\
        (db.vertex.cs3==c[2]) &\
        (db.vertex.cs2==c[1]) &\
        (db.vertex.cs1==c[0])
        if vs != None:
            q = q & (db.vertex.vs == vs[i])
        q = q & (db.vertex.netID == db.nets.id) & db.nets.id.belongs(nets)
        #if i != 0:
        #    q = q & (db.vertex.netID == db.nets.id) & db.nets.id.belongs(nets)
        #else:
        #    q = q & (db.vertex.netID == db.nets.id) # & (db.nets.id)
        nets = db(q).select(db.nets.id,distinct=True)
    nets = db(db.nets.id.belongs(nets)).select(db.nets.ALL)
    return nets

