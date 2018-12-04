#!/usr/bin/env python
# -*- coding: utf-8 -*-

import xmlrpclib
import molsys
import logging

logger = logging.getLogger("mofplus")

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
                            return molsys.topo.fromString(lines)
                        else:
                            return molsys.mol.fromString(lines)
                if binary == False:
                    if dtype == "orients":
                        f=open(str(args[1])+'.orients', 'w')
                    else:
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


def nolocal(func):
    def inner(self,*args,**kwargs):
        if self.local == True:
            logger.error("Method %s not allowed for local DB usage" % func.__name__)
            return
        else:
            return func(*args, **kwargs)
    return inner


