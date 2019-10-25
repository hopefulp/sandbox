#!/home/noische/Enthought/Canopy_64bit/User/bin/python

#import lammpstrj as lt
from lammps.trj import *
import memory_footprint

def myhandler(obj):
    assert isinstance(obj, Trj)
    yield obj.trj_file
    yield obj.data_file
    yield obj.timesteps
    yield obj._dump_style
    yield obj._dump_keywords
    yield obj._dump_keywords_r
    yield obj.natoms
    yield obj.nheader
    yield obj.coord
    yield obj.xlo
    yield obj.xhi
    yield obj.ylo
    yield obj.yhi
    yield obj.zlo
    yield obj.zhi
    yield obj.pbc
    yield obj._is_loaded
    yield obj._is_dumped

mytrj = Trj('temp')
mytrj.load(force=True)
mytrj.dump()
print memory_footprint.getsizeof(mytrj.coord)
print memory_footprint.getsizeof(mytrj)
print memory_footprint.getsizeof(mytrj, {Trj: myhandler})


