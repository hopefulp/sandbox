"""

Timer Module stolen from RSDFT5 and modified to fit with pydlpoly

"""

import time
import os
from mpi4py import MPI
pid = os.getpid()

is_master = MPI.COMM_WORLD.Get_rank()

class timer:
    '''
    Implements a Timer class that can be switched on and off and which can
    have an arbitrary number of labels to timer
    
    just call the report method to wirte out the timing info
    '''

    def __init__(self, name):
        ''' make a timer with name
        '''
        self.name = name
        self.on = 0
        self.current = None
        self.last_time = 0.0
        self.timers = {}
        return
        
    def switch_on(self, t):
        if self.on: raise "timer is already on"
        if not self.timers.has_key(t):
            self.timers[t] = (0.0, 0)
        self.current = t
        self.on   = 1
        self.last_time = time.time()
        return
        
    def switch_off(self):
        _t = time.time()-self.last_time
        if not self.on: raise "timer is not on"
        _ot = self.timers[self.current]
        self.timers[self.current] = (_ot[0]+_t, _ot[1]+1)
        self.on = 0
        return
        
    def switch_to(self, t):
        _t = time.time()-self.last_time
        if not self.on: raise "timer is not on"
        _ot = self.timers[self.current]
        self.timers[self.current] = (_ot[0]+_t, _ot[1]+1)
        if not self.timers.has_key(t):
            self.timers[t] = (0.0, 0)
        self.current = t
        self.last_time = time.time()
        return

    def report(self, ident=""):
        if is_master == 0:
            _st = 0.0
            for t in self.timers.keys():
                _st = _st + self.timers[t][0]
            print("\n%-5stimer %-40s total time %10.3f sec [wall]" % (ident, self.name, _st))
            print("%-5s%-50s %10s %12s %8s" % (ident, "name", "total", "per call", "# calls"))
            _tl = self.timers.keys()
            _tl.sort()
            for t in _tl:
                _t = self.timers[t]
                if _t[1] == 0:
                    _percall= 0.0
                else:
                    _percall = _t[0]/_t[1]
                print("%-5s%-50s %10.3f %12.5f %8d" % (ident, t, _t[0], _percall, _t[1]))
            print ("\n")
        return

    def zero(self):
        if self.on: raise "cant zero when timer is on"
        for t in self.timers.keys():
            self.timers[t] = (0.0, 0)
        return

    def reset(self):
        if self.on: raise "cant reset when timer is on"
        self.current = None
        self.last_time = 0.0
        self.timers = {}
        return
            
        
