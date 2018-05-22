import os
import sys
import numpy as np

class combine_fit:

    def __init__(self, objectives, weights = None):
        self.objectives = objectives
        if weights == None:
            self.weights = np.zeros((len(objectives)), dtype = 'float64')
            self.weights[:] = 1.0/len(objectives)
        else:
            assert len(self.objectives) == len(weights)
            self.weights = np.array(weights)
        return

    def initialize(self, pd, ref):
        for i in range(len(self.objectives)):
            self.objectives[i].initialize(pd, ref)
        self.msds = np.zeros((len(self.objectives)), dtype = 'float64')
        self.cycle = 0
        self.fdiagnostics = open('combine_fit.punch', 'w')
        return

    def __call__(self):
        self.a_msds = []
        for i in range(len(self.objectives)):
            self.msds[i], a_msd = self.objectives[i].__call__()
            #self.a_msds.append(a_msd)
            self.a_msds += a_msd
        self.msd = np.sum(self.msds * self.weights)
        self.fdiagnostics.write('%s %6.6f\n' % (self.cycle, self.msd))
        self.cycle += 1
        return self.msd, self.a_msds
