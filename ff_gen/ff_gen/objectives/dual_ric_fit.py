import os
import sys
import numpy as np
import ric_fit
import force_ric_fit

class dual_ric_fit:

    def __init__(self, wgt_ric_fit = 0.5):
        assert wgt_ric_fit <= 1.0
        assert wgt_ric_fit >= 0.0
        self.wgt_ric_fit = wgt_ric_fit
        self.wgt_force_ric_fit = 1 - self.wgt_ric_fit 
        return

    def initialize(self, pd, ref):
        self.ric_fit = ric_fit.ric_fit()
        self.ric_fit.initialize(pd, ref)
        self.force_ric_fit = force_ric_fit.force_ric_fit()
        self.force_ric_fit.initialize(pd, ref)
        self.cycle = 0
        self.fdiagnostics = open('dual_ric_fit.punch', 'w')
        return

    def __call__(self):
        msd_ric = self.ric_fit()
        msd_force = self.force_ric_fit()
        self.msd = (self.wgt_ric_fit * msd_ric) + (self.wgt_force_ric_fit * msd_force)
        self.fdiagnostics.write("%s %6.6f\n" % (self.cycle, self.msd))
        self.cycle += 1
        return self.msd


