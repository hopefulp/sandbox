#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import shutil
import unittest
import numpy as np

from ff_gen.main import main
from ff_gen.optimizer import cma

class PeriodicTest(unittest.TestCase):

#    def __init__(self):
#        super(PeriodicTest, self).__init__()

    def setUp(self):
        self.testdir = os.getcwd()
        self.ff = main('ZIF-FF', use_range = True, out = 'blurp.out')
        self.ff.init_objective('ric_fit', 'ZIF8', 'primary', offline = True, fpar='ZIF-FF',
                refsys="ZIF8", objkwargs={'wstress':1.0,"wtor":0.0, 'catchneg':False})

    def TestSingleFitness(self):
        self.ff.setup()
        fitness = self.ff.__call__(self.ff.initial_guess())
        np.testing.assert_allclose(fitness, 1.01949905257e-06)
        self.ff.use_range=False
        fitness = self.ff.__call__(self.ff.initial_guess())
        np.testing.assert_allclose(fitness, 1.01949905257e-06)

    def TestSingleOpt(self):
        self.ff.setup()
        es = cma(self.ff, maxiter = 2, initials = self.ff.initial_guess(), sigma = 0.4, 
                cmaopts={'popsize':16,'seed':42})
        es()
        fitness = self.ff.__call__(es.result)
        np.testing.assert_allclose(fitness, 1.943837267321726e-02, atol = 5e-7)
        return

    def TestMultiFitness(self):
        self.ff.init_objective('ric_fit', 'ZIF6', 'primary', offline = True, fpar='ZIF-FF',
                refsys="ZIF8", objkwargs={'wstress':1.0,"wtor":1.0, 'catchneg':False})
        self.ff.setup()
        fitness = self.ff.__call__(self.ff.initial_guess())
        print fitness
        np.testing.assert_allclose(fitness, 0.0122700830156)
        self.ff.use_range=False
        fitness = self.ff.__call__(self.ff.initial_guess())
        np.testing.assert_allclose(fitness, 0.0122700830156)

    def tearDown(self):
        os.chdir(self.testdir)
        shutil.rmtree('ZIF-FF')
        for i in ['blurp.out', 'molsys.log', 'FFgen.log']: 
            if os.path.isfile(i): os.remove(i)
        return


if __name__=="__main__":
    suite = unittest.TestSuite()
    suite.addTest(PeriodicTest('TestSingleFitness'))
    suite.addTest(PeriodicTest('TestMultiFitness'))
#    suite.addTest(PeriodicTest('TestSingleOpt'))
    unittest.TextTestRunner(verbosity=3).run(suite)
