#! /usb/bin/env python

import unittest as ut
import numpy as np
import numdifftools as nd

from ric import Group, Dihedral

class DihedralTests(ut.TestCase):

  def setUp(self):

    self.grp1 = Group([1])
    self.grp2 = Group([2,3])
    self.grp3 = Group([3,4,5])
    self.grp4 = Group([1,5])

    self.dihed1 = Dihedral([1], [2,3], [4])
    self.dihed2 = Dihedral([1,2,3,4,5], [6,7], [8,9,10,11])
    self.dihed3 = Dihedral([self.grp1], [self.grp2,self.grp3], [self.grp4])
    self.dihed4 = Dihedral([self.grp1], [2,self.grp3], [6])

  def test_evaluate(self):

    for coords, ref in [([[ 0,-1, 0],[0,0,0],[0,0,1],[1,0,1]], -90),
                        ([[ 1,-1, 0],[0,0,0],[0,0,1],[1,0,1]], -45),
                        ([[ 1, 0, 0],[0,0,0],[0,0,1],[1,0,1]],   0),
                        ([[ 1, 0,-1],[0,0,0],[0,0,1],[1,0,2]],   0),
                        ([[ 1, 1, 0],[0,0,0],[0,0,1],[1,0,1]],  45),
                        ([[ 2, 2, 0],[0,0,0],[0,0,1],[3,0,1]],  45),
                        ([[ 0, 1, 0],[0,0,0],[0,0,1],[1,0,1]],  90),
                        ([[-1, 0, 0],[0,0,0],[0,0,1],[1,0,1]], 180)]:
      coords = np.array(coords, dtype=np.float64)
      res = np.degrees(self.dihed1.evaluate(coords))
      self.assertAlmostEqual(ref, res)

  def test_evaluate_ivals(self):

    coords = np.array([[ 1., 0., 0.],  # 1:    0
                       [ 1., 1.,-1.],  # 2:   45
                       [ 0., 1.,-2.],  # 3:   90
                       [-1., 1.,-3.],  # 4:  135
                       [-1., 0.,-4.],  # 5:  180
                       [ 0., 0., 0.],  # central
                       [ 0., 0., 1.],  # central
                       [ 1., 0., 1.],  # 1:    0
                       [ 1.,-1., 2.],  # 2:  -45
                       [ 0.,-1., 3.],  # 3:  -90
                       [-1.,-1., 4.]]) # 4: -135

    for ivals, ref in [([1,1],   0),([1,2],  45),([1,3],  90),([1,4], 135),
                       ([2,1],  45),([2,2],  90),([2,3], 135),([2,4], 180),
                       ([3,1],  90),([3,2], 135),([3,3], 180),([3,4],-135),
                       ([4,1], 135),([4,2], 180),([4,3],-135),([4,4], -90),
                       ([5,1], 180),([5,2],-135),([5,3], -90),([5,4], -45)]:
      res = np.degrees(self.dihed2.evaluate(coords, ivals=ivals))
      self.assertAlmostEqual(ref, res, delta=1e-5)

  def test_project(self):

    dihed = self.dihed1
    for _ in range(10):
      coords = 20*np.random.random((4,3)) - 10
      grad = nd.Gradient(lambda x: dihed.evaluate(x.reshape((-1,3))))
      ref = grad(coords.ravel())
      res = dihed.project(coords)
      self.assertTrue(np.allclose(ref, res))

    dihed = self.dihed3
    for _ in range(10):
      coords = 20*np.random.random((5,3)) - 10
      grad = nd.Gradient(lambda x: dihed.evaluate(x.reshape((-1,3))))
      ref = grad(coords.ravel())
      res = dihed.project(coords)
      self.assertTrue(np.allclose(ref, res))

    dihed = self.dihed4
    for _ in range(10):
      coords = 20*np.random.random((6,3)) - 10
      grad = nd.Gradient(lambda x: dihed.evaluate(x.reshape((-1,3))))
      ref = grad(coords.ravel())
      res = dihed.project(coords)
      self.assertTrue(np.allclose(ref, res))

if __name__ == '__main__':
  ut.main()

