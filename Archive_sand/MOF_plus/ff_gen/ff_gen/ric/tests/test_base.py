#! /usb/bin/env python

import unittest as ut
import numpy as np

from ric import Group
from ric.base import _Base

class BaseTests(ut.TestCase):

  def setUp(self):

    self.grp1 = Group([1,2,3])
    self.grp2 = Group([3,4])

    self.base = _Base([self.grp1,5,self.grp2,2])

  def test_points(self):
    self.assertListEqual(self.base.points, [self.grp1,5,self.grp2,2])

  def test_indices(self):
    self.assertListEqual(list(self.base.indices), range(5))

  def test_bmat_indices(self):
    self.assertListEqual(list(self.base.bmat_indices), range(15))

  def test_get_points(self):
    for _ in range(10):
      coords = 20*np.random.random((5,3)) - 10
      ref = [np.mean(coords[self.grp1.indices], axis=0),
             coords[4],
             np.mean(coords[self.grp2.indices], axis=0),
             coords[1]]
      res = self.base.get_points(coords)
      self.assertTrue(np.allclose(ref, res))

  def test_get_vectors(self):
    for _ in range(10):
      coords = 20*np.random.random((5,3)) - 10
      pnts = self.base.get_points(coords)
      ref = [pnts[0] - pnts[1],
             pnts[2] - pnts[0],
             pnts[3] - pnts[1]]
      res = self.base.get_vectors([[1,0],[0,2],[1,3]], coords)
      self.assertTrue(np.allclose(ref, res))

  def test_project_points(self):
    coords = np.zeros((5,3), dtype=np.float64) # Dummy coords
    for _ in range(1):
      vecs = 20*np.random.random((4,3)) - 10
      ref = np.zeros((5,3), dtype=np.float64)
      ref[(0,1,2),:] += vecs[0,:]/3
      ref[ 4     ,:] += vecs[1,:]
      ref[(2,3)  ,:] += vecs[2,:]/2
      ref[ 1     ,:] += vecs[3,:]
      ref = ref.ravel()
      res = self.base.project_points(vecs, coords)
      self.assertTrue(np.allclose(ref, res))


if __name__ == '__main__':
  ut.main()

