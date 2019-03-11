#! /usb/bin/env python

import unittest as ut
import numpy as np
import numdifftools as nd

from ric import Group, Distance

class DistanceTests(ut.TestCase):

  def setUp(self):

    self.grp1 = Group([1,2,4])
    self.dist = Distance([self.grp1,3])

  def test_evaluate(self):

    for _ in range(10):
      coords = 20*np.random.random((4,3)) - 10
      ref = coords[(0,1,3),:].mean(axis=0) - coords[2,:]
      ref = np.sqrt(np.sum(ref**2))
      res = self.dist.evaluate(coords)
      self.assertTrue(np.allclose(ref, res))

  def test_project(self):

    for _ in range(10):
      coords = 20*np.random.random((4,3)) - 10
      grad = nd.Gradient(lambda x: self.dist.evaluate(x.reshape((-1,3))))
      ref = grad(coords.ravel())
      res = self.dist.project(coords)
      self.assertTrue(np.allclose(ref, res))


if __name__ == '__main__':
  ut.main()

