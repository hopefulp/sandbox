#! /usb/bin/env python

import unittest as ut
import numpy as np

from ric import Group

class GroupTests(ut.TestCase):

  def setUp(self):
    self.grp = Group([1,3,2,5])

  def test_len(self):
    self.assertEqual(len(self.grp), 4)

  def test_indices(self):
    self.assertListEqual(list(self.grp.indices), [0,1,2,4])

  def test_bmat_indices(self):
    self.assertListEqual(list(self.grp.bmat_indices), [0,1,2,3,4,5,6,7,8,12,13,14])

  def test_evaluate(self):
    coords = np.array([[  1,  2,  3], # 1
                       [  4,  5,  6], # 2
                       [  7,  8,  9], # 3
                       [ 10, 11, 12], # 4 -- unused
                       [ 13, 14, 15], # 5
                       [ 16, 17, 18], # 6 -- unused
                      ], dtype=np.float64)
    ref = np.mean(coords[self.grp.indices], axis=0)
    res = self.grp.evaluate(coords)
    self.assertTrue(np.allclose(ref, res))

  def test_project(self):
    vec = np.array([4, 0, -2], dtype=np.float64)
    ref = [1, 0, -.5]*4
    res = self.grp.project(vec)
    self.assertTrue(np.allclose(ref, res))


if __name__ == '__main__':
  ut.main()

