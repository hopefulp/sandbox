#! /usb/bin/env python

import numpy        as np
import numdifftools as nd
import unittest     as ut

from ric import RedIntCoords as RIC

class InPlaneBendTests(ut.TestCase):

  def test_val_0D_1(self):
    """Values without PBC"""

    masses = np.array([1.]*3)

    ric = RIC()
    ric.add_in_bend([1,2,3])
    ric.setup(masses)

    for coords, ref in [([[ 1, 0, 0],[ 0, 0, 0],[ 1, 0, 0]],  0),
                        ([[ 1, 0, 0],[ 0, 0, 0],[ 1, 1, 0]], 45),
                        ([[ 1, 0, 0],[ 0, 0, 0],[ 0, 1, 0]], 90),
                        ([[ 1, 0, 0],[ 0, 0, 0],[-1, 1, 0]],135),
                        ([[ 1, 0, 0],[ 0, 0, 0],[-1, 0, 0]],180),
                        ([[ 1, 1, 1],[ 0, 0, 0],[-1,-1,-1]],180)]:

      ric.construct_b_matrix(None,np.array(coords,dtype=np.float64))
      res = ric.get_val_in_bends()[0]
      self.assertAlmostEqual(res,np.radians(ref))

  def test_grad_0D_1(self):

    masses = np.array([1.]*3)

    ric = RIC()
    ric.add_in_bend([1,2,3])
    ric.setup(masses)

    def func(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_in_bends())
      return res
    grad = nd.Gradient(func,step_nom=[0.01]*9)

    for coords in [[[ 1, 0, 0],[ 0, 0, 0],[ 1, 1, 0]],
                   [[ 1, 0, 0],[ 0, 0, 0],[ 0, 1, 0]],
                   [[ 1, 0, 0],[ 0, 0, 0],[-1, 1, 0]],
                   [[ 0, 1, 2],[ 5, 4, 3],[ 6, 7, 8]]]:
      coords = np.array(coords,dtype=np.float64)
      ref = grad(coords.flatten())
      res = ric.construct_b_matrix(None,coords)[0,:]
      self.assertLess(np.max(np.abs(ref-res)),1.e-10)

  def test_grad_0D_rand(self):

    masses = np.array([1.]*3)

    ric = RIC()
    ric.add_in_bend([1,2,3])
    ric.setup(masses)

    def func(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_in_bends())
      return res
    grad = nd.Gradient(func,step_nom=[0.01]*9)

    for coords in 2*np.random.random((10,3,3))-1:
      coords = np.array(coords,dtype=np.float64)
      ref = grad(coords.flatten())
      res = ric.construct_b_matrix(None,coords)[0,:]
      self.assertLess(np.max(np.abs(ref-res)),1.e-10)


if __name__ == '__main__':
  ut.main()

