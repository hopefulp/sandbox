#! /usb/bin/env python

import numpy        as np
import numdifftools as nd
import unittest     as ut

from ric import RedIntCoords as RIC

class OutBendTests(ut.TestCase):

  def test_values_1(self):

    masses = np.array([1.]*4)
    hmat   = np.array([[0.,0.,0.]]*3)

    for coords, ref in [
      ([[0,0,0],[0,0, 1],[-1,-1,0],[1,-1,0]],   0),
      ([[0,0,0],[0,1, 1],[-1,-1,0],[1,-1,0]],  45),
      ([[0,0,0],[0,1, 0],[-1,-1,0],[1,-1,0]],  90),
      ([[0,0,0],[0,1,-1],[-1,-1,0],[1,-1,0]], 135),
      ([[0,0,0],[0,0,-1],[-1,-1,0],[1,-1,0]], 180),
                       ]:
      ric = RIC()
      ric.add_out_bend([1,2,3,4])
      ric.setup(masses)
      ric.construct_b_matrix(hmat,np.array(coords,dtype=np.float64))
      res = ric.get_val_out_bends()[0]
      self.assertAlmostEqual(res,np.radians(ref))
      del ric

  def test_grad_0D(self):

    masses = np.array([1.]*4)

    ric = RIC()
    ric.add_out_bend([1,2,3,4])
    ric.setup(masses)

    def func(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_out_bends())
      return res
    grad = nd.Gradient(func,step_nom=[0.01]*12)

    for coords in [
      [[ 0 ,0, 0],[0,1, 0],[-1,-1,0],[1,-1,0]],
      [[ 0 ,0, 0],[0,1, 0],[-2,-2,0],[1,-1,0]],
      [[ 0 ,0, 0],[0,1, 1],[-1,-1,0],[1,-1,0]],
      [[ 0 ,0, 0],[0,1, 0],[-1,-1,0],[1,-1,1]],
      [[ 0 ,0, 0],[0,1, 1],[-2,-2,1],[2,-2,1]],
      [[.1,.1,.9],[0,1, 0],[-1,-1,0],[1,-1,0]],
                  ]:
      coords = np.array(coords,dtype=np.float64)
      ref = grad(coords.flatten())
      res = ric.construct_b_matrix(None,coords)[0,:]
      self.assertLess(np.max(np.abs(ref-res)),1.e-10)

  def test_grad_0D_rand(self):

    masses = np.array([1.]*4)

    ric = RIC()
    ric.add_out_bend([1,2,3,4])
    ric.setup(masses)

    def func(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_out_bends())
      return res
    grad = nd.Gradient(func,step_nom=[0.01]*12)

    for coords in 2*np.random.random((10,4,3))-1:
      coords = np.array(coords,dtype=np.float64)
      ref = grad(coords.flatten())
      res = ric.construct_b_matrix(None,coords)[0,:]
      self.assertLess(np.max(np.abs(ref-res)),1.e-9)

if __name__ == '__main__':
  ut.main()

