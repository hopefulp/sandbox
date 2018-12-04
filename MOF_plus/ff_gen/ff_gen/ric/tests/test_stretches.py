#! /usb/bin/env python

import numpy        as np
import numdifftools as nd
import unittest     as ut

from ric import RedIntCoords as RIC

class StretchTests(ut.TestCase):

  def test_val_0D_1(self):
    """Values without PBC"""

    masses = np.array([1.]*2)

    ric = RIC()
    ric.add_stretch([1,2])
    ric.setup(masses)

    for coords, ref in [([[-1, 0, 0],[ 1, 0, 0]],            2),
                        ([[ 0,-1, 0],[ 0, 1, 0]],            2),
                        ([[ 0, 0,-1],[ 0, 0, 1]],            2),
                        ([[ 1, 1, 1],[-1,-1,-1]], 2*np.sqrt(3))]:

      ric.construct_b_matrix(None,np.array(coords,dtype=np.float64))
      res = ric.get_val_stretches()[0]
      self.assertAlmostEqual(res,ref)

  def test_grad_0D_1(self):

    masses = np.array([1.]*2)

    ric = RIC()
    ric.add_stretch([1,2])
    ric.setup(masses)

    def func(x):
      ric.construct_b_matrix(None,np.reshape(x,(2,3)))
      return ric.get_val_stretches()[0]
    grad = nd.Gradient(func,step_nom=[0.01]*6)

    for coords in [[[-1, 0, 0],[ 1, 0, 0]],
                   [[ 0,-1, 0],[ 0, 1, 0]],
                   [[ 0, 0,-1],[ 0, 0, 1]],
                   [[ 2, 2, 2],[-1,-1,-1]]]:
      coords = np.array(coords,dtype=np.float64)
      ref = grad(coords.flatten())
      res = ric.construct_b_matrix(None,coords)[0,:]
      self.assertLess(np.max(np.abs(ref-res)),1.e-12)

  def test_val_3D_1(self):
    """Values with PBC in an orthogonal box"""

    masses = np.array([1.]*2)
    hmat   = np.array([[5.,0.,0.],
                       [0.,6.,0.],
                       [0.,0.,7.]])

    ric = RIC()
    ric.add_stretch([1,2])
    ric.setup(masses)

    for coords, ref in [([[-1, 0, 0],[ 1, 0, 0]],            2),
                        ([[ 0,-1, 0],[ 0, 1, 0]],            2),
                        ([[ 0, 0,-1],[ 0, 0, 1]],            2),
                        ([[ 1, 1, 1],[-1,-1,-1]], 2*np.sqrt(3))]:

      ric.construct_b_matrix(hmat,np.array(coords,dtype=np.float64))
      res = ric.get_val_stretches()[0]
      self.assertAlmostEqual(res,ref)

  def test_grad_3D_1(self):

    masses = np.array([1.]*2)
    hmat   = np.array([[5.,0.,0.],
                       [0.,6.,0.],
                       [0.,0.,7.]])

    ric = RIC()
    ric.add_stretch([1,2])
    ric.setup(masses)

    def func(x):
      ric.construct_b_matrix(hmat,np.reshape(x,(2,3)))
      return ric.get_val_stretches()[0]
    grad = nd.Gradient(func,step_nom=[0.01]*12)

    for coords in [[[-1, 0, 0],[ 1, 0, 0]],
                   [[ 0,-1, 0],[ 0, 1, 0]],
                   [[ 0, 0,-1],[ 0, 0, 1]],
                   [[ 1, 1, 1],[-1,-1,-1]]]:
      coords = np.array(coords,dtype=np.float64)
      ref = grad(coords.flatten())
      res = ric.construct_b_matrix(hmat,coords)[0,:]
      self.assertLess(np.max(np.abs(ref-res)),1.e-12)

if __name__ == '__main__':
  ut.main()

