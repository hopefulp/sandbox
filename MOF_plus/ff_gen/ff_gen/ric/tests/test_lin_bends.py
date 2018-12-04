#! /usb/bin/env python

import numpy        as np
import numdifftools as nd
import unittest     as ut

from ric import RedIntCoords as RIC

class LinearBendTests(ut.TestCase):

  def test_val_0D_xy(self):
    """Values without PBC"""

    masses = np.array([1.]*3)

    ric = RIC()
    ric.add_lin_bend([1,2,3],'xy')
    ric.setup(masses)

    for coords, ref in [
      ([[ 0, 0,-1],[ 0, 0, 0],[ 0, 0, 1]], [  0,   0]), # Ref
      ([[-1,-1,-1],[ 0, 0, 0],[ 1, 1, 1]], [  0,   0]), # Diagonal
      ([[ 0, 0, 0],[ 1, 1, 1],[ 2, 2, 2]], [  0,   0]), # Offset
      ([[ 0, 0,-3],[ 0, 0, 0],[ 0, 0, 5]], [  0,   0]), # Longer
      ([[ 0, 0,-1],[ 0, 0, 0],[-1, 0, 1]], [  0,  45]), # -x
      ([[ 0, 0,-1],[ 0, 0, 0],[ 1, 0, 1]], [  0, -45]), # +x
      ([[ 0, 0,-1],[ 0, 0, 0],[ 0,-1, 1]], [-45,   0]), # -y
      ([[ 0, 0,-1],[ 0, 0, 0],[ 0, 1, 1]], [ 45,   0]), # +y
                       ]:

      ric.construct_b_matrix(None,np.array(coords,dtype=np.float64))
      res = ric.get_val_lin_bends()
      self.assertAlmostEqual(res[0],np.radians(ref[0]))
      self.assertAlmostEqual(res[1],np.radians(ref[1]))

  def test_val_0D_xz(self):
    """Values without PBC"""

    masses = np.array([1.]*3)

    ric = RIC()
    ric.add_lin_bend([1,2,3],'xz')
    ric.setup(masses)

    for coords, ref in [
      ([[ 0,-1, 0],[ 0, 0, 0],[ 0, 1, 0]], [  0,   0]), # Ref
      ([[-1,-1,-1],[ 0, 0, 0],[ 1, 1, 1]], [  0,   0]), # Diagonal
      ([[ 0, 0, 0],[ 1, 1, 1],[ 2, 2, 2]], [  0,   0]), # Offset
      ([[ 0,-3, 0],[ 0, 0, 0],[ 0, 5, 0]], [  0,   0]), # Longer
      ([[ 0,-1, 0],[ 0, 0, 0],[-1, 1, 0]], [  0, -45]), # -x
      ([[ 0,-1, 0],[ 0, 0, 0],[ 1, 1, 0]], [  0,  45]), # +x
      ([[ 0,-1, 0],[ 0, 0, 0],[ 0, 1,-1]], [ 45,   0]), # -z
      ([[ 0,-1, 0],[ 0, 0, 0],[ 0, 1, 1]], [-45,   0]), # +z
                       ]:

      ric.construct_b_matrix(None,np.array(coords,dtype=np.float64))
      res = ric.get_val_lin_bends()
      self.assertAlmostEqual(res[0],np.radians(ref[0]))
      self.assertAlmostEqual(res[1],np.radians(ref[1]))

  def test_val_0D_yz(self):
    """Values without PBC"""

    masses = np.array([1.]*3)

    ric = RIC()
    ric.add_lin_bend([1,2,3],'yz')
    ric.setup(masses)

    for coords, ref in [
      ([[-1, 0, 0],[ 0, 0, 0],[ 1, 0, 0]], [  0,   0]), # Ref
      ([[-1,-1,-1],[ 0, 0, 0],[ 1, 1, 1]], [  0,   0]), # Diagonal
      ([[ 0, 0, 0],[ 1, 1, 1],[ 2, 2, 2]], [  0,   0]), # Offset
      ([[-3, 0, 0],[ 0, 0, 0],[ 5, 0, 0]], [  0,   0]), # Longer
      ([[-1, 0, 0],[ 0, 0, 0],[ 1,-1, 0]], [  0,  45]), # -y
      ([[-1, 0, 0],[ 0, 0, 0],[ 1, 1, 0]], [  0, -45]), # +y
      ([[-1, 0, 0],[ 0, 0, 0],[ 1, 0,-1]], [-45,   0]), # -z
      ([[-1, 0, 0],[ 0, 0, 0],[ 1, 0, 1]], [ 45,   0]), # +z
                       ]:

      ric.construct_b_matrix(None,np.array(coords,dtype=np.float64))
      res = ric.get_val_lin_bends()
      self.assertAlmostEqual(res[0],np.radians(ref[0]))
      self.assertAlmostEqual(res[1],np.radians(ref[1]))


  def test_val_0D_ind(self):
    """Values without PBC"""

    masses = np.array([1.]*4)

    ric = RIC()
    ric.add_lin_bend([1,2,3],4)
    ric.setup(masses)

    for coords, ref in [
      ([[-1, 0, 0],[ 0, 0, 0],[ 1, 0, 0], [ 0, 1, 0]], [  0,  0]), # x-y
      ([[-1, 0, 0],[ 0, 0, 0],[ 1, 0, 0], [ 0, 0, 1]], [  0,  0]), # x-z
      ([[-1, 0, 0],[ 0, 0, 0],[ 1, 0, 0], [ 0, 1, 1]], [  0,  0]), # x-yz
      ([[-1, 0, 0],[ 0, 0, 0],[ 1, 0, 0], [ 1, 1, 1]], [  0,  0]), # x-xyz
      ([[ 0,-1, 0],[ 0, 0, 0],[ 0, 1, 0], [ 2, 0, 0]], [  0,  0]), # y-x
      ([[ 0,-1, 0],[ 0, 0, 0],[ 0, 1, 0], [ 0, 0, 2]], [  0,  0]), # y-z
      ([[ 0,-1, 0],[ 0, 0, 0],[ 0, 1, 0], [ 0, 2, 2]], [  0,  0]), # y-yz
      ([[ 0,-1, 0],[ 0, 0, 0],[ 0, 1, 0], [ 2, 2, 2]], [  0,  0]), # y-xyz
      ([[ 0, 0,-1],[ 0, 0, 0],[ 0, 0, 1], [ 3, 0, 0]], [  0,  0]), # z-x
      ([[ 0, 0,-1],[ 0, 0, 0],[ 0, 0, 1], [ 0, 3, 0]], [  0,  0]), # z-y
      ([[ 0, 0,-1],[ 0, 0, 0],[ 0, 0, 1], [ 3, 3, 0]], [  0,  0]), # z-xy
      ([[ 0, 0,-1],[ 0, 0, 0],[ 0, 0, 1], [ 3, 3, 3]], [  0,  0]), # z-xyz
                       ]:

      ric.construct_b_matrix(None,np.array(coords,dtype=np.float64))
      res = ric.get_val_lin_bends()
      #print coords
      #print ref
      #print res
      self.assertAlmostEqual(res[0],np.radians(ref[0]))
      self.assertAlmostEqual(res[1],np.radians(ref[1]))

  def test_grad_0D_xy(self):

    masses = np.array([1.]*3)

    ric = RIC()
    ric.add_lin_bend([1,2,3],'xy')
    ric.setup(masses)

    def func(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res
    jac = nd.Jacobian(func,step_nom=[0.01]*9) # This is buggy!

    def func0(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res[0]
    def func1(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res[1]
    grad0 = nd.Gradient(func0,step_nom=[0.01]*9)
    grad1 = nd.Gradient(func1,step_nom=[0.01]*9)

    for coords in [
      [[ 0, 0,-1],[ 0, 0, 0],[ 0, 0, 1]], # Ref
      [[-1,-1,-1],[ 0, 0, 0],[ 1, 1, 1]], # Diagonal
      [[ 0, 0,-3],[ 0, 0, 0],[ 0, 0, 5]], # Longer
      [[ 0, 0,-1],[ 0, 0, 0],[-1, 0, 1]], # -x
      [[ 0, 0,-1],[ 0, 0, 0],[ 1, 0, 1]], # +x
      [[ 0, 0,-1],[ 0, 0, 0],[ 0,-1, 1]], # -y
      [[ 0, 0,-1],[ 0, 0, 0],[ 0, 1, 1]], # +y
                  ]:
      coords = np.array(coords,dtype=np.float64)
      ref  =  jac(coords.flatten())
      ref0 = grad0(coords.flatten())
      ref1 = grad1(coords.flatten())
      res = ric.construct_b_matrix(None,coords)
      #print coords
      #print ref
      #print ref0
      #print ref1
      #print res
      #self.assertLess(np.max(np.abs(ref-res)),1.e-8)
      self.assertLess(np.max(np.abs(ref0-res[0])),1.e-8)
      self.assertLess(np.max(np.abs(ref1-res[1])),1.e-8)

  def test_grad_0D_xz(self):

    masses = np.array([1.]*3)

    ric = RIC()
    ric.add_lin_bend([1,2,3],'xz')
    ric.setup(masses)

    def func(x):
      assert x.dtype == np.float64
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      assert res.dtype == np.float64
      return res
    jac = nd.Jacobian(func,step_nom=[0.01]*9) # This is buggy!

    def func0(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res[0]
    def func1(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res[1]
    grad0 = nd.Gradient(func0,step_nom=[0.01]*9)
    grad1 = nd.Gradient(func1,step_nom=[0.01]*9)

    for coords in [
      [[ 0,-1, 0],[ 0, 0, 0],[ 0, 1, 0]], # Ref
      [[-1,-1,-1],[ 0, 0, 0],[ 1, 1, 1]], # Diagonal
      [[ 0, 0, 0],[ 1, 1, 1],[ 2, 2, 2]], # Offset
      [[ 0,-3, 0],[ 0, 0, 0],[ 0, 5, 0]], # Longer
      [[ 0,-1, 0],[ 0, 0, 0],[-1, 1, 0]], # -x
      [[ 0,-1, 0],[ 0, 0, 0],[ 1, 1, 0]], # +x
      [[ 0,-1, 0],[ 0, 0, 0],[ 0, 1,-1]], # -z
      [[ 0,-1, 0],[ 0, 0, 0],[ 0, 1, 1]], # +z
                  ]:
      coords = np.array(coords,dtype=np.float64)
      #ref  =  jac(coords.flatten())
      ref0 = grad0(coords.flatten())
      ref1 = grad1(coords.flatten())
      res = ric.construct_b_matrix(None,coords)
      #print coords
      #print ref
      #print ref0
      #print ref1
      #print res
      #self.assertLess(np.max(np.abs(ref-res)),1.e-8)
      self.assertLess(np.max(np.abs(ref0-res[0])),1.e-8)
      self.assertLess(np.max(np.abs(ref1-res[1])),1.e-8)

  def test_grad_0D_yz(self):

    masses = np.array([1.]*3)

    ric = RIC()
    ric.add_lin_bend([1,2,3],'yz')
    ric.setup(masses)

    def func(x):
      assert x.dtype == np.float64
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      assert res.dtype == np.float64
      return res
    jac = nd.Jacobian(func,step_nom=[0.01]*9) # This is buggy!

    def func0(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res[0]
    def func1(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res[1]
    grad0 = nd.Gradient(func0,step_nom=[0.01]*9)
    grad1 = nd.Gradient(func1,step_nom=[0.01]*9)

    for coords in [
      [[-1, 0, 0],[ 0, 0, 0],[ 1, 0, 0]], # Ref
      [[-1,-1,-1],[ 0, 0, 0],[ 1, 1, 1]], # Diagonal
      [[ 0, 0, 0],[ 1, 1, 1],[ 2, 2, 2]], # Offset
      [[-3, 0, 0],[ 0, 0, 0],[ 5, 0, 0]], # Longer
      [[-1, 0, 0],[ 0, 0, 0],[ 1,-1, 0]], # -y
      [[-1, 0, 0],[ 0, 0, 0],[ 1, 1, 0]], # +y
      [[-1, 0, 0],[ 0, 0, 0],[ 1, 0,-1]], # -z
      [[-1, 0, 0],[ 0, 0, 0],[ 1, 0, 1]], # +z
                  ]:
      coords = np.array(coords,dtype=np.float64)
      #ref  =  jac(coords.flatten())
      ref0 = grad0(coords.flatten())
      ref1 = grad1(coords.flatten())
      res = ric.construct_b_matrix(None,coords)
      #print coords
      #print ref
      #print ref0
      #print ref1
      #print res
      #self.assertLess(np.max(np.abs(ref-res)),1.e-8)
      self.assertLess(np.max(np.abs(ref0-res[0])),1.e-8)
      self.assertLess(np.max(np.abs(ref1-res[1])),1.e-8)

  def test_grad_0D_ind(self):

    masses = np.array([1.]*4)

    ric = RIC()
    ric.add_lin_bend([1,2,3], 4)
    ric.setup(masses)

    def func(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res
    jac = nd.Jacobian(func,step_nom=[0.01]*12) # This is buggy!

    def func0(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res[0]
    def func1(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res[1]
    grad0 = nd.Gradient(func0,step_nom=[0.01]*12)
    grad1 = nd.Gradient(func1,step_nom=[0.01]*12)

    for coords in [
      [[-1, 0, 0],[.1, 0, 0],[ 1, 0, 0], [ 0, 1, 0]], # x-y
      [[-1, 0, 0],[.1, 0, 0],[ 1, 0, 0], [ 0, 0, 1]], # x-z
      [[-1, 0, 0],[.1, 0, 0],[ 1, 0, 0], [ 0, 1, 1]], # x-yz
      [[-1, 0, 0],[.1, 0, 0],[ 1, 0, 0], [ 1, 1, 1]], # x-xyz
      [[ 0,-1, 0],[ 0,.1, 0],[ 0, 1, 0], [ 2, 0, 0]], # y-x
      [[ 0,-1, 0],[ 0,.1, 0],[ 0, 1, 0], [ 0, 0, 2]], # y-z
      [[ 0,-1, 0],[ 0,.1, 0],[ 0, 1, 0], [ 0, 2, 2]], # y-yz
      [[ 0,-1, 0],[ 0,.1, 0],[ 0, 1, 0], [ 2, 2, 2]], # y-xyz
      [[ 0, 0,-1],[ 0, 0,.1],[ 0, 0, 1], [ 3, 0, 0]], # z-x
      [[ 0, 0,-1],[ 0, 0,.1],[ 0, 0, 1], [ 0, 3, 0]], # z-y
      [[ 0, 0,-1],[ 0, 0,.1],[ 0, 0, 1], [ 3, 3, 0]], # z-xy
      [[ 0, 0,-1],[ 0, 0,.1],[ 0, 0, 1], [ 3, 3, 3]], # z-xyz
                  ]:
      coords = np.array(coords,dtype=np.float64)
      res = ric.construct_b_matrix(None,coords)

      # HACK: disable the axes updates
      _inds = np.copy(ric._ric.ric_lin_bend_inds)
      ric._ric.ric_lin_bend_inds[:] = 0

      #ref  =   jac(coords.flatten())
      ref0 = grad0(coords.flatten())
      ref1 = grad1(coords.flatten())

      # HACK: enbale the axes updates
      ric._ric.ric_lin_bend_inds[:] = _inds

      #print coords
      #print ref
      #print ref0
      #print ref1
      #print res
      #self.assertLess(np.max(np.abs(ref-res)),1.e-8)
      self.assertLess(np.max(np.abs(ref0-res[0])),1.e-8)
      self.assertLess(np.max(np.abs(ref1-res[1])),1.e-8)

  def test_grad_0D_ind_rand(self):

    masses = np.array([1.]*4)

    ric = RIC()
    ric.add_lin_bend([1,2,3], 4)
    ric.setup(masses)

    def func(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res
    jac = nd.Jacobian(func,step_nom=[0.01]*12) # This is buggy!

    def func0(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res[0]
    def func1(x):
      ric.construct_b_matrix(None,np.reshape(x,(-1,3)))
      res = np.copy(ric.get_val_lin_bends())
      return res[1]
    grad0 = nd.Gradient(func0,step_nom=[0.01]*12)
    grad1 = nd.Gradient(func1,step_nom=[0.01]*12)

    for coords in 2*np.random.random((10,4,3))-1:
      coords = np.array(coords,dtype=np.float64)
      res  = ric.construct_b_matrix(None,coords)
      vals = np.degrees(ric.get_val_lin_bends())
      #if np.any(np.abs(vals) > 45): continue

      # HACK: disable the axes updates
      _inds = np.copy(ric._ric.ric_lin_bend_inds)
      ric._ric.ric_lin_bend_inds[:] = 0

      #ref  =   jac(coords.flatten())
      ref0 = grad0(coords.flatten())
      ref1 = grad1(coords.flatten())

      # HACK: enbale the axes updates
      ric._ric.ric_lin_bend_inds = _inds

      #print coords
      #print vals
      #print ref
      #print ref0
      #print ref1
      #print res
      #self.assertLess(np.max(np.abs(ref-res)),1.e-8)
      self.assertLess(np.max(np.abs(ref0-res[0])),1.e-8)
      self.assertLess(np.max(np.abs(ref1-res[1])),1.e-8)


if __name__ == '__main__':
  ut.main()

