#! /usb/bin/env python

import numpy        as np
import numpy.linalg as nl
import unittest     as ut

from ric import RedIntCoords as RIC

class H2OTests(ut.TestCase):

  def setUp(self):

    self.masses = np.array([1.,16.,1.])
    self.hmat   = np.array([[0.,0.,0.]]*3)
    self.coords = np.array([[0.,-1.,0.],
                            [0., 0.,1.],
                            [0., 1.,0.]])
    self.hess   = np.random.randn(9,9)
    self.hess   = .5*(self.hess + np.transpose(self.hess)) # Make it symmetric

  def test_1(self):

    ric = RIC()
    ric.add_stretch([1,2])
    ric.add_stretch([2,3])
    ric.add_in_bend([2,1,3])
    ric.add_eckart()
    ric.setup(self.masses)

    bmat = ric.construct_b_matrix(self.hmat,self.coords)

    self.assertEqual(ric.get_val_stretches()[0],
                     np.sqrt(np.sum((self.coords[0,:]-self.coords[1,:])**2)))

    bmat_inv, rank = ric.invert_b_matrix()
    bmat_inv_ = nl.pinv(bmat)
    self.assertLess(np.max(np.abs(bmat_inv - bmat_inv_)), 2.e-13)

    hess = ric.project_hessian(self.hess)
    hess_ = np.dot(bmat_inv.T,np.dot(self.hess,bmat_inv))
    self.assertLess(np.max(np.abs(hess - hess_)),3.e-14)

  def test_2(self):
    """Multiple repeated calls"""

    ric = RIC()
    ric.add_stretch([1,2])
    ric.add_stretch([2,3])
    ric.add_in_bend([2,1,3])
    ric.add_eckart()
    ric.setup(self.masses)

    bmat = ric.construct_b_matrix(None,self.coords)
    b1 = np.copy(bmat)
    bmat = ric.construct_b_matrix(None,self.coords)
    b2 = np.copy(bmat)
    self.assertEqual(np.max(np.abs(b1-b2)),0)

    bmat_inv, rank = ric.invert_b_matrix()
    i1 = np.copy(bmat_inv)
    bmat = ric.construct_b_matrix(None,self.coords)
    b3 = np.copy(bmat)
    bmat_inv, rank = ric.invert_b_matrix()
    i2 = np.copy(bmat_inv)
    self.assertEqual(np.max(np.abs(b1-b3)),0)
    self.assertEqual(np.max(np.abs(i1-i2)),0)

    hess = ric.project_hessian(self.hess)
    h1 = np.copy(hess)
    bmat = ric.construct_b_matrix(None,self.coords)
    b4 = np.copy(bmat)
    bmat_inv, rank = ric.invert_b_matrix()
    i3 = np.copy(bmat_inv)
    hess = ric.project_hessian(self.hess)
    h2 = np.copy(hess)
    self.assertEqual(np.max(np.abs(b1-b4)),0)
    self.assertEqual(np.max(np.abs(i1-i3)),0)
    self.assertEqual(np.max(np.abs(h1-h2)),0)


if __name__ == '__main__':
  ut.main()

