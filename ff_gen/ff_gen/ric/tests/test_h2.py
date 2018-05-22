#! /usb/bin/env python

import numpy        as np
import numpy.linalg as nl
import unittest     as ut

from ric import RedIntCoords as RIC

class H2Tests(ut.TestCase):

  def setUp(self):

    self.masses = np.array([1.,1.])
    self.hmat   = np.array([[0.,0.,0.]]*3)
    self.coords = np.array([[1.,2.,3.],
                            [4.,5.,6.]])
    self.hess   = np.random.randn(6,6)
    self.hess   = .5*(self.hess + np.transpose(self.hess)) # Make it symmetric

  def test_1(self):

    ric = RIC()
    ric.add_stretch([1,2])
    ric.add_eckart()
    ric.setup(self.masses)

    bmat = ric.construct_b_matrix(self.hmat,self.coords)

    self.assertEqual(ric.get_val_stretches()[0],
                     np.sqrt(np.sum((self.coords[0,:]-self.coords[1,:])**2)))

    bmat_inv, rank = ric.invert_b_matrix()
    bmat_inv_ = nl.pinv(bmat)
    self.assertLess(np.max(np.abs(bmat_inv -  bmat_inv_)), 1.e-14)

    hess = ric.project_hessian(self.hess)
    hess_ = np.dot(bmat_inv.T,np.dot(self.hess,bmat_inv))
    self.assertLess(np.max(np.abs(hess - hess_)),1.e-14)


if __name__ == '__main__':
  ut.main()

