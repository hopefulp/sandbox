#! /usb/bin/env python

import numpy        as np
import numpy.linalg as nl
import unittest     as ut

from ric import RedIntCoords as RIC
from ric import Group

class H2CCCH2Tests(ut.TestCase):

  def setUp(self):

    self.masses = np.array([1,1,12,12,12,1,1], dtype=np.float64)
    self.coords = np.array([[-2,-1, 0], # 1
                            [-2, 1, 0], # 2
                            [-1, 0, 0], # 3
                            [ 0, 0, 0], # 4
                            [ 1, 0, 0], # 5
                            [ 2,-1, 0], # 6
                            [ 2, 1, 0], # 7
                           ], dtype=np.float64)
    self.hess   = np.random.random((self.coords.size,self.coords.size))
    self.hess   = .5*(self.hess + np.transpose(self.hess)) # Make it symmetric

  def test_1(self):

    ric = RIC()
    ric.add_stretch([1,3])
    ric.add_stretch([2,3])
    ric.add_stretch([3,4])
    ric.add_stretch([4,5])
    ric.add_stretch([5,6])
    ric.add_stretch([5,7])
    ric.add_in_bend([1,3,2])
    ric.add_in_bend([1,3,4])
    ric.add_in_bend([2,3,4])
    ric.add_in_bend([4,5,6])
    ric.add_in_bend([4,5,7])
    ric.add_in_bend([6,5,7])
    ric.add_lin_bend([3,4,5], 'yz')
    ric.add_out_bend([3,1,2,4])
    ric.add_out_bend([3,2,4,1])
    ric.add_out_bend([3,4,1,2])
    ric.add_out_bend([5,4,6,7])
    ric.add_out_bend([5,6,7,4])
    ric.add_out_bend([5,7,4,6])
    ric.add_torsion([1,2,0,0,0,3,5,6,7,0,0,0])
    ric.add_eckart()
    ric.setup(self.masses)

    bmat = ric.construct_b_matrix(None, self.coords)

    bmat_inv, rank = ric.invert_b_matrix()
    bmat_inv_ = nl.pinv(bmat)
    self.assertEqual(rank, self.coords.size)
    self.assertLess(np.max(np.abs(bmat_inv - bmat_inv_)), 7.e-14)

    hess = ric.project_hessian(self.hess)
    hess_ = np.dot(bmat_inv.T,np.dot(self.hess,bmat_inv))
    self.assertLess(np.max(np.abs(hess - hess_)),3.e-13)

  def test_2(self):

    g1, g2 = Group([1,2,3]), Group([5,6,7])

    ric = RIC()
    ric.add_stretch([1,3])
    ric.add_stretch([2,3])
    ric.add_distance([g1,4])
    ric.add_distance([4,g2])
    ric.add_stretch([5,6])
    ric.add_stretch([5,7])
    ric.add_in_bend([1,3,2])
    ric.add_in_bend([1,3,4])
    ric.add_in_bend([2,3,4])
    ric.add_in_bend([4,5,6])
    ric.add_in_bend([4,5,7])
    ric.add_in_bend([6,5,7])
    ric.add_lin_bend([3,4,5], 'yz')
    ric.add_out_bend([3,1,2,4])
    ric.add_out_bend([3,2,4,1])
    ric.add_out_bend([3,4,1,2])
    ric.add_out_bend([5,4,6,7])
    ric.add_out_bend([5,6,7,4])
    ric.add_out_bend([5,7,4,6])
    ric.add_dihedral([1,2], [g1,g2], [6,7])
    ric.add_eckart()
    ric.setup(self.masses)

    bmat = ric.construct_b_matrix(None, self.coords)

    bmat_inv, rank = ric.invert_b_matrix()
    bmat_inv_ = nl.pinv(bmat)
    self.assertEqual(rank, self.coords.size)
    self.assertLess(np.max(np.abs(bmat_inv - bmat_inv_)), 1.e-13)

    hess = ric.project_hessian(self.hess)
    hess_ = np.dot(bmat_inv.T,np.dot(self.hess,bmat_inv))
    self.assertLess(np.max(np.abs(hess - hess_)),1.e-13)


if __name__ == '__main__':
  ut.main()

