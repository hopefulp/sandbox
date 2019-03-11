#! /usb/bin/env python

import numpy    as np
import unittest as ut

from ric import RedIntCoords as RIC

class BasicTests(ut.TestCase):

  def test_1(self):
    """Mutiple instances safeguard"""

    ric = RIC()
    del ric
    ric = RIC()
    try:
      ric = RIC()
    except AssertionError:
      pass
    del ric

  def test_2(self):

    masses = np.array(range(1,5),dtype=np.float64)

    ric = RIC()
    ric.add_stretch([1,2])
    ric.add_stretch([3,4])
    ric.setup(masses)

    self.assertListEqual(list(ric._ric.ric_def_stretches[:,0]),[1,2])
    self.assertListEqual(list(ric._ric.ric_def_stretches[:,1]),[3,4])
    self.assertListEqual(list(ric._ric.ric_ibr_stretches),[1,2])

    self.assertEqual(ric.num_ric     ,2)
    self.assertEqual(ric.num_stretch ,2)
    self.assertEqual(ric.num_in_bend ,0)
    self.assertEqual(ric.num_out_bend,0)
    self.assertEqual(ric.num_lin_bend,0)
    self.assertEqual(ric.num_torsion ,0)
    self.assertEqual(ric.num_eckart  ,0)

  def test_3(self):

    masses = np.array(range(1,6),dtype=np.float64)

    ric = RIC()
    ric.add_in_bend([1,2,3])
    ric.add_in_bend([3,4,5])
    ric.setup(masses)

    self.assertListEqual(list(ric._ric.ric_def_in_bends[:,0]),[1,2,3])
    self.assertListEqual(list(ric._ric.ric_def_in_bends[:,1]),[3,4,5])
    self.assertListEqual(list(ric._ric.ric_ibr_in_bends),[1,2])

    self.assertEqual(ric.num_ric     ,2)
    self.assertEqual(ric.num_stretch ,0)
    self.assertEqual(ric.num_in_bend ,2)
    self.assertEqual(ric.num_out_bend,0)
    self.assertEqual(ric.num_lin_bend,0)
    self.assertEqual(ric.num_torsion ,0)
    self.assertEqual(ric.num_eckart  ,0)

  def test_4(self):

    masses = np.array(range(1,7),dtype=np.float64)

    ric = RIC()
    ric.add_out_bend([1,2,3,4])
    ric.add_out_bend([3,4,5,6])
    ric.setup(masses)

    self.assertListEqual(list(ric._ric.ric_def_out_bends[:,0]),[1,2,3,4])
    self.assertListEqual(list(ric._ric.ric_def_out_bends[:,1]),[3,4,5,6])
    self.assertListEqual(list(ric._ric.ric_ibr_out_bends),[1,2])

    self.assertEqual(ric.num_ric     ,2)
    self.assertEqual(ric.num_stretch ,0)
    self.assertEqual(ric.num_in_bend ,0)
    self.assertEqual(ric.num_out_bend,2)
    self.assertEqual(ric.num_lin_bend,0)
    self.assertEqual(ric.num_torsion ,0)
    self.assertEqual(ric.num_eckart  ,0)

  def test_5(self):

    masses = np.array(range(1,14),dtype=np.float64)

    ric = RIC()
    ric.add_lin_bend([1,2,3],'xy')
    ric.add_lin_bend([4,5,6],'xz')
    ric.add_lin_bend([7,8,9],'yz')
    ric.add_lin_bend([10,11,12],13)
    ric.setup(masses)

    self.assertListEqual(list(ric._ric.ric_def_lin_bends[:,0]),[1,2,3])
    self.assertListEqual(list(ric._ric.ric_def_lin_bends[:,1]),[1,2,3])
    self.assertListEqual(list(ric._ric.ric_def_lin_bends[:,2]),[4,5,6])
    self.assertListEqual(list(ric._ric.ric_def_lin_bends[:,3]),[4,5,6])
    self.assertListEqual(list(ric._ric.ric_def_lin_bends[:,4]),[7,8,9])
    self.assertListEqual(list(ric._ric.ric_def_lin_bends[:,5]),[7,8,9])
    self.assertListEqual(list(ric._ric.ric_def_lin_bends[:,6]),[10,11,12])
    self.assertListEqual(list(ric._ric.ric_def_lin_bends[:,7]),[10,11,12])
    self.assertListEqual(list(ric._ric.ric_ibr_lin_bends),[1,2,3,4,5,6,7,8])
    self.assertListEqual(list(ric._ric.ric_lin_bend_inds),[0,0,0,0,0,0,13,-13])
    self.assertListEqual(list(ric._ric.ric_lin_bend_axes[:,0]),[1,0,0])
    self.assertListEqual(list(ric._ric.ric_lin_bend_axes[:,1]),[0,1,0])
    self.assertListEqual(list(ric._ric.ric_lin_bend_axes[:,2]),[1,0,0])
    self.assertListEqual(list(ric._ric.ric_lin_bend_axes[:,3]),[0,0,1])
    self.assertListEqual(list(ric._ric.ric_lin_bend_axes[:,4]),[0,1,0])
    self.assertListEqual(list(ric._ric.ric_lin_bend_axes[:,5]),[0,0,1])
    self.assertListEqual(list(ric._ric.ric_lin_bend_axes[:,6]),[0,0,0])
    self.assertListEqual(list(ric._ric.ric_lin_bend_axes[:,7]),[0,0,0])

    self.assertEqual(ric.num_ric     ,8)
    self.assertEqual(ric.num_stretch ,0)
    self.assertEqual(ric.num_in_bend ,0)
    self.assertEqual(ric.num_out_bend,0)
    self.assertEqual(ric.num_lin_bend,8)
    self.assertEqual(ric.num_torsion ,0)
    self.assertEqual(ric.num_eckart  ,0)

  def test_6(self):

    masses = np.array(range(1,10),dtype=np.float64)

    ric = RIC()
    ric.add_torsion([1,0,0,0,0,2,3,4,0,0,0,0])
    ric.add_torsion([3,4,0,0,0,5,6,7,8,9,0,0],ivals=[2,3])
    ric.setup(masses)

    self.assertListEqual(list(ric._ric.ric_def_torsions[:,0]),[1,0,0,0,0,2,3,4,0,0,0,0])
    self.assertListEqual(list(ric._ric.ric_def_torsions[:,1]),[3,4,0,0,0,5,6,7,8,9,0,0])
    self.assertListEqual(list(ric._ric.ric_torsion_ivals[:,0]),[1,1])
    self.assertListEqual(list(ric._ric.ric_torsion_ivals[:,1]),[2,3])
    self.assertListEqual(list(ric._ric.ric_ibr_torsions),[1,2])

    self.assertEqual(ric.num_ric     ,2)
    self.assertEqual(ric.num_stretch ,0)
    self.assertEqual(ric.num_in_bend ,0)
    self.assertEqual(ric.num_out_bend,0)
    self.assertEqual(ric.num_lin_bend,0)
    self.assertEqual(ric.num_torsion ,2)
    self.assertEqual(ric.num_eckart  ,0)

  def test_7(self):

    masses = np.array(range(1,2),dtype=np.float64)

    ric = RIC()
    ric.add_eckart(rots=[False,0,'a'])
    ric.setup(masses)

    self.assertListEqual(list(ric._ric.ric_def_eckart_rots),[3])
    self.assertListEqual(list(ric._ric.ric_ibr_eckart_rots),[4])

    self.assertEqual(ric.num_ric     ,4)
    self.assertEqual(ric.num_stretch ,0)
    self.assertEqual(ric.num_in_bend ,0)
    self.assertEqual(ric.num_out_bend,0)
    self.assertEqual(ric.num_lin_bend,0)
    self.assertEqual(ric.num_torsion ,0)
    self.assertEqual(ric.num_eckart  ,4)

  def test_8(self):

    masses = np.array(range(1,10),dtype=np.float64)

    ric = RIC()
    ric.add_stretch([2,3])
    ric.add_in_bend([3,4,5])
    ric.add_out_bend([4,5,6,7])
    ric.add_lin_bend([5,6,7],'xy')
    ric.add_torsion([1,0,0,0,0,2,3,4,5,6,7,8])
    ric.add_torsion([3,4,0,0,0,5,6,7,8,9,0,0])
    ric.add_eckart(rots=[True,False,True])
    ric.setup(masses)

    self.assertTrue(np.all(ric._ric.atomic_masses == masses))
    self.assertListEqual(list(ric._ric.cart_coords.shape),[3,9])
    self.assertListEqual(list(ric._ric.cart_hessian.shape),[3*9,3*9])
    self.assertListEqual(list(ric._ric.ric_ibr_stretches),[1])
    self.assertListEqual(list(ric._ric.ric_ibr_in_bends ),[2])
    self.assertListEqual(list(ric._ric.ric_ibr_out_bends),[3])
    self.assertListEqual(list(ric._ric.ric_ibr_lin_bends),[4,5])
    self.assertListEqual(list(ric._ric.ric_ibr_torsions) ,[6,7])
    self.assertListEqual(list(ric._ric.ric_ibr_eckart_trans),[8,9,10])
    self.assertListEqual(list(ric._ric.ric_ibr_eckart_rots) ,[11,12])
    self.assertListEqual(list(ric._ric.bmat.shape),[3*9,12])
    self.assertEqual(ric._ric.ric_val_stretches.size,1)
    self.assertEqual(ric._ric.ric_val_in_bends.size,1)
    self.assertEqual(ric._ric.ric_val_out_bends.size,1)
    self.assertEqual(ric._ric.ric_val_lin_bends.size,2)
    self.assertListEqual(list(ric._ric.ric_hessian.shape),[12,12])

    self.assertEqual(ric.num_ric     ,12)
    self.assertEqual(ric.num_stretch ,1)
    self.assertEqual(ric.num_in_bend ,1)
    self.assertEqual(ric.num_out_bend,1)
    self.assertEqual(ric.num_lin_bend,2)
    self.assertEqual(ric.num_torsion ,2)
    self.assertEqual(ric.num_eckart  ,5)


if __name__ == '__main__':
  ut.main()

