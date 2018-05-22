#! /usb/bin/env python

import numpy    as np
import unittest as ut

from ric import RedIntCoords as RIC

class TorsionTests(ut.TestCase):

  def test_values_1(self):

    masses = np.array([1.]*4)
    hmat   = np.array([[0.,0.,0.]]*3)

    for coords, ref in [([[ 0,-1, 0],[0,0,0],[0,0,1],[1,0,1]], -90),
                        ([[ 1,-1, 0],[0,0,0],[0,0,1],[1,0,1]], -45),
                        ([[ 1, 0, 0],[0,0,0],[0,0,1],[1,0,1]],   0),
                        ([[ 1, 0,-1],[0,0,0],[0,0,1],[1,0,2]],   0),
                        ([[ 1, 1, 0],[0,0,0],[0,0,1],[1,0,1]],  45),
                        ([[ 2, 2, 0],[0,0,0],[0,0,1],[3,0,1]],  45),
                        ([[ 0, 1, 0],[0,0,0],[0,0,1],[1,0,1]],  90),
                        ([[-1, 0, 0],[0,0,0],[0,0,1],[1,0,1]], 180)]:
      ric = RIC()
      ric.add_torsion([1,2,3,4])
      ric.setup(masses)
      ric.construct_b_matrix(hmat,np.array(coords,dtype=np.float64))
      res = ric.get_val_torsions()[0]
      self.assertAlmostEqual(res,np.radians(ref))
      del ric

  def test_values_2(self):

    masses = np.array([1.]*11)
    hmat   = np.array([[0.,0.,0.]]*3)
    coords = np.array([[ 1., 0., 0.],  # 1:    0
                       [ 1., 1.,-1.],  # 2:   45
                       [ 0., 1.,-2.],  # 3:   90
                       [-1., 1.,-3.],  # 4:  135
                       [-1., 0.,-4.],  # 5:  180
                       [ 0., 0., 0.],  # central
                       [ 0., 0., 1.],  # central
                       [ 1., 0., 1.],  # 1:    0
                       [ 1.,-1., 2.],  # 2:  -45
                       [ 0.,-1., 3.],  # 3:  -90
                       [-1.,-1., 4.]]) # 4: -135

    for ivals, ref in [([1,1],   0),([1,2],  45),([1,3],  90),([1,4], 135),
                       ([2,1],  45),([2,2],  90),([2,3], 135),([2,4], 180),
                       ([3,1],  90),([3,2], 135),([3,3], 180),([3,4],-135),
                       ([4,1], 135),([4,2], 180),([4,3],-135),([4,4], -90),
                       ([5,1], 180),([5,2],-135),([5,3], -90),([5,4], -45)]:
      ric = RIC()
      ric.add_torsion([1,2,3,4,5,6,7,8,9,10,11,0],ivals=ivals)
      ric.setup(masses)
      ric.construct_b_matrix(hmat,coords)
      res = ric.get_val_torsions()[0]
      self.assertAlmostEqual(res,np.radians(ref))
      del ric

if __name__ == '__main__':
  ut.main()

