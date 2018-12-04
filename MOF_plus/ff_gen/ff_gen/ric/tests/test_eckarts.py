#! /usb/bin/env python

import numpy        as np
import numdifftools as nd
import unittest     as ut

from ric import RedIntCoords as RIC

class EckartTests(ut.TestCase):

  def test_val_0D_1(self):
    """Values without PBC"""

    masses = np.array([1.,2.,1.])

    ric = RIC()
    ric.add_eckart(trans=[True]*3,rots=[False]*3)
    ric.setup(masses)

    for coords, ref in [([[ 0, 0, 0],[ 0, 0, 0],[ 0, 0, 0]],[0,0,0]),
                        ([[-1,-1,-1],[ 0, 0, 0],[ 1, 1, 1]],[0,0,0]),
                        ([[-2, 0, 0],[ 2, 0, 0],[-2, 0, 0]],[0,0,0])]:

      ric.construct_b_matrix(None,np.array(coords,dtype=np.float64))
      res = ric.get_val_eckarts()
      self.assertListEqual(list(res),list(ref))

  def test_grad_0D_1(self):

    masses = np.array([1.,2.,1.])

    ric = RIC()
    ric.add_eckart(trans=[True]*3,rots=[False]*3)
    ric.setup(masses)

    def func(x):
      ric.construct_b_matrix(None,np.reshape(x,(3,3)))
      return np.copy(ric.get_val_eckarts())
    jac = nd.Jacobian(func,step_nom=[0.01]*9)

    for coords in [[[ 0, 0, 0],[ 0, 0, 0],[ 0, 0, 0]],
                   [[-1,-1,-1],[ 0, 0, 0],[ 1, 1, 1]],
                   [[-2, 0, 0],[ 2, 0, 0],[-2, 0, 0]]]:
      coords = np.array(coords,dtype=np.float64)
      ref = jac(coords.flatten())
      res = ric.construct_b_matrix(None,coords)
      #print ref
      #print res
      self.assertLess(np.max(np.abs(ref-res)),1.e-12)


if __name__ == '__main__':
  ut.main()

