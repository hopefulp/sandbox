#! /usb/bin/env python

import unittest as ut
import numpy as np
import numpy.linalg as npl
import numdifftools as nd

from ric import Group, Angle

class AngleTests(ut.TestCase):

  def setUp(self):

    self.grp1 = Group([1])
    self.grp2 = Group([2,3])
    self.grp3 = Group([3,4,5])
    self.grp4 = Group([1,5])

    self.angl1 = Angle([1,2,3], axis=None)
    self.angl2 = Angle([1,2,3], axis='x')
    self.angl3 = Angle([1,2,3], axis='y')
    self.angl4 = Angle([1,2,3], axis='z')
    self.angl5 = Angle([1,2,3], axis=np.array([1,2,3],dtype=np.float64))
    self.angl6 = Angle([1,2,3], axis=4, inplane=True)
    self.angl7 = Angle([1,2,3], axis=4, inplane=False)

    self.angl8 = Angle([self.grp1,self.grp2,self.grp3], axis=None)

  def test_get_axis(self):

    coords = 20*np.random.random((4,3)) - 10

    a = self.angl1.get_axis(coords)
    self.assertAlmostEqual(1, npl.norm(a))

    ref = [1,0,0]
    res = self.angl2.get_axis(coords)
    self.assertTrue(np.allclose(ref, res))

    ref = [0,1,0]
    res = self.angl3.get_axis(coords)
    self.assertTrue(np.allclose(ref, res))

    ref = [0,0,1]
    res = self.angl4.get_axis(coords)
    self.assertTrue(np.allclose(ref, res))

    ref  = [1,2,3]
    ref /= npl.norm(ref)
    res = self.angl5.get_axis(coords)
    self.assertTrue(np.allclose(ref, res))

    a1 = self.angl6.get_axis(coords)
    a2 = self.angl7.get_axis(coords)
    self.assertAlmostEqual(0, np.sum(a1*a2))

  def test_evaluate(self):

    angl = self.angl1
    for coords, ref in [([[ 1, 0, 0],[ 0, 0, 0],[ 1, 0, 0]],   0),
                        ([[ 1, 0, 0],[ 0, 0, 0],[ 1, 1, 0]],  45),
                        ([[ 1, 0, 0],[ 0, 0, 0],[ 0, 1, 0]],  90),
                        ([[ 1, 0, 0],[ 0, 0, 0],[-1, 1, 0]], 135),
                        ([[ 1, 0, 0],[ 0, 0, 0],[-1, 0, 0]], 180),
                        ([[ 2, 2, 2],[ 0, 0, 0],[-1,-1,-1]], 180)]:

      res = angl.evaluate(np.array(coords))
      self.assertAlmostEqual(res,np.radians(ref))

    angl = self.angl4
    for coords, ref in [([[ 1, 0, 0],[ 0, 0, 0],[ 1, 1, 0]], -135),
                        ([[ 1, 0, 0],[ 0, 0, 0],[ 0, 1, 0]],  -90),
                        ([[ 1, 0, 0],[ 0, 0, 0],[-1, 1, 0]],  -45),
                        ([[ 1, 0, 0],[ 0, 0, 0],[-1, 0, 0]],    0),
                        ([[ 1, 0, 0],[ 0, 0, 0],[-1,-1, 0]],   45),
                        ([[ 1, 0, 0],[ 0, 0, 0],[ 0,-1, 0]],   90),
                        ([[ 1, 0, 0],[ 0, 0, 0],[ 1,-1, 0]],  135),
                        ([[ 1, 0, 0],[ 0, 0, 0],[ 1, 0, 0]],  180)]:

      res = angl.evaluate(np.array(coords))
      self.assertAlmostEqual(res,np.radians(ref))

  def test_project(self):

    angl = self.angl1
    for _ in range(1):
      coords = 20*np.random.random((3,3)) - 10
      grad = nd.Gradient(lambda x: angl.evaluate(x.reshape((-1,3))))
      ref = grad(coords.ravel())
      res = angl.project(coords)
      self.assertTrue(np.allclose(ref, res))

    angl = self.angl2
    for _ in range(1):
      coords = 20*np.random.random((3,3)) - 10
      grad = nd.Gradient(lambda x: angl.evaluate(x.reshape((-1,3))))
      ref = grad(coords.ravel())
      res = angl.project(coords)
      self.assertTrue(np.allclose(ref, res))

    angl = self.angl3
    for _ in range(1):
      coords = 20*np.random.random((3,3)) - 10
      grad = nd.Gradient(lambda x: angl.evaluate(x.reshape((-1,3))))
      ref = grad(coords.ravel())
      res = angl.project(coords)
      self.assertTrue(np.allclose(ref, res))

    angl = self.angl4
    for _ in range(1):
      coords = 20*np.random.random((3,3)) - 10
      grad = nd.Gradient(lambda x: angl.evaluate(x.reshape((-1,3))))
      ref = grad(coords.ravel())
      res = angl.project(coords)
      self.assertTrue(np.allclose(ref, res))

    angl = self.angl5
    for _ in range(1):
      coords = 20*np.random.random((3,3)) - 10
      grad = nd.Gradient(lambda x: angl.evaluate(x.reshape((-1,3))))
      ref = grad(coords.ravel())
      res = angl.project(coords)
      self.assertTrue(np.allclose(ref, res))

    angl = self.angl6
    for _ in range(1):
      coords = 20*np.random.random((4,3)) - 10
      angl_ = Angle(angl.points, axis=angl.get_axis(coords))
      grad = nd.Gradient(lambda x: angl_.evaluate(x.reshape((-1,3))))
      ref = grad(coords.ravel())
      res = angl.project(coords)
      self.assertTrue(np.allclose(ref, res))

    angl = self.angl7
    for _ in range(1):
      coords = 20*np.random.random((4,3)) - 10
      angl_ = Angle(angl.points, axis=angl.get_axis(coords))
      grad = nd.Gradient(lambda x: angl_.evaluate(x.reshape((-1,3))))
      ref = grad(coords.ravel())
      res = angl.project(coords)
      self.assertTrue(np.allclose(ref, res))

    angl = self.angl8
    for _ in range(1):
      coords = 20*np.random.random((5,3)) - 10
      grad = nd.Gradient(lambda x: angl.evaluate(x.reshape((-1,3))))
      ref = grad(coords.ravel())
      res = angl.project(coords)
      self.assertTrue(np.allclose(ref, res))

if __name__ == '__main__':
  ut.main()

