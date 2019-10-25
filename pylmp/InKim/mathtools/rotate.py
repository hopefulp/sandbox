import numpy as np
from numpy import sin, cos

def rotate_matrix(tx, ty, tz):
	"""
    Input: three angles (phi, theta, rho)
	Returns the rotation matrix.
	Note that tx, ty, tz are radians.
	"""

	Rx = np.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
	Ry = np.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0, cos(ty)]])
	Rz = np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0,0,1]])

	return np.dot(Rx, np.dot(Ry, Rz))


