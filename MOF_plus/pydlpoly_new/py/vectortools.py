# -*- coding: utf-8 -*-

from numpy import *


def random_quat():
    rand = random.random(3)
    r1 = sqrt(1.0 - rand[0])
    r2 = sqrt(rand[0])
    pi2 = pi * 2.0
    t1 = pi2 * rand[1]
    t2 = pi2 * rand[2]
    return array([cos(t2)*r2, sin(t1)*r1, cos(t1)*r1, sin(t2)*r2])

def quat_to_mat(quat):
    q = array(quat, copy=True)
    n = dot(q, q)
    if n < 1.0e-15:
        return identity(3)
    q *= math.sqrt(2.0 / n)
    q = outer(q, q)
    return array([
        [1.0-q[2, 2]-q[3, 3],     q[1, 2]-q[3, 0],     q[1, 3]+q[2, 0]],
        [    q[1, 2]+q[3, 0], 1.0-q[1, 1]-q[3, 3],     q[2, 3]-q[1, 0]],
        [    q[1, 3]-q[2, 0],     q[2, 3]+q[1, 0], 1.0-q[1, 1]-q[2, 2]]])

def apply_mat(m,v):
    return dot(v,m)

def rotate_random(v):
    return apply_mat(quat_to_mat(random_quat()),v)
    
    