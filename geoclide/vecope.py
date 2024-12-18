#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide.basic import Vector, Point, Normal
import math
import numpy as np


def dot(a, b):
    if (isinstance(a, Vector) or isinstance(a, Normal)) and \
       (isinstance(b, Vector) or isinstance(b, Normal)):
        return (a.x*b.x + a.y*b.y + a.z*b.z)
    else:
        raise NameError('The dot parameters have to be Vector or Normal classes')


def cross(a, b):
    if (isinstance(a, Vector) and isinstance(b, Vector)) or \
    (isinstance(a, Vector) and isinstance(b, Normal)) or \
    (isinstance(a, Normal) and isinstance(b, Vector)):
        return Vector((a.y*b.z)-(a.z*b.y), (a.z*b.x)-(a.x*b.z), (a.x*b.y)-(a.y*b.x))
    elif isinstance(a, Normal) and isinstance(b, Normal):
        raise NameError('Only 1 Normal is tolerate not 2')
    else:
        raise NameError('The cross parameters must be Vector or Normal')
    

def normalize(v):
    if isinstance(v, Vector) or isinstance(v, Normal):
        return v / v.length()
    else:
        raise NameError('Normalize argument have to be Vector or Normal class')


def coordinateSystem(v1):
    """
    Create an orthogonal coordinate system from 1 vector (v1)
    """
    if (abs(v1.x) > abs(v1.y)):
        invLen = 1/ math.sqrt(v1.x*v1.x + v1.z*v1.z)
        v2 = Vector(-v1.z*invLen, 0, v1.x*invLen)
    else:
        invLen = 1/ math.sqrt(v1.y*v1.y + v1.z*v1.z)
        v2 = Vector(0, v1.z*invLen, -v1.y*invLen)
    v3 = cross(v1, v2)
    return v2, v3


def distance(p1, p2):
    if isinstance(p1, Point) and isinstance(p2, Point):
        return (p1 - p2).length()
    else:
        raise NameError('The distance parameters have to be Point classes)')


def faceForward(a, b):
    """
    Flip the Vector/Normal a if the Vector/Normal b is in the opposite direction.
    For exemple, it can be useful to flip a surface normal so that it lies in the
    same hemisphere as a given vector.
    Args : Vector or Normal a, b
    Output : Possibly fliped Vector or Normal a
    """
    if (isinstance(a, Vector) or isinstance(a, Normal)) and \
    (isinstance(b, Vector) or isinstance(b, Normal)):
        return (a*-1) if (dot(a, b) < 0) else a
    else:
        raise NameError('FaceForward args have to be Vector or Normal classes')


def rotation3D(theta, u):
    """
    rotation matrix of an angle theta in degree around unit vector u
    """
    # Rodrigues rotation formula
    ct = math.cos(math.radians(theta))
    st = math.sqrt(1-ct*ct)
    matA = np.zeros((3,3))
    for k in range(3) : matA[k,k] = 1.
    matB = np.zeros((3,3))
    matB[0,1] = -u.z
    matB[0,2] =  u.y
    matB[1,0] =  u.z
    matB[1,2] = -u.x
    matB[2,0] = -u.y
    matB[2,1] =  u.x
    matC = matB.dot(matB)
    return matA + matB*st + matC*(1-ct)


def maxDimension(v):
    if isinstance(v, Vector):
        val = 0
        if (v.x > v.y):
            if(v.x > v.z):
                val = 0
            else:
                val = 2
        else:
            if(v.y > v.z):
                val = 1
            else:
                val = 2
        return int(val)
    else:
        raise NameError('v argument must be a Vector')


def permute(v, xx, yy, zz):
    if isinstance(v, Vector):
        return Vector(v[xx], v[yy], v[zz])
    else:
        raise NameError('v argument must be a Vector')