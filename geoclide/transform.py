#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide import Vector, Point, Normal, Ray, BBox
from geoclide import normalize
import numpy as np
from numpy.linalg import inv
import math


class Transform(object):
    '''
    Tool to perform tranlation(s) and/or rotation(s) to objects

    Parameters
    ----------
    m : Transform | np.ndarray, optional
    mInv : np.ndarray, optional

    Exemples
    --------
    >>> import geoclide as gc
    >>> t1 = gc.Transform()
    >>> t1
    m=
    array(
    [[1. 0. 0. 0.]
    [0. 1. 0. 0.]
    [0. 0. 1. 0.]
    [0. 0. 0. 1.]] )
    mInv=
    array(
    [[1. 0. 0. 0.]
    [0. 1. 0. 0.]
    [0. 0. 1. 0.]
    [0. 0. 0. 1.]] )
    '''

    def __init__(self, m = None, mInv = None):
        if (isinstance(m, Transform)):
            self.m = m.m
            self.mInv = m.mInv
        elif (m is None and mInv is None):
            self.m = np.identity(4)
            self.mInv = self.m.copy()
        elif (isinstance(m, np.ndarray) and mInv is None):
            if (m.shape != (4,4)):
                raise ValueError("The m parameter must be an np.array of shape (4,4)")
            self.m = m
            self.mInv = inv(m)
        elif (m is None and isinstance(mInv, np.ndarray)):
            if (mInv.shape != (4,4)):
                raise ValueError("The mInv parameter must be an np.array of shape (4,4)")
            self.m = inv(m)
            self.mInv = mInv
        elif (isinstance(m, np.ndarray) and isinstance(mInv, np.ndarray)):
            if (m.shape != (4,4) or mInv.shape != (4,4)):
                raise ValueError("The matrix shape of m and mInv must be (4,4)")
            self.m = m
            self.mInv = mInv
        else:
            raise ValueError("Wrong parameter value(s) for Transform")

    def __eq__(self, t):
        if isinstance(t, Transform):
            self.m = t.m
            self.mInv = t.mInv
        else:
            raise ValueError("Equality with a Transfor must be only with another Transfor")
        
    def __mul__(self, T): 
        if (isinstance(T, Transform)):
            return Transform(np.dot(self.m, T.m), np.dot(T.mInv, self.mInv))
        else:
            raise NameError('A transform can be multiplied only by another Transform')
    
    def __getitem__(self, c):
        """"
        Apply the transformations to: Point, Vector, Normal, Ray or BBox
        """
        if isinstance(c, Point):
            xp = self.m[0,0]*c.x + self.m[0,1]*c.y + self.m[0,2]*c.z + self.m[0,3]
            yp = self.m[1,0]*c.x + self.m[1,1]*c.y + self.m[1,2]*c.z + self.m[1,3]
            zp = self.m[2,0]*c.x + self.m[2,1]*c.y + self.m[2,2]*c.z + self.m[2,3]
            wp = self.m[3,0]*c.x + self.m[3,1]*c.y + self.m[3,2]*c.z + self.m[3,3]
            if (wp == 1):
                return Point(xp, yp, zp)
            else: 
                return Point(xp, yp, zp)/wp
        elif isinstance(c, Vector):
            xv = self.m[0,0]*c.x + self.m[0,1]*c.y + self.m[0,2]*c.z
            yv = self.m[1,0]*c.x + self.m[1,1]*c.y + self.m[1,2]*c.z
            zv = self.m[2,0]*c.x + self.m[2,1]*c.y + self.m[2,2]*c.z
            return Vector(xv, yv, zv)
        elif isinstance(c, Normal):
            xn = self.mInv[0,0]*c.x + self.mInv[1,0]*c.y + self.mInv[2,0]*c.z
            yn = self.mInv[0,1]*c.x + self.mInv[1,1]*c.y + self.mInv[2,1]*c.z
            zn = self.mInv[0,2]*c.x + self.mInv[1,2]*c.y + self.mInv[2,2]*c.z
            return Normal(xn, yn, zn)
        elif isinstance(c, Ray):
            R = Ray(c.o, c.d)
            R.o = self[R.o]
            R.d = self[R.d]
            return R
        elif isinstance(c, BBox):
            b = BBox()
            p0 = self[c.p0]
            v0 = self[c.p1-c.p0]
            v1 = self[c.p3-c.p0]
            v2 = self[c.p4-c.p0]
            b = b.union(p0)
            b = b.union(p0+v0)
            b = b.union(p0+(v0+v1))
            b = b.union(p0+v1)
            b = b.union(p0+v2)
            b = b.union(p0+(v0+v2))
            b = b.union(p0+(v0+v1+v2))
            b = b.union(p0+(v1+v2))
            return b
        else:
            raise NameError('Unknown type for transformations')

    def __str__(self):
        print("m=\n", self.m, "\nmInv=\n", self.mInv)
        return ""
    
    def __repr__(self):
        print("m=\narray(\n", self.m, ")\nmInv=\narray(\n",self.mInv,")")
        return ""

    def inverse(self, t):
        return Transform(t.mInv, t.m)

    def isIdentity(self):
        return (self.m[0,0] == 1) and (self.m[0,1] == 0) and (self.m[0,2] == 0) and \
            (self.m[0,3] == 0) and (self.m[1,0] == 0) and (self.m[1,1] == 1) and \
            (self.m[1,2] == 0) and (self.m[1,3] == 0) and (self.m[2,0] == 0) and \
            (self.m[2,1] == 0) and (self.m[2,2] == 1) and (self.m[2,3] == 0) and \
            (self.m[3,0] == 0) and (self.m[3,1] == 0) and (self.m[3,2] == 0) and \
            (self.m[3,3] == 1)

    def translate(self, v):
        """
        Apply translate to initial transformation
        
        Parameters
        ----------
        v : Vector
            The vector used for the transformation

        Returns
        -------
        t : Transform
            The product of the initial transformation and the translate transformation

        examples
        --------
        >>> import geoclide as gc
        >>> t = Transform()
        >>> t = t.translate(gc.Vector(5.,0.,0.))
        >>> t
        m=
        array(
        [[1. 0. 0. 5.]
        [0. 1. 0. 0.]
        [0. 0. 1. 0.]
        [0. 0. 0. 1.]] )
        mInv=
        array(
        [[ 1.  0.  0. -5.]
        [ 0.  1.  0.  0.]
        [ 0.  0.  1.  0.]
        [ 0.  0.  0.  1.]] )
        """
        t = get_translate_tf(v)
        return self*t

    def scale(self, x, y, z):
        """
        Apply scale to initial transformation

        Parameters
        ----------
        x : float
            The scale factor to apply (x axis)
        y : float
            The scale factor to apply (y axis)
        y : float
            The scale factor to apply (z axis)

        Returns
        -------
        t : Transform
            The product of the initial transformation and the scale transformation
        """
        t = get_scale_tf(x,y,z)
        return self*t

    def rotateX(self, angle):
        """
        Apply rotateX to initial transformation

        Parameters
        ----------
        angle : float
            The angle in degrees for the rotation around the x axis

        Returns
        -------
        t : Transform
            The product of the initial transformation and the rotateX transformation
        """
        t = get_rotateX_tf(angle)
        return self*t

    def rotateY(self, angle):
        """
        Apply rotateY to initial transformation

        Parameters
        ----------
        angle : float
            The angle in degrees for the rotation around the y axis

        Returns
        -------
        t : Transform
            The product of the initial transformation and the rotateY transformation
        """
        t = get_rotateY_tf(angle)
        return self*t

    def rotateZ(self, angle):
        """
        Apply rotateZ to initial transformation

        Parameters
        ----------
        v : Vector
            The angle in degrees for the rotation around the Z axis

        Returns
        -------
        t : Transform
            The product of the initial transformation and the rotateZ transformation
        """
        t = get_rotateZ_tf(angle)
        return self*t

    def rotate(self, angle, axis):
        """
        Apply rotate to initial transformation

        Parameters
        ----------
        angle : float
            The angle in degrees for the rotation
        axis : Vector | Normal
            The rotation is performed arount the parameter axis 

        Returns
        -------
        t : Transform
            The product of the initial transformation and the rotate transformation
        """
        t = get_rotate_tf(angle, axis)
        return self*t


def get_translate_tf(v):
    """
    Get the translate Transform

    Parameters
    ----------
    v : Vector
        The vector used for the transformation

    Returns
    -------
    t : Transform
        The translate transformation

    examples
    --------
    >>> import geoclide as gc
    >>> t = gc.get_translate_tf(gc.Vector(5.,0.,0.))
    >>> t
    m=
    array(
    [[1. 0. 0. 5.]
    [0. 1. 0. 0.]
    [0. 0. 1. 0.]
    [0. 0. 0. 1.]] )
    mInv=
    array(
    [[ 1.  0.  0. -5.]
    [ 0.  1.  0. -0.]
    [ 0.  0.  1. -0.]
    [ 0.  0.  0.  1.]] )
    """
    mat = np.identity(4)
    mat[0,3] = v.x
    mat[1,3] = v.y
    mat[2,3] = v.z
    matInv = np.identity(4)
    matInv[0,3] = (v.x)*-1
    matInv[1,3] = (v.y)*-1
    matInv[2,3] = (v.z)*-1
    return Transform(mat, matInv)

def get_scale_tf(x, y, z):
    """
    Get the scale Transform

    Parameters
    ----------
    x : float
        The scale factor to apply (x axis)
    y : float
        The scale factor to apply (y axis)
    y : float
        The scale factor to apply (z axis)

    Returns
    -------
    t : Transform
        The scale transformation
    """
    mat = np.identity(4)
    mat[0,0] = x
    mat[1,1] = y
    mat[2,2] = z
    matInv = np.identity(4)
    matInv[0,0] = 1./x
    matInv[1,1] = 1./y
    matInv[2,2] = 1./z
    return Transform(mat, matInv)

def get_rotateX_tf(angle):
    """
    Get the rotateX Transform

    Parameters
    ----------
    angle : float
        The angle in degrees for the rotation around the x axis

    Returns
    -------
    t : Transform
        The rotateX transformation
    """
    sin_t = math.sin(angle*(np.pi / 180.))
    cos_t = math.cos(angle*(np.pi / 180.))
    myM = np.identity(4)
    myM[1,1] = cos_t
    myM[1,2] = -1.*sin_t
    myM[2,1] = sin_t
    myM[2,2] = cos_t
    return Transform(myM, np.transpose(myM))

def get_rotateY_tf(angle):
    """
    Get the rotateY Transform

    Parameters
    ----------
    angle : float
        The angle in degrees for the rotation around the y axis

    Returns
    -------
    t : Transform
        The rotateY transformation
    """
    sin_t = math.sin(angle*(np.pi / 180.))
    cos_t = math.cos(angle*(np.pi / 180.))
    mat = np.identity(4)
    mat[0,0] = cos_t
    mat[2,0] = -1.*sin_t
    mat[0,2] = sin_t
    mat[2,2] = cos_t
    return Transform(mat, np.transpose(mat))

def get_rotateZ_tf(angle):
    """
    Get the rotateZ Transform

    Parameters
    ----------
    v : Vector
        The angle in degrees for the rotation around the Z axis

    Returns
    -------
    t : Transform
        The rotateZ transformation
    """
    sin_t = math.sin(angle*(np.pi / 180.))
    cos_t = math.cos(angle*(np.pi / 180.))
    mat = np.identity(4)
    mat[0,0] = cos_t
    mat[0,1] = -1.*sin_t
    mat[1,0] = sin_t
    mat[1,1] = cos_t
    return Transform(mat, np.transpose(mat))

def get_rotate_tf(angle, axis):
    """
    Get the rotate Transform

    Parameters
    ----------
    angle : float
        The angle in degrees for the rotation
    axis : Vector | Normal
        The rotation is performed around the Vector/Normal axis 

    Returns
    -------
    t : Transform
        The rotate transformation
    """
    a = Vector(normalize(axis))
    s = math.sin(angle*(np.pi / 180.))
    c = math.cos(angle*(np.pi / 180.))
    mat = np.identity(4)

    mat[0,0] = a.x*a.x+(1-a.x*a.x)*c
    mat[0,1] = a.x*a.y*(1-c)-a.z*s
    mat[0,2] = a.x*a.z*(1-c)+a.y*s

    mat[1,0] = a.x*a.y*(1-c)+a.z*s
    mat[1,1] = a.y*a.y+(1-a.y*a.y)*c
    mat[1,2] = a.y*a.z*(1-c)-a.x*s

    mat[2,0] = a.x*a.z*(1-c)-a.y*s
    mat[2,1] = a.y*a.z*(1-c)+a.x*s
    mat[2,2] = a.z*a.z+(1-a.z*a.z)*c
    return Transform(mat, np.transpose(mat))


