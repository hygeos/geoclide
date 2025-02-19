#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide.basic import Vector, Point, Normal, Ray, BBox
from geoclide.vecope import normalize
import numpy as np
from numpy.linalg import inv
import math


class Transform(object):
    '''
    Tool to perform tranlation(s) and/or rotation(s) to objects

    Parameters
    ----------
    m : Transform | 2-D ndarray, optional
        The transformation matrix
    mInv : Transform | 2-D ndarray, optional
        The inverse transformation matrix

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
        if (not isinstance(t, Transform)):
            raise ValueError("Equality with a Transform must be only with another Transform")
        
        self.m = t.m
        self.mInv = t.mInv
        
    def __mul__(self, t): 
        if (not isinstance(t, Transform)):
            raise ValueError('A transform can be multiplied only by another Transform')
        
        return Transform(np.dot(self.m, t.m), np.dot(t.mInv, self.mInv))
    
    def __getitem__(self, c):
        """"
        Apply the transformations

        Parameters
        ----------
        c : Vector | Point | Normal | Ray | BBox
            The Vector/Point/Normal/Ray/BBox to which the transformation is applied
        
        Results
        -------
        out : Vector | Point | Normal | Ray | BBox
            The Vector/Point/Normal/Ray/BBox after the transformation

        Examples
        --------
        >>> import geoclide as gc
        >>> t = gc.get_translate_tf(gc.Vector(5., 5., 5.))
        >>> p = gc.Point(0., 0., 0.)
        >>> t[p]
        Point(5.0, 5.0, 5.0)
        """
        if isinstance(c, Vector):
            xv = self.m[0,0]*c.x + self.m[0,1]*c.y + self.m[0,2]*c.z
            yv = self.m[1,0]*c.x + self.m[1,1]*c.y + self.m[1,2]*c.z
            zv = self.m[2,0]*c.x + self.m[2,1]*c.y + self.m[2,2]*c.z
            return Vector(xv, yv, zv)
        elif isinstance(c, Point):
            xp = self.m[0,0]*c.x + self.m[0,1]*c.y + self.m[0,2]*c.z + self.m[0,3]
            yp = self.m[1,0]*c.x + self.m[1,1]*c.y + self.m[1,2]*c.z + self.m[1,3]
            zp = self.m[2,0]*c.x + self.m[2,1]*c.y + self.m[2,2]*c.z + self.m[2,3]
            wp = self.m[3,0]*c.x + self.m[3,1]*c.y + self.m[3,2]*c.z + self.m[3,3]
            if ((not isinstance(wp, np.ndarray) and wp == 1) or
                (isinstance(wp, np.ndarray) and np.all(wp == 1)) ):
                return Point(xp, yp, zp)
            else: 
                return Point(xp, yp, zp)/wp
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
            raise ValueError('Unknown type for transformations')

    def __str__(self):
        print("m=\n", self.m, "\nmInv=\n", self.mInv)
        return ""
    
    def __repr__(self):
        print("m=\narray(\n", self.m, ")\nmInv=\narray(\n",self.mInv,")")
        return ""

    def inverse(self):
        """
        Inverse the initial transformation

        Parameters
        ----------
        t : Transform
            The transformation to be inversed

        Returns
        -------
        out : Transform
            The inversed transformation
        """
        return get_inverse_tf(self)

    def is_identity(self):
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


def get_inverse_tf(t):
    """
    Get the inverse transformation

    Parameters
    ----------
    t : Transform
        The transformation to be inversed

    Returns
    -------
    out : Transform
        The inversed transformation
    """
    return Transform(t.mInv, t.m)
    

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
    if (not isinstance(v, Vector)):
        raise ValueError("The parameter v must be a Vector")
    
    m = np.identity(4)
    m[0,3] = v.x
    m[1,3] = v.y
    m[2,3] = v.z
    mInv = np.identity(4)
    mInv[0,3] = (v.x)*-1
    mInv[1,3] = (v.y)*-1
    mInv[2,3] = (v.z)*-1
    return Transform(m, mInv)


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
    if ( (not np.isscalar(x)) or
         (not np.isscalar(y)) or
         (not np.isscalar(z)) ):
        raise ValueError("The parameters x, y and z must be all scalars")
    
    m = np.identity(4)
    m[0,0] = x
    m[1,1] = y
    m[2,2] = z
    mInv = np.identity(4)
    mInv[0,0] = 1./x
    mInv[1,1] = 1./y
    mInv[2,2] = 1./z
    return Transform(m, mInv)


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
    if (not np.isscalar(angle)):
        raise ValueError("The parameter angle must be a scalar")
    
    sin_t = math.sin(angle*(math.pi / 180.))
    cos_t = math.cos(angle*(math.pi / 180.))
    m = np.identity(4)
    m[1,1] = cos_t
    m[1,2] = -1.*sin_t
    m[2,1] = sin_t
    m[2,2] = cos_t
    return Transform(m, np.transpose(m))


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
    if (not np.isscalar(angle)):
        raise ValueError("The parameter angle must be a scalar")
    
    sin_t = math.sin(angle*(math.pi / 180.))
    cos_t = math.cos(angle*(math.pi / 180.))
    m = np.identity(4)
    m[0,0] = cos_t
    m[2,0] = -1.*sin_t
    m[0,2] = sin_t
    m[2,2] = cos_t
    return Transform(m, np.transpose(m))


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
    if (not np.isscalar(angle)):
        raise ValueError("The parameter angle must be a scalar")
    
    sin_t = math.sin(angle*(math.pi / 180.))
    cos_t = math.cos(angle*(math.pi / 180.))
    m = np.identity(4)
    m[0,0] = cos_t
    m[0,1] = -1.*sin_t
    m[1,0] = sin_t
    m[1,1] = cos_t
    return Transform(m, np.transpose(m))


def get_rotate_tf(angle, axis):
    """
    Get the rotate Transform around a given axis

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
    if (not np.isscalar(angle)):
        raise ValueError("The parameter angle must be a scalar")
    if ( (not isinstance(axis, Vector)) and
         (not isinstance(axis, Normal)) ):
        raise ValueError("The parameter axis must be a Vector or a Normal")

    a = Vector(normalize(axis))
    s = math.sin(angle*(math.pi / 180.))
    c = math.cos(angle*(math.pi / 180.))
    m = np.identity(4)

    m[0,0] = a.x*a.x+(1-a.x*a.x)*c
    m[0,1] = a.x*a.y*(1-c)-a.z*s
    m[0,2] = a.x*a.z*(1-c)+a.y*s

    m[1,0] = a.x*a.y*(1-c)+a.z*s
    m[1,1] = a.y*a.y+(1-a.y*a.y)*c
    m[1,2] = a.y*a.z*(1-c)-a.x*s

    m[2,0] = a.x*a.z*(1-c)-a.y*s
    m[2,1] = a.y*a.z*(1-c)+a.x*s
    m[2,2] = a.z*a.z+(1-a.z*a.z)*c
    return Transform(m, np.transpose(m))


