#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide.basic import Vector, Point, Normal, Ray, BBox
from geoclide.vecope import normalize
import numpy as np
from numpy.linalg import inv
import math
import warnings


class Transform(object):
    '''
    Tool to perform tranlation(s) and/or rotation(s) to objects

    Parameters
    ----------
    m : Transform | 2-D ndarray | 3-D ndarray, optional
        The matrix of the transformation(s)
    mInv : Transform | 2-D ndarray | 3-D ndarray, optional
        The matrix of the inverse transformation(s)

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
            if ( (len(m.shape) == 2 and m.shape != (4,4)) or
                 (len(m.shape) == 3 and m.shape[1] != 4 and m.shape[2] != 4)):
                raise ValueError("The m parameter must be an np.array of shape (4,4) or (nT,4,4)")
            self.m = m
            self.mInv = inv(m)
        elif (m is None and isinstance(mInv, np.ndarray)):
            if ( (len(mInv.shape) == 2 and mInv.shape != (4,4)) or
                 (len(mInv.shape) == 3 and mInv.shape[1] != 4 and mInv.shape[2] != 4)):
                raise ValueError("The mInv parameter must be an np.array of shape (4,4) or (nT,4,4)")
            self.m = inv(mInv)
            self.mInv = mInv
        elif (isinstance(m, np.ndarray) and isinstance(mInv, np.ndarray)):
            if ( (len(m.shape) == 2 and m.shape != (4,4)) or
                 (len(m.shape) == 3 and m.shape[1] != 4 and m.shape[2] != 4)):
                raise ValueError("The m parameter must be an np.array of shape (4,4) or (nT,4,4)")
            if ( (len(mInv.shape) == 2 and mInv.shape != (4,4)) or
                 (len(mInv.shape) == 3 and mInv.shape[1] != 4 and mInv.shape[2] != 4)):
                raise ValueError("The mInv parameter must be an np.array of shape (4,4) or (nT,4,4)")
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
    
    def __call__(self, c, calc_diag=False):
        """
        Apply the transformations

        Parameters
        ----------
        c : Vector | Point | Normal | Ray | BBox
            The Vector/Point/Normal/Ray/BBox to which the transformation is applied
        diag_calc : bool, optional
            Perform diagonal calculations between c(i) and tranformation(i). The number of 
            transformations must be equal to the number of vectors / points / ... 
        
        Returns
        -------
        out : Vector | Point | Normal | Ray | BBox | 1-D array
            The Vector/Point/Normal/Ray/BBox after the transformation, or in case 
            several transformations are given return a 1-D ndarray where dtype is equal 
            to one of the previously mentionned classes.

        Examples
        --------
        >>> import geoclide as gc
        >>> t = gc.get_translate_tf(gc.Vector(5., 5., 5.))
        >>> p = gc.Point(0., 0., 0.)
        >>> t[p]
        Point(5.0, 5.0, 5.0)
        """
        if len(self.m.shape) == 3: # Case with several transformations in Transform
            is_vector = isinstance(c, Vector)
            is_point = isinstance(c, Point)
            is_normal = isinstance(c, Normal)
            if is_vector or is_point or is_normal:
                nT = self.m.shape[0]
                if is_vector or is_point: mat = np.moveaxis(self.m, 0,2)
                else: mat = np.moveaxis(self.mInv, 0,2) # if is_normal
                is_c_arr = isinstance(c.x, np.ndarray)
                key_bis = np.arange(nT)
                if is_c_arr and not calc_diag:
                    mat = mat[:,:,np.newaxis,:]
                    x = c.x[:,np.newaxis]
                    y = c.y[:,np.newaxis]
                    z = c.z[:,np.newaxis]
                    keys =  [(slice(None), k) for k in key_bis]
                else: # if diag_calc = True or if not is_c_arr
                    x = c.x
                    y = c.y
                    z = c.z
                    keys = key_bis
            if is_vector:
                xv = mat[0,0]*x + mat[0,1]*y + mat[0,2]*z
                yv = mat[1,0]*x + mat[1,1]*y + mat[1,2]*z
                zv = mat[2,0]*x + mat[2,1]*y + mat[2,2]*z
                vectors = np.empty(nT, dtype=Vector)
                for iv in range (0, nT):
                    vectors[iv] = Vector(xv[keys[iv]], yv[keys[iv]], zv[keys[iv]])
                return vectors
            elif is_point:
                xp = mat[0,0]*x + mat[0,1]*y + mat[0,2]*z + mat[0,3]
                yp = mat[1,0]*x + mat[1,1]*y + mat[1,2]*z + mat[1,3]
                zp = mat[2,0]*x + mat[2,1]*y + mat[2,2]*z + mat[2,3]
                wp = mat[3,0]*x + mat[3,1]*y + mat[3,2]*z + mat[3,3]
                points = np.empty(nT, dtype=Point)
                for ip in range (0, nT):
                    if ((not isinstance(wp[keys[ip]], np.ndarray) and wp[keys[ip]] == 1) or
                        (isinstance(wp[keys[ip]], np.ndarray) and np.all(wp[keys[ip]] == 1)) ):
                        points[ip] = Point(xp[keys[ip]], yp[keys[ip]], zp[keys[ip]])
                    else: 
                        points[ip] = Point(xp[keys[ip]], yp[keys[ip]], zp[keys[ip]])/wp[keys[ip]]
                return points
            elif is_normal:
                xn = mat[0,0]*x + mat[1,0]*y + mat[2,0]*z
                yn = mat[0,1]*x + mat[1,1]*y + mat[2,1]*z
                zn = mat[0,2]*x + mat[1,2]*y + mat[2,2]*z
                normals = np.empty(nT, dtype=Vector)
                for inorm in range (0, nT):
                    normals[inorm] = Normal(xv[keys[inorm]], yv[keys[inorm]], zv[keys[inorm]])
                return normals
            elif isinstance(c, Ray):
                origins = self(c.o, calc_diag)
                directions = self(c.d, calc_diag)
                nT = self.m.shape[0]
                rays = np.empty(nT, dtype=Ray)
                for ir in range (0, nT):
                    rays[ir] = Ray(origins[ir], directions[ir], mint=c.mint, maxt=c.maxt)
                return rays
            elif isinstance(c, BBox):
                p0 = self(c.p0, calc_diag)
                v0 = self(c.p1-c.p0, calc_diag)
                v1 = self(c.p3-c.p0, calc_diag)
                v2 = self(c.p4-c.p0, calc_diag)
                nT = self.m.shape[0]
                bboxes = np.empty(nT, dtype=BBox)
                for ib in range (0, nT):
                    b = BBox()
                    b = b.union(p0[ib])
                    b = b.union(p0[ib]+v0[ib])
                    b = b.union(p0[ib]+(v0[ib]+v1[ib]))
                    b = b.union(p0[ib]+v1[ib])
                    b = b.union(p0[ib]+v2[ib])
                    b = b.union(p0[ib]+(v0[ib]+v2[ib]))
                    b = b.union(p0[ib]+(v0[ib]+v1[ib]+v2[ib]))
                    b = b.union(p0[ib]+(v1[ib]+v2[ib]))
                    bboxes[ib] = b
                return bboxes
            else:
                raise ValueError('Unknown type for transformations')
        else:
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
                R.o = self(R.o)
                R.d = self(R.d)
                return R
            elif isinstance(c, BBox):
                b = BBox()
                p0 = self(c.p0)
                v0 = self(c.p1-c.p0)
                v1 = self(c.p3-c.p0)
                v2 = self(c.p4-c.p0)
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
    
    def __getitem__(self, c, diag_calc):
        """
        Apply the transformations

        Parameters
        ----------
        c : Vector | Point | Normal | Ray | BBox
            The Vector/Point/Normal/Ray/BBox to which the transformation is applied
        diag_calc : bool, optional
            Perform diagonal calculations between c(i) and tranformation(i). The number of 
            transformations must be equal to the number of vectors / points / ... 
        
        Returns
        -------
        out : Vector | Point | Normal | Ray | BBox | 1-D array
            The Vector/Point/Normal/Ray/BBox after the transformation, or in case 
            several transformations are given return a 1-D ndarray where dtype is equal 
            to one of the previously mentionned classes.

        Examples
        --------
        >>> import geoclide as gc
        >>> t = gc.get_translate_tf(gc.Vector(5., 5., 5.))
        >>> p = gc.Point(0., 0., 0.)
        >>> t[p]
        Point(5.0, 5.0, 5.0)
        """
        warnings.simplefilter('always', DeprecationWarning)
        warn_message = "\nApplying the transformation through square brackets is deprecated\n" + \
            "as of version 2.1.0 and will be no more possible in the future.\n" + \
            "Please use parenthesis instead."
        warnings.warn(warn_message, DeprecationWarning, stacklevel=1)
        return self(c, diag_calc)

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
        The vector(s) used for the transformation

    Returns
    -------
    t : Transform
        The translate transformation(s)

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
    if isinstance(v.x, np.ndarray):
        nc = len(v.x)
        m = np.tile(np.identity(4, dtype=np.float64), (nc,1)).reshape(nc,4,4)
        mInv = m.copy()
        m[:,0,3] = v.x
        m[:,1,3] = v.y
        m[:,2,3] = v.z
        mInv[:,0,3] = -v.x
        mInv[:,1,3] = -v.y
        mInv[:,2,3] = -v.z
    else:
        m = np.identity(4)
        mInv = m.copy()
        m[0,3] = v.x
        m[1,3] = v.y
        m[2,3] = v.z
        mInv[0,3] = -v.x
        mInv[1,3] = -v.y
        mInv[2,3] = -v.z
    return Transform(m, mInv)


def get_scale_tf(x, y, z):
    """
    Get the scale Transform

    Parameters
    ----------
    x : float | 1-D ndarray
        The scale factor(s) to apply (x axis)
    y : float | 1-D ndarray
        The scale factor(s) to apply (y axis)
    y : float | 1-D ndarray
        The scale factor(s) to apply (z axis)

    Returns
    -------
    t : Transform
        The scale transformation(s)

    Notes
    -----
    x, y and z must be all scalars or all 1-D arrays with the same size
    """
    is_x_arr = isinstance(x, np.ndarray)
    is_y_arr = isinstance(y, np.ndarray)
    is_z_arr = isinstance(z, np.ndarray)
    
    if (is_x_arr and is_y_arr and is_z_arr):
        sx = len(x)
        sy = len(y)
        sz = len(z)
        if (sx != sy or sx != sz):
            raise ValueError("The parameters x, y and z must have the same size")
        nc = sx
        m = np.tile(np.identity(4, dtype=np.float64), (nc,1)).reshape(nc,4,4)
        mInv = m.copy()
        m[:,0,0] = x
        m[:,1,1] = y
        m[:,2,2] = z
        mInv[:,0,0] = 1./x
        mInv[:,1,1] = 1./y
        mInv[:,2,2] = 1./z
    elif(not is_x_arr and not is_y_arr and not is_z_arr):
        m = np.identity(4)
        mInv = m.copy()
        m[0,0] = x
        m[1,1] = y
        m[2,2] = z
        mInv[0,0] = 1./x
        mInv[1,1] = 1./y
        mInv[2,2] = 1./z
    else:
        raise ValueError("The parameters x, y and z must be all scalars " + \
                         "or all 1-D ndarrays")
    return Transform(m, mInv)


def get_rotateX_tf(angle):
    """
    Get the rotateX Transform

    Parameters
    ----------
    angle : float | 1-D ndarray
        The angle in degrees for the rotation around the x axis

    Returns
    -------
    t : Transform
        The rotateX transformation(s)
    """
    is_ang_arr = isinstance(angle, np.ndarray)

    if (is_ang_arr):
        sin_t = np.sin(angle*(math.pi / 180.))
        cos_t = np.cos(angle*(math.pi / 180.))
        nc = len(angle)
        m = np.tile(np.identity(4, dtype=np.float64), (nc,1)).reshape(nc,4,4)
        m[:,1,1] = cos_t
        m[:,1,2] = -1.*sin_t
        m[:,2,1] = sin_t
        m[:,2,2] = cos_t
        return Transform(m, np.transpose(m, axes=(0,2,1)))
    else:
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
    angle : float | 1-D ndarray
        The angle(s) in degrees for the rotation around the y axis

    Returns
    -------
    t : Transform
        The rotateY transformation(s)
    """
    is_ang_arr = isinstance(angle, np.ndarray)

    if (is_ang_arr):
        sin_t = np.sin(angle*(math.pi / 180.))
        cos_t = np.cos(angle*(math.pi / 180.))
        nc = len(angle)
        m = np.tile(np.identity(4, dtype=np.float64), (nc,1)).reshape(nc,4,4)
        m[:,0,0] = cos_t
        m[:,2,0] = -1.*sin_t
        m[:,0,2] = sin_t
        m[:,2,2] = cos_t
        return Transform(m, np.transpose(m, axes=(0,2,1)))
    else:
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
    angle : float | 1-D ndarray
        The angle(s) in degrees for the rotation around the Z axis

    Returns
    -------
    t : Transform
        The rotateZ transformation(s)
    """
    is_ang_arr = isinstance(angle, np.ndarray)

    if (is_ang_arr):
        sin_t = np.sin(angle*(math.pi / 180.))
        cos_t = np.cos(angle*(math.pi / 180.))
        nc = len(angle)
        m = np.tile(np.identity(4, dtype=np.float64), (nc,1)).reshape(nc,4,4)
        m[:,0,0] = cos_t
        m[:,0,1] = -sin_t
        m[:,1,0] = sin_t
        m[:,1,1] = cos_t
        return Transform(m, np.transpose(m, axes=(0,2,1)))
    else:
        sin_t = math.sin(angle*(math.pi / 180.))
        cos_t = math.cos(angle*(math.pi / 180.))
        m = np.identity(4)
        m[0,0] = cos_t
        m[0,1] = -sin_t
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


