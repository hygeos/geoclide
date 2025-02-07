#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide.basic import Vector, Point, Normal
from geoclide.vecope import normalize, cross, face_forward
import numpy as np
from geoclide.transform import Transform


class Shape(object):
    '''
    The parent class of the different shape classes
    '''
    indShape = 0
    def __init__(self, ObjectToWorld, WorldToObject):
        if (not isinstance(ObjectToWorld, Transform) or not isinstance(WorldToObject, Transform)):
            raise ValueError('The parameters oTw and wTo must be both Transform')
        self.oTw = ObjectToWorld
        self.wTo = WorldToObject


class DifferentialGeometry(object):
    '''
    The general parametric description of most of the shapes

    - A point p is described by a function depending to variables u and v such that p=f(u,v)
    - From u and v we can get two directions (partial derivative of p) parallel to the shape surface

    Parameters
    ----------
    p : Point, optional
        The concerned position, at the surface of the given shape
    dpdu : Vector, optional
        The surface partial derivative of p with respect to u
    dpdv : vector, optional
        The surface partial derivative of p with respect to v
    u : float, optional
        The u coordinate of the parametric representation
    v : float, optional
        The v coordinate of the parametric representation
    ray_dir : float, optional
        The direction of the ray hitting the shape
    shape : Sphere | Triangle | ...
        The shape used
    '''
    def __init__(self, p=None, dpdu=None, dpdv=None, u=None, v=None, ray_dir=None, shape = None):
        if ( (p is None) and (dpdu is None) and (dpdv is None) and (u is None) and
             (v is None) and (ray_dir is None) and (shape is None) ):
            self.p = Point(0., 0., 0.)
            self.dpdu = Vector(0., 0., 0.)
            self.dpdv = Vector(0., 0., 0.)
            self.n = Normal(0., 0., 0.)
            self.u = 0.
            self.v = 0.
            self.ray_dir = Vector(0., 0., 0.)
            self.shape = None
        elif ( isinstance(p, Point)        and
               isinstance(dpdu, Vector)    and
               isinstance(dpdv, Vector)    and
               np.isscalar(u)              and
               np.isscalar(v)              and
               isinstance(ray_dir, Vector)):
            self.p = p
            self.dpdu = dpdu
            self.dpdv = dpdv
            self.ray_dir = ray_dir
            self.n = face_forward(Normal(normalize(cross(dpdu, dpdv))), -ray_dir)
            self.u = u
            self.v = v
            self.shape = shape
        else:
            raise ValueError('Problem with parameter(s)')
