#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide import Vector, Point, Normal
from geoclide import normalize, cross
import numpy as np


class Shape(object):
    '''
    The parent class of the different shape classes
    '''
    indShape = 0
    def __init__(self, ObjectToWorld, WorldToObject):
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
    shape : Sphere | Triangle | ...
        The shape used
    '''
    def __init__(self, p=None, dpdu=None, dpdv=None, u=None, v=None, shape = None):
        if ( (p is None) and (dpdu is None) and (dpdv is None) and (u is None) and (v is None) and (shape is None) ):
            self.p = Point(0., 0., 0.)
            self.dpdu = Vector(0., 0., 0.)
            self.dpdv = Vector(0., 0., 0.)
            self.n = Normal(0., 0., 0.)
            self.u = 0.
            self.v = 0.
            self.shape = None
        elif ( isinstance(p, Point)     and
               isinstance(dpdu, Vector) and
               isinstance(dpdv, Vector) and
               np.isscalar(u)           and
               np.isscalar(v) ):
            self.p = p
            self.dpdu = dpdu
            self.dpdv = dpdv
            self.n = Normal(normalize(cross(self.dpdu, self.dpdv)))
            self.u = u
            self.v = v
            self.shape = shape
        else:
            raise ValueError('Problem with parameter(s)')
