#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide.basic import Vector
from geoclide.vecope import normalize
from geoclide.transform import get_rotateY_tf, get_rotateZ_tf


def ang2vec(theta, phi, vec_view='zenith'):
    """
    Convert a direction described by 2 angles into a direction described by a vector

    - coordinate system convention:
    >>>  y
    >>>  ^   x : right; y : front; z : top
    >>>  |
    >>> z X -- > x

    Parameters
    ----------
    Theta : float
        The polar angle in degrees, starting at z+ in the zx plane and going 
        in the trigonometric direction around the y axis

    Phi : float
        The azimuthal angle in degrees, starting at x+ in the xy plane and going in 
        the trigonométric direction around the z axis

    vec_view : str, optional
        Two choices (concerning intial direction at theta=phi=0): 'zenith' (i.e. pointing above) or 
        'bellow' (i.e. pointing bellow)

    Returns
    -------
    v : Vector
        The direction described by a vector
    """
    if (vec_view == "zenith"): # initial vector is facing zenith (pointing above)
        v = Vector(0., 0., 1.)
    elif (vec_view == "nadir"): # initial vector is facing nadir (pointing bellow)
        v = Vector(0., 0., -1.)
    else:
        raise ValueError("The value of vec_view parameter must be: 'zenith' or 'nadir")
    
    v = get_rotateY_tf(theta)[v]
    v = get_rotateZ_tf(phi)[v]
    v = normalize(v)
    
    return v
