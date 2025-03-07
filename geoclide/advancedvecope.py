#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide.basic import Vector
from geoclide.vecope import normalize
from geoclide.transform import get_rotateY_tf, get_rotateZ_tf
from geoclide.mathope import clamp
import numpy as np
import math


def ang2vec(theta, phi, vec_view='zenith'):
    """
    Convert a direction described by 2 angles into a direction described by a vector

    - coordinate system convention:
    >>>   y
    >>>   ^   x : right; y : front; z : top
    >>>   |
    >>> z X -- > x

    Parameters
    ----------
    theta : float
        The polar angle in degrees, starting at z+ in the zx plane and going 
        in the trigonometric direction around the y axis

    phi : float
        The azimuthal angle in degrees, starting at x+ in the xy plane and going in 
        the trigonométric direction around the z axis

    vec_view : str, optional
        Two choices (concerning intial direction at theta=phi=0): 'zenith' (i.e. pointing above) or 
        'bellow' (i.e. pointing bellow)

    Returns
    -------
    v : Vector
        The direction described by a vector
    
    Examples
    --------
    >>> import geoclide as gc
    >>> th = 30.
    >>> ph = 0.
    >>> v1 = gc.ang2vec(theta=th, phi=ph, vec_view='zenith')
    >>> v1
    Vector(0.49999999999999994, 0.0, 0.8660254037844387)
    >>> v2 = gc.ang2vec(theta=th, phi=ph, vec_view='nadir')
    >>> v2
    Vector(-0.49999999999999994, 0.0, -0.8660254037844387)
    """
    if (vec_view == "zenith"): # initial vector is facing zenith (pointing above)
        v = Vector(0., 0., 1.)
    elif (vec_view == "nadir"): # initial vector is facing nadir (pointing bellow)
        v = Vector(0., 0., -1.)
    else:
        raise ValueError("The value of vec_view parameter must be: 'zenith' or 'nadir")
    
    v = get_rotateY_tf(theta)(v)
    v = get_rotateZ_tf(phi)(v)
    v = normalize(v)
    
    return v


def vec2ang(v, vec_view='zenith'):
    """
    Convert a direction described by a vector into a direction described by 2 angles

    - coordinate system convention:
    >>>   y
    >>>   ^   x : right; y : front; z : top
    >>>   |
    >>> z X -- > x

    Parameters
    ----------
    v : Vector
        The direction described by a vector

    vec_view : str, optional
        Two choices (concerning intial direction at theta=phi=0): 'zenith' (i.e. pointing above) or 
        'nadir' (i.e. pointing bellow)

    Returns
    -------
    theta : float
        The polar angle in degrees, starting at z+ in the zx plane and going 
        in the trigonometric direction around the y axis

    phi : float
        The azimuthal angle in degrees, starting at x+ in the xy plane and going in 
        the trigonométric direction around the z axis

    Examples
    --------
    >>> import geoclide as gc
    >>> th = 30.
    >>> ph = 0.
    >>> v1 = gc.ang2vec(theta=th, phi=ph, vec_view='zenith')
    >>> v1
    Vector(0.49999999999999994, 0.0, 0.8660254037844387)
    >>> theta, phi = gc.vec2ang(v1, vec_view='zenith')
    >>> theta, phi
    (29.999999999999993, 0.0)
    >>> v2 = gc.ang2vec(theta=th, phi=ph, vec_view='nadir')
    >>> v2
    Vector(-0.49999999999999994, 0.0, -0.8660254037844387)
    >>> theta, phi = gc.vec2ang(v1, vec_view='nadir')
    >>> theta, phi
    (29.999999999999993, 0.0)
    """
    if (not isinstance(v, Vector)):
        raise ValueError('The parameter v must be a vector')
    
    if (vec_view == "zenith"): # initial vector is facing zenith (pointing above)
        pass
    elif (vec_view == "nadir"): # initial vector is facing nadir (pointing bellow)
        v = -v # by doing that, we can keep a v_ini facing upwards
    else:
        raise ValueError("The value of vec_view parameter must be: 'zenith' or 'nadir")
    
    v = normalize(v) # ensure v is normalized
    v_ini = Vector(0., 0., 1.)

    # In case v = v_ini -> no rotations
    if (np.all(np.isclose(v.to_numpy()-v_ini.to_numpy(), 0., 0., 1e-14))):
        return 0., 0.
    
    for icase in range (1, 6):
        if icase == 1:
            roty_rad = math.acos(v.z)
            if (v.x == 0 and roty_rad == 0): cosphi = 0.
            else: cosphi = clamp(v.x/math.sin(roty_rad), -1., 1.)
            rotz_rad = math.acos(cosphi)
        elif(icase == 2):
            roty_rad = math.acos(v.z)
            if (v.x == 0 and roty_rad == 0): cosphi = 0.
            else: cosphi = clamp(v.x/math.sin(roty_rad), -1., 1.)
            rotz_rad = -math.acos(cosphi)
        elif(icase == 3):
            roty_rad = -math.acos(v.z)
            if (v.x == 0 and roty_rad == 0): cosphi = 0.
            else: cosphi = clamp(v.x/math.sin(roty_rad), -1., 1.)
            rotz_rad = math.acos(cosphi)
        elif(icase == 4):
            roty_rad = -math.acos(v.z)
            if (v.x == 0 and roty_rad == 0): cosphi = 0.
            else: cosphi = clamp(v.x/math.sin(roty_rad), -1., 1.)
            rotz_rad = -math.acos(cosphi)
        else:
            raise NameError('No rotation has been found!')
        
        theta = math.degrees(roty_rad)
        phi = math.degrees(rotz_rad)
        rotzy = get_rotateZ_tf(phi)*get_rotateY_tf(theta)
        v_ini_rotated = normalize(rotzy(v_ini))

        if (np.all(np.isclose(v.to_numpy()-v_ini_rotated.to_numpy(), 0., 0., 1e-14))):
            break
    
    return theta, phi
