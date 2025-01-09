#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math


def clamp(val, val_min, val_max):
    """
    Clamps val into the range [val_min, val_max]

    Parameters
    ----------
    val : float
        The scalar to be clamped
    val_min : float
        The minumum value
    val_max : The maximum value

    Results
    -------
    out : float
        The result of the clamp

    Examples
    --------
    >>> import geoclide as gc
    >>> gc.clamp(4, val_min=5, val_max=11)
    5
    """
    if (not np.isscalar(val)     or 
        not np.isscalar(val_min) or
        not np.isscalar(val_max) ):
        raise ValueError('The parameters must be all scalars')
    
    return val_min if (val < val_min) else (val_max if (val > val_max) else val)


def swap(a, b):
    """
    Swaps the values of two scalar variables

    Parameters
    ----------
    a : float
        The first variable to swap
    b : The second variable to swap

    Results
    -------
    a : float
        The a variable with the value of the variable b
    b : float
        The b variable with the value of the variable a

    Examples
    --------
    >>> import geoclide as gc
    >>> a = 5
    >>> a
    5
    >>> b = 12
    >>> b
    12
    >>> a, b = gc.swap(a, b)
    >>> a
    12
    >>> b
    5
    """
    if (not np.isscalar(a) or not np.isscalar(b) ):
        raise ValueError('The parameters must be all scalars')
    
    return b, a


def quadratic(a, b, c):
    """
    Resolve the quadratic polynomial: ax**2 + bx + c

    - where x is the quadratic polynomial variable and a, b and c the coefficients

    Parameters
    ----------
    a : float
        The first coefficient of the quadratic polynomial
    b : float
        The second coefficient of the quadratic polynomial
    c : float
        The third coefficient of the quadratic polynomial

    Results
    -------
    b : bool
        If the quadratic can be solved -> True, else False
    x0 : float | None
        The first solution
    x1 : float | None 
        The second solution

    Notes
    -----
    If There are 2 solutions x0 < x1. And if there is only one solution x0 = x1.

    Examples
    --------
    >>> import geoclide as gc
    >>> a = 2
    >>> b = -5
    >>> c = 0
    >>> gc.quadratic(a, b, c)
    (True, 0.0, 2.5)
    """
    if (not np.isscalar(a) or not np.isscalar(b) or not np.isscalar(c)):
        raise ValueError('The parameters must be all scalars')
    
    #  Find quadratic discriminant
    discrim = (b * b) - (4 * a * c)

    if (discrim < 0): return False, None, None

    rootDiscrim = math.sqrt(discrim)

    # Compute quadratic xi values
    if (b < 0): q = -0.5 * (b - rootDiscrim)
    else: q = -0.5 * (b + rootDiscrim)

    if (a != 0): x0 = q / a
    else: x0 = c / q

    x1 = c / q

    if (x0 > x1): x0, x1 = swap(x0, x1)

    return True, x0, x1