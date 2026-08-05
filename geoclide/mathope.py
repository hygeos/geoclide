"""
Basic mathematical utilities.

This module gathers the scalar/array helper functions used by the
geoclide shapes: clamping of a value into a range, solving of
quadratic equations, and the gamma terms bounding the
floating-point rounding errors (in simple and double precision)
used by the ray-shape intersection tests.
"""

from __future__ import annotations

import math
from typing import overload

import numpy as np


def clamp(val: float, val_min: float, val_max: float) -> float:
    """
    Clamps val into the range [val_min, val_max]

    Parameters
    ----------
    val : float
        The scalar to be clamped
    val_min : float
        The minimum value
    val_max : float
        The maximum value

    Returns
    -------
    float
        The result of the clamp

    Examples
    --------
    >>> import geoclide as gc
    >>> gc.clamp(4, val_min=5, val_max=11)
    5
    """
    if not all(np.isscalar(v) for v in (val, val_min, val_max)):
        raise ValueError("The parameters must be all scalars")

    return val_min if val < val_min else (val_max if val > val_max else val)


@overload
def quadratic(
    a: np.ndarray, b: np.ndarray, c: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]: ...


@overload
def quadratic(
    a: float, b: float, c: float
) -> tuple[bool, float | None, float | None]: ...


def quadratic(
    a: float | np.ndarray,
    b: float | np.ndarray,
    c: float | np.ndarray,
) -> tuple[
    bool | np.ndarray,
    float | np.ndarray | None,
    float | np.ndarray | None,
]:
    """
    Resolve the quadratic polynomial: ax**2 + bx + c

    - where x is the quadratic polynomial variable and a, b and c the
      coefficients

    Parameters
    ----------
    a : float or ndarray
        The first coefficient(s) of the quadratic polynomial. In
        case of an ndarray, it must be 1-D
    b : float or ndarray
        The second coefficient(s) of the quadratic polynomial. In
        case of an ndarray, it must be 1-D
    c : float or ndarray
        The third coefficient(s) of the quadratic polynomial. In
        case of an ndarray, it must be 1-D

    Returns
    -------
    b : bool or ndarray
        If the quadratic can be solved return True, else False.
        In case of an ndarray, it is a 1-D ndarray of booleans
    x0 : float or None or ndarray
        The first solution(s). In case of an ndarray, it is 1-D
    x1 : float or None or ndarray
        The second solution(s). In case of an ndarray, it is 1-D

    Notes
    -----
    If There are 2 solutions x0 < x1. And if there is only one
    solution x0 = x1.

    Examples
    --------
    >>> import geoclide as gc
    >>> a = 2
    >>> b = -5
    >>> c = 0
    >>> gc.quadratic(a, b, c)
    (True, 0.0, 2.5)
    """
    if isinstance(a, np.ndarray):
        # an ndarray a implies ndarrays b and c
        assert isinstance(b, np.ndarray)
        assert isinstance(c, np.ndarray)
        with np.errstate(divide="ignore", invalid="ignore"):
            # Find quadratic discriminant
            discrim = (b * b) - (4 * a * c)

            c1 = discrim < 0
            is_solution = np.logical_not(c1)
            root_discrim = np.sqrt(discrim)

            # Compute quadratic xi values
            q = np.where(
                b < 0,
                -0.5 * (b - root_discrim),
                -0.5 * (b + root_discrim),
            )

            x1 = c / q
            x0 = np.where(a != 0, q / a, x1)

            # keep the smallest solution in x0
            c4 = x0 > x1
            x0, x1 = np.where(c4, x1, x0), np.where(c4, x0, x1)

            x0[c1] = None
            x1[c1] = None

        return is_solution, x0, x1
    else:
        # Find quadratic discriminant
        discrim = (b * b) - (4 * a * c)

        if discrim < 0:
            return False, None, None

        root_discrim = math.sqrt(discrim)

        # Compute quadratic xi values
        if b < 0:
            q = -0.5 * (b - root_discrim)
        else:
            q = -0.5 * (b + root_discrim)

        if a != 0:
            x0 = q / a
        else:
            x0 = c / q

        x1 = c / q

        if x0 > x1:
            x0, x1 = x1, x0

        return True, x0, x1


@overload
def gamma_f32(n: float) -> float: ...


@overload
def gamma_f32(n: np.ndarray) -> np.ndarray: ...


def gamma_f32(n: float | np.ndarray) -> float | np.floating | np.ndarray:
    """
    :meta private:

    Gamma function from pbrt v3
    """
    epsi = np.finfo(np.float32).eps * 0.5
    return (n * epsi) / (1 - n * epsi)


@overload
def gamma_f64(n: float) -> float: ...


@overload
def gamma_f64(n: np.ndarray) -> np.ndarray: ...


def gamma_f64(n: float | np.ndarray) -> float | np.floating | np.ndarray:
    """
    :meta private:

    Gamma function from pbrt v3 but in double precision
    """
    epsi = np.finfo(np.float64).eps * 0.5
    return (n * epsi) / (1 - n * epsi)
