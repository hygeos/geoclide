"""
Geometric transformations of the geoclide objects.

This module implements the Transform class, which encapsulates a
4x4 transformation matrix along with its inverse. A transform is
applied by calling it directly on a Vector, Point, Normal, Ray or
BBox, returning an object of the same nature, and two transforms
can be combined by multiplication. Helper functions create the
common transformations: translation, scale, and rotations around
the x, y or z axis or around an arbitrary axis.
"""

from __future__ import annotations

import math
import warnings
from typing import Any, cast, overload

import numpy as np
from numpy.linalg import inv

from geoclide.basic import BBox, Normal, Point, Ray, Vector
from geoclide.vecope import normalize


class Transform:
    """
    Represents 3D geometric transformation(s) using a 4x4 matrix or
    ntx4x4 matrix, where nt is the number of transformations

    It allows translation, rotation and scalling. It can be applied to
    vectors, points, normals and rays

    Parameters
    ----------
    m : Transform or ndarray, optional
        The matrix of the transformation(s), of shape (4, 4) or
        (nt, 4, 4) for a set of nt transformations
    m_inv : Transform or ndarray, optional
        The inverse matrix of the transformation(s), of shape
        (4, 4) or (nt, 4, 4) for a set of nt transformations

    Examples
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
    m_inv=
    array(
    [[1. 0. 0. 0.]
    [0. 1. 0. 0.]
    [0. 0. 1. 0.]
    [0. 0. 0. 1.]] )
    """

    def __init__(
        self,
        m: Transform | np.ndarray | None = None,
        m_inv: np.ndarray | None = None,
    ):
        if isinstance(m, Transform):
            self.m = m.m
            self.m_inv = m.m_inv
        elif m is None and m_inv is None:
            self.m = np.identity(4)
            self.m_inv = self.m.copy()
        elif isinstance(m, np.ndarray) and m_inv is None:
            if (len(m.shape) == 2 and m.shape != (4, 4)) or (
                len(m.shape) == 3 and m.shape[1] != 4 and m.shape[2] != 4
            ):
                raise ValueError(
                    "The m parameter must be an np.array of shape "
                    "(4,4) or (nt,4,4)"
                )
            self.m = m
            self.m_inv = inv(m)
        elif m is None and isinstance(m_inv, np.ndarray):
            if (len(m_inv.shape) == 2 and m_inv.shape != (4, 4)) or (
                len(m_inv.shape) == 3
                and m_inv.shape[1] != 4
                and m_inv.shape[2] != 4
            ):
                raise ValueError(
                    "The m_inv parameter must be an np.array of shape "
                    "(4,4) or (nt,4,4)"
                )
            self.m = inv(m_inv)
            self.m_inv = m_inv
        elif isinstance(m, np.ndarray) and isinstance(m_inv, np.ndarray):
            if (len(m.shape) == 2 and m.shape != (4, 4)) or (
                len(m.shape) == 3 and m.shape[1] != 4 and m.shape[2] != 4
            ):
                raise ValueError(
                    "The m parameter must be an np.array of shape "
                    "(4,4) or (nt,4,4)"
                )
            if (len(m_inv.shape) == 2 and m_inv.shape != (4, 4)) or (
                len(m_inv.shape) == 3
                and m_inv.shape[1] != 4
                and m_inv.shape[2] != 4
            ):
                raise ValueError(
                    "The m_inv parameter must be an np.array of shape "
                    "(4,4) or (nt,4,4)"
                )
            self.m = m
            self.m_inv = m_inv
        else:
            raise ValueError("Wrong parameter value(s) for Transform")

    def __eq__(self, t):
        if not isinstance(t, Transform):
            raise ValueError(
                "Equality with a Transform must be only with another Transform"
            )

        return bool(np.all(self.m == t.m) and np.all(self.m_inv == t.m_inv))

    def __mul__(self, t: Transform) -> Transform:
        if not isinstance(t, Transform):
            raise ValueError(
                "A transform can be multiplied only by another Transform"
            )

        return Transform(self.m @ t.m, t.m_inv @ self.m_inv)

    @overload
    def __call__(
        self, c: Vector, diag_calc: bool = ..., flatten: bool = ...
    ) -> Vector: ...

    @overload
    def __call__(
        self, c: Point, diag_calc: bool = ..., flatten: bool = ...
    ) -> Point: ...

    @overload
    def __call__(
        self, c: Normal, diag_calc: bool = ..., flatten: bool = ...
    ) -> Normal: ...

    @overload
    def __call__(
        self, c: Ray, diag_calc: bool = ..., flatten: bool = ...
    ) -> Ray: ...

    @overload
    def __call__(
        self, c: BBox, diag_calc: bool = ..., flatten: bool = ...
    ) -> BBox: ...

    def __call__(
        self,
        c: Vector | Point | Normal | Ray | BBox,
        diag_calc: bool = False,
        flatten: bool = False,
    ) -> Vector | Point | Normal | Ray | BBox | np.ndarray:
        """
        Apply the transformations

        Parameters
        ----------
        c : Vector or Point or Normal or Ray or BBox
            The vector(s)/point(s)/normal(s)/ray(s)/bounding box(es) to
            which the transformation is applied
        diag_calc : bool, optional
            Perform diagonal calculations between c(i) and
            tranformation(i). The number of transformations must be
            equal to the number of vectors/points/ ...

        Returns
        -------
        Vector or Point or Normal or Ray or BBox or ndarray
            The vector(s)/point(s)/normal(s)/ray(s)/bounding box(es)
            after the application of the transformation(s). In case of
            several transformations, it returns a 1-D ndarray of dtype
            equals to the c parameter type, but if flatten is True
            returns directly an object of same type as the c parameter.

        Examples
        --------
        >>> import geoclide as gc
        >>> t = gc.get_translate_tf(gc.Vector(5., 5., 5.))
        >>> p = gc.Point(0., 0., 0.)
        >>> t[p]
        Point(5.0, 5.0, 5.0)
        """
        # Case with several transformations in Transform
        if len(self.m.shape) == 3:
            is_vector = isinstance(c, Vector)
            is_point = isinstance(c, Point)
            is_normal = isinstance(c, Normal)
            use_flatten = False
            # default bindings, reassigned in the branch just below
            nt = 0
            mat = self.m
            x: Any = None
            y: Any = None
            z: Any = None
            keys: Any = None
            if isinstance(c, (Vector, Point, Normal)):
                nt = self.m.shape[0]
                if is_vector or is_point:
                    mat = np.moveaxis(self.m, 0, 2)
                else:  # if is_normal
                    mat = np.moveaxis(self.m_inv, 0, 2)
                cx, cy, cz = c.x, c.y, c.z
                key_bis = np.arange(nt)
                if isinstance(cx, np.ndarray) and not diag_calc:
                    mat = mat[:, :, np.newaxis, :]
                    x = cx[:, np.newaxis]
                    y = cast(np.ndarray, cy)[:, np.newaxis]
                    z = cast(np.ndarray, cz)[:, np.newaxis]
                    if flatten:
                        keys = (slice(None), key_bis)
                        use_flatten = True
                    else:
                        keys = [(slice(None), k) for k in key_bis]
                else:  # if diag_calc = True or if c is not an array
                    x = cx
                    y = cy
                    z = cz
                    keys = key_bis
            if is_vector:
                xv = mat[0, 0] * x + mat[0, 1] * y + mat[0, 2] * z
                yv = mat[1, 0] * x + mat[1, 1] * y + mat[1, 2] * z
                zv = mat[2, 0] * x + mat[2, 1] * y + mat[2, 2] * z
                if flatten:
                    if use_flatten:
                        vectors = Vector(
                            xv[keys].flatten("F"),
                            yv[keys].flatten("F"),
                            zv[keys].flatten("F"),
                        )
                    else:
                        vectors = Vector(xv[keys], yv[keys], zv[keys])
                else:
                    vectors = np.empty(nt, dtype=Vector)
                    for iv in range(0, nt):
                        vectors[iv] = Vector(
                            xv[keys[iv]], yv[keys[iv]], zv[keys[iv]]
                        )
                return vectors
            elif is_point:
                xp = mat[0, 0] * x + mat[0, 1] * y + mat[0, 2] * z + mat[0, 3]
                yp = mat[1, 0] * x + mat[1, 1] * y + mat[1, 2] * z + mat[1, 3]
                zp = mat[2, 0] * x + mat[2, 1] * y + mat[2, 2] * z + mat[2, 3]
                wp = mat[3, 0] * x + mat[3, 1] * y + mat[3, 2] * z + mat[3, 3]
                if flatten:
                    if use_flatten:
                        points = Point(
                            xp[keys].flatten("F"),
                            yp[keys].flatten("F"),
                            zp[keys].flatten("F"),
                        )
                        points /= wp[keys].flatten("F")
                    else:
                        points = Point(xp[keys], yp[keys], zp[keys]) / wp[keys]
                else:
                    points = np.empty(nt, dtype=Point)
                    for ip in range(0, nt):
                        if (
                            not isinstance(wp[keys[ip]], np.ndarray)
                            and wp[keys[ip]] == 1
                        ) or (
                            isinstance(wp[keys[ip]], np.ndarray)
                            and np.all(wp[keys[ip]] == 1)
                        ):
                            points[ip] = Point(
                                xp[keys[ip]], yp[keys[ip]], zp[keys[ip]]
                            )
                        else:
                            points[ip] = (
                                Point(
                                    xp[keys[ip]],
                                    yp[keys[ip]],
                                    zp[keys[ip]],
                                )
                                / wp[keys[ip]]
                            )
                return points
            elif is_normal:
                xn = mat[0, 0] * x + mat[1, 0] * y + mat[2, 0] * z
                yn = mat[0, 1] * x + mat[1, 1] * y + mat[2, 1] * z
                zn = mat[0, 2] * x + mat[1, 2] * y + mat[2, 2] * z
                if flatten:
                    if use_flatten:
                        normals = Normal(
                            xn[keys].flatten("F"),
                            yn[keys].flatten("F"),
                            zn[keys].flatten("F"),
                        )
                    else:
                        normals = Normal(xn[keys], yn[keys], zn[keys])
                else:
                    normals = np.empty(nt, dtype=Vector)
                    for inorm in range(0, nt):
                        normals[inorm] = Normal(
                            xn[keys[inorm]],
                            yn[keys[inorm]],
                            zn[keys[inorm]],
                        )
                return normals
            elif isinstance(c, Ray):
                origins = self(c.o, diag_calc, flatten)
                directions = self(c.d, diag_calc, flatten)
                if flatten:
                    rays = Ray(origins, directions, mint=c.mint, maxt=c.maxt)
                else:
                    # without flatten the recursive calls give back
                    # object ndarrays instead of Point/Vector
                    origins_arr = cast(np.ndarray, origins)
                    directions_arr = cast(np.ndarray, directions)
                    nt = self.m.shape[0]
                    rays = np.empty(nt, dtype=Ray)
                    for ir in range(0, nt):
                        rays[ir] = Ray(
                            origins_arr[ir],
                            directions_arr[ir],
                            mint=c.mint,
                            maxt=c.maxt,
                        )
                return rays
            elif isinstance(c, BBox):
                p0 = self(c.p0, diag_calc, flatten)
                v0 = self(c.p1 - c.p0, diag_calc, flatten)
                v1 = self(c.p3 - c.p0, diag_calc, flatten)
                v2 = self(c.p4 - c.p0, diag_calc, flatten)
                if flatten:
                    b = BBox()
                    b = b.union(p0)
                    b = b.union(p0 + v0)
                    b = b.union(p0 + (v0 + v1))
                    b = b.union(p0 + v1)
                    b = b.union(p0 + v2)
                    b = b.union(p0 + (v0 + v2))
                    b = b.union(p0 + (v0 + v1 + v2))
                    b = b.union(p0 + (v1 + v2))
                    bboxes = b
                else:
                    # without flatten the recursive calls give back
                    # object ndarrays instead of Point/Vector
                    p0a = cast(np.ndarray, p0)
                    v0a = cast(np.ndarray, v0)
                    v1a = cast(np.ndarray, v1)
                    v2a = cast(np.ndarray, v2)
                    nt = self.m.shape[0]
                    bboxes = np.empty(nt, dtype=BBox)
                    for ib in range(0, nt):
                        b = BBox()
                        b = b.union(p0a[ib])
                        b = b.union(p0a[ib] + v0a[ib])
                        b = b.union(p0a[ib] + (v0a[ib] + v1a[ib]))
                        b = b.union(p0a[ib] + v1a[ib])
                        b = b.union(p0a[ib] + v2a[ib])
                        b = b.union(p0a[ib] + (v0a[ib] + v2a[ib]))
                        b = b.union(p0a[ib] + (v0a[ib] + v1a[ib] + v2a[ib]))
                        b = b.union(p0a[ib] + (v1a[ib] + v2a[ib]))
                        bboxes[ib] = b
                return bboxes
            else:
                raise ValueError("Unknown type for transformations")
        else:
            if isinstance(c, Vector):
                xv = (
                    self.m[0, 0] * c.x
                    + self.m[0, 1] * c.y
                    + self.m[0, 2] * c.z
                )
                yv = (
                    self.m[1, 0] * c.x
                    + self.m[1, 1] * c.y
                    + self.m[1, 2] * c.z
                )
                zv = (
                    self.m[2, 0] * c.x
                    + self.m[2, 1] * c.y
                    + self.m[2, 2] * c.z
                )
                return Vector(xv, yv, zv)
            elif isinstance(c, Point):
                xp = (
                    self.m[0, 0] * c.x
                    + self.m[0, 1] * c.y
                    + self.m[0, 2] * c.z
                    + self.m[0, 3]
                )
                yp = (
                    self.m[1, 0] * c.x
                    + self.m[1, 1] * c.y
                    + self.m[1, 2] * c.z
                    + self.m[1, 3]
                )
                zp = (
                    self.m[2, 0] * c.x
                    + self.m[2, 1] * c.y
                    + self.m[2, 2] * c.z
                    + self.m[2, 3]
                )
                wp = (
                    self.m[3, 0] * c.x
                    + self.m[3, 1] * c.y
                    + self.m[3, 2] * c.z
                    + self.m[3, 3]
                )
                if (not isinstance(wp, np.ndarray) and wp == 1) or (
                    isinstance(wp, np.ndarray) and np.all(wp == 1)
                ):
                    return Point(xp, yp, zp)
                else:
                    return Point(xp, yp, zp) / wp
            elif isinstance(c, Normal):
                xn = (
                    self.m_inv[0, 0] * c.x
                    + self.m_inv[1, 0] * c.y
                    + self.m_inv[2, 0] * c.z
                )
                yn = (
                    self.m_inv[0, 1] * c.x
                    + self.m_inv[1, 1] * c.y
                    + self.m_inv[2, 1] * c.z
                )
                zn = (
                    self.m_inv[0, 2] * c.x
                    + self.m_inv[1, 2] * c.y
                    + self.m_inv[2, 2] * c.z
                )
                return Normal(xn, yn, zn)
            elif isinstance(c, Ray):
                return Ray(self(c.o), self(c.d), mint=c.mint, maxt=c.maxt)
            elif isinstance(c, BBox):
                b = BBox()
                p0 = self(c.p0)
                v0 = self(c.p1 - c.p0)
                v1 = self(c.p3 - c.p0)
                v2 = self(c.p4 - c.p0)
                b = b.union(p0)
                b = b.union(p0 + v0)
                b = b.union(p0 + (v0 + v1))
                b = b.union(p0 + v1)
                b = b.union(p0 + v2)
                b = b.union(p0 + (v0 + v2))
                b = b.union(p0 + (v0 + v1 + v2))
                b = b.union(p0 + (v1 + v2))
                return b
            else:
                raise ValueError("Unknown type for transformations")

    def __getitem__(
        self,
        c: Vector | Point | Normal | Ray | BBox,
        diag_calc: bool = False,
    ) -> Vector | Point | Normal | Ray | BBox | np.ndarray:
        """
        Apply the transformations

        Parameters
        ----------
        c : Vector or Point or Normal or Ray or BBox
            The Vector/Point/Normal/Ray/BBox to which the
            transformation is applied
        diag_calc : bool, optional
            Perform diagonal calculations between c(i) and
            tranformation(i). The number of transformations must be
            equal to the number of vectors / points / ...

        Returns
        -------
        Vector or Point or Normal or Ray or BBox or ndarray
            The Vector/Point/Normal/Ray/BBox after the transformation,
            or in case several transformations are given return a 1-D
            ndarray where dtype is equal to one of the previously
            mentioned classes.

        Examples
        --------
        >>> import geoclide as gc
        >>> t = gc.get_translate_tf(gc.Vector(5., 5., 5.))
        >>> p = gc.Point(0., 0., 0.)
        >>> t[p]
        Point(5.0, 5.0, 5.0)
        """
        warnings.simplefilter("always", DeprecationWarning)
        warn_message = (
            "\nApplying the transformation through square brackets is "
            "deprecated\nas of version 2.1.0 and will be no more "
            "possible in the future.\nPlease use parenthesis instead."
        )
        warnings.warn(warn_message, DeprecationWarning, stacklevel=1)
        return self(c, diag_calc)

    def __str__(self) -> str:
        print("m=\n", self.m, "\nm_inv=\n", self.m_inv)
        return ""

    def __repr__(self) -> str:
        print("m=\narray(\n", self.m, ")\nm_inv=\narray(\n", self.m_inv, ")")
        return ""

    def inverse(self) -> Transform:
        """
        Inverse the transformation(s) matrix

        Parameters
        ----------
        t : Transform
            The transformation(s) to be inversed

        Returns
        -------
        Transform
            The inversed transformation(s)
        """
        return get_inverse_tf(self)

    def is_identity(self) -> bool:
        return (
            (self.m[0, 0] == 1)
            and (self.m[0, 1] == 0)
            and (self.m[0, 2] == 0)
            and (self.m[0, 3] == 0)
            and (self.m[1, 0] == 0)
            and (self.m[1, 1] == 1)
            and (self.m[1, 2] == 0)
            and (self.m[1, 3] == 0)
            and (self.m[2, 0] == 0)
            and (self.m[2, 1] == 0)
            and (self.m[2, 2] == 1)
            and (self.m[2, 3] == 0)
            and (self.m[3, 0] == 0)
            and (self.m[3, 1] == 0)
            and (self.m[3, 2] == 0)
            and (self.m[3, 3] == 1)
        )

    def translate(self, v: Vector) -> Transform:
        """
        Update the self transformation(s) by adding a translate
        transformation(s)

        Parameters
        ----------
        v : Vector
            The vector(s) used for the transformation(s)

        Returns
        -------
        Transform
            The product of the self transformation(s) and the translate
            transformation(s)

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
        m_inv=
        array(
        [[ 1.  0.  0. -5.]
        [ 0.  1.  0.  0.]
        [ 0.  0.  1.  0.]
        [ 0.  0.  0.  1.]] )
        """
        t = get_translate_tf(v)
        return self * t

    def scale(self, v: Vector) -> Transform:
        """
        Update the self transformation(s) by adding a scale
        transformation(s)

        Parameters
        ----------
        v : Vector
            The vector(s) used for scale transformation(s)

        Returns
        -------
        Transform
            The product of the self transformation(s) and the scale
            transformation(s) matrices
        """
        t = get_scale_tf(v)
        return self * t

    def rotate_x(self, angle: float | np.ndarray) -> Transform:
        """
        Update the self transformation(s) by adding a rotate_x
        transformation(s)

        Parameters
        ----------
        angle : float or ndarray
            The angle(s) in degrees for the rotation(s) around the x
            axis. In case of an ndarray, it must be 1-D

        Returns
        -------
        Transform
            The product of the self transformation(s) and the rotate_x
            transformation(s) matrices
        """
        t = get_rotate_x_tf(angle)
        return self * t

    def rotate_y(self, angle: float | np.ndarray) -> Transform:
        """
        Update the self transformation(s) by adding a rotate_y
        transformation(s)

        Parameters
        ----------
        angle : float or ndarray
            The angle(s) in degrees for the rotation(s) around the y
            axis. In case of an ndarray, it must be 1-D

        Returns
        -------
        Transform
            The product of the self transformation(s) and the rotate_y
            transformation(s) matrices
        """
        t = get_rotate_y_tf(angle)
        return self * t

    def rotate_z(self, angle: float | np.ndarray) -> Transform:
        """
        Update the self transformation(s) by adding a rotate_z
        transformation(s)

        Parameters
        ----------
        angle : float or ndarray
            The angle(s) in degrees for the rotation(s) around the Z
            axis. In case of an ndarray, it must be 1-D

        Returns
        -------
        Transform
            The product of the initial transformation(s) and the
            rotate_z transformation(s) matrices
        """
        t = get_rotate_z_tf(angle)
        return self * t

    def rotate(
        self,
        angle: float | np.ndarray,
        axis: Vector | Normal,
        diag_calc: bool = False,
    ) -> Transform:
        """
        Update the self transformation(s) by adding a rotate
        transformation(s)

        .. warning::
            The angle parameter can be a 1-D array only if axis
            parameter is a Vector/Normal with scalar x, y, z
            components, or if the parameter diag_calc=True

        Parameters
        ----------
        angle : float or ndarray
            The angle(s) in degrees for the rotation(s). In case
            of an ndarray, it must be 1-D
        axis : Vector or Normal
            The rotation(s) is/are performed around the
            vector(s)/normal(s) axis/axes
        diag_calc : bool, optional
                Perform diagonal calculations in case angle is a 1-D
                ndarray and axis is a Vector/Normal with 1-D ndarray
                x, y, z components. Use angle(i) with axis(i) to
                calculate transformation(i)

        Returns
        -------
        Transform
            The product of the self transformation(s) and the rotate
            transformation(s) matrices
        """
        t = get_rotate_tf(angle, axis, diag_calc=diag_calc)
        return self * t


def get_inverse_tf(t: Transform) -> Transform:
    """
    Get the inverse transformation(s)

    Parameters
    ----------
    t : Transform
        The transformation(s) to be inversed

    Returns
    -------
    Transform
        The inversed transformation(s)
    """
    return Transform(t.m_inv, t.m)


def get_translate_tf(v: Vector) -> Transform:
    """
    Get the translate transformation(s)

    Parameters
    ----------
    v : Vector
        The vector(s) used for the translate transformation(s)

    Returns
    -------
    Transform
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
    m_inv=
    array(
    [[ 1.  0.  0. -5.]
    [ 0.  1.  0. -0.]
    [ 0.  0.  1. -0.]
    [ 0.  0.  0.  1.]] )
    """
    if not isinstance(v, Vector):
        raise ValueError("The parameter v must be a Vector")
    if isinstance(v.x, np.ndarray):
        nc = len(v.x)
        m = np.tile(np.identity(4, dtype=np.float64), (nc, 1)).reshape(
            nc, 4, 4
        )
        m_inv = m.copy()
        m[:, 0, 3] = v.x
        m[:, 1, 3] = v.y
        m[:, 2, 3] = v.z
        m_inv[:, 0, 3] = -v.x
        m_inv[:, 1, 3] = -v.y
        m_inv[:, 2, 3] = -v.z
    else:
        m = np.identity(4)
        m_inv = m.copy()
        m[0, 3] = v.x
        m[1, 3] = v.y
        m[2, 3] = v.z
        m_inv[0, 3] = -v.x
        m_inv[1, 3] = -v.y
        m_inv[2, 3] = -v.z
    return Transform(m, m_inv)


def get_scale_tf(v: Vector) -> Transform:
    """
    Get the scale transformation(s)

    Parameters
    ----------
    v : Vector
        The vector(s) used for scale transformation(s)

    Returns
    -------
    Transform
        The scale transformation(s)
    """
    if not isinstance(v, Vector):
        raise ValueError("The parameter v must be a Vector")

    if isinstance(v.x, np.ndarray):
        nc = len(v.x)
        m = np.tile(np.identity(4, dtype=np.float64), (nc, 1)).reshape(
            nc, 4, 4
        )
        m_inv = m.copy()
        m[:, 0, 0] = v.x
        m[:, 1, 1] = v.y
        m[:, 2, 2] = v.z
        m_inv[:, 0, 0] = 1.0 / v.x
        m_inv[:, 1, 1] = 1.0 / v.y
        m_inv[:, 2, 2] = 1.0 / v.z
    else:
        m = np.identity(4)
        m_inv = m.copy()
        m[0, 0] = v.x
        m[1, 1] = v.y
        m[2, 2] = v.z
        m_inv[0, 0] = 1.0 / v.x
        m_inv[1, 1] = 1.0 / v.y
        m_inv[2, 2] = 1.0 / v.z
    return Transform(m, m_inv)


def get_rotate_x_tf(angle: float | np.ndarray) -> Transform:
    """
    Get the rotate_x transformation(s)

    Parameters
    ----------
    angle : float or ndarray
        The angle(s) in degrees for the rotation(s) around the x
        axis. In case of an ndarray, it must be 1-D

    Returns
    -------
    Transform
        The rotate_x transformation(s)
    """
    is_ang_arr = isinstance(angle, np.ndarray)

    if is_ang_arr:
        sin_t = np.sin(angle * (math.pi / 180.0))
        cos_t = np.cos(angle * (math.pi / 180.0))
        nc = len(angle)
        m = np.tile(np.identity(4, dtype=np.float64), (nc, 1)).reshape(
            nc, 4, 4
        )
        m[:, 1, 1] = cos_t
        m[:, 1, 2] = -1.0 * sin_t
        m[:, 2, 1] = sin_t
        m[:, 2, 2] = cos_t
        return Transform(m, np.transpose(m, axes=(0, 2, 1)))
    else:
        sin_t = math.sin(angle * (math.pi / 180.0))
        cos_t = math.cos(angle * (math.pi / 180.0))
        m = np.identity(4)
        m[1, 1] = cos_t
        m[1, 2] = -1.0 * sin_t
        m[2, 1] = sin_t
        m[2, 2] = cos_t
        return Transform(m, np.transpose(m))


def get_rotate_y_tf(angle: float | np.ndarray) -> Transform:
    """
    Get the rotate_y transformation(s)

    Parameters
    ----------
    angle : float or ndarray
        The angle(s) in degrees for the rotation(s) around the y
        axis. In case of an ndarray, it must be 1-D

    Returns
    -------
    Transform
        The rotate_y transformation(s)
    """
    is_ang_arr = isinstance(angle, np.ndarray)

    if is_ang_arr:
        sin_t = np.sin(angle * (math.pi / 180.0))
        cos_t = np.cos(angle * (math.pi / 180.0))
        nc = len(angle)
        m = np.tile(np.identity(4, dtype=np.float64), (nc, 1)).reshape(
            nc, 4, 4
        )
        m[:, 0, 0] = cos_t
        m[:, 2, 0] = -1.0 * sin_t
        m[:, 0, 2] = sin_t
        m[:, 2, 2] = cos_t
        return Transform(m, np.transpose(m, axes=(0, 2, 1)))
    else:
        sin_t = math.sin(angle * (math.pi / 180.0))
        cos_t = math.cos(angle * (math.pi / 180.0))
        m = np.identity(4)
        m[0, 0] = cos_t
        m[2, 0] = -1.0 * sin_t
        m[0, 2] = sin_t
        m[2, 2] = cos_t
        return Transform(m, np.transpose(m))


def get_rotate_z_tf(angle: float | np.ndarray) -> Transform:
    """
    Get the rotate_z transformation(s)

    Parameters
    ----------
    angle : float or ndarray
        The angle(s) in degrees for the rotation(s) around the Z
        axis. In case of an ndarray, it must be 1-D

    Returns
    -------
    Transform
        The rotate_z transformation(s)
    """
    is_ang_arr = isinstance(angle, np.ndarray)

    if is_ang_arr:
        sin_t = np.sin(angle * (math.pi / 180.0))
        cos_t = np.cos(angle * (math.pi / 180.0))
        nc = len(angle)
        m = np.tile(np.identity(4, dtype=np.float64), (nc, 1)).reshape(
            nc, 4, 4
        )
        m[:, 0, 0] = cos_t
        m[:, 0, 1] = -sin_t
        m[:, 1, 0] = sin_t
        m[:, 1, 1] = cos_t
        return Transform(m, np.transpose(m, axes=(0, 2, 1)))
    else:
        sin_t = math.sin(angle * (math.pi / 180.0))
        cos_t = math.cos(angle * (math.pi / 180.0))
        m = np.identity(4)
        m[0, 0] = cos_t
        m[0, 1] = -sin_t
        m[1, 0] = sin_t
        m[1, 1] = cos_t
        return Transform(m, np.transpose(m))


def get_rotate_tf(
    angle: float | np.ndarray,
    axis: Vector | Normal,
    diag_calc: bool = False,
) -> Transform:
    """
    Get the rotate transformation(s) around a given axis/axes

    .. warning::
            The angle parameter can be a 1-D array only if axis
            parameter is a Vector/Normal with scalar x, y, z
            components, or if the parameter diag_calc=True

    Parameters
    ----------
    angle : float or ndarray
        The angle(s) in degrees for the rotation(s). In case of
        an ndarray, it must be 1-D
    axis : Vector or Normal
        The rotation(s) is/are performed around the
        vector(s)/normal(s) axis/axes
    diag_calc : bool, optional
            Perform diagonal calculations in case angle is a 1-D
            ndarray and axis is a Vector/Normal with 1-D ndarray x, y,
            z components. Use angle(i) with axis(i) to calculate
            transformation(i)

    Returns
    -------
    Transform
        The rotate transformation(s)
    """
    if not isinstance(axis, (Vector, Normal)):
        raise ValueError("The parameter axis must be a Vector or a Normal")

    is_ang_arr = isinstance(angle, np.ndarray)
    is_axis_arr = isinstance(axis.x, np.ndarray)

    if is_ang_arr and is_axis_arr and not diag_calc:
        raise ValueError(
            "1-D array for angle parameter is allowed only if axis "
            "parameter is a Vector/Normal with scalar components, or "
            "if diag_calc=True"
        )

    a = Vector(normalize(axis))
    if is_ang_arr or is_axis_arr:
        nc = 1
        if is_ang_arr:
            nc = max(nc, len(angle))
            s = np.sin(angle * (math.pi / 180.0))
            c = np.cos(angle * (math.pi / 180.0))
        else:
            s = math.sin(angle * (math.pi / 180.0))
            c = math.cos(angle * (math.pi / 180.0))

        if is_axis_arr:
            nc = max(nc, len(cast(np.ndarray, axis.x)))
        m = np.tile(np.identity(4, dtype=np.float64), (nc, 1)).reshape(
            nc, 4, 4
        )

        m[:, 0, 0] = a.x * a.x + (1 - a.x * a.x) * c
        m[:, 0, 1] = a.x * a.y * (1 - c) - a.z * s
        m[:, 0, 2] = a.x * a.z * (1 - c) + a.y * s

        m[:, 1, 0] = a.x * a.y * (1 - c) + a.z * s
        m[:, 1, 1] = a.y * a.y + (1 - a.y * a.y) * c
        m[:, 1, 2] = a.y * a.z * (1 - c) - a.x * s

        m[:, 2, 0] = a.x * a.z * (1 - c) - a.y * s
        m[:, 2, 1] = a.y * a.z * (1 - c) + a.x * s
        m[:, 2, 2] = a.z * a.z + (1 - a.z * a.z) * c
        return Transform(m, np.transpose(m, axes=(0, 2, 1)))
    else:
        s = math.sin(angle * (math.pi / 180.0))
        c = math.cos(angle * (math.pi / 180.0))
        m = np.identity(4)

        m[0, 0] = a.x * a.x + (1 - a.x * a.x) * c
        m[0, 1] = a.x * a.y * (1 - c) - a.z * s
        m[0, 2] = a.x * a.z * (1 - c) + a.y * s

        m[1, 0] = a.x * a.y * (1 - c) + a.z * s
        m[1, 1] = a.y * a.y + (1 - a.y * a.y) * c
        m[1, 2] = a.y * a.z * (1 - c) - a.x * s

        m[2, 0] = a.x * a.z * (1 - c) - a.y * s
        m[2, 1] = a.y * a.z * (1 - c) + a.x * s
        m[2, 2] = a.z * a.z + (1 - a.z * a.z) * c
        return Transform(m, np.transpose(m))
