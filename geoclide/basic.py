"""
Basic geometric objects used across geoclide.

This module implements the elementary objects on which the whole
package is built: vectors, points, normals, rays and axis-aligned
bounding boxes. All of them can represent either a single element
or a set of elements (their components being ndarrays), allowing
vectorized calculations.

Key Classes
-----------
Vector
    A direction in the three-dimensional space, with the common
    operators (addition, subtraction, scaling, ...).
Point
    A position in the three-dimensional space. The subtraction of
    two points gives a vector.
Normal
    A vector perpendicular to a surface at a particular position.
    It is not necessarily normalized and is transformed
    differently from a vector.
Ray
    A semi-infinite line described by an origin (Point), a
    direction (Vector) and the parametric range [mint, maxt].
BBox
    An axis-aligned bounding box described by its pmin and pmax
    corner points, supporting union operations and ray
    intersection tests.
"""

from __future__ import annotations

import math
import warnings
from datetime import datetime
from typing import Literal, TypeVar, cast, overload

import numpy as np
import xarray as xr

from geoclide.constants import GAMMA3_F64, VERSION


def _xyz_arrays(
    vpn: Vector | Point | Normal,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    :meta private:

    Give the x, y and z components, narrowed to ndarrays
    """
    return (
        cast(np.ndarray, vpn.x),
        cast(np.ndarray, vpn.y),
        cast(np.ndarray, vpn.z),
    )


_VPN = TypeVar("_VPN", "Vector", "Point", "Normal")

FLOAT64 = np.dtype(np.float64)


def _new_vpn(
    cls: type[_VPN],
    x: float | np.ndarray,
    y: float | np.ndarray,
    z: float | np.ndarray,
) -> _VPN:
    """
    :meta private:

    Create a Vector/Point/Normal without checking the components

    It is used by the operators, where the components are already
    the floats or the float64 ndarrays returned by numpy.
    """
    vpn = object.__new__(cls)
    vpn.x = x
    vpn.y = y
    vpn.z = z
    return vpn


def _init_xyz(
    x: float | np.ndarray | Vector | Point | Normal | None,
    y: float | np.ndarray | None,
    z: float | np.ndarray | None,
    copy: bool = False,
) -> tuple[
    float | np.ndarray, float | np.ndarray, float | np.ndarray
]:
    """
    :meta private:

    Give the x, y and z components of a Vector, Point or Normal

    It gathers the constructor logic of the 3 classes. Unless copy
    is True, the components given as 3 ndarrays are not copied (see
    the notes of the classes).
    """
    if isinstance(x, np.ndarray):
        if isinstance(y, np.ndarray) and isinstance(z, np.ndarray):
            if copy:
                return (
                    np.array(x, dtype=np.float64),
                    np.array(y, dtype=np.float64),
                    np.array(z, dtype=np.float64),
                )
            if (
                x.dtype is FLOAT64
                and y.dtype is FLOAT64
                and z.dtype is FLOAT64
            ):
                return x, y, z
            return (
                np.asarray(x, dtype=np.float64),
                np.asarray(y, dtype=np.float64),
                np.asarray(z, dtype=np.float64),
            )
        if y is None and z is None:
            if x.ndim == 1 and len(x) == 3:
                return float(x[0]), float(x[1]), float(x[2])
            if x.ndim == 2 and x.shape[1] == 3:
                # the columns are copied to get contiguous components
                return (
                    x[:, 0].astype(np.float64),
                    x[:, 1].astype(np.float64),
                    x[:, 2].astype(np.float64),
                )
        raise ValueError("Wrong parameter value(s)")
    # the exact type is checked, a numpy float is converted below to
    # keep python floats as components
    if type(x) is float and type(y) is float and type(z) is float:
        return x, y, z
    if x is None and y is None and z is None:
        return 0.0, 0.0, 0.0
    if isinstance(x, (Vector, Point, Normal)):
        if copy and isinstance(x.x, np.ndarray):
            xa, ya, za = _xyz_arrays(x)
            return xa.copy(), ya.copy(), za.copy()
        return x.x, x.y, x.z
    if np.isscalar(x) and np.isscalar(y) and np.isscalar(z):
        return (
            float(cast(float, x)),
            float(cast(float, y)),
            float(cast(float, z)),
        )
    raise ValueError("Wrong parameter value(s)")


class Vector:
    """
    Parameters
    ----------
    x : float or ndarray or Point or Vector or Normal, optional
        The x component(s) of the vector (see notes)
    y : float or ndarray, optional
        The y component(s) of the vector. In case of an ndarray,
        it must be 1-D
    z : float or ndarray, optional
        The z component(s) of the vector. In case of an ndarray,
        it must be 1-D
    copy : bool, optional
        If True the given ndarrays are copied, else they are used
        as they are (see notes)

    Notes
    -----
    - if the parameter x is a 1-D ndarray of size 3 and y and z are
      None, the values of x, y and z will be equal to respectively
      x[0], x[1], and x[2]
    - if the parameter x is a 2-D ndarray of shape (n,3) and y and z
      are None, the values of x, y and z will be equal to respectively
      x[:,0], x[:,1], and x[:,2]
    - if the parameter x is a Point, Vector or Normal, it will
      circumvent the y and z parameters and take the components of the
      Point/Vector/Normal for x, y and z values
    - the x, y and z ndarrays given as parameters are not copied,
      modifying them afterwards modifies the vector. Use copy=True to
      get a vector with its own components

    Examples
    --------
    >>> import geoclide as gc
    >>> v1 = gc.Vector(0.,0.,1.)
    >>> v1
    Vector(0,0,1)
    """

    __array_priority__ = 1
    fmt = ".8f"

    def __init__(
        self,
        x: float | np.ndarray | Vector | Point | Normal | None = None,
        y: float | np.ndarray | None = None,
        z: float | np.ndarray | None = None,
        copy: bool = False,
    ):
        self.x, self.y, self.z = _init_xyz(x, y, z, copy)

    def __eq__(self, v2):
        if isinstance(v2, Vector):
            if isinstance(self.x, np.ndarray) or isinstance(v2.x, np.ndarray):
                return np.logical_and.reduce(
                    (self.x == v2.x, self.y == v2.y, self.z == v2.z)
                )
            else:
                return (
                    (self.x == v2.x) and (self.y == v2.y) and (self.z == v2.z)
                )
        else:
            raise ValueError(
                "Equality with a Vector must be only with another Vector"
            )

    def __add__(self, v2: Vector) -> Vector:
        if isinstance(v2, Vector):
            return _new_vpn(
                Vector, self.x + v2.x, self.y + v2.y, self.z + v2.z
            )
        else:
            raise ValueError(
                "Addition with a Vector must be only with another Vector"
            )

    def __sub__(self, v2: Vector) -> Vector:
        if isinstance(v2, Vector):
            return _new_vpn(
                Vector, self.x - v2.x, self.y - v2.y, self.z - v2.z
            )
        else:
            raise ValueError(
                "Substraction with a Vector must be only with another Vector"
            )

    def __truediv__(self, sca: float | np.ndarray) -> Vector:
        div = 1.0 / sca
        return _new_vpn(
            Vector, self.x * div, self.y * div, self.z * div
        )

    def __mul__(self, sca: float | np.ndarray) -> Vector:
        return _new_vpn(
            Vector, sca * self.x, sca * self.y, sca * self.z
        )

    def __rmul__(self, sca: float | np.ndarray) -> Vector:
        return _new_vpn(
            Vector, sca * self.x, sca * self.y, sca * self.z
        )

    def __neg__(self) -> Vector:
        return _new_vpn(Vector, -self.x, -self.y, -self.z)

    def __getitem__(self, ind: int) -> float | np.ndarray:
        if not isinstance(ind, (int, np.integer)):
            raise IndexError("Only an integer is a valid index")
        if ind == 0:
            return self.x
        elif ind == 1:
            return self.y
        elif ind == 2:
            return self.z
        else:
            raise IndexError(f"Index {ind} is out of range")

    def __str__(self) -> str:
        return print_basic(self)

    def __repr__(self) -> str:
        return print_basic(self, self.__class__.__name__)

    def length_squared(self) -> float | np.ndarray:
        return self.x * self.x + self.y * self.y + self.z * self.z

    def length(self) -> float | np.ndarray:
        if isinstance(self.x, np.ndarray):
            return np.sqrt(self.length_squared())
        else:
            return math.sqrt(self.length_squared())

    def to_numpy(self) -> np.ndarray:
        if isinstance(self.x, np.ndarray):
            return np.array([self.x, self.y, self.z], dtype=np.float64).T
        else:
            return np.array([self.x, self.y, self.z], dtype=np.float64)


class Point:
    """
    Parameters
    ----------
    x : float or ndarray or Point or Vector or Normal, optional
        The x component(s) of the point (see notes)
    y : float or ndarray, optional
        The y component(s) of the point. In case of an ndarray,
        it must be 1-D
    z : float or ndarray, optional
        The z component(s) of the point. In case of an ndarray,
        it must be 1-D
    copy : bool, optional
        If True the given ndarrays are copied, else they are used
        as they are (see notes)

    Notes
    -----
    - if the parameter x is a 1-D ndarray of size 3 and y and z are
      None, the values of x, y and z will be equal to respectively
      x[0], x[1], and x[2]
    - if the parameter x is a 2-D ndarray of shape (n,3) and y and z
      are None, the values of x, y and z will be equal to respectively
      x[:,0], x[:,1], and x[:,2]
    - if the parameter x is a Point, Vector or Normal, it will
      circumvent the y and z parameters and take the components of the
      Point/Vector/Normal for x, y and z values
    - the x, y and z ndarrays given as parameters are not copied,
      modifying them afterwards modifies the point. Use copy=True to
      get a point with its own components

    Examples
    --------
    >>> import geoclide as gc
    >>> p1 = gc.Point(0.,0.,1.)
    >>> p1
    Point(0,0,1)
    """

    __array_priority__ = 1
    fmt = ".8f"

    def __init__(
        self,
        x: float | np.ndarray | Vector | Point | Normal | None = None,
        y: float | np.ndarray | None = None,
        z: float | np.ndarray | None = None,
        copy: bool = False,
    ):
        self.x, self.y, self.z = _init_xyz(x, y, z, copy)

    def __eq__(self, p2):
        if isinstance(p2, Point):
            if isinstance(self.x, np.ndarray) or isinstance(p2.x, np.ndarray):
                return np.logical_and.reduce(
                    (self.x == p2.x, self.y == p2.y, self.z == p2.z)
                )
            else:
                return (
                    (self.x == p2.x) and (self.y == p2.y) and (self.z == p2.z)
                )
        else:
            raise ValueError(
                "Equality with a Point must be only with another Point"
            )

    def __add__(self, v: Vector | Point) -> Point:
        if isinstance(v, (Vector, Point)):
            return _new_vpn(
                Point, self.x + v.x, self.y + v.y, self.z + v.z
            )
        else:
            raise ValueError(
                "Addition with a Point must be only with a Vector or"
                " (exceptionally tolerated) another Point"
            )

    @overload
    def __sub__(self, vp2: Vector) -> Point: ...

    @overload
    def __sub__(self, vp2: Point) -> Vector: ...

    def __sub__(self, vp2: Vector | Point) -> Point | Vector:
        if isinstance(vp2, Vector):
            return _new_vpn(
                Point, self.x - vp2.x, self.y - vp2.y, self.z - vp2.z
            )
        elif isinstance(vp2, Point):
            return _new_vpn(
                Vector, self.x - vp2.x, self.y - vp2.y, self.z - vp2.z
            )
        else:
            raise ValueError(
                "Substraction with a Point must be with another Point "
                "or a Vector"
            )

    def __truediv__(self, sca: float | np.ndarray) -> Point:
        div = 1.0 / sca
        return _new_vpn(
            Point, self.x * div, self.y * div, self.z * div
        )

    def __mul__(self, sca: float | np.ndarray) -> Point:
        return _new_vpn(
            Point, sca * self.x, sca * self.y, sca * self.z
        )

    def __rmul__(self, sca: float | np.ndarray) -> Point:
        return _new_vpn(
            Point, sca * self.x, sca * self.y, sca * self.z
        )

    def __neg__(self) -> Point:
        return _new_vpn(Point, -self.x, -self.y, -self.z)

    def __getitem__(self, ind: int) -> float | np.ndarray:
        if not isinstance(ind, (int, np.integer)):
            raise IndexError("Only an integer is a valid index")
        if ind == 0:
            return self.x
        elif ind == 1:
            return self.y
        elif ind == 2:
            return self.z
        else:
            raise IndexError(f"Index {ind} is out of range")

    def __str__(self) -> str:
        return print_basic(self)

    def __repr__(self) -> str:
        return print_basic(self, self.__class__.__name__)

    def to_numpy(self) -> np.ndarray:
        if isinstance(self.x, np.ndarray):
            return np.array([self.x, self.y, self.z], dtype=np.float64).T
        else:
            return np.array([self.x, self.y, self.z], dtype=np.float64)


class Normal:
    """
    Parameters
    ----------
    x : float or ndarray or Point or Vector or Normal, optional
        The x component(s) of the normal (see notes)
    y : float or ndarray, optional
        The y component(s) of the normal. In case of an ndarray,
        it must be 1-D
    z : float or ndarray, optional
        The z component(s) of the normal. In case of an ndarray,
        it must be 1-D
    copy : bool, optional
        If True the given ndarrays are copied, else they are used
        as they are (see notes)

    Notes
    -----
    - if the parameter x is a 1-D ndarray of size 3 and y and z are
      None, the values of x, y and z will be equal to respectively
      x[0], x[1], and x[2]
    - if the parameter x is a 2-D ndarray of shape (n,3) and y and z
      are None, the values of x, y and z will be equal to respectively
      x[:,0], x[:,1], and x[:,2]
    - if the parameter x is a Point, Vector or Normal, it will
      circumvent the y and z parameters and take the components of the
      Point/Vector/Normal for x, y and z values
    - the x, y and z ndarrays given as parameters are not copied,
      modifying them afterwards modifies the normal. Use copy=True to
      get a normal with its own components

    Examples
    --------
    >>> import geoclide as gc
    >>> n1 = gc.Normal(0.,0.,1.)
    >>> n1
    Normal(0,0,1)
    """

    __array_priority__ = 1
    fmt = ".8f"

    def __init__(
        self,
        x: float | np.ndarray | Vector | Point | Normal | None = None,
        y: float | np.ndarray | None = None,
        z: float | np.ndarray | None = None,
        copy: bool = False,
    ):
        self.x, self.y, self.z = _init_xyz(x, y, z, copy)

    def __eq__(self, n2):
        if isinstance(n2, Normal):
            if isinstance(self.x, np.ndarray) or isinstance(n2.x, np.ndarray):
                return np.logical_and.reduce(
                    (self.x == n2.x, self.y == n2.y, self.z == n2.z)
                )
            else:
                return (
                    (self.x == n2.x) and (self.y == n2.y) and (self.z == n2.z)
                )
        else:
            raise ValueError(
                "Equality with a Normal must be only with another Normal"
            )

    def __add__(self, n2: Normal) -> Normal:
        if isinstance(n2, Normal):
            return _new_vpn(
                Normal, self.x + n2.x, self.y + n2.y, self.z + n2.z
            )
        else:
            raise ValueError(
                "Addition with a Normal must be only with another Normal"
            )

    def __sub__(self, n2: Normal) -> Normal:
        if isinstance(n2, Normal):
            return _new_vpn(
                Normal, self.x - n2.x, self.y - n2.y, self.z - n2.z
            )
        else:
            raise ValueError(
                "Substraction with a Normal must be only with another Normal"
            )

    def __truediv__(self, sca: float | np.ndarray) -> Normal:
        div = 1.0 / sca
        return _new_vpn(
            Normal, self.x * div, self.y * div, self.z * div
        )

    def __mul__(self, sca: float | np.ndarray) -> Normal:
        return _new_vpn(
            Normal, sca * self.x, sca * self.y, sca * self.z
        )

    def __rmul__(self, sca: float | np.ndarray) -> Normal:
        return _new_vpn(
            Normal, sca * self.x, sca * self.y, sca * self.z
        )

    def __neg__(self) -> Normal:
        return _new_vpn(Normal, -self.x, -self.y, -self.z)

    def __getitem__(self, ind: int) -> float | np.ndarray:
        if not isinstance(ind, (int, np.integer)):
            raise IndexError("Only an integer is a valid index")
        if ind == 0:
            return self.x
        elif ind == 1:
            return self.y
        elif ind == 2:
            return self.z
        else:
            raise IndexError(f"Index {ind} is out of range")

    def __str__(self) -> str:
        return print_basic(self)

    def __repr__(self) -> str:
        return print_basic(self, self.__class__.__name__)

    def length_squared(self) -> float | np.ndarray:
        return self.x * self.x + self.y * self.y + self.z * self.z

    def length(self) -> float | np.ndarray:
        if isinstance(self.x, np.ndarray):
            return np.sqrt(self.length_squared())
        else:
            return math.sqrt(self.length_squared())

    def to_numpy(self) -> np.ndarray:
        if isinstance(self.x, np.ndarray):
            return np.array([self.x, self.y, self.z], dtype=np.float64).T
        else:
            return np.array([self.x, self.y, self.z], dtype=np.float64)


class Ray:
    """
    Definition of ray:

    r(t) = o + t*d, where:

    - o is/are the origin point(s) of the ray(s)
    - d is/are the direction(s) of the ray(s)
    - t belongs to stricly positive real numbers

    Parameters
    ----------
    o : Point or Ray
        Origin point(s) of the ray(s). If the o parameter is a Ray ->
        circumvent all the parameters by the ray attributs
    d : Vector
        Direction(s) of the ray(s)
    mint : float, optional
        The minimum t value
    maxt : float, optional
        The maximum t value

    Examples
    --------
    >>> import geoclide as gc
    >>> o = gc.Point(0., 50., 2.)
    >>> d = gc.Vector(0.,0.,1.)
    >>> r = gc.Ray(o, d, mint=20, maxt=100)
    >>> r
    r(t) = (0.0, 50.0, 2.0) + t*(0.0, 0.0, 1.0) with t ∈ [20,100[
    """

    def __init__(
        self,
        o: Point | Ray,
        d: Vector | None = None,
        mint: float = 0,
        maxt: float = float("inf"),
    ):
        if isinstance(o, Ray):
            self.o = o.o
            self.d = o.d
            self.mint = float(o.mint)
            self.maxt = float(o.maxt)
        else:
            if not isinstance(o, Point):
                raise ValueError("The parameter o must be a Point or a Ray")
            if not isinstance(d, Vector):
                raise ValueError("The parameter d must only be a Vector")
            if not all(np.isscalar(v) for v in (mint, maxt)):
                raise ValueError(
                    "The parameters mint and maxt must be both scalars"
                )
            if mint > maxt:
                raise ValueError("maxt must be greater than mint")
            self.o = o
            self.d = d
            self.mint = mint
            self.maxt = maxt

    def __call__(self, t: float | np.ndarray) -> Point:
        """
        Solve ray(s) equation(s)

        Parameters
        ----------
        t : float or ndarray
            The t rays(s) values(s). The value(s) must lie between mint
            and maxt. In case of an ndarray, it must be 1-D

        Returns
        -------
        Point
            The result(s) of the equation r(t) = o + t*d

        Examples
        --------
        >>> import geoclide as gc
        >>> o = gc.Point(0., 0., 0.)
        >>> d = gc.Vector(1., 0., 0.)
        >>> r = gc.Ray(o, d)
        >>> t = 10.
        >>> r(t)
        Point(10., 0., 0.)
        """
        if (
            isinstance(t, np.ndarray)
            and np.any(np.logical_or(t < self.mint, t > self.maxt))
        ) or (
            not isinstance(t, np.ndarray) and (t < self.mint or t > self.maxt)
        ):
            raise ValueError(
                f"The value {t} is out of bounds. It must be between "
                f"{self.mint} and {self.maxt}"
            )
        else:
            return self.o + self.d * t

    def __getitem__(self, t: float | np.ndarray) -> Point:
        warnings.simplefilter("always", DeprecationWarning)
        warn_message = (
            "\nThe use of square brackets is deprecated as of version "
            "2.1.0 and will be\nno more possible in the future. "
            "Please use parenthesis instead."
        )
        warnings.warn(warn_message, DeprecationWarning, stacklevel=1)
        return self(t)

    def __str__(self) -> str:
        if not isinstance(self.o.x, np.ndarray):
            return (
                f"({self.o.x}, {self.o.y}, {self.o.z}) + "
                f"t*({self.d.x}, {self.d.y}, {self.d.z})"
                f" with t ∈ [{self.mint},{self.maxt}["
            )
        else:
            ox, oy, oz = _xyz_arrays(self.o)
            dx, dy, dz = _xyz_arrays(self.d)
            nrays = len(ox)
            mint = np.full(nrays, self.mint, dtype=np.float64)
            maxt = np.full(nrays, self.maxt, dtype=np.float64)
            output = ""
            if nrays <= 100:
                for ir in range(0, nrays):
                    output += (
                        f"({ox[ir]}, {oy[ir]}, "
                        f"{oz[ir]}) + "
                        f"t{ir}*({dx[ir]}, {dy[ir]}, "
                        f"{dz[ir]})"
                        f" with t{ir} ∈ [{mint[ir]},{maxt[ir]}["
                    )
                    if ir < nrays - 1:
                        output += "\n"
            else:
                for ir in range(0, 97):
                    output += (
                        f"({ox[ir]}, {oy[ir]}, "
                        f"{oz[ir]}) + "
                        f"t{ir}*({dx[ir]}, {dy[ir]}, "
                        f"{dz[ir]})"
                        f" with t{ir} ∈ [{mint[ir]},{maxt[ir]}["
                    )
                    output += "\n"
                output += "       ...\n"
                for ir in range(nrays - 3, nrays):
                    output += (
                        f"({ox[ir]}, {oy[ir]}, "
                        f"{oz[ir]}) + "
                        f"t{ir}*({dx[ir]}, {dy[ir]}, "
                        f"{dz[ir]})"
                        f" with t{ir} ∈ [{mint[ir]},{maxt[ir]}["
                    )
                    if ir < nrays - 1:
                        output += "\n"
            return output

    def __repr__(self) -> str:
        if not isinstance(self.o.x, np.ndarray):
            return (
                f"r(t) = ({self.o.x}, {self.o.y}, {self.o.z}) + "
                f"t*({self.d.x}, {self.d.y}, {self.d.z})"
                f" with t ∈ [{self.mint},{self.maxt}["
            )
        else:
            ox, oy, oz = _xyz_arrays(self.o)
            dx, dy, dz = _xyz_arrays(self.d)
            nrays = len(ox)
            mint = np.full(nrays, self.mint, dtype=np.float64)
            maxt = np.full(nrays, self.maxt, dtype=np.float64)
            output = ""
            if nrays <= 100:
                for ir in range(0, nrays):
                    output += (
                        f"r(t{ir}) = ({ox[ir]}, {oy[ir]}, "
                        f"{oz[ir]}) + "
                        f"t{ir}*({dx[ir]}, {dy[ir]}, "
                        f"{dz[ir]})"
                        f" with t{ir} ∈ [{mint[ir]},{maxt[ir]}["
                    )
                    if ir < nrays - 1:
                        output += "\n"
            else:
                for ir in range(0, 97):
                    output += (
                        f"r(t{ir}) = ({ox[ir]}, {oy[ir]}, "
                        f"{oz[ir]}) + "
                        f"t{ir}*({dx[ir]}, {dy[ir]}, "
                        f"{dz[ir]})"
                        f" with t{ir} ∈ [{mint[ir]},{maxt[ir]}["
                    )
                    output += "\n"
                output += "       ...\n"
                for ir in range(nrays - 3, nrays):
                    output += (
                        f"r(t{ir}) = ({ox[ir]}, {oy[ir]}, "
                        f"{oz[ir]}) + "
                        f"t{ir}*({dx[ir]}, {dy[ir]}, "
                        f"{dz[ir]})"
                        f" with t{ir} ∈ [{mint[ir]},{maxt[ir]}["
                    )
                    if ir < nrays - 1:
                        output += "\n"
            return output


class BBox:
    """
    Bounding Box

    Parameters
    ----------
    p1 : Point, optional
        Frist point(s) to use to create the bounding box(es)
    p2 : Point, optional
        Second point(s) to use to create the bounding box(es)

    Examples
    --------
    >>> import geoclide as gc
    >>> p1 = gc.Point(0., 0., 0.)
    >>> p2 = gc.Point(1., 1., 1.)
    >>> b1 = gc.BBox(p1, p2)
    >>> b1
    pmin=Point(0.0, 0.0, 0.0), pmax=Point(1.0, 1.0, 1.0)
    """

    def __init__(self, p1: Point | None = None, p2: Point | None = None):
        if isinstance(p1, Point) and isinstance(p2, Point):
            if isinstance(p1.x, np.ndarray) or isinstance(p2.x, np.ndarray):
                self.pmin = Point(
                    np.minimum(p1.x, p2.x),
                    np.minimum(p1.y, p2.y),
                    np.minimum(p1.z, p2.z),
                )
                self.pmax = Point(
                    np.maximum(p1.x, p2.x),
                    np.maximum(p1.y, p2.y),
                    np.maximum(p1.z, p2.z),
                )
            else:
                self.pmin = Point(
                    min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z)
                )
                self.pmax = Point(
                    max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z)
                )
        elif p1 is None and p2 is None:
            self.pmin = Point(float("inf"), float("inf"), float("inf"))
            self.pmax = Point(float("-inf"), float("-inf"), float("-inf"))
        elif isinstance(p1, Point) and p2 is None:
            self.pmin = p1
            self.pmax = p1
        elif p1 is None and isinstance(p2, Point):
            self.pmin = p2
            self.pmax = p2
        else:
            raise ValueError("The only parameters accepted are Point objects")

        # The 8 vertices of the BBox
        # - p0=pmin, then next 3 points are in the XY plane at z=pmin.z
        #   the order being anti-clockwise
        # - next 4 points are in the XY plane at z=pmax.z, starting
        #   with point p4 just above p0, so p6=pmax
        self.p0 = Point(self.pmin.x, self.pmin.y, self.pmin.z)
        self.p1 = Point(self.pmax.x, self.pmin.y, self.pmin.z)
        self.p2 = Point(self.pmax.x, self.pmax.y, self.pmin.z)
        self.p3 = Point(self.pmin.x, self.pmax.y, self.pmin.z)
        self.p4 = Point(self.pmin.x, self.pmin.y, self.pmax.z)
        self.p5 = Point(self.pmax.x, self.pmin.y, self.pmax.z)
        self.p6 = Point(self.pmax.x, self.pmax.y, self.pmax.z)
        self.p7 = Point(self.pmin.x, self.pmax.y, self.pmax.z)
        self.vertices = [
            self.p0,
            self.p1,
            self.p2,
            self.p3,
            self.p4,
            self.p5,
            self.p6,
            self.p7,
        ]

    def __str__(self) -> str:
        return (
            f"pmin=({self.pmin.x}, {self.pmin.y}, {self.pmin.z}), "
            f"pmax=({self.pmax.x}, {self.pmax.y}, {self.pmax.z})"
        )

    def __repr__(self) -> str:
        return (
            f"pmin=Point({self.pmin.x}, {self.pmin.y}, {self.pmin.z}), "
            f"pmax=Point({self.pmax.x}, {self.pmax.y}, {self.pmax.z})"
        )

    def union(self, b: Point | BBox) -> BBox:
        """
        Union with a point/set of points or a bounding box/set of
        bounding boxes

        Parameters
        ----------
        b : Point or BBox
            The point(s) or bounding box(es) to use for the union

        Returns
        -------
        BBox
            The new bounding box(es) after the union

        Examples
        --------
        >>> import geoclide as gc
        >>> p1 = gc.Point(0., 0., 0.)
        >>> p2 = gc.Point(1., 1., 1.)
        >>> p3 = gc.Point(1., 1., 3.)
        >>> b1 = gc.BBox(p1, p2)
        >>> b1
        pmin=Point(0.0, 0.0, 0.0), pmax=Point(1.0, 1.0, 1.0)
        >>> b2 = b1.union(p3)
        >>> b2
        pmin=Point(0.0, 0.0, 0.0), pmax=Point(1.0, 1.0, 3.0)
        """
        pmin = Point()
        pmax = Point()
        if isinstance(b, Point):
            if isinstance(b.x, np.ndarray) or isinstance(
                self.p0.x, np.ndarray
            ):
                pmin.x = np.minimum(self.pmin.x, b.x)
                pmin.y = np.minimum(self.pmin.y, b.y)
                pmin.z = np.minimum(self.pmin.z, b.z)
                pmax.x = np.maximum(self.pmax.x, b.x)
                pmax.y = np.maximum(self.pmax.y, b.y)
                pmax.z = np.maximum(self.pmax.z, b.z)
            else:
                pmin.x = min(self.pmin.x, b.x)
                pmin.y = min(self.pmin.y, b.y)
                pmin.z = min(self.pmin.z, b.z)
                pmax.x = max(self.pmax.x, b.x)
                pmax.y = max(self.pmax.y, b.y)
                pmax.z = max(self.pmax.z, b.z)
        elif isinstance(b, BBox):
            if isinstance(b.p0.x, np.ndarray) or isinstance(
                self.p0.x, np.ndarray
            ):
                pmin.x = np.minimum(self.pmin.x, b.pmin.x)
                pmin.y = np.minimum(self.pmin.y, b.pmin.y)
                pmin.z = np.minimum(self.pmin.z, b.pmin.z)
                pmax.x = np.maximum(self.pmax.x, b.pmax.x)
                pmax.y = np.maximum(self.pmax.y, b.pmax.y)
                pmax.z = np.maximum(self.pmax.z, b.pmax.z)
            else:
                pmin.x = min(self.pmin.x, b.pmin.x)
                pmin.y = min(self.pmin.y, b.pmin.y)
                pmin.z = min(self.pmin.z, b.pmin.z)
                pmax.x = max(self.pmax.x, b.pmax.x)
                pmax.y = max(self.pmax.y, b.pmax.y)
                pmax.z = max(self.pmax.z, b.pmax.z)
        else:
            raise ValueError("The union must be with another BBox or Point")

        return BBox(pmin, pmax)

    def is_inside(self, p: Point) -> bool | np.ndarray:
        """
        Test if point(s) p is/are included in the bounding box(es)
        """
        if isinstance(self.p0.x, np.ndarray) or isinstance(p.x, np.ndarray):
            conds = cast(
                "tuple[np.ndarray, ...]",
                (
                    (p.x >= self.pmin.x),
                    (p.x <= self.pmax.x),
                    (p.y >= self.pmin.y),
                    (p.y <= self.pmax.y),
                    (p.z >= self.pmin.z),
                    (p.z <= self.pmax.z),
                ),
            )
            return np.logical_and.reduce(conds)
        else:
            return (
                (p.x >= self.pmin.x)
                and (p.x <= self.pmax.x)
                and (p.y >= self.pmin.y)
                and (p.y <= self.pmax.y)
                and (p.z >= self.pmin.z)
                and (p.z <= self.pmax.z)
            )

    def is_intersection(
        self, r: Ray, diag_calc: bool = False
    ) -> bool | np.ndarray:
        """
        Test if a ray/rays intersect(s) the bounding box(es)

        Parameters
        ----------
        r : Ray
            The ray(s) to use for the intersection test(s)
        diag_calc : bool
            Perform diagonal calculations in case of multiple bounding
            boxes and rays, the output is a 1-D array instead of a 2-D
            array where out[i] is calculated using r(i) and bbox(i). The
            same size for the BBox and Ray objects is required.

        Returns
        -------
        bool or ndarray
            If there is at least 1 intersection returns True, else
            False. In case of an ndarray, it is an ndarray of
            booleans, 1-D, or 2-D for a set of rays and a set of
            bounding boxes.

        Examples
        --------
        >>> import geoclide as gc
        >>> p1 = gc.Point(0., 0., 0.)
        >>> p2 = gc.Point(1., 1., 1.)
        >>> b1 = gc.BBox(p1, p2)
        pmin=Point(0.0, 0.0, 0.0), pmax=Point(1.0, 1.0, 1.0)
        >>> p3 = gc.Point(0.5, 0.5, 0.1)
        >>> v1 = gc.Vector(0., 0., 1.)
        >>> r1 = gc.Ray(p3, v1)
        >>> r1
        r(t) = (0.5, 0.5, 0.1) + t*(0.0, 0.0, 1.0) with t ∈ [0,inf[
        >>> b1.is_intersection(r1)
        True
        """
        t0, t1, is_intersection = self.intersect(
            r, diag_calc=diag_calc, ds_output=False
        )
        return is_intersection

    @overload
    def intersect(
        self,
        r: Ray,
        diag_calc: bool = ...,
        *,
        ds_output: Literal[True] = ...,
    ) -> xr.Dataset: ...

    @overload
    def intersect(
        self,
        r: Ray,
        diag_calc: bool = ...,
        *,
        ds_output: Literal[False],
    ) -> tuple[
        float | np.ndarray, float | np.ndarray, bool | np.ndarray
    ]: ...

    @overload
    def intersect(
        self, r: Ray, diag_calc: bool = ..., *, ds_output: bool
    ) -> xr.Dataset | tuple: ...

    def intersect(
        self, r: Ray, diag_calc: bool = False, ds_output: bool = True
    ) -> xr.Dataset | tuple:
        """
        Test if a ray/rays intersect(s) the bounding box(es)

        There are 3 possibilities:

        - no intersection
        - only 1 intersection (case of ray located initially inside the
          BBox)
        - 2 intersections

        Parameters
        ----------
        r : Ray
            The ray(s) to use for the intersection test(s)
        diag_calc : bool, optional
            Perform diagonal calculations in case of multiple bounding
            boxes and rays, the output is a 1-D array instead of a 2-D
            array where out[i] is calculated using r(i) and bbox(i). The
            same size for the BBox and Ray objects is required.
        ds_output : bool, optional
            If True the output is a dataset, else returns a tuple with
            intersection information variables

        Returns
        -------
        Dataset or tuple
            Xarray dataset containing the intersection information
            if ds_output is True (see the get_bbox_intersect_dataset
            function for its variables), else a tuple. Form of the
            tuple:

            * t0 : None or float or ndarray
                -> The t ray variable of the first intersection. In case
                of only 1 intersection it represents nothing. An
                ndarray is 1-D, or 2-D for a set of rays and a set
                of bounding boxes.
            * t1 : None or float or ndarray
                -> The t ray variable of the second intersection. In
                case of only 1 intersection, t1 becomes the t ray
                variable of the first intersection. An ndarray is
                1-D, or 2-D for a set of rays and a set of
                bounding boxes.
            * is_intersection : bool or ndarray
                -> If there is at least 1 intersection return True, else
                False. An ndarray is an ndarray of booleans, 1-D,
                or 2-D for a set of rays and a set of bounding
                boxes.

        Examples
        --------
        >>> import geoclide as gc
        >>> p1 = gc.Point(0., 0., 0.)
        >>> p2 = gc.Point(1., 1., 1.)
        >>> b1 = gc.BBox(p1, p2)
        pmin=Point(0.0, 0.0, 0.0), pmax=Point(1.0, 1.0, 1.0)
        >>> p3 = gc.Point(0.5, 0.5, 0.1)
        >>> v1 = gc.Vector(0., 0., 1.)
        >>> r1 = gc.Ray(p3, v1)
        >>> r1
        r(t) = (0.5, 0.5, 0.1) + t*(0.0, 0.0, 1.0) with t ∈ [0,inf[
        >>> t0, t1, is_intersection = b1.intersect(r1, ds_output=False)
        >>> t0, t1, is_intersection
        (0.0, 0.9, True)
        >>> r1[t1]
        Point(0.5, 0.5, 1.0)
        """
        if not isinstance(r, Ray):
            raise ValueError("The given parameter must be a Ray")
        is_r_arr = isinstance(r.o.x, np.ndarray)
        is_bbox_arr = isinstance(self.pmin.x, np.ndarray)
        if is_r_arr and is_bbox_arr and not diag_calc:
            with np.errstate(divide="ignore", invalid="ignore"):
                b_size = len(cast(np.ndarray, self.pmin.x))
                r_size = len(cast(np.ndarray, r.o.x))
                t0 = np.zeros((b_size, r_size), dtype=np.float64)
                t1 = np.full((b_size, r_size), r.maxt, dtype=np.float64)
                is_intersection = np.full((b_size, r_size), True)
                for i in range(3):
                    rdi = cast(np.ndarray, r.d[i])
                    roi = cast(np.ndarray, r.o[i])
                    pmini = cast(np.ndarray, self.pmin[i])
                    pmaxi = cast(np.ndarray, self.pmax[i])
                    inv_ray_dir = np.where(
                        rdi != 0, 1.0 / rdi, math.inf
                    )
                    t_near = (
                        pmini[:, None] - roi[None, :]
                    ) * inv_ray_dir
                    t_far = (
                        pmaxi[:, None] - roi[None, :]
                    ) * inv_ray_dir
                    c1 = t_near > t_far
                    t_near, t_far = (
                        np.where(c1, t_far, t_near),
                        np.where(c1, t_near, t_far),
                    )
                    t_far *= 1 + 2 * GAMMA3_F64
                    c2 = np.logical_and(t_near > t0, is_intersection)
                    c3 = np.logical_and(t_far < t1, is_intersection)
                    t0 = np.where(c2, t_near, t0)
                    t1 = np.where(c3, t_far, t1)
                    c4 = t0 > t1
                    is_intersection = np.logical_and(
                        is_intersection, np.logical_not(c4)
                    )
                    t0[c4] = 0.0
                    t1[c4] = 0.0
            if ds_output:
                return get_bbox_intersect_dataset(
                    self, r, t0, t1, is_intersection
                )
            else:
                return t0, t1, is_intersection
        elif is_r_arr or is_bbox_arr:
            with np.errstate(divide="ignore", invalid="ignore"):
                size = 1
                if is_bbox_arr:
                    size = max(size, len(cast(np.ndarray, self.pmin.x)))
                if is_r_arr:
                    size = max(size, len(cast(np.ndarray, r.o.x)))
                t0 = np.zeros(size, dtype=np.float64)
                t1 = np.full(size, r.maxt, dtype=np.float64)
                is_intersection = np.full(size, True)
                for i in range(3):
                    rdi = r.d[i]
                    inv_ray_dir: float | np.ndarray
                    if isinstance(rdi, np.ndarray):
                        inv_ray_dir = np.where(
                            rdi != 0, 1.0 / rdi, math.inf
                        )
                    elif rdi != 0:
                        inv_ray_dir = 1.0 / rdi
                    else:
                        inv_ray_dir = math.inf
                    t_near = cast(
                        np.ndarray,
                        (self.pmin[i] - r.o[i]) * inv_ray_dir,
                    )
                    t_far = cast(
                        np.ndarray,
                        (self.pmax[i] - r.o[i]) * inv_ray_dir,
                    )
                    c1 = t_near > t_far
                    t_near, t_far = (
                        np.where(c1, t_far, t_near),
                        np.where(c1, t_near, t_far),
                    )
                    t_far *= 1 + 2 * GAMMA3_F64
                    c2 = np.logical_and(t_near > t0, is_intersection)
                    c3 = np.logical_and(t_far < t1, is_intersection)
                    t0 = np.where(c2, t_near, t0)
                    t1 = np.where(c3, t_far, t1)
                    c4 = t0 > t1
                    is_intersection = np.logical_and(
                        is_intersection, np.logical_not(c4)
                    )
                    t0[c4] = 0.0
                    t1[c4] = 0.0
            if ds_output:
                return get_bbox_intersect_dataset(
                    self, r, t0, t1, is_intersection
                )
            else:
                return t0, t1, is_intersection
        else:
            t0 = 0.0
            t1 = r.maxt
            for i in range(3):
                rdi = cast(float, r.d[i])
                roi = cast(float, r.o[i])
                pmini = cast(float, self.pmin[i])
                pmaxi = cast(float, self.pmax[i])
                if rdi != 0:
                    inv_ray_dir = 1.0 / rdi
                else:
                    inv_ray_dir = math.inf
                t_near = (pmini - roi) * inv_ray_dir
                t_far = (pmaxi - roi) * inv_ray_dir
                if t_near > t_far:
                    t_near, t_far = t_far, t_near
                t_far *= 1 + 2 * GAMMA3_F64
                t0 = t_near if t_near > t0 else t0
                t1 = t_far if t_far < t1 else t1
                if t0 > t1:
                    if ds_output:
                        return get_bbox_intersect_dataset(
                            self, r, 0.0, 0.0, False
                        )
                    else:
                        return 0.0, 0.0, False
            if ds_output:
                return get_bbox_intersect_dataset(self, r, t0, t1, True)
            else:
                return t0, t1, True

    def common_vertices(self, b: BBox) -> np.ndarray:
        """
        Get a list of boolean checking which vertices (self) are common
        to the bounding box(es) b

        Parameters
        ----------
        b : BBox
            The secondary bounding box(es)

        Returns
        -------
        ndarray
            Returns an array of boolean values indicating if the
            bounding box(es) vertices are common to the secondary b
            bounding box(es) vertices

        Examples
        --------
        >>> import geoclide as gc
        >>> b0 = gc.BBox(gc.Point(0., 0., 0.), gc.Point(1., 1., 1.))
        >>> b1 = gc.BBox(gc.Point(1., 0., 0.), gc.Point(2., 1., 1.))
        >>> b0.common_vertices(b1)
        array([False,  True,  True, False, False,  True,  True, False])
        >>> b1.common_vertices(b0)
        array([ True, False, False,  True,  True, False, False,  True])
        """
        return get_common_vertices(self, b)

    def common_face(
        self, b: BBox, fill_value: int | float | None = None
    ) -> int | float | np.ndarray | None:
        """
        Get the face index/indices which is/are common with one of the
        face(s) of bounding box(es) b2

        The convention of index from face 0 to 5, for +X,-X,+Y,-Y,+Z,-Z:

        >>>    |F2|                     |+Y|
        >>> |F1|F4|F0|F5|  where ->  |-X|+Z|+X|-Z|
        >>>    |F3|                     |-Y|

        `More information <https://en.wikipedia.org/wiki/Cube_mapping>`_

        Parameters
        ----------
        b : BBox
            The secondary bounding box(es)
        fill_value : integer, optional
            In case there is no common face(s) returns fill_value

        Returns
        -------
        int or ndarray
            Returns the index/indices of the common face(s) or
            fill_value. In case of an ndarray, it is 1-D

        Examples
        --------
        >>> import geoclide as gc
        >>> b0 = gc.BBox(gc.Point(0., 0., 0.), gc.Point(1., 1., 1.))
        >>> b1 = gc.BBox(gc.Point(1., 0., 0.), gc.Point(2., 1., 1.))
        >>> gc.get_common_face(b1, b2)
        0
        >>> gc.get_common_face(b2, b1)
        1
        """
        return get_common_face(self, b, fill_value=fill_value)


def get_common_vertices(b1: BBox, b2: BBox) -> np.ndarray:
    """
    Check which vertices of bounding box(es) b1 are common to the
    vectices of bounding box(es) b2

    Parameters
    ----------
    b1 : BBox
        The principal bounding box(es)
    b2 : BBox
        The secondary bounding box(es)

    Returns
    -------
    ndarray
        Returns an array of boolean values indicating whether the
        principal bounding box(es) b1 vertices are common to the
        secondary bounding box(es) b2 vertices. It is 1-D, or 2-D
        in case of a set of bounding boxes.

    Examples
    --------
    >>> import geoclide as gc
    >>> b0 = gc.BBox(gc.Point(0., 0., 0.), gc.Point(1., 1., 1.))
    >>> b1 = gc.BBox(gc.Point(1., 0., 0.), gc.Point(2., 1., 1.))
    >>> gc.get_common_vertices(b1, b2)
    array([False,  True,  True, False, False,  True,  True, False])
    >>> gc.get_common_vertices(b1, b2)
    array([ True, False, False,  True,  True, False, False,  True])
    """
    if not isinstance(b1, BBox) or not isinstance(b2, BBox):
        raise ValueError("The parameters b1 and b2 must be both BBox objects")

    size = 1
    if isinstance(b1.p0.x, np.ndarray):
        size = max(len(b1.p0.x), size)
    if isinstance(b2.p0.x, np.ndarray):
        size = max(len(b2.p0.x), size)

    if size > 1:
        res = np.full((size, 8), False, dtype=bool)
        for i in range(0, 8):
            res[:, i] = np.logical_or.reduce(
                (
                    b1.vertices[i] == b2.vertices[0],
                    b1.vertices[i] == b2.vertices[1],
                    b1.vertices[i] == b2.vertices[2],
                    b1.vertices[i] == b2.vertices[3],
                    b1.vertices[i] == b2.vertices[4],
                    b1.vertices[i] == b2.vertices[5],
                    b1.vertices[i] == b2.vertices[6],
                    b1.vertices[i] == b2.vertices[7],
                )
            )
        return res
    else:
        return np.array(list(map(lambda x: x in b2.vertices, b1.vertices)))


def get_common_face(
    b1: BBox, b2: BBox, fill_value: int | float | None = None
) -> int | float | np.ndarray | None:
    """

    Get the face index/indices of the bounding box(es) b1 which is/are
    common to the bounding box(es) b2

    The convention of index from face 0 to 5, for +X,-X,+Y,-Y,+Z,-Z:

    >>>    |F2|                     |+Y|
    >>> |F1|F4|F0|F5|  where ->  |-X|+Z|+X|-Z|
    >>>    |F3|                     |-Y|

    `More information <https://en.wikipedia.org/wiki/Cube_mapping>`_

    Parameters
    ----------
    b1 : BBox
        The principal bounding box(es)
    b2 : BBox
        The secondary bounding box(es)
    fill_value : integer, optional
            In case there is no common face(s) returns fill_value

    Returns
    -------
    int or ndarray
        Returns the index/indices of the common face(s) or
        fill_value. In case of an ndarray, it is 1-D

    Examples
    --------
    >>> import geoclide as gc
    >>> b0 = gc.BBox(gc.Point(0., 0., 0.), gc.Point(1., 1., 1.))
    >>> b1 = gc.BBox(gc.Point(1., 0., 0.), gc.Point(2., 1., 1.))
    >>> gc.get_common_face(b1, b2)
    0
    >>> gc.get_common_face(b2, b1)
    1
    """
    ok = get_common_vertices(b1, b2)
    if len(ok.shape) > 1:
        size1 = ok.shape[0]
        cond = ok.sum(axis=1) == 4
        res = np.full(size1, fill_value, dtype=np.float64)
        if any(cond):
            n = np.zeros((size1, 8), dtype=np.int32)
            for i in range(0, 8):
                n[:, i] = i
            size2 = cond.sum()
            res_bis = res[cond]
            n = n[cond][ok[cond]].reshape(size2, 4)
            m1 = np.repeat(np.array([[1, 2, 5, 6]]), size2, axis=0)
            m2 = np.repeat(np.array([[0, 3, 4, 7]]), size2, axis=0)
            m3 = np.repeat(np.array([[2, 3, 6, 7]]), size2, axis=0)
            m4 = np.repeat(np.array([[0, 1, 4, 5]]), size2, axis=0)
            m5 = np.repeat(np.array([[4, 5, 6, 7]]), size2, axis=0)
            m6 = np.repeat(np.array([[0, 1, 2, 3]]), size2, axis=0)
            res_bis[np.all(n == m1, axis=1)] = 0
            res_bis[np.all(n == m2, axis=1)] = 1
            res_bis[np.all(n == m3, axis=1)] = 2
            res_bis[np.all(n == m4, axis=1)] = 3
            res_bis[np.all(n == m5, axis=1)] = 4
            res_bis[np.all(n == m6, axis=1)] = 5
            res[cond] = res_bis
        return res
    else:
        if ok.sum() == 4:
            n = np.arange(8)[ok]
            if np.array_equal(n, np.array([1, 2, 5, 6])):
                return 0
            elif np.array_equal(n, np.array([0, 3, 4, 7])):
                return 1
            elif np.array_equal(n, np.array([2, 3, 6, 7])):
                return 2
            elif np.array_equal(n, np.array([0, 1, 4, 5])):
                return 3
            elif np.array_equal(n, np.array([4, 5, 6, 7])):
                return 4
            elif np.array_equal(n, np.array([0, 1, 2, 3])):
                return 5
            else:
                return fill_value
        else:
            return fill_value


def get_bbox_intersect_dataset(
    bbox: BBox,
    r: Ray,
    t0: float | np.ndarray | None = None,
    t1: float | np.ndarray | None = None,
    is_intersection: bool | np.ndarray = False,
    diag_calc: bool = False,
) -> xr.Dataset:
    """
    Create dataset containing the intersection test information

    - The intersect method return of BBox class gives the t0, t1 and
      is_intersection inputs of this function

    Parameters
    ----------
    bbox : BBox
        The bounding box(es) used for the intersection test
    r : Ray
        The ray(s) used for the intersection test
    t0 : float or ndarray
        The t ray variable of the first intersection. In case of
        an ndarray, it is 1-D, or 2-D for a set of rays and a set
        of bounding boxes
    t1 : float or ndarray
        The t ray variable of the second intersection. In case of
        an ndarray, it is 1-D, or 2-D for a set of rays and a set
        of bounding boxes
    is_intersection : bool or ndarray, optional
        If there is an intersection returns True, else False. In
        case of an ndarray, it is an ndarray of booleans, 1-D, or
        2-D for a set of rays and a set of bounding boxes
    diag_calc : bool, optional
        This indicates whether diagonal calculations have been
        performed

    Returns
    -------
    Dataset
        Xarray dataset containing the intersection information.

        Key variables included:

        - **o**: The origin(s) of the ray(s) [xyz]
        - **d**: The direction(s) of the ray(s) [xyz]
        - **mint**: The mint attribute of the ray(s)
        - **maxt**: The maxt attribute of the ray(s)
        - **is_intersection**: If there is an intersection ->
          True, else False
        - **thit**: The t ray variable(s) of the intersection
          point(s)
        - **phit**: The intersection point(s) [xyz]

        In case of a set of rays and/or a set of bounding boxes,
        the variables get an extra nrays/nobj dimension.
    """
    is_r_arr = isinstance(r.o.x, np.ndarray)
    is_bbox_arr = isinstance(bbox.pmin.x, np.ndarray)

    ds = xr.Dataset(coords={"xyz": np.arange(3)})

    # bind defaults, reassigned below when is_r_arr is True
    nrays = 0
    ro = np.empty(0)
    rd = np.empty(0)
    mint = np.empty(0)
    maxt = np.empty(0)
    thit = np.empty(0)

    if is_r_arr:
        nrays = len(cast(np.ndarray, r.o.x))
        ro = r.o.to_numpy()
        rd = r.d.to_numpy()
        ds["o"] = xr.DataArray(ro, dims=["nrays", "xyz"])
        ds["d"] = xr.DataArray(rd, dims=["nrays", "xyz"])
        mint = np.full(nrays, r.mint, dtype=np.float64)
        maxt = np.full(nrays, r.maxt, dtype=np.float64)
        ds["mint"] = xr.DataArray(mint, dims=["nrays"])
        ds["maxt"] = xr.DataArray(maxt, dims=["nrays"])
    else:
        ds["o"] = xr.DataArray(r.o.to_numpy(), dims=["xyz"])
        ds["d"] = xr.DataArray(r.d.to_numpy(), dims=["xyz"])
        ds["mint"] = xr.DataArray(r.mint)
        ds["maxt"] = xr.DataArray(r.maxt)

    if is_r_arr or is_bbox_arr:
        t0_arr = cast(np.ndarray, t0)
        t1_arr = cast(np.ndarray, t1)
        c1 = t0_arr > 0
        not_c1 = np.logical_not(c1)
        thit = np.full(t0_arr.shape, np.nan, dtype=np.float64)
        c2 = np.logical_and(is_intersection, c1)
        c3 = np.logical_and(is_intersection, not_c1)
        if np.any(c2):
            thit[c2] = t0_arr[c2]
        if np.any(c3):
            thit[c3] = t1_arr[c3]

    if is_r_arr and is_bbox_arr and not diag_calc:
        nobj = len(cast(np.ndarray, bbox.p0.x))
        ds.attrs.update({"nobj": nobj, "nrays": nrays})
        ds["is_intersection"] = xr.DataArray(
            is_intersection, dims=["nobj", "nrays"]
        )
        ds["thit"] = xr.DataArray(thit, dims=["nobj", "nrays"])
        phit = np.zeros((nobj, nrays, 3), dtype=np.float64)
        for ir in range(0, nrays):
            ri = Ray(Point(ro[ir, :]), Vector(rd[ir, :]), mint[ir], maxt[ir])
            phit[:, ir, :] = ri[thit[:, ir]].to_numpy()
        ds["phit"] = xr.DataArray(phit, dims=["nobj", "nrays", "xyz"])
    elif is_r_arr or is_bbox_arr:
        if diag_calc:
            dim_name = "ndiag"
            size = nrays
            ds.attrs.update({"nobj": size, "nrays": size})
        elif is_r_arr:
            dim_name = "nrays"
            size = nrays
        else:
            dim_name = "nobj"
            size = len(cast(np.ndarray, bbox.p0.x))
        ds.attrs.update({dim_name: size})
        phit = r(thit).to_numpy()
        ds["is_intersection"] = xr.DataArray(is_intersection, dims=[dim_name])
        ds["thit"] = xr.DataArray(thit, dims=[dim_name])
        ds["phit"] = xr.DataArray(phit, dims=[dim_name, "xyz"])
    else:
        if t0 is None:
            thit = None
        elif t0 > 0:
            thit = t0
        else:
            thit = t1
        phit = r(cast(float, thit)).to_numpy()
        ds["is_intersection"] = xr.DataArray(is_intersection)
        ds["thit"] = xr.DataArray(thit)
        ds["phit"] = xr.DataArray(phit, dims=["xyz"])

    ds["o"].attrs = {
        "type": "Point",
        "description": "the x, y and z components of the ray point",
    }
    ds["d"].attrs = {
        "type": "Vector",
        "description": "the x, y and z components of the ray vector",
    }
    ds["mint"].attrs = {"description": "the mint attribut of the ray"}
    ds["maxt"].attrs = {"description": "the maxt attribut of the ray"}
    ds["is_intersection"].attrs = {
        "description": "this variable tells if there is an intersection "
        "between the ray and the shape"
    }
    ds["thit"].attrs = {
        "description": "the t ray factor for the intersection point "
        "calculation"
    }
    ds["phit"].attrs = {
        "type": "Point",
        "description": "the x, y and z components of the intersection point",
    }
    ds.attrs = {"shape": bbox.__class__.__name__}
    date = datetime.now().strftime("%Y-%m-%d")
    ds.attrs.update({"date": date, "version": VERSION})
    return ds


def print_basic(basic: Vector | Point | Normal, name: str = "") -> str:
    """
    :meta private:

    Parameters
    ----------
    basic : Vector or Point or Normal
        The basic object
    name : str, optional
        The str name to show at the start

    Returns
    -------
    str
        The return for the method __repr__ or __str___
    """
    if not isinstance(basic.x, np.ndarray):
        return f"{name}({basic.x}, {basic.y}, {basic.z})"
    else:
        size_name = len(name)
        if size_name > 0:
            first_space = f"{'':{size_name + 2}}"
        else:
            first_space = "  "
        bx, by, bz = _xyz_arrays(basic)
        ncomponents = len(bx)
        values = basic.to_numpy()
        space = np.empty_like(values, dtype=str)
        space[values >= 0] = " "
        space[values < 0] = ""
        fmt = basic.fmt
        output = ""
        if ncomponents <= 100:
            for i in range(0, ncomponents):
                if i == 0:
                    output += (
                        f"{name}([[{space[i, 0]}{bx[i]:{fmt}}, "
                        f"{space[i, 1]}{by[i]:{fmt}}, "
                        f"{space[i, 2]}{bz[i]:{fmt}}],\n"
                    )
                elif i == ncomponents - 1:
                    output += (
                        f"{first_space}[{space[i, 0]}{bx[i]:{fmt}}, "
                        f"{space[i, 1]}{by[i]:{fmt}}, "
                        f"{space[i, 2]}{bz[i]:{fmt}}]])"
                    )
                else:
                    output += (
                        f"{first_space}[{space[i, 0]}{bx[i]:{fmt}}, "
                        f"{space[i, 1]}{by[i]:{fmt}}, "
                        f"{space[i, 2]}{bz[i]:{fmt}}],\n"
                    )
        else:
            for i in range(0, 97):
                if i == 0:
                    output += (
                        f"{name}([[{space[i, 0]}{bx[i]:{fmt}}, "
                        f"{space[i, 1]}{by[i]:{fmt}}, "
                        f"{space[i, 2]}{bz[i]:{fmt}}],\n"
                    )
                else:
                    output += (
                        f"{first_space}[{space[i, 0]}{bx[i]:{fmt}}, "
                        f"{space[i, 1]}{by[i]:{fmt}}, "
                        f"{space[i, 2]}{bz[i]:{fmt}}],\n"
                    )
            output += "       ...\n"
            for i in range(ncomponents - 3, ncomponents):
                if i == ncomponents - 1:
                    output += (
                        f"{first_space}[{space[i, 0]}{bx[i]:{fmt}}, "
                        f"{space[i, 1]}{by[i]:{fmt}}, "
                        f"{space[i, 2]}{bz[i]:{fmt}}]])"
                    )
                else:
                    output += (
                        f"{first_space}[{space[i, 0]}{bx[i]:{fmt}}, "
                        f"{space[i, 1]}{by[i]:{fmt}}, "
                        f"{space[i, 2]}{bz[i]:{fmt}}],\n"
                    )
        return output
