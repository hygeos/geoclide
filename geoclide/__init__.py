"""
Geoclide
========

A python package for geometric calculations in the
three-dimensional Euclidean space.

Provides
  1. Basic geometric objects: vectors, points, normals, rays and
     bounding boxes
  2. Geometric transformations: translations, scales and
     rotations
  3. Ray intersection tests with shapes: spheres, spheroids,
     disks, triangles and triangle meshes

How to use the documentation
----------------------------
Documentation is available in two forms: docstrings provided with
the code, and a standalone reference guide, available from `the
geoclide homepage <https://hygeos.github.io/geoclide/>`_.

Code snippets are indicated by three greater-than signs::

    >>> import geoclide as gc
    >>> v = gc.Vector(0., 0., 1.)

Use the built-in ``help`` function to view a function's
docstring::

    >>> help(gc.calc_intersection)
"""

from geoclide.advancedvecope import ang2vec, vec2ang
from geoclide.basic import (
    BBox,
    Normal,
    Point,
    Ray,
    Vector,
    get_common_face,
    get_common_vertices,
)
from geoclide.constants import (
    GAMMA2_F32,
    GAMMA2_F64,
    GAMMA3_F32,
    GAMMA3_F64,
    GAMMA5_F32,
    GAMMA5_F64,
    TWO_PI,
    VERSION,
)
from geoclide.intersection import calc_intersection
from geoclide.mathope import clamp, gamma_f32, gamma_f64, quadratic
from geoclide.quadrics import Disk, Sphere, Spheroid
from geoclide.transform import (
    Transform,
    get_inverse_tf,
    get_rotate_tf,
    get_rotate_x_tf,
    get_rotate_y_tf,
    get_rotate_z_tf,
    get_scale_tf,
    get_translate_tf,
)
from geoclide.trianglemesh import (
    Triangle,
    TriangleMesh,
    create_disk_trianglemesh,
    create_sphere_trianglemesh,
    read_trianglemesh,
)
from geoclide.vecope import (
    coordinate_system,
    cross,
    distance,
    dot,
    face_forward,
    normalize,
    permute,
    vabs,
    vargmax,
    vargmin,
    vmax,
    vmin,
)

__all__ = [
    "GAMMA2_F32",
    "GAMMA2_F64",
    "GAMMA3_F32",
    "GAMMA3_F64",
    "GAMMA5_F32",
    "GAMMA5_F64",
    "TWO_PI",
    "VERSION",
    "BBox",
    "Disk",
    "Normal",
    "Point",
    "Ray",
    "Sphere",
    "Spheroid",
    "Transform",
    "Triangle",
    "TriangleMesh",
    "Vector",
    "ang2vec",
    "calc_intersection",
    "clamp",
    "coordinate_system",
    "create_disk_trianglemesh",
    "create_sphere_trianglemesh",
    "cross",
    "distance",
    "dot",
    "face_forward",
    "gamma_f32",
    "gamma_f64",
    "get_common_face",
    "get_common_vertices",
    "get_inverse_tf",
    "get_rotate_tf",
    "get_rotate_x_tf",
    "get_rotate_y_tf",
    "get_rotate_z_tf",
    "get_scale_tf",
    "get_translate_tf",
    "normalize",
    "permute",
    "quadratic",
    "read_trianglemesh",
    "vabs",
    "vargmax",
    "vargmin",
    "vec2ang",
    "vmax",
    "vmin",
]
