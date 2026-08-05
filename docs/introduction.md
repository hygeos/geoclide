A python package for geometric calculations in the three-dimensional Euclidean space

Mustapha Moulana  
[HYGEOS website](https://hygeos.com/en/)

-----------------------------------------

# Features
- Basic geometric objects: vectors, points, normals, rays and bounding boxes
- Geometric transformations: translations, scales and rotations
- Ray intersection tests with shapes: spheres, spheroids, disks, triangles
  and triangle meshes
- Vectorized calculations with numpy: sets of objects and/or sets of rays
  can be processed at once
- Intersection results returned as xarray datasets, gathering the
  intersection information and the shape attributes
- Visualization of the quadrics and triangle meshes, and reading/writing of
  triangle meshes (netcdf4, and formats supported by trimesh as stl, ply, ...)

# Quickstart
Perform an intersection test between a ray and a sphere, and get the
intersection information as an xarray dataset:
```python
>>> import geoclide as gc
>>> sphere = gc.Sphere(radius=1.)
>>> ray = gc.Ray(o=gc.Point(-2., 0., 0.8), d=gc.Vector(1., 0., 0.))
>>> ds = gc.calc_intersection(sphere, ray)
>>> ds['thit'].values, ds['phit'].values
(array(1.4), array([-0.6,  0. ,  0.8]))
```
More complete walkthroughs are available in the [Examples](examples.rst)
section.

# Documentation contents

## Examples
Four example notebooks, from the first steps to more applied use cases:

- [Some Basics](01_basic_example.ipynb) — create points, vectors and rays,
  build a small triangle mesh with a transformation and perform a first
  intersection test
- [Examples for remote sensing applications](02_remote_sensing_example.ipynb)
  — compute a satellite position from viewing angles (flat and spherical
  earth) and the pixel directions of a satellite camera
- [How to create and visualize quadrics](03_quadrics_example.ipynb) —
  construct and plot disks, annuli, spheres, partial spheres and spheroids
- [Acceleration with numpy](04_numpy_acceleration_example.ipynb) —
  vectorized intersection tests with sets of rays and/or sets of objects,
  including the diagonal calculation mode

## Geoclide package
The [API reference](geoclide.rst), organized by module:

- [basic](geoclide.basic.rst) — the Vector, Point, Normal, Ray and BBox
  classes
- [vecope](geoclide.vecope.rst) and
  [advancedvecope](geoclide.advancedvecope.rst) — operations on vectors,
  points and normals, and conversions between angles and directions
- [transform](geoclide.transform.rst) — the Transform class and the
  translation, scale and rotation transformations
- [quadrics](geoclide.quadrics.rst) — the Sphere, Spheroid and Disk shapes
- [trianglemesh](geoclide.trianglemesh.rst) — the Triangle and TriangleMesh
  shapes
- [intersection](geoclide.intersection.rst) — the calc_intersection
  function, standardized intersection tests for any shape
- [shapes](geoclide.shapes.rst) — the Shape base class and the intersection
  dataset builder
- [mathope](geoclide.mathope.rst) — basic mathematical utilities

## Releases
The [Releases](changelog_link.rst) section lists the versions of geoclide
and the changes they introduced.

# Installation
The installation can be performed using one of the following commands:
```shell
$ conda install -c conda-forge geoclide
```
```shell
$ pip install geoclide
```
```shell
$ pip install git+https://github.com/hygeos/geoclide.git
```

# Testing
Run the command `pytest geoclide/tests/ -s -v` to check that everything is running correctly.
