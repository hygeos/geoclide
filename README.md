<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/hygeos/geoclide/refs/heads/main/geoclide/img/geoclide-logo-dark.png">
  <img alt="geoclide" src="https://raw.githubusercontent.com/hygeos/geoclide/refs/heads/main/geoclide/img/geoclide-logo-light.png" width="450">
</picture>
</p>

------------------------------------------------

<p align="center">
<a href="https://pypi.python.org/pypi/geoclide"><img alt="pypi" src="https://img.shields.io/pypi/v/geoclide.svg"></a>
<a href="https://anaconda.org/conda-forge/geoclide"><img alt="conda-forge" src="https://img.shields.io/conda/vn/conda-forge/geoclide.svg"></a>
<a href="https://github.com/hygeos/geoclide"><img alt="github" src="https://img.shields.io/github/v/tag/hygeos/geoclide?label=github&amp;color=blue"></a>
<a href="https://pepy.tech/project/geoclide"><img alt="downloads" src="https://static.pepy.tech/badge/geoclide"></a>
</p>

<p align="center">
<a href="https://github.com/hygeos/geoclide/actions/workflows/tests.yml"><img alt="tests" src="https://github.com/hygeos/geoclide/actions/workflows/tests.yml/badge.svg?branch=main"></a>
<a href="https://hygeos.github.io/geoclide/"><img alt="docs" src="https://github.com/hygeos/geoclide/actions/workflows/docs_github_pages.yml/badge.svg?branch=main"></a>
<a href="https://github.com/hygeos/geoclide/blob/main/LICENSE.txt"><img alt="license" src="https://img.shields.io/badge/license-Apache%202.0-blue.svg"></a>
<a href="https://pixi.sh"><img alt="pixi" src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json"></a>
<a href="https://github.com/astral-sh/ruff"><img alt="ruff" src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json"></a>
</p>

A python package for geometric calculations in the three-dimensional Euclidean space

Mustapha Moulana  
[HYGEOS website](https://hygeos.com/en/)  
[Documentation](https://hygeos.github.io/geoclide/)

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

Create a triangle mesh, place it in the scene with a transformation, then
intersect it with the same ray:
```python
>>> import numpy as np
>>> vertices = np.array([[-5., -5., 0.], [5., -5., 0.],
...                      [-5., 5., 0.], [5., 5., 0.]])
>>> faces = np.array([[0, 1, 2], [2, 3, 1]])
>>> translate = gc.get_translate_tf(gc.Vector(2.5, 0., 0.))
>>> rotate = gc.get_rotate_y_tf(-90.)
>>> mesh = gc.TriangleMesh(vertices, faces, otw=translate*rotate)
>>> ds = gc.calc_intersection(mesh, ray)
>>> ds['is_intersection'].values, ds['thit'].values
(array(True), array(4.5))
>>> ds['phit'].values
array([2.5, 0. , 0.8])
```

# Documentation
The complete documentation is available at
[hygeos.github.io/geoclide](https://hygeos.github.io/geoclide/). It
includes example notebooks (basics, remote sensing applications, quadrics
visualization and numpy acceleration) and the full API reference. The
docstrings are also available from the built-in `help` function, e.g.
`help(gc.calc_intersection)`.

# License
Geoclide is licensed under the Apache License 2.0, see
[LICENSE.txt](https://github.com/hygeos/geoclide/blob/main/LICENSE.txt).
