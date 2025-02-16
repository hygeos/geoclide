<p align="center">
<img src="geoclide/img/geoclide_logo.png" width="300">
</p>

------------------------------------------------

[![image](https://img.shields.io/pypi/v/geoclide.svg)](https://pypi.python.org/pypi/geoclide)
[![image](https://img.shields.io/conda/vn/conda-forge/geoclide.svg)](https://anaconda.org/conda-forge/geoclide)
[![image](https://pepy.tech/badge/geoclide)](https://pepy.tech/project/geoclide)

The python package for geometric calculations in the three-dimentional Euclidian space

Mustapha Moulana  
[HYGEOS](https://hygeos.com/en/)

-----------------------------------------

## Installation
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

## Testing
Run the command `pytest geoclide/tests/ -s -v` to check that everything is running correctly.

## Available classes/functions
<details>
  <summary>Click here</summary>

  | Class/Function | Type | Description |
  | -------------- | ---- | ----------- |
  | `Vector`| Class | vector with x, y and z components |
  | `Point` | Class | point with x, y and z components |
  | `Normal` | Class | normal with x, y and z components |
  | `Ray` | Class | the ray: r(t) = o + t*d, with 'o' a Point, 'd' a vector and t ∈ [0,inf[ |
  | `BBox` | Class | bounding box |
  | `Sphere` | Class | sphere object. It can be a partial sphere|
  | `Spheroid` | Class | spheroid object (oblate or prolate) |
  | `Disk` | Class | disk object. It can be a partial disk or an annulus/partial annulus |
  | `Triangle` | Class | triangle object |
  | `TriangleMesh` | Class | triangle mesh object |
  | `Transform` | Class | transformation to translate and/or rotate every objects except a BBox |
  | `calc_intersection` | Function | intersection test between a shape and a ray and returns dataset |
  | `get_common_vertices` | Function | gives the vertices of BBox b1 which are common to another BBox b2 |
  | `get_common_face` | Function | same as `get_common_vertices` but with faces |
  | `dot` | Function | dot product (only vector or normal) |
  | `cross` | Function | cross product (only vector or normal) |
  | `normalize` | Function | normalize a vector/normal |
  | `coordinate_system` | Function | from a vector v1 compute vectors v2 and v3 such that v1, v2 and v3 are unit vectors of an orthogonal coordinate system |
  | `distance` | Function | compute the distance between 2 points |
  | `face_forward` | Function | ensure a vector/normal is in the same hemipherical direction than another given vector/normal |
  | `vmax` | Function | largest component value of the vector/point/normal |
  | `vmin` | Function | smallest component value of the vector/point/normal |
  | `vargmax` | Function | index of the vector/point/normal component with the largest value |
  | `vargmin` | Function | index of the vector/point/normal component with the smallest value |
  | `vabs` | Function | absolute value of each components of the vector/point/normal |
  | `permute` | Function | permutes the vector/point/normal values according to the given indices |
  | `clamp` | Function | clamp a value into the range [val_min, val_max] |
  | `quadratic` | Function | resolve the quadratic polynomial: ax**2 + bx + c |
  | `gamma_f32` | Function | gamma function from pbrt v3 |
  | `gamma_f64` | Function | gamma function from pbrt v3 but in double precision |
  | `get_inverse_tf` | Function | get the inverse transform from a another transform |
  | `get_translate_tf` | Function | get the translate transfrom from a given vector |
  | `get_scale_tf` | Function | get scale transform giving factors in x, y and z |
  | `get_rotateX_tf` | Function | get the rotate (around x axis) transform from scalar in degrees |
  | `get_rotateY_tf` | Function | get the rotate (around y axis) transform from scalar in degrees |
  | `get_rotateZ_tf` | Function | get the rotate (around z axis) transform from scalar in degrees |
  | `get_rotate_tf` | Function | get the rotate transform around a given vector/normal |
  | `ang2vec` | Function | convert a direction described by 2 angles into a direction described by a vector |
  | `vec2ang` | Function | convert a direction described by a vector into a direction described by 2 angles |
  | `create_sphere_trianglemesh` | Function | create a sphere / partial sphere triangleMesh |
</details>


## Examples
<details>
  <summary>Basic exemple</summary>

  ```python
  >>> import geoclide as gc
  >>> import numpy as np
  >>> # Some basics
  >>> p1 = gc.Point(0., 0., 0.) # create a point
  >>> v1 = gc.normalize(gc.Vector(0.5, 0.5, 0.1)) # create a vector and normalize it
  >>> v1
  Vector(0.4082482904638631, 0.4082482904638631, 0.8164965809277261)
  >>> # With a point and a vector we can create a ray
  >>> r1 = gc.Ray(o=p1, d=v1)
  >>> r1
  r(t) = (0.0, 0.0, 0.0) + t*(0.4082482904638631, 0.4082482904638631, 0.8164965809277261) with t ∈ [0,inf[
  >>> # Let's create a triangle mesh with 2 triangles
  >>> # We have 4 vertices
  >>> p0 = gc.Point(-5, -5, 0.)
  >>> p1 = gc.Point(5, -5, 0.)
  >>> p2 = gc.Point(-5, 5, 0.)
  >>> p3 = gc.Point(5, 5, 0.)
  >>> v = np.array([p0, p1, p2, p3], dtype=gc.Point)
  >>> # Get the vertices indices of each triangle
  >>> vid_t0 = np.array([0, 1, 2]) # the vertices indices of triangle 0
  >>> vid_t1 = np.array([2, 3, 1]) # the vertices indices of triangle 1
  >>> vi = np.concatenate((vid_t0, vid_t1)) # regroup everything
  >>> # Here if we create the triangle mesh, it would be a square of dimension 10*10
  >>> # centered at origin (0.,0.,0.) and parallel to the xy plane
  >>> # We can create a transformation to translate and rotate it
  >>> translate = gc.get_translate_tf(gc.Vector(2.5, 0., 0.)) # translation of 2.5 in x axis
  >>> rotate = gc.get_rotateY_tf(-90.) # rotation of -90 degrees around the y axis
  >>> oTw = translate*rotate # object to world transformation to apply to the triangle mesh
  >>> tri_mesh = gc.TriangleMesh(vi, v, oTw=oTw) # create the triangle mesh
  >>> ds = gc.calc_intersection(tri_mesh, r1) # see if the ray r1 intersect the triangle mesh
  >>> ds
  <xarray.Dataset> Size: 801B
  Dimensions:          (xyz: 3, nvertices: 4, ntriangles: 2, p0p1p2: 3, dim_0: 4,
                        dim_1: 4)
  Coordinates:
    * xyz              (xyz) int64 24B 0 1 2
  Dimensions without coordinates: nvertices, ntriangles, p0p1p2, dim_0, dim_1
  Data variables: (12/14)
      is_intersection  bool 1B True
      o                (xyz) float64 24B 0.0 0.0 0.0
      d                (xyz) float64 24B 0.7001 0.7001 0.14
      mint             int64 8B 0
      maxt             float64 8B inf
      v                (nvertices, xyz) float64 96B -5.0 -5.0 0.0 ... 5.0 5.0 0.0
      ...               ...
      wTo_mInv         (dim_0, dim_1) float64 128B 6.123e-17 0.0 -1.0 ... 0.0 1.0
      oTw_m            (dim_0, dim_1) float64 128B 6.123e-17 0.0 -1.0 ... 0.0 1.0
      oTw_mInv         (dim_0, dim_1) float64 128B 6.123e-17 0.0 1.0 ... 0.0 1.0
      thit             float64 8B 3.571
      phit             (xyz) float64 24B 2.5 2.5 0.5
      nhit             (xyz) float64 24B 1.0 0.0 -8.882e-17
  Attributes:
      shape:       TriangleMesh
      ntriangles:  2
      nvertices:   4
      date:        2025-01-27
      version:     1.0.0
  >>> # Here there is intersection, see more detail on intersection point phit
  >>> ds['phit']
  <xarray.DataArray 'phit' (xyz: 3)> Size: 24B
  array([2.5, 2.5, 0.5])
  Coordinates:
    * xyz      (xyz) int64 24B 0 1 2
  Attributes:
      type:         Point
      description:  the x, y and z components of the intersection point
  >>> # We can convert it into a point object
  >>> phit = gc.Point(ds['phit'].values)
  >>> phit
  Point(2.5, 2.5, 0.5)
  ```
</details>

<details>
  <summary>Example for remote sensing applications</summary>


  ```python
  import geoclide as gc
  import math

  # Find satellite x an y positions knowing its altitude and its viewing zenith and azimuth angles
  vza = 45. # viewing zenith angle in degrees
  vaa = 45. # viewing azimuth angle in degrees
  sat_altitude = 700.  # satellite altitude in kilometers
  origin = gc.Point(0., 0., 0.) # origin is the viewer seeing the satellite
  # The vaa start from north going clockwise.
  # Let's assume that in our coordinate system the x axis is in the north direction
  # Then theta (zenith) angle = vza and phi (azimuth) angle = -vaa
  theta = vza
  phi = -vaa

  # Get the vector from ground to the satellite
  dir_to_sat = gc.ang2vec(theta=theta, phi=phi)
  ray = gc.Ray(o=origin, d=dir_to_sat) # create the ray, starting from origin going in dir_to_sat direction

  # Here without considering the sphericity of the earth
  b1 = gc.BBox(p1=gc.Point(-math.inf, -math.inf, 0.), p2=gc.Point(math.inf, math.inf, sat_altitude))
  ds_pp = gc.calc_intersection(b1, ray) # return an xarray dataset

  # Here with the consideration of the sphericity of the earth
  earth_radius = 6378. # the equatorial earth radius in kilometers
  oTw = gc.get_translate_tf(gc.Vector(0., 0., -earth_radius))
  sphere_sat_alti = gc.Sphere(radius=earth_radius+sat_altitude, oTw=oTw)  # apply oTw to move the sphere center to earth center
  ds_sp = gc.calc_intersection(sphere_sat_alti, ray) # return an xarray dataset

  print ("Satellite position (pp case) :", ds_pp['phit'].values)
  print ("Satellite position (sp case) ", ds_sp['phit'].values)
  ```
</details>