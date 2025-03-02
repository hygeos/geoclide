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
  | `create_disk_trianglemesh` | Function | create a disk / partial disk / annulus / partial annulus triangleMesh |
  | `read_trianglemesh` | Function | read mesh file (gcnc, stl, obj, ...) and convert it to a TriangleMesh class object |

</details>


## Examples

### Basic example
<details>
  <summary>Click here</summary>

  Create a point and a Vector
  ```python
  >>> import geoclide as gc
  >>> import numpy as np
  >>> 
  >>> p1 = gc.Point(0., 0., 0.) # create a point
  >>> v1 = gc.normalize(gc.Vector(0.5, 0.5, 0.1)) # create a vector and normalize it
  >>> p1, v1
  (Point(0.0, 0.0, 0.0), Vector(0.7001400420140049, 0.7001400420140049, 0.140028008402801))
  >>> v1.length()
  1.0
  ```

  With a point and a vector we can create a ray
  ```python
  >>> r1 = gc.Ray(o=p1, d=v1)
  >>> r1
  (0.0, 0.0, 0.0) + t*(0.7001400420140049, 0.7001400420140049, 0.140028008402801) with t ∈ [0,inf[
  ```

  Create a simple triangle mesh composed of 2 triangles
  ```python
  >>> v0 = np.array([-5, -5, 0.])
  >>> v1 = np.array([5, -5, 0.])
  >>> v2 = np.array([-5, 5, 0.])
  >>> v3 = np.array([5, 5, 0.])
  >>> vertices = np.array([v0, v1, v2, v3])
  >>> f0 = np.array([0, 1, 2]) # the vertices indices of triangle 0 / face 0
  >>> f1 = np.array([2, 3, 1]) # the vertices indices of triangle 1 / face 1
  >>> faces = np.array([f0, f1])
  >>> # We can create a transformation to translate and rotate it
  >>> translate = gc.get_translate_tf(gc.Vector(2.5, 0., 0.)) # translation of 2.5 in x axis
  >>> rotate = gc.get_rotateY_tf(-90.) # rotation of -90 degrees around the y axis
  >>> oTw = translate*rotate # object to world transformation to apply to the triangle mesh
  >>> tri_mesh = gc.TriangleMesh(vertices, faces, oTw=oTw) # create the triangle mesh
  >>> ds = gc.calc_intersection(tri_mesh, r1) # see if the ray r1 intersect the triangle mesh
  >>> ds
  <xarray.Dataset> Size: 865B
  Dimensions:          (xyz: 3, nvertices: 4, ntriangles: 2, p0p1p2: 3, dim_0: 4,
                        dim_1: 4)
  Coordinates:
    * xyz              (xyz) int64 24B 0 1 2
  Dimensions without coordinates: nvertices, ntriangles, p0p1p2, dim_0, dim_1
  Data variables: (12/18)
      o                (xyz) float64 24B 0.0 0.0 0.0
      d                (xyz) float64 24B 0.7001 0.7001 0.14
      mint             int64 8B 0
      maxt             float64 8B inf
      is_intersection  bool 1B True
      thit             float64 8B 3.571
      ...               ...
      vertices         (nvertices, xyz) float64 96B -5.0 -5.0 0.0 ... 5.0 5.0 0.0
      faces            (ntriangles, p0p1p2) int64 48B 0 1 2 2 3 1
      wTo_m            (dim_0, dim_1) float64 128B 6.123e-17 0.0 1.0 ... 0.0 1.0
      wTo_mInv         (dim_0, dim_1) float64 128B 6.123e-17 0.0 -1.0 ... 0.0 1.0
      oTw_m            (dim_0, dim_1) float64 128B 6.123e-17 0.0 -1.0 ... 0.0 1.0
      oTw_mInv         (dim_0, dim_1) float64 128B 6.123e-17 0.0 1.0 ... 0.0 1.0
  Attributes:
      shape:       TriangleMesh
      date:        2025-03-02
      version:     2.0.0
      ntriangles:  2
      nvertices:   4
  ```

  We can access the details of each variable. Bellow an example with the variable 'phit'.
  ```python
  >>> ds['phit']
  <xarray.DataArray 'phit' (xyz: 3)> Size: 24B
  array([2.5, 2.5, 0.5])
  Coordinates:
    * xyz      (xyz) int64 24B 0 1 2
  Attributes:
      type:         Point
      description:  the x, y and z components of the intersection point
  ```
</details>

### Example for remote sensing applications
<details>
  <summary>Click here</summary>


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
  ```bash
  Satellite position (pp case) : [ 494.97474683 -494.97474683  700.        ]
  Satellite position (sp case)  [ 472.61058011 -472.61058011  668.37229212]
  ```
</details>


### How to create quadrics (disk, sphere and spheroid)

<details>
  <summary>Click here</summary>

  #### disk, annulus and partial annulus

  ```python

  >>> import geoclide as gc
  >>> disk = gc.Disk(radius=1.)
  >>> disk.plot(color='green', edgecolor='k')
  >>> annulus = gc.Disk(radius=1., inner_radius=0.5)
  >>> annulus.plot(color='green', edgecolor='k')
  >>> partial_annulus = gc.Disk(radius=1., inner_radius=0.5, phimax=270)
  >>> partial_annulus.plot(color='green', edgecolor='k')
  ```

  <p align="center">
  <img src="geoclide/img/disk.png" width="250">
  <img src="geoclide/img/annulus.png" width="250">
  <img src="geoclide/img/partial_annulus.png" width="250">
  </p>


  #### sphere and partial spheres

  ```python
  >>> import geoclide as gc
  >>> sphere = gc.Sphere(radius=1.)
  >>> sphere.plot(color='blue', edgecolor='k')
  >>> partial_sphere1 = gc.Sphere(radius=1., zmax=0.5)
  >>> partial_sphere1.plot(color='blue', edgecolor='k')
  >>> partial_sphere2 = gc.Sphere(radius=1., zmax=0.5, phimax=180.)
  >>> partial_sphere2.plot(color='blue', edgecolor='k')
  ```
  <p align="center">
  <img src="geoclide/img/sphere.png" width="250">
  <img src="geoclide/img/sphere_partial1.png" width="250">
  <img src="geoclide/img/sphere_partial2.png" width="250">
  </p>


  #### spheroid (prolate and oblate)

  ```python
  >>> import geoclide as gc
  >>> prolate = gc.Spheroid(radius_xy=1, radius_z=3)
  >>> prolate.plot(color='red', edgecolor='k')
  >>> oblate = gc.Spheroid(radius_xy=1, radius_z=0.8)
  >>> oblate.plot(color='cyan', edgecolor='k')
  ```
  <p align="center">
  <img src="geoclide/img/prolate.png" width="250">
  <img src="geoclide/img/oblate.png" width="250">
  </p>
</details>

### Accelerate the calculations using numpy ndarray (only since geoclide 2.0.0)

#### Bounding BBox - Ray intersection test (multiples bboxes and 1 ray)
<details>
  <summary>Click here</summary>

  Here we create 1000000 bounding boxes and 1 ray
  ```python
  import numpy as np
  import geoclide as gc
  from time import process_time
  nx = 100
  ny = 100
  nz = 100
  x = np.linspace(0., nx-1, nx, np.float64)
  y = np.linspace(0., ny-1, ny, np.float64)
  z = np.linspace(0., nz-1, nz, np.float64)
  x_, y_, z_ = np.meshgrid(x,y,z, indexing='ij')
  pmin_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
  x = np.linspace(1., nx, nx, np.float64)
  y = np.linspace(1., ny, ny, np.float64)
  z = np.linspace(1., nz, nz, np.float64)
  x_, y_, z_ = np.meshgrid(x,y,z, indexing='ij')
  pmax_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
  r0 = gc.Ray(gc.Point(-2., 0., 0.25), gc.normalize(gc.Vector(1, 0., 0.5)))

  # Test intersection tests using a loop
  start = process_time()
  nboxes = pmin_arr.shape[0]
  t0_ = np.zeros(nboxes, dtype=np.float64)
  t1_ = np.zeros_like(t0_)
  is_int_ = np.full(nboxes, False, dtype=bool)
  for ib in range (0, nboxes):
      bi = gc.BBox(gc.Point(pmin_arr[ib,:]), gc.Point(pmax_arr[ib,:]))
      t0_[ib], t1_[ib], is_int_[ib] = bi.intersect(r0, ds_output=False)
  end = process_time()
  print("elapsed time (s) using loop: ", end - start)

  #Test intersection tests using ndarray calculations
  start = process_time()
  pmin = gc.Point(pmin_arr)
  pmax = gc.Point(pmax_arr)
  b_set = gc.BBox(pmin, pmax)
  t0, t1, is_int1 = b_set.intersect(r0, ds_output=False)
  end = process_time()
  print("elapsed time (s) using numpy: ", end - start)
  ```
  ``` bash
  elapsed time (s) using loop:  7.527952158
  elapsed time (s) using numpy:  0.06523970699999992
  ```
  
  In this example, we are approximately 115 times faster by using numpy ndarray calculations.
  </details>

  #### Bounding BBox - Ray intersection test (multiples bboxes and multiple rays)
  <details>
  <summary>Click here</summary>

  We create 10000 bounding boxes and 10000 rays
  ```python
  import numpy as np
  import geoclide as gc
  from time import process_time
  nx = 100
  ny = 100
  nz = 1
  x = np.linspace(0., nx-1, nx, np.float64)
  y = np.linspace(0., ny-1, ny, np.float64)
  z = np.linspace(0., nz-1, nz, np.float64)
  x_, y_, z_ = np.meshgrid(x,y,z, indexing='ij')
  pmin_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
  x = np.linspace(1., nx, nx, np.float64)
  y = np.linspace(1., ny, ny, np.float64)
  z = np.linspace(1., nz, nz, np.float64)
  x_, y_, z_ = np.meshgrid(x,y,z, indexing='ij')
  pmax_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
  nboxes = pmin_arr.shape[0]
  x_, y_, z_ = np.meshgrid(np.linspace(0.5, nx-0.5, nx, np.float64),
                          np.linspace(0.5, ny-0.5, ny, np.float64),
                          nz+1, indexing='ij')

  o_set_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
  nrays = o_set_arr.shape[0]
  d_set_arr = np.zeros_like(o_set_arr)
  d_set_arr[:,0] = 0.
  d_set_arr[:,1] = 0.
  d_set_arr[:,2] = -1.
  o_set = gc.Point(o_set_arr)
  d_set = gc.Vector(d_set_arr)

  # === Case 1: for each ray, perform intersection test with all the bounding boxes
  # The tests using loops
  start = process_time()
  t0_ = np.zeros((nboxes, nrays), dtype=np.float64)
  t1_ = np.zeros_like(t0_)
  is_int_ = np.full((nboxes,nrays), False, dtype=bool)
  list_rays = []
  for ir in range(0, nrays):
    list_rays.append(gc.Ray(gc.Point(o_set_arr[ir,:]),
                            gc.normalize(gc.Vector(d_set_arr[ir,:]))))
  for ib in range (0, nboxes):
    bi = gc.BBox(gc.Point(pmin_arr[ib,:]), gc.Point(pmax_arr[ib,:]))
    for ir in range(0, nrays):
        t0_[ib,ir], t1_[ib,ir], is_int_[ib,ir] = bi.intersect(list_rays[ir], ds_output=False)
  end = process_time()
  print("case 1 - elapsed time (s) using loops:", end-start)

  # The tests using numpy calculations
  start = process_time()
  r_set = gc.Ray(o_set, d_set)
  pmin = gc.Point(pmin_arr)
  pmax = gc.Point(pmax_arr)
  b_set = gc.BBox(pmin, pmax)
  t0, t1, is_int1 = b_set.intersect(r_set, ds_output=False)
  end = process_time()
  time_fast = end-start
  print("case 1 - elapsed time (s) using numpy:", end-start)

  # === Case 2: perform intersection test only between ray(i) and bbox(i) i.e. diagonal calculations
  # The tests using lo
  start = process_time()
  t0_ = np.zeros((nboxes), dtype=np.float64)
  t1_ = np.zeros_like(t0_)
  is_int_ = np.full((nboxes), False, dtype=bool)
  list_rays = []
  for ib in range(0, nboxes):
      bi = gc.BBox(gc.Point(pmin_arr[ib,:]), gc.Point(pmax_arr[ib,:]))
      ri = gc.Ray(gc.Point(o_set_arr[ib,:]), gc.Vector(d_set_arr[ib,:]))
      t0_[ib], t1_[ib], is_int_[ib] = bi.intersect(ri, ds_output=False)
  end = process_time()
  print("case 2 - elapsed time (s) using loops:", end-start)

  # The tests using numpy calculations
  start = process_time()
  r_set = gc.Ray(o_set, d_set)
  pmin = gc.Point(pmin_arr)
  pmax = gc.Point(pmax_arr)
  b_set = gc.BBox(pmin, pmax)
  t0, t1, is_int1 = b_set.intersect(r_set, diag_calc=True, ds_output=False)
  end = process_time()
  print("case 2 - elapsed time (s) using numpy:", end-start)
  ```
  ``` bash
  case 1 - elapsed time (s) using loop: 201.16944833099998
  case 1 - elapsed time (s) using numpy: 3.383698623000015
  case 2 - elapsed time (s) using loop: 0.10501369499999669
  case 2 - elapsed time (s) using numpy: 0.04252339900000379
  ```
</details>

