
# GEOCLIDE CHANGELOG


## v3.0.3
Release date: 07-06-2025

* Fix bug in `vec2ang` function with float32 x, y, z vector components
  - Also add acc parameter to control the tolerance of numerical errors


## v3.0.2
Release date: 06-06-2025

* Fix bug in `vec2ang` function


## v3.0.1
Release date: 13-03-2025

* Several doctring corrections

* Clearer docsting

* Improvement of the sphinx docs

* A mistake has been found concerning the parameter name for diagonal calculations 
  only for the the Transform `__call__` method and the function `ang2vec` 
  - The name calc_diag has been replaced by the correct one i.e. diag_calc

* Add missing diag_calc parameter to Transform `rotate` method

* Correct the Transform `scale` method
  - the method uses the `get_scale_tf` function, which takes only one parameter
    of type Vector instead of three


## v3.0.0
Release date: 12-03-2025

This release includes only a breacking change with the function `get_scale_tf`.
The Transform and Ray `__call__` method is now recommended in replacement 
of the method `__getitem__`. The method `__getitem__` is still here but has a
depracated message

* Improve the `__str__` and `__repr__` methods in case of multiple vectors, points, 
  normals or rays

* Replace the Transform and Ray classes method `__getitem__` by the method `__call__`

* Multiple transformations is now possible with the Transform class
  - Allow 1-D calculations with several transformations and 1 point, vectors, ...
  - Allow 2-D calculations with several transformations and several points, vectors, ...
  - Allow diagonal calculations (case 2-D with option diag_calc=True)

* Add sphinx documentation

* Function `get_scale_tf` takes as input a Vector variable instead of 3 scalars

* Add new option flatten to Tranform `__call__` method
  - Instead of having a numpy array of Vector objects, the output can be only 1 Vector 
    object where a flatten operation is made to x, y, z components of all vectors
  - Avoid loop operation then reduce the computational time

* Function `vec2ang` now accepts Vector object with multiple vectors

* Function `ang2vec` now accepts 1d ndarray theta and phi angles


## v2.0.2
Release date: 05-03-2025

* Several corrections in TriangleMesh `intersect`, `is_intersection_t` and 
  `is_intersection` methods

* Transformations were not taken into account when using the plot method, 
  now it's corrected

* Addition of more tests


## v2.0.1
Release date: 03-03-2025

* Correction of a bug occuring while using a version of numpy < 2.0.0
  - np.atan2 and np.acos replaced by np.arctan2 and np.arccos

* README image links updated

* Correct file format and add test files in MANIFEST file

* Ignore certain unwanted warnings


## v2.0.0
Release date: 02-03-2025

This release includes some breaking changes like the input of TriangleMesh where
the variable names are different, but also the input (2d ndarray). Other breaks are
for example changes of certain method outputs.

* Acceleration of geometric calculations thanks to numpy multi-dimentional array operations
    - Point, Vector and Normal classes can now have 1-D ndarray x, y, z components
    - Functions in vecope.py have been adapted to allow Point, Vector and Normal
      with 1-D ndarray x, y, z components
    - Readaptation of the Ray and BBox classes with the new Point, Vector and Normal
      classes
    - BBox `intersect` method allow 1-D, 2-D and 2-D diagonal calculations
    - Readaptation of the quadric classes enabling intersection tests with a class
      Ray containing several rays
    - Readaptation of Triangle and TriangleMesh classes. Allow 1-D, 2-D and 2-D
      diagonal calculations
    - ...

* The 3D objects `intersect` method output is now an xarray Dataset

* Several new functionalities
    - New method `to_triangle_mesh` to convert a quadric to a TriangleMesh object
    - New method `to_plot` for TriangleMesh, and also for all quadrics
    - New TriangleMesh method `write` to save the triangle mesh in nectcdf4 and also in
    - several other formats as stl, ply, ...
    - The `read_trianglemesh` function allows to read triangle mesh files (gcnc, stl, ply,
      etc.) and return a TriangleMesh object
    - ...

* A lot of tests have been added, more documentations, cleaning, ...


## v1.2.2
Release date: 13-02-2025

* Correction of a bug in TriangleMesh constructor
    - The bug occurred when we set the vi parameter as a 2d array.

* Add tests detecting the corrected bug

* Add docstring to area method for all objects


## v1.2.1
Release date: 13-02-2025

* Correction in `vec2ang` function
    - small correction in the docstring
    - correct a bug appearing when the parameter v is not a normalized vector
    - add the missing default value of the parameter vec_view

* Add tests enabling to detect the corrected bug


## v1.2.0
Release date: 10-02-2025

* Add new functions
    - `ang2vec`
    - `vec2ang`

* Add the Disk class
    - It can be a disk, a partial disk, an annulus or a partial annulus

* Add new method to all 3d objects `is_intersection_t`
    - the method is faster than `intersect` method, but gives only thit

* Move tests folder into geoclide folder
* And some improvements


## v1.1.1
Release date: 07-02-2025

* Correction of a bug with the normal orientation at the intersection with a 3d object


## v1.1.0
Release date: 05-02-2025

* Add the Spheroid class (oblate or prolate object)
* Small optimization of sphere intersection test


## v1.0.0

Release date: 2025-01-27

First public release.