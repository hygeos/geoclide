# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.5
#   kernelspec:
#     display_name: geoclide-py313
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Acceleration with numpy

# %%
import os
import sys

sys.path.insert(0, os.path.abspath(".."))

from time import process_time

import numpy as np

import geoclide as gc

# %% [markdown]
# ## BBox - Ray intersection test, multiples bboxes and 1 ray

# %%
# Here we create 100000 bounding boxes and 1 ray
nx = 100
ny = 100
nz = 10
x = np.linspace(0., nx-1, nx, dtype=np.float64)
y = np.linspace(0., ny-1, ny, dtype=np.float64)
z = np.linspace(0., nz-1, nz, dtype=np.float64)
x_, y_, z_ = np.meshgrid(x,y,z, indexing='ij')
pmin_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
x = np.linspace(1., nx, nx, dtype=np.float64)
y = np.linspace(1., ny, ny, dtype=np.float64)
z = np.linspace(1., nz, nz, dtype=np.float64)
x_, y_, z_ = np.meshgrid(x,y,z, indexing='ij')
pmax_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
r0 = gc.Ray(gc.Point(-2., 0., 0.25), gc.normalize(gc.Vector(1, 0., 0.5)))

# %% [markdown]
# ### Intersection tests using a loop

# %%
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

# %% [markdown]
# ### Intersection tests using numpy

# %%
start = process_time()
pmin = gc.Point(pmin_arr)
pmax = gc.Point(pmax_arr)
b_set = gc.BBox(pmin, pmax)
t0, t1, is_int1 = b_set.intersect(r0, ds_output=False)
end = process_time()
print("elapsed time (s) using numpy: ", end - start)

# %% [markdown]
# ## BBox - Ray intersection test, multiples bboxes and multiple rays

# %% [markdown]
# ### Case 1: test each ray against all the bounding boxes

# %%
# We create 400 bounding boxes and 400 rays
nx = 20
ny = 20
nz = 1
x = np.linspace(0., nx-1, nx, dtype=np.float64)
y = np.linspace(0., ny-1, ny, dtype=np.float64)
z = np.linspace(0., nz-1, nz, dtype=np.float64)
x_, y_, z_ = np.meshgrid(x,y,z, indexing='ij')
pmin_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
x = np.linspace(1., nx, nx, dtype=np.float64)
y = np.linspace(1., ny, ny, dtype=np.float64)
z = np.linspace(1., nz, nz, dtype=np.float64)
x_, y_, z_ = np.meshgrid(x,y,z, indexing='ij')
pmax_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
nboxes = pmin_arr.shape[0]
x_, y_, z_ = np.meshgrid(np.linspace(0.5, nx-0.5, nx, dtype=np.float64),
                        np.linspace(0.5, ny-0.5, ny, dtype=np.float64),
                        nz+1, indexing='ij')

o_set_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
nrays = o_set_arr.shape[0]
d_set_arr = np.zeros_like(o_set_arr)
d_set_arr[:,0] = 0.
d_set_arr[:,1] = 0.
d_set_arr[:,2] = -1.
o_set = gc.Point(o_set_arr)
d_set = gc.Vector(d_set_arr)

# %% [markdown]
# #### The tests using loops

# %%
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
      t0_[ib,ir], t1_[ib,ir], is_int_[ib,ir] = bi.intersect(
          list_rays[ir], ds_output=False)
end = process_time()
print("case 1 - elapsed time (s) using loops:", end-start)

# %% [markdown]
# #### The tests using numpy calculations

# %%
start = process_time()
r_set = gc.Ray(o_set, d_set)
pmin = gc.Point(pmin_arr)
pmax = gc.Point(pmax_arr)
b_set = gc.BBox(pmin, pmax)
t0, t1, is_int1 = b_set.intersect(r_set, ds_output=False)
end = process_time()
time_fast = end-start
print("case 1 - elapsed time (s) using numpy:", end-start)

# %% [markdown]
# ### Case 2: diagonal calculations, only ray(i) with bbox(i)

# %%
# We create 40000 bounding boxes and 40000 rays
nx = 200
ny = 200
nz = 1
x = np.linspace(0., nx-1, nx, dtype=np.float64)
y = np.linspace(0., ny-1, ny, dtype=np.float64)
z = np.linspace(0., nz-1, nz, dtype=np.float64)
x_, y_, z_ = np.meshgrid(x,y,z, indexing='ij')
pmin_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
x = np.linspace(1., nx, nx, dtype=np.float64)
y = np.linspace(1., ny, ny, dtype=np.float64)
z = np.linspace(1., nz, nz, dtype=np.float64)
x_, y_, z_ = np.meshgrid(x,y,z, indexing='ij')
pmax_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
nboxes = pmin_arr.shape[0]
x_, y_, z_ = np.meshgrid(np.linspace(0.5, nx-0.5, nx, dtype=np.float64),
                        np.linspace(0.5, ny-0.5, ny, dtype=np.float64),
                        nz+1, indexing='ij')

o_set_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
nrays = o_set_arr.shape[0]
d_set_arr = np.zeros_like(o_set_arr)
d_set_arr[:,0] = 0.
d_set_arr[:,1] = 0.
d_set_arr[:,2] = -1.
o_set = gc.Point(o_set_arr)
d_set = gc.Vector(d_set_arr)

# %% [markdown]
# #### The tests using loop

# %%
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

# %% [markdown]
# #### The tests using numpy calculations

# %%
start = process_time()
r_set = gc.Ray(o_set, d_set)
pmin = gc.Point(pmin_arr)
pmax = gc.Point(pmax_arr)
b_set = gc.BBox(pmin, pmax)
t0, t1, is_int1 = b_set.intersect(r_set, diag_calc=True, ds_output=False)
end = process_time()
print("case 2 - elapsed time (s) using numpy:", end-start)
