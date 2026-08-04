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
# # Some Basics

# %%
import os
import sys

sys.path.insert(0, os.path.abspath(".."))

import numpy as np

import geoclide as gc

# %% [markdown]
# ## Create a point and a vector

# %%
p1 = gc.Point(0., 0., 0.) # create a point
# create a vector and normalize it
v1 = gc.normalize(gc.Vector(0.5, 0.5, 0.1))
p1, v1

# %%
v1.length()

# %% [markdown]
# ## Create a ray from the created point and vector

# %%
r1 = gc.Ray(o=p1, d=v1)
r1

# %%
r1(8)

# %% [markdown]
# ## Create a simple triangle mesh composed of 2 triangles

# %%
v0 = np.array([-5, -5, 0.])
v1 = np.array([5, -5, 0.])
v2 = np.array([-5, 5, 0.])
v3 = np.array([5, 5, 0.])
vertices = np.array([v0, v1, v2, v3])
f0 = np.array([0, 1, 2]) # the vertices indices of triangle 0 / face 0
f1 = np.array([2, 3, 1]) # the vertices indices of triangle 1 / face 1
faces = np.array([f0, f1])
# We can create a transformation to translate and rotate it
# translation of 2.5 in x axis
translate = gc.get_translate_tf(gc.Vector(2.5, 0., 0.))
# rotation of -90 degrees around the y axis
rotate = gc.get_rotate_y_tf(-90.)
# object to world transformation to apply to the triangle mesh
otw = translate*rotate
# create the triangle mesh
tri_mesh = gc.TriangleMesh(vertices, faces, otw=otw)
# see if the ray r1 intersect the triangle mesh
ds = gc.calc_intersection(tri_mesh, r1)
ds

# %%
ds['phit']
