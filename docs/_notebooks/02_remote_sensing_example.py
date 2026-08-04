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
# # Examples for remote sensing applications

# %%
import os, sys
sys.path.insert(0, os.path.abspath(".."))

import geoclide as gc
import math
import numpy as np

# %% [markdown]
# ## Find the x and y components of the satellite position knowing its altitude and its viewing zenith and azimuth angles

# %%
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

# %% [markdown]
# ### Here without considering the sphericity of the earth

# %%
b1 = gc.BBox(p1=gc.Point(-math.inf, -math.inf, 0.), p2=gc.Point(math.inf, math.inf, sat_altitude))
ds_pp = gc.calc_intersection(b1, ray) # return an xarray dataset
ds_pp['phit'].values

# %% [markdown]
# ### Here with the consideration of the sphericity of the earth

# %%
earth_radius = 6378. # the equatorial earth radius in kilometers
otw = gc.get_translate_tf(gc.Vector(0., 0., -earth_radius))
sphere_sat_alti = gc.Sphere(radius=earth_radius+sat_altitude, otw=otw)  # apply otw to move the sphere center to earth center
ds_sp = gc.calc_intersection(sphere_sat_alti, ray) # return an xarray dataset
ds_sp['phit'].values

# %% [markdown]
# ## Satellite camera directions (3MI example)

# %% [markdown]
# ### Compute all the pixel directions

# %%
# Number of camera pixels in x and y
nx = 509
ny = 255

# Satellite altitude (in km)
z_sat = 830.

# Nadir resolution (in km2)
nad_res = 4.

# swath (in km2)
swath = 2200

# Field of view of 1 pixel and of the camera
pixel_fov = np.rad2deg(np.arctan(nad_res/z_sat))
full_fov = np.rad2deg(np.arctan(0.5*swath/z_sat)) - pixel_fov

focal_length = ((nx-1)*0.5) / ((np.tan(np.radians(full_fov))))
focal_pos  = gc.Point(x=0., y=0., z=focal_length)

x_ = -(nx-1)*0.5 + np.arange(nx)
y_ = -(ny-1)*0.5 + np.arange(ny)
x, y = np.meshgrid(x_, y_)
x = x.flatten()
y = y.flatten()
z = np.zeros_like(x)
id_pixels = gc.Point(x,y,z) # id = (0, 0, 0) corresponds to the pixel at the center of the camera
dir_pixels = gc.normalize(id_pixels - focal_pos)

# %%
# directions of first y pixels (x, y, z components)
dir_pixels.to_numpy().reshape(ny,nx,3)[0,:]

# %% [markdown]
# ### Select only pixels that view a specific box zone

# %%
# box of size 50 km2 in x and y, and 10 km in z
box = gc.BBox(p1=gc.Point(-25., -25., 0.), p2=gc.Point(25.,25.,10.))

# satellite position/pixel positions. We duplicate to get same size as the number of pixels
sat_pos = gc.Point(np.zeros_like(dir_pixels.x), np.zeros_like(dir_pixels.x),
                   np.full_like(dir_pixels.x, z_sat))
# create th rays
r_sat = gc.Ray(sat_pos, dir_pixels)
is_intersection = box.is_intersection(r_sat)


# %%
# pixels directions and x and y id of pixels
dir_pixels.to_numpy()[is_intersection][0:10], id_pixels.to_numpy()[is_intersection][0:10]
