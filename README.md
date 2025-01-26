# Geoclide

The python package for geometric calculations in the three-dimentional Euclidian space

Mustapha Moulana  
[HYGEOS](https://hygeos.com/en/)

-----------------------------------------

## Installation
The installation can be performed using the following command:
```shell
$ pip install git+https://github.com/hygeos/gatiab.git
```

## Testing
Run the command `pytest tests/ -s -v` to check that everything is running correctly.

## Examples
```python
import numpy as np
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

# Find the direction from ground to the satellite
dir_to_sat = gc.Vector(0., 0., 1.)  # start facing nadir
dir_to_sat = gc.get_rotateY_tf(theta)[dir_to_sat] # perform a rotation around y axis to consider vza
dir_to_sat = gc.get_rotateZ_tf(phi)[dir_to_sat]   # then a rotation around z axis to consider vaa
ray = gc.Ray(o=origin, d=dir_to_sat) # create the ray, starting from origin going in dir_to_sat direction

# Here without considering the sphericity of the earth
b1 = gc.BBox(p1=gc.Point(-math.inf, -math.inf, 0.), p2=gc.Point(math.inf, math.inf, sat_altitude))
ds_pp = gc.calc_intersection(b1, ray) # return an xarray dataset

# Here with the consideration of the sphericity of the earth
earth_radius = 6378. # the equatorial earth radius in kilometers
oTw = gc.get_translate_tf(gc.Vector(0., 0., -earth_radius))
sphere_sat_alti = gc.Sphere(radius=earth_radius+sat_altitude, oTw=oTw)  # apply oTw to move the sphere center to earth center
ds_sp = gc.calc_intersection(sphere_sat_alti, ray) # return an xarray dataset

print ("Satellete position (pp case) :", ds_pp['phit'].values)
print ("Satellete position (sp case) ", ds_sp['phit'].values)
```