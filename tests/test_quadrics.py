#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
ROOTPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
import sys
sys.path.insert(0, ROOTPATH)
import geoclide as gc


def test_sphere():
    # Take README exemple
    vza = 45.
    vaa = 45.
    sat_altitude = 700.
    origin = gc.Point(0., 0., 0.)
    theta = vza
    phi = -vaa

    dir_to_sat = gc.Vector(0., 0., 1.)
    dir_to_sat = gc.get_rotateY_tf(theta)[dir_to_sat]
    dir_to_sat = gc.get_rotateZ_tf(phi)[dir_to_sat]
    ray = gc.Ray(o=origin, d=dir_to_sat)

    earth_radius = 6378. 
    oTw = gc.get_translate_tf(gc.Vector(0., 0., -earth_radius))
    sphere_sat_alti = gc.Sphere(radius=earth_radius+sat_altitude, oTw=oTw)  # apply oTw to move the sphere center to earth center
    ds_sp = gc.calc_intersection(sphere_sat_alti, ray)

    p = ds_sp['phit'].values

    assert (np.isclose(p[0], 472.61058011386376, 0., 1e-15))
    assert (np.isclose(p[1], -472.61058011386365, 0., 1e-15))
    assert (np.isclose(p[2], 668.3722921180424, 0., 1e-15))