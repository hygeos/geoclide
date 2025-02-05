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


def test_spheroid():
    oblate = gc.Spheroid(radius_xy=3., radius_z=1.5)
    prolate = gc.Spheroid(radius_xy=1.5, radius_z=3.)
    r1 = gc.Ray(o=gc.Point(2.5, 0., 10.), d=(gc.Vector(0., 0., -1.)))
    r2 = gc.Ray(o=gc.Point(10., 0., 2.5), d=(gc.Vector(-1., 0., 0.)))
    
    ds = gc.calc_intersection(oblate, r1)
    p = ds['phit'].values
    assert (p[0] == 2.5)
    assert (p[1] == 0.)
    assert (np.isclose(p[2], 0.8291561975888655, 0., 1e-15))
    n = ds['nhit'].values
    assert (np.isclose(n[0], 0.6019292654288356, 0., 1e-15))
    assert (n[1] == 0.)
    assert (np.isclose(n[2], 0.7985494095046983, 0., 1e-15))

    ds = gc.calc_intersection(prolate, r2)
    p = ds['phit'].values
    assert (np.isclose(p[0], 0.8291561975888655, 0., 1e-15))
    assert (p[1] == 0.)
    assert (p[2] == 2.5)
    n = ds['nhit'].values
    assert (np.isclose(n[0], 0.7985494095046906, 0., 1e-15))
    assert (n[1] == 0.)
    assert (np.isclose(n[2], 0.6019292654288461, 0., 1e-15))


