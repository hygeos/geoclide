#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
ROOTPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
import sys
sys.path.insert(0, ROOTPATH)
import geoclide as gc
import numpy as np


def test_triangle_intersection():
    p0 = gc.Point(-1., -1., 0.)
    p1 = gc.Point(1., -1., 0.)
    p2 = gc.Point(1., 1., 0.)

    tri = gc.Triangle(p0, p1, p2)
    ray = gc.Ray(o=gc.Point(0., 0., 1.), d=gc.normalize(gc.Vector(0.999,0.999,-1.)))

    thit0, dg0, is_int0 = tri.intersect_v2(ray)
    thit1, dg1, is_int1 = tri.intersect_v3(ray)

    assert (is_int0 is True), 'Problem with v2 intersection test'
    assert (is_int1 is True), 'Problem with v3 intersection test'
    assert (np.isclose(1.73089629960896, thit0, 0., 1e-14)), 'Problem with v2 intersection test'
    assert (np.isclose(1.73089629960896, thit1, 0., 1e-14)), 'Problem with v3 intersection test'
    assert (dg0.n == gc.Normal(0., 0., 1.)), 'Problem with v2 intersection test'
    assert (dg1.n == gc.Normal(0., 0., 1.)), 'Problem with v3 intersection test'

    p3 = gc.Point(0.999, 0.999, 0.)

    assert (np.isclose(dg0.p.x, p3.x, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dg0.p.y, p3.y, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dg0.p.z, p3.z, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dg1.p.x, p3.x, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(dg1.p.y, p3.y, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(dg1.p.z, p3.z, 0., 1e-15)), 'Problem with v3 triangle intersection test'

    # Bellow the ray cannot reach the triangle
    ray = gc.Ray(o=gc.Point(0., 0., 1.), d=gc.normalize(gc.Vector(0.999,0.999,-1.)), maxt=1.7)

    thit0, dg0, is_int0 = tri.intersect_v2(ray)
    thit1, dg1, is_int1 = tri.intersect_v3(ray)
    
    assert (is_int0 is False), 'Problem with v2 intersection test'
    assert (is_int1 is False), 'Problem with v3 intersection test'

