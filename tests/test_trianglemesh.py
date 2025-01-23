#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
ROOTPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
import sys
sys.path.insert(0, ROOTPATH)
import geoclide as gc
import numpy as np
import math


def test_triangle_intersection():
    p0 = gc.Point(-1., -1., 0.)
    p1 = gc.Point(1., -1., 0.)
    p2 = gc.Point(1., 1., 0.)

    tri = gc.Triangle(p0, p1, p2)
    ray = gc.Ray(o=gc.Point(0., 0., 1.), d=gc.normalize(gc.Vector(0.999,0.999,-1.)))

    thitv2, dgv2, is_intv2 = tri.intersect_v2(ray)
    thitv3, dgv3, is_intv3 = tri.intersect_v3(ray)

    assert (is_intv2 is True), 'Problem with v2 intersection test'
    assert (is_intv3 is True), 'Problem with v3 intersection test'
    assert (np.isclose(1.73089629960896, thitv2, 0., 1e-14)), 'Problem with v2 intersection test'
    assert (np.isclose(1.73089629960896, thitv3, 0., 1e-14)), 'Problem with v3 intersection test'
    assert (dgv2.n == gc.Normal(0., 0., 1.)), 'Problem with v2 intersection test'
    assert (dgv3.n == gc.Normal(0., 0., 1.)), 'Problem with v3 intersection test'

    p3 = gc.Point(0.999, 0.999, 0.)

    assert (np.isclose(dgv2.p.x, p3.x, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dgv2.p.y, p3.y, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dgv2.p.z, p3.z, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dgv3.p.x, p3.x, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(dgv3.p.y, p3.y, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(dgv3.p.z, p3.z, 0., 1e-15)), 'Problem with v3 triangle intersection test'

    # Bellow the ray cannot reach the triangle
    ray = gc.Ray(o=gc.Point(0., 0., 1.), d=gc.normalize(gc.Vector(0.999,0.999,-1.)), maxt=1.7)

    thitv2, dgv2, is_intv2 = tri.intersect_v2(ray)
    thitv3, dgv3, is_intv3 = tri.intersect_v3(ray)
    
    assert (is_intv2 is False), 'Problem with v2 intersection test'
    assert (is_intv3 is False), 'Problem with v3 intersection test'

def test_triangle_transform():
    p0 = gc.Point(-0.5, 0.5, 0.)
    p1 = gc.Point(0.5, 0.5, 0.)
    p2 = gc.Point(0.5, -0.5, 0.)

    wTo = gc.get_translate_tf(gc.Vector(10., 0., 5.)) * gc.get_rotateY_tf(45.)
    oTw = wTo.inverse()

    tri = gc.Triangle(p0, p1, p2, wTo=wTo, oTw=oTw)
    ray = gc.Ray(o=gc.Point(0., 0., 4.8), d=gc.normalize(gc.Vector(1.,0.,0.)))

    thitv2, dgv2, is_intv2 = tri.intersect(ray, method='v2')
    thitv3, dgv3, is_intv3 = tri.intersect(ray, method='v3')

    assert (is_intv2 is True), 'Problem with v2 intersection test'
    assert (is_intv3 is True), 'Problem with v3 intersection test'
    assert (np.isclose(10.2, thitv2, 0., 1e-14)), 'Problem with v2 intersection test'
    assert (np.isclose(10.2, thitv3, 0., 1e-14)), 'Problem with v3 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., dgv2.n.x, 0., 1e-13)), \
        'Problem with v2 intersection test'
    assert (dgv2.n.y == 0.), 'Problem with v2 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., dgv2.n.z, 0., 1e-13)), \
        'Problem with v2 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., dgv3.n.x, 0., 1e-13)), \
        'Problem with v3 intersection test'
    assert (dgv3.n.y == 0.), 'Problem with v3 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., dgv3.n.z, 0., 1e-13)), \
        'Problem with v3 intersection test'

    p3 = gc.Point(10.2, 0., 4.8)

    assert (np.isclose(dgv2.p.x, p3.x, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dgv2.p.y, p3.y, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dgv2.p.z, p3.z, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dgv3.p.x, p3.x, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(dgv3.p.y, p3.y, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(dgv3.p.z, p3.z, 0., 1e-15)), 'Problem with v3 triangle intersection test'

    assert (tri.is_intersection(ray, method='v2')), 'Problem with v2 is_intersection test'
    assert (tri.is_intersection(ray, method='v3')), 'Problem with v3 is_intersection test'


def test_triangle_mesh():
    # list of vertices
    v = np.array([ [-0.5, -0.5, 0.],               # v0
                [0.5, -0.5, 0.],                   # v1
                [-0.5, 0.5, 0.],                   # v2
                [0.5, 0.5, 0.]], dtype=np.float64) # v3

    vi = np.array([0, 1, 2,                # vertices index of T0
                2, 3, 1], dtype=np.int32)  # vertices index of T1

    wTo = gc.get_translate_tf(gc.Vector(10., 0., 5.)) * gc.get_rotateY_tf(45.)
    oTw = wTo.inverse()

    tri_mesh = gc.TriangleMesh(vi=vi, v=v, wTo=wTo, oTw=oTw)
    ray = gc.Ray(o=gc.Point(0., 0., 4.8), d=gc.normalize(gc.Vector(1.,0.,0.)))

    thitv2, dgv2, is_intv2 = tri_mesh.intersect(ray, method='v2')
    thitv3, dgv3, is_intv3 = tri_mesh.intersect(ray, method='v3')

    assert (is_intv2 is True), 'Problem with v2 intersection test'
    assert (is_intv3 is True), 'Problem with v3 intersection test'
    assert (np.isclose(10.2, thitv2, 0., 1e-14)), 'Problem with v2 intersection test'
    assert (np.isclose(10.2, thitv3, 0., 1e-14)), 'Problem with v3 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., dgv2.n.x, 0., 1e-13)), \
        'Problem with v2 intersection test'
    assert (dgv2.n.y == 0.), 'Problem with v2 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., dgv2.n.z, 0., 1e-13)), \
        'Problem with v2 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., dgv3.n.x, 0., 1e-13)), \
        'Problem with v3 intersection test'
    assert (dgv3.n.y == 0.), 'Problem with v3 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., dgv3.n.z, 0., 1e-13)), \
        'Problem with v3 intersection test'

    p3 = gc.Point(10.2, 0., 4.8)

    assert (np.isclose(dgv2.p.x, p3.x, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dgv2.p.y, p3.y, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dgv2.p.z, p3.z, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dgv3.p.x, p3.x, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(dgv3.p.y, p3.y, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(dgv3.p.z, p3.z, 0., 1e-15)), 'Problem with v3 triangle intersection test'

    assert (tri_mesh.is_intersection(ray, method='v2')), 'Problem with v2 is_intersection test'
    assert (tri_mesh.is_intersection(ray, method='v3')), 'Problem with v3 is_intersection test'