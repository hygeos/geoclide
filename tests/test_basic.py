#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
import math
import os
ROOTPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
import sys
sys.path.insert(0, ROOTPATH)
import geoclide as gc

P1 = [np.array([0.,0.,0.]), np.array([1.,2.,3.])]
V1 = [np.array([0.,0.,1.]), np.array([1.,1.,1.])]
N1 = [np.array([0.,0.,1.]), np.array([1.,1.,1.])]


@pytest.mark.parametrize('v_arr', V1)
def test_vector(v_arr):
    v1 = gc.Vector(v_arr[0], v_arr[1], v_arr[2])
    v2 = gc.Vector(v_arr)
    assert (v1 == v2)
    assert (np.all(v1.to_numpy() == v_arr))
    assert (v1[0] == v_arr[0])
    assert (v1[1] == v_arr[1])
    assert (v1[2] == v_arr[2])
    assert (v1[0] == v1.x)
    assert (v1[1] == v1.y)
    assert (v1[2] == v1.z)


@pytest.mark.parametrize('p_arr', P1)
def test_point(p_arr):
    p1 = gc.Point(p_arr[0], p_arr[1], p_arr[2])
    p2 = gc.Point(p_arr)
    assert (p1 == p2)
    assert (np.all(p1.to_numpy() == p_arr))
    assert (p1[0] == p_arr[0])
    assert (p1[1] == p_arr[1])
    assert (p1[2] == p_arr[2])
    assert (p1[0] == p1.x)
    assert (p1[1] == p1.y)
    assert (p1[2] == p1.z)


@pytest.mark.parametrize('n_arr', N1)
def test_normal(n_arr):
    n1 = gc.Normal(n_arr[0], n_arr[1], n_arr[2])
    n2 = gc.Normal(n_arr)
    assert (n1 == n2)
    assert (np.all(n1.to_numpy() == n_arr))
    assert (n1[0] == n_arr[0])
    assert (n1[1] == n_arr[1])
    assert (n1[2] == n_arr[2])
    assert (n1[0] == n1.x)
    assert (n1[1] == n1.y)
    assert (n1[2] == n1.z)


@pytest.mark.parametrize('v_arr', V1)
def test_ope_vector(v_arr):
    v1 = gc.Vector(v_arr)
    v2 = gc.Vector(1., 3., 0.5)
    assert (v1+v2 == gc.Vector(v_arr+v2.to_numpy()))
    assert (v1-v2 == gc.Vector(v_arr-v2.to_numpy()))
    assert (v1*2 == gc.Vector(v_arr*2.))
    assert (v1/2 == gc.Vector(v_arr/2.))
    v3 = v1+v2
    assert (v3.length() == math.sqrt(v3[0]**2+v3[1]**2+v3[2]**2))


@pytest.mark.parametrize('p_arr', P1)
def test_ope_point(p_arr):
    p1 = gc.Point(p_arr)
    v1 = gc.Vector(1., 3., 0.5)
    p2 = gc.Point(1.,1.,1.)
    assert (p1+v1 == gc.Point(p_arr+v1.to_numpy()))
    assert (p1-p2 == gc.Vector(p_arr-p2.to_numpy()))
    assert (p1-v1 == gc.Point(p_arr-v1.to_numpy()))
    assert (p1*2 == gc.Point(p_arr*2.))
    assert (p1/2 == gc.Point(p_arr/2.))


@pytest.mark.parametrize('n_arr', N1)
def test_ope_normal(n_arr):
    n1 = gc.Normal(n_arr)
    n2 = gc.Normal(1., 3., 0.5)
    assert (n1+n2 == gc.Normal(n_arr+n2.to_numpy()))
    assert (n1-n2 == gc.Normal(n_arr-n2.to_numpy()))
    assert (n1*2 == gc.Normal(n_arr*2.))
    assert (n1/2 == gc.Normal(n_arr/2.))


@pytest.mark.parametrize('p_arr', P1)
@pytest.mark.parametrize('v_arr', V1)
def test_ray(p_arr, v_arr):
    p1 = gc.Point(p_arr)
    v1 = gc.Vector(v_arr)
    r1 = gc.Ray(p1,v1)
    assert (r1.o == p1)
    assert (r1.d == v1)
    assert (r1[10.] == gc.Point(p_arr + 10.*v_arr))
    assert (r1.mint == 0)
    assert (r1.maxt == float("inf"))
    r2 = gc.Ray(p1,v1,0.5, 20.)
    assert (r2.mint == 0.5)
    assert (r2.maxt == 20.)
    r3 = gc.Ray(r2)
    assert (r3.o == r2.o)
    assert (r3.d == r2.d)
    assert (r3.mint == r2.mint)
    assert (r3.maxt == r2.maxt)


def test_bbox():
    # test attributs
    p1 = gc.Point(1., 1., 1.)
    p2 = gc.Point(2., 2., 3.)
    b1 = gc.BBox(p1, p2)
    minx = min(p1.x, p2.x)
    miny = min(p1.y, p2.y)
    minz = min(p1.z, p2.z)
    maxx = max(p1.x, p2.x)
    maxy = max(p1.y, p2.y)
    maxz = max(p1.z, p2.z)
    pmin = gc.Point(minx, miny, minz)
    pmax = gc.Point(maxx, maxy, maxz)
    assert (b1.pmin == pmin)
    assert (b1.pmax == pmax)
    p = [pmin]
    p.append(gc.Point(pmax.x,pmin.y,pmin.z))
    p.append(gc.Point(pmax.x,pmax.y,pmin.z))
    p.append(gc.Point(pmin.x,pmax.y,pmin.z))
    p.append(gc.Point(pmin.x,pmin.y,pmax.z))
    p.append(gc.Point(pmax.x,pmin.y,pmax.z))
    p.append(pmax)
    p.append(gc.Point(pmin.x,pmax.y,pmax.z))
    assert (b1.p0 == p[0])
    assert (b1.p1 == p[1])
    assert (b1.p2 == p[2])
    assert (b1.p3 == p[3])
    assert (b1.p4 == p[4])
    assert (b1.p5 == p[5])
    assert (b1.p6 == p[6])
    assert (b1.p7 == p[7])
    for i in range (len(b1.vertices)):
        assert (b1.vertices[i] == p[i])

    # test union method
    p3 = gc.Point(3., 3., 3.)
    b2 = b1.union(p3)
    assert (b2.pmin.x == min(b1.pmin.x, p3.x))
    assert (b2.pmin.y == min(b1.pmin.y, p3.y))
    assert (b2.pmin.z == min(b1.pmin.z, p3.z))
    assert (b2.pmax.x == max(b1.pmax.x, p3.x))
    assert (b2.pmax.y == max(b1.pmax.y, p3.y))
    assert (b2.pmax.z == max(b1.pmax.z, p3.z))
    b3 = b1.union(b2)
    assert (b3.pmin.x == min(b1.pmin.x, b2.pmin.x))
    assert (b3.pmin.y == min(b1.pmin.y, b2.pmin.y))
    assert (b3.pmin.z == min(b1.pmin.z, b2.pmin.z))
    assert (b3.pmax.x == max(b1.pmax.x, b2.pmax.x))
    assert (b3.pmax.y == max(b1.pmax.y, b2.pmax.y))
    assert (b3.pmax.z == max(b1.pmax.z, b2.pmax.z))

    # test is_inside method
    pIn = gc.Point(1.5, 1.5, 1.5)
    pOut = gc.Point(0., 0., 0.)
    assert (b3.is_inside(pIn))
    assert (not b3.is_inside(pOut))

    # test common_vertices method annd function get_common_vertices
    bc1 = gc.BBox(gc.Point(0., 0., 0.), gc.Point(2.5, 2.5, 2.5))
    bc2 = gc.BBox(gc.Point(2.5, 0., 0.), gc.Point(5., 2.5, 2.5))
    # The vertices of bc1 in common with the bc2 vertices are p1, p2, p5 and p6: 
    assert (np.all(bc1.common_vertices(bc2) == np.array([False, True, True, False, False, True, True, False])))
    assert (np.all(gc.get_common_vertices(bc1, bc2) == np.array([False, True, True, False, False, True, True, False])))
    # The vertices of bc2 in common with the bc1 vertices are p0, p3, p4 and p7: 
    assert (np.all(bc2.common_vertices(bc1) == np.array([True, False, False, True, True, False, False, True])))
    assert (np.all(gc.get_common_vertices(bc2, bc1) == np.array([True, False, False, True, True, False, False, True])))

    # test common_face method
    b_ref = gc.BBox(gc.Point(0., 0., 0.), gc.Point(1., 1., 1.))
    b_f0 = gc.BBox(gc.Point(1., 0., 0.), gc.Point(2., 1., 1.))
    b_f1 = gc.BBox(gc.Point(-1., 0., 0.), gc.Point(0., 1., 1.))
    b_f2 = gc.BBox(gc.Point(0., 1., 0.), gc.Point(1., 2., 1.))
    b_f3 = gc.BBox(gc.Point(0., -1., 0.), gc.Point(1., 0., 1.))
    b_f4 = gc.BBox(gc.Point(0., 0., 1.), gc.Point(1., 1., 2.))
    b_f5 = gc.BBox(gc.Point(0., 0., -1.), gc.Point(1., 1., 0.))
    assert (b_ref.common_face(b_f0) == 0)
    assert (b_ref.common_face(b_f1) == 1)
    assert (b_ref.common_face(b_f2) == 2)
    assert (b_ref.common_face(b_f3) == 3)
    assert (b_ref.common_face(b_f4) == 4)
    assert (b_ref.common_face(b_f5) == 5)
    assert (gc.get_common_face(b_ref, b_f0) == 0)
    assert (gc.get_common_face(b_ref, b_f1) == 1)
    assert (gc.get_common_face(b_ref, b_f2) == 2)
    assert (gc.get_common_face(b_ref, b_f3) == 3)
    assert (gc.get_common_face(b_ref, b_f4) == 4)
    assert (gc.get_common_face(b_ref, b_f5) == 5)

    # test method intersect_p
    b_int = gc.BBox(gc.Point(0., 0., 0.), gc.Point(1., 1., 1.))
    r_int_1 = gc.Ray(gc.Point(-0.5, 0.5, 0.5), gc.Vector(1.,0.,0.))
    r_int_2 = gc.Ray(gc.Point(0.5, 0.5, 0.5), gc.Vector(1.,0.,1.))
    r_int_3 = gc.Ray(gc.Point(0.5, 0.5, 2.), gc.Vector(0.,0.,1.))
    # case where the ray origin is outside and where the ray intersect 2 times the BBox
    t0, t1, is_int = b_int.intersect_p(r_int_1)
    assert (is_int)
    assert (r_int_1[t0] == gc.Point(0.0, 0.5, 0.5))
    assert (r_int_1[t1] == gc.Point(1.0, 0.5, 0.5))
    # case where the ray origin is inside and where the ray intersect 1 times the BBox
    t0, t1, is_int = b_int.intersect_p(r_int_2)
    assert (is_int)
    assert (r_int_2[t1] == gc.Point(1.0, 0.5, 1.0))
    # case where the ray origin is outsie and where the ray does not intersect with the BBox
    t0, t1, is_int = b_int.intersect_p(r_int_3)
    assert (not is_int)