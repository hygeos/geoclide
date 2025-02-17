#!/usr/bin/env python
# -*- coding: utf-8 -*-

import geoclide as gc
import numpy as np
import math


def test_triangle_intersection():
    p0 = gc.Point(-1., -1., 0.)
    p1 = gc.Point(1., -1., 0.)
    p2 = gc.Point(1., 1., 0.)

    tri = gc.Triangle(p0, p1, p2)
    assert (tri.area() == 2.)

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

    oTw = gc.get_translate_tf(gc.Vector(10., 0., 5.)) * gc.get_rotateY_tf(45.)
    tri = gc.Triangle(p0, p1, p2, oTw=oTw)
    assert (tri.area() == 0.5)

    ray = gc.Ray(o=gc.Point(0., 0., 4.8), d=gc.normalize(gc.Vector(1.,0.,0.)))

    thitv2, dgv2, is_intv2 = tri.intersect(ray, method='v2')
    thitv3, dgv3, is_intv3 = tri.intersect(ray, method='v3')

    assert (is_intv2), 'Problem with v2 intersection test'
    assert (is_intv3), 'Problem with v3 intersection test'
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

    thitv2_, is_intv2_ = tri.is_intersection_t(ray, method='v2')
    thitv3_, is_intv3_ = tri.is_intersection_t(ray, method='v3')
    assert(is_intv2_)
    assert(is_intv3_)
    assert (thitv2 == thitv2_)
    assert (thitv3 == thitv3_)

def test_triangle_mesh():
    # list of vertices
    vertices = np.array([ [-0.5, -0.5, 0.],                  # v0
                          [0.5, -0.5, 0.],                   # v1
                          [-0.5, 0.5, 0.],                   # v2
                          [0.5, 0.5, 0.]], dtype=np.float64) # v3

    faces = np.array([[0, 1, 2],                   # vertices index of T0
                      [2, 3, 1]], dtype=np.int32)  # vertices index of T1

    oTw = gc.get_translate_tf(gc.Vector(10., 0., 5.)) * gc.get_rotateY_tf(45.)
    tri_mesh = gc.TriangleMesh(vertices=vertices, faces=faces, oTw=oTw)
    assert (tri_mesh.area() == 1.)

    ray = gc.Ray(o=gc.Point(0., 0., 4.8), d=gc.normalize(gc.Vector(1.,0.,0.)))

    ds_v2 = gc.calc_intersection(tri_mesh, ray, method='v2')
    ds_v3 = gc.calc_intersection(tri_mesh, ray, method='v3')
    thitv2 = ds_v2['thit'].values
    nhitv2 = gc.Normal(ds_v2['nhit'].values)
    is_intv2 = bool(ds_v2['is_intersection'].values)
    thitv3 = ds_v3['thit'].values
    nhitv3 = gc.Normal(ds_v3['nhit'].values)
    is_intv3 = bool(ds_v3['is_intersection'].values)

    assert (is_intv2 is True), 'Problem with v2 intersection test'
    assert (is_intv3 is True), 'Problem with v3 intersection test'
    assert (np.isclose(10.2, thitv2, 0., 1e-14)), 'Problem with v2 intersection test'
    assert (np.isclose(10.2, thitv3, 0., 1e-14)), 'Problem with v3 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., nhitv2.x, 0., 1e-13)), \
        'Problem with v2 intersection test'
    assert (nhitv2.y == 0.), 'Problem with v2 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., nhitv2.z, 0., 1e-13)), \
        'Problem with v2 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., nhitv3.x, 0., 1e-13)), \
        'Problem with v3 intersection test'
    assert (nhitv3.y == 0.), 'Problem with v3 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., nhitv3.z, 0., 1e-13)), \
        'Problem with v3 intersection test'

    p3 = gc.Point(10.2, 0., 4.8)
    phitv2 = gc.Point(ds_v2['phit'].values)
    phitv3 = gc.Point(ds_v3['phit'].values)

    assert (np.isclose(phitv2.x, p3.x, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(phitv2.y, p3.y, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(phitv2.z, p3.z, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(phitv3.x, p3.x, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(phitv3.y, p3.y, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(phitv3.z, p3.z, 0., 1e-15)), 'Problem with v3 triangle intersection test'

    assert (tri_mesh.is_intersection(ray, method='v2')), 'Problem with v2 is_intersection test'
    assert (tri_mesh.is_intersection(ray, method='v3')), 'Problem with v3 is_intersection test'


def test_read_gcnc_trianglemesh():
    msh_read = gc.read_gcnc_trianglemesh('geoclide/tests/data/sphere_r1p5_resth18_resph36.gcnc')
    msh = gc.Sphere(1.5).to_trianglemesh(reso_theta=18, reso_phi=36)
    assert (np.all(np.isclose(msh_read.vertices, msh.vertices, 0., 1e-15)))
    assert (np.all(msh_read.faces == msh.faces))

    msh_read = gc.read_gcnc_trianglemesh('geoclide/tests/data/sphere_r1p5_resth18_resph18_phimax180.gcnc')
    msh = gc.Sphere(1.5, phi_max=180.).to_trianglemesh(reso_theta=18, reso_phi=18)
    assert (np.all(np.isclose(msh_read.vertices, msh.vertices, 0., 1e-15)))
    assert (np.all(msh_read.faces == msh.faces))

    msh_read = gc.read_gcnc_trianglemesh('geoclide/tests/data/oblate_rxy1p5_rz1p2_resth18_resph36.gcnc')
    msh = gc.Spheroid(1.5, 1.2).to_trianglemesh(reso_theta=18, reso_phi=36)
    assert (np.all(np.isclose(msh_read.vertices, msh.vertices, 0., 1e-15)))
    assert (np.all(msh_read.faces == msh.faces))

    msh_read = gc.read_gcnc_trianglemesh('geoclide/tests/data/prolate_rxy1p5_rz3_resth18_resph36.gcnc')
    msh = gc.Spheroid(1.5, 3.).to_trianglemesh(reso_theta=18, reso_phi=36)
    assert (np.all(np.isclose(msh_read.vertices, msh.vertices, 0., 1e-15)))
    assert (np.all(msh_read.faces == msh.faces))

