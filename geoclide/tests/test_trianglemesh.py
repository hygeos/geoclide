#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
import geoclide as gc
import os
from geoclide.shapes import get_intersect_dataset

ROOTPATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))


def test_triangle_intersection():
    p0 = gc.Point(-1., -1., 0.)
    p1 = gc.Point(1., -1., 0.)
    p2 = gc.Point(1., 1., 0.)

    tri = gc.Triangle(p0, p1, p2)
    assert (tri.area() == 2.)

    ray = gc.Ray(o=gc.Point(0., 0., 1.), d=gc.normalize(gc.Vector(0.999,0.999,-1.)))

    ds_v2 = get_intersect_dataset(*tri.intersect_v2(ray))
    thitv3, dgv3, is_intv3 = tri.intersect_v3(ray)

    assert (ds_v2['is_intersection'].item()), 'Problem with v2 intersection test'
    assert (is_intv3), 'Problem with v3 intersection test'
    assert (np.isclose(1.73089629960896, ds_v2['thit'].item(), 0., 1e-14)), 'Problem with v2 intersection test'
    assert (np.isclose(1.73089629960896, thitv3, 0., 1e-14)), 'Problem with v3 intersection test'
    assert (gc.Normal(ds_v2['nhit'].values) == gc.Normal(0., 0., 1.)), 'Problem with v2 intersection test'
    assert (dgv3.n == gc.Normal(0., 0., 1.)), 'Problem with v3 intersection test'

    p3 = gc.Point(0.999, 0.999, 0.)

    assert (np.isclose(ds_v2['phit'].values[0], p3.x, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(ds_v2['phit'].values[1], p3.y, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(ds_v2['phit'].values[2], p3.z, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dgv3.p.x, p3.x, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(dgv3.p.y, p3.y, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(dgv3.p.z, p3.z, 0., 1e-15)), 'Problem with v3 triangle intersection test'

    # Bellow the ray cannot reach the triangle
    ray = gc.Ray(o=gc.Point(0., 0., 1.), d=gc.normalize(gc.Vector(0.999,0.999,-1.)), maxt=1.7)

    ds_v2 = get_intersect_dataset(*tri.intersect_v2(ray))
    thitv3, dgv3, is_intv3 = tri.intersect_v3(ray)
    
    assert (ds_v2['is_intersection'].item() is False), 'Problem with v2 intersection test'
    assert (is_intv3 is False), 'Problem with v3 intersection test'

def test_triangle_transform():
    p0 = gc.Point(-0.5, 0.5, 0.)
    p1 = gc.Point(0.5, 0.5, 0.)
    p2 = gc.Point(0.5, -0.5, 0.)

    oTw = gc.get_translate_tf(gc.Vector(10., 0., 5.)) * gc.get_rotateY_tf(45.)
    tri = gc.Triangle(p0, p1, p2, oTw=oTw)
    assert (tri.area() == 0.5)

    ray = gc.Ray(o=gc.Point(0., 0., 4.8), d=gc.normalize(gc.Vector(1.,0.,0.)))

    ds_v2 = get_intersect_dataset(*tri.intersect(ray, method='v2'))
    thitv3, dgv3, is_intv3 = tri.intersect(ray, method='v3')

    assert (ds_v2['is_intersection'].item()), 'Problem with v2 intersection test'
    assert (is_intv3), 'Problem with v3 intersection test'
    assert (np.isclose(10.2, ds_v2['thit'].item(), 0., 1e-14)), 'Problem with v2 intersection test'
    assert (np.isclose(10.2, thitv3, 0., 1e-14)), 'Problem with v3 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., ds_v2['nhit'].values[0], 0., 1e-13)), \
        'Problem with v2 intersection test'
    assert (ds_v2['nhit'].values[1] == 0.), 'Problem with v2 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., ds_v2['nhit'].values[2], 0., 1e-13)), \
        'Problem with v2 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., dgv3.n.x, 0., 1e-13)), \
        'Problem with v3 intersection test'
    assert (dgv3.n.y == 0.), 'Problem with v3 intersection test'
    assert (np.isclose(-math.sqrt(2.)/2., dgv3.n.z, 0., 1e-13)), \
        'Problem with v3 intersection test'

    p3 = gc.Point(10.2, 0., 4.8)

    assert (np.isclose(ds_v2['phit'].values[0], p3.x, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(ds_v2['phit'].values[1], p3.y, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(ds_v2['phit'].values[2], p3.z, 0., 1e-15)), 'Problem with v2 triangle intersection test'
    assert (np.isclose(dgv3.p.x, p3.x, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(dgv3.p.y, p3.y, 0., 1e-15)), 'Problem with v3 triangle intersection test'
    assert (np.isclose(dgv3.p.z, p3.z, 0., 1e-15)), 'Problem with v3 triangle intersection test'

    assert (tri.is_intersection(ray, method='v2')), 'Problem with v2 is_intersection test'
    assert (tri.is_intersection(ray, method='v3')), 'Problem with v3 is_intersection test'

    thitv2_, is_intv2_ = tri.is_intersection_t(ray, method='v2')
    thitv3_, is_intv3_ = tri.is_intersection_t(ray, method='v3')
    assert (is_intv2_)
    assert (is_intv3_)
    assert (ds_v2['thit'].item() == thitv2_)
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
    msh_read = gc.read_trianglemesh(ROOTPATH + '/geoclide/tests/data/sphere_r1p5_resth18_resph36.gcnc')
    msh = gc.Sphere(1.5).to_trianglemesh(reso_theta=18, reso_phi=36)
    assert (np.all(np.isclose(msh_read.vertices, msh.vertices, 0., 1e-15)))
    assert (np.all(msh_read.faces == msh.faces))

    msh_read = gc.read_trianglemesh(ROOTPATH + '/geoclide/tests/data/sphere_r1p5_resth18_resph18_phimax180.gcnc')
    msh = gc.Sphere(1.5, phi_max=180.).to_trianglemesh(reso_theta=18, reso_phi=18)
    assert (np.all(np.isclose(msh_read.vertices, msh.vertices, 0., 1e-15)))
    assert (np.all(msh_read.faces == msh.faces))

    msh_read = gc.read_trianglemesh(ROOTPATH + '/geoclide/tests/data/oblate_rxy1p5_rz1p2_resth18_resph36.gcnc')
    msh = gc.Spheroid(1.5, 1.2).to_trianglemesh(reso_theta=18, reso_phi=36)
    assert (np.all(np.isclose(msh_read.vertices, msh.vertices, 0., 1e-15)))
    assert (np.all(msh_read.faces == msh.faces))

    msh_read = gc.read_trianglemesh(ROOTPATH + '/geoclide/tests/data/prolate_rxy1p5_rz3_resth18_resph36.gcnc')
    msh = gc.Spheroid(1.5, 3.).to_trianglemesh(reso_theta=18, reso_phi=36)
    assert (np.all(np.isclose(msh_read.vertices, msh.vertices, 0., 1e-15)))
    assert (np.all(msh_read.faces == msh.faces))

    msh_read = gc.read_trianglemesh(ROOTPATH + '/geoclide/tests/data/disk_r1_reso36.gcnc')
    msh = gc.Disk(1.).to_trianglemesh(reso=36)
    assert (np.all(np.isclose(msh_read.vertices, msh.vertices, 0., 1e-15)))
    assert (np.all(msh_read.faces == msh.faces))

    msh_read = gc.read_trianglemesh(ROOTPATH + '/geoclide/tests/data/disk_r1_zh5_reso36.gcnc')
    msh = gc.Disk(1., z_height=5.).to_trianglemesh(reso=36)
    assert (np.all(np.isclose(msh_read.vertices, msh.vertices, 0., 1e-15)))
    assert (np.all(msh_read.faces == msh.faces))

    msh_read = gc.read_trianglemesh(ROOTPATH + '/geoclide/tests/data/annulus_r1_ir0p5_reso36.gcnc')
    msh = gc.Disk(1., inner_radius=0.5).to_trianglemesh(reso=36)
    assert (np.all(np.isclose(msh_read.vertices, msh.vertices, 0., 1e-15)))
    assert (np.all(msh_read.faces == msh.faces))

    msh_read = gc.read_trianglemesh(ROOTPATH + '/geoclide/tests/data/annulus_partial_r1_ir0p5_reso27_phimax270.gcnc')
    msh = gc.Disk(1., inner_radius=0.5, phi_max=270).to_trianglemesh(reso=27)
    assert (np.all(np.isclose(msh_read.vertices, msh.vertices, 0., 1e-15)))
    assert (np.all(msh_read.faces == msh.faces))


def test_trianglemesh_to_dataset():
    msh1 = gc.Sphere(1.5).to_trianglemesh(reso_theta=18, reso_phi=36)
    msh2 = msh1.to_dataset()
    assert (np.all(np.isclose(msh1.vertices, msh2['vertices'].values[0,:,:], 0., 1e-15)))
    assert (np.all(msh1.faces == msh2['faces'].values[0,:,:]))


def test_trianglemesh_fast_tests():
    msh = gc.Spheroid(radius_xy=1., radius_z=3).to_trianglemesh(reso_theta=45, reso_phi=90)
    r0 = gc.Ray(gc.Point(5., -0.2, 2.5), gc.Vector(-1., 0., 0.))
    ds_v2 = get_intersect_dataset(*msh.intersect(r0, 'v2', fast_test=False))
    ds_v2f = get_intersect_dataset(*msh.intersect(r0, 'v2', fast_test=True))
    t_v3f, dg_v3f, is_int_v3f = msh.intersect(r0, 'v3', fast_test=True)
    assert (ds_v2['is_intersection'].item() == ds_v2f['is_intersection'].item())
    assert (ds_v2['is_intersection'].item() == is_int_v3f)
    assert (np.isclose(ds_v2['thit'].item(), ds_v2f['thit'].item(), 0., 1e-15))
    assert (np.isclose(ds_v2['thit'].item(), t_v3f, 0., 1e-15))
    assert (np.isclose(ds_v2['phit'].values[0], ds_v2f['phit'].values[0], 0., 1e-15))
    assert (np.isclose(ds_v2['phit'].values[1], ds_v2f['phit'].values[1], 0., 1e-15))
    assert (np.isclose(ds_v2['phit'].values[2], ds_v2f['phit'].values[2], 0., 1e-15))
    assert (np.isclose(ds_v2['phit'].values[0], dg_v3f.p.x, 0., 1e-15))
    assert (np.isclose(ds_v2['phit'].values[1], dg_v3f.p.y, 0., 1e-15))
    assert (np.isclose(ds_v2['phit'].values[2], dg_v3f.p.z, 0., 1e-15))
    assert (np.isclose(ds_v2['dpdu'].values[0], ds_v2f['dpdu'].values[0], 0., 1e-15))
    assert (np.isclose(ds_v2['dpdu'].values[1], ds_v2f['dpdu'].values[1], 0., 1e-15))
    assert (np.isclose(ds_v2['dpdu'].values[2], ds_v2f['dpdu'].values[2], 0., 1e-15))
    assert (np.isclose(ds_v2['dpdu'].values[0], dg_v3f.dpdu.x, 0., 1e-15))
    assert (np.isclose(ds_v2['dpdu'].values[1], dg_v3f.dpdu.y, 0., 1e-15))
    assert (np.isclose(ds_v2['dpdu'].values[2], dg_v3f.dpdu.z, 0., 1e-15))
    assert (np.isclose(ds_v2['dpdv'].values[0], ds_v2f['dpdv'].values[0], 0., 1e-15))
    assert (np.isclose(ds_v2['dpdv'].values[1], ds_v2f['dpdv'].values[1], 0., 1e-15))
    assert (np.isclose(ds_v2['dpdv'].values[2], ds_v2f['dpdv'].values[2], 0., 1e-15))
    assert (np.isclose(ds_v2['dpdv'].values[0], dg_v3f.dpdv.x, 0., 1e-15))
    assert (np.isclose(ds_v2['dpdv'].values[1], dg_v3f.dpdv.y, 0., 1e-15))
    assert (np.isclose(ds_v2['dpdv'].values[2], dg_v3f.dpdv.z, 0., 1e-15))
    assert (np.isclose(ds_v2['u'].item(), ds_v2f['u'].item(), 0., 1e-15))
    assert (np.isclose(ds_v2['u'].item(), dg_v3f.u, 0., 1e-15))
    assert (np.isclose(ds_v2['v'].item(), ds_v2['v'].item(), 0., 1e-15))
    assert (np.isclose(ds_v2['v'].item(), dg_v3f.v, 0., 1e-15))


def test_triangle_2d_arr1():
    msh = gc.Sphere(1.).to_trianglemesh(reso_theta=5, reso_phi=5)
    x_, y_, z_ = np.meshgrid(np.linspace(-0.4, 0.4, 4, np.float64),
                            np.linspace(-0.4, 0.4, 4, np.float64),
                            2., indexing='ij')
    o_set_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
    nrays = o_set_arr.shape[0]
    nobj = msh.ntriangles
    d_set_arr = np.zeros_like(o_set_arr)
    d_set_arr[:,0] = 0.
    d_set_arr[:,1] = 0.
    d_set_arr[:,2] = -1.
    o_set = gc.Point(o_set_arr)
    d_set = gc.Vector(d_set_arr)
    
    p0 = gc.Point(msh.vertices[msh.faces[:,0],:])
    p1 = gc.Point(msh.vertices[msh.faces[:,1],:])
    p2 = gc.Point(msh.vertices[msh.faces[:,2],:])
    triangles = gc.Triangle(p0, p1, p2)
    r_set = gc.Ray(o_set, d_set)
    res = triangles.intersect_v2(r_set, diag_calc=False)
    ds = get_intersect_dataset(*res)

    is_int_2d = np.full((nobj, nrays), True, dtype=bool)
    t_2d = np.zeros((nobj, nrays), dtype=np.float64)
    p_2d = np.zeros((nobj, nrays,3), dtype=np.float64)
    n_2d = np.zeros_like(p_2d)
    dpdu_2d = np.zeros_like(p_2d)
    dpdv_2d = np.zeros_like(p_2d)
    u_2d = np.zeros_like(t_2d)
    v_2d = np.zeros_like(t_2d)
    list_rays = []
    for ir in range (0, nrays):
        list_rays.append(gc.Ray(gc.Point(o_set_arr[ir,:]), gc.Vector(d_set_arr[ir,:])))

    for itri in range(0, nobj):
        p0 = gc.Point(msh.vertices[msh.faces[itri,0],:])
        p1 = gc.Point(msh.vertices[msh.faces[itri,1],:])
        p2 = gc.Point(msh.vertices[msh.faces[itri,2],:])
        triangle = gc.Triangle(p0, p1, p2)
        for ir in range (0, nrays):
            res_sca = triangle.intersect_v2(list_rays[ir])
            if res_sca[2] is not None:
                ds_sca = get_intersect_dataset(*res_sca)
                is_int_2d[itri,ir] = ds_sca['is_intersection'].values
                t_2d[itri,ir] = ds_sca['thit'].values
                p_2d[itri,ir,:] = ds_sca['phit'].values
                n_2d[itri,ir,:] = ds_sca['nhit'].values
                dpdu_2d[itri,ir,:] = ds_sca['dpdu'].values
                dpdv_2d[itri,ir,:] = ds_sca['dpdv'].values
                u_2d[itri,ir] = ds_sca['u'].values
                v_2d[itri,ir] = ds_sca['v'].values
            else:
                is_int_2d[itri,ir] = res_sca[3]
                t_2d[itri,ir] = res_sca[2]
                p_2d[itri,ir,:] = None
                n_2d[itri,ir,:] = None
                dpdu_2d[itri,ir,:] = None
                dpdv_2d[itri,ir,:] = None
                u_2d[itri,ir] = None
                v_2d[itri,ir] = None

    assert (np.array_equal(ds['thit'].values, t_2d, equal_nan=True))
    assert (np.array_equal(ds['phit'].values, p_2d, equal_nan=True))
    assert (np.array_equal(ds['nhit'].values, n_2d, equal_nan=True))
    assert (np.array_equal(ds['u'].values, u_2d, equal_nan=True))
    assert (np.array_equal(ds['v'].values, v_2d, equal_nan=True))
    assert (np.array_equal(ds['dpdu'].values, dpdu_2d, equal_nan=True))
    assert (np.array_equal(ds['dpdv'].values, dpdv_2d, equal_nan=True))


def test_triangle_2d_arr2():
    msh = gc.Sphere(1.).to_trianglemesh(reso_theta=4, reso_phi=4)
    x_, y_, z_ = np.meshgrid(np.linspace(-0.4, 0.4, 5, np.float64),
                            np.linspace(-0.4, 0.4, 5, np.float64),
                            2., indexing='ij')
    o_set_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
    nrays = o_set_arr.shape[0]
    nobj = msh.ntriangles
    d_set_arr = np.zeros_like(o_set_arr)
    d_set_arr[:,0] = 0.
    d_set_arr[:,1] = 0.
    d_set_arr[:,2] = -1.
    o_set = gc.Point(o_set_arr)
    d_set = gc.Vector(d_set_arr)
    
    p0 = gc.Point(msh.vertices[msh.faces[:,0],:])
    p1 = gc.Point(msh.vertices[msh.faces[:,1],:])
    p2 = gc.Point(msh.vertices[msh.faces[:,2],:])
    triangles = gc.Triangle(p0, p1, p2)
    r_set = gc.Ray(o_set, d_set)
    res = triangles.intersect_v2(r_set, diag_calc=False)
    ds = get_intersect_dataset(*res)

    is_int_2d = np.full((nobj, nrays), True, dtype=bool)
    t_2d = np.zeros((nobj, nrays), dtype=np.float64)
    p_2d = np.zeros((nobj, nrays,3), dtype=np.float64)
    n_2d = np.zeros_like(p_2d)
    dpdu_2d = np.zeros_like(p_2d)
    dpdv_2d = np.zeros_like(p_2d)
    u_2d = np.zeros_like(t_2d)
    v_2d = np.zeros_like(t_2d)
    list_rays = []
    for ir in range (0, nrays):
        list_rays.append(gc.Ray(gc.Point(o_set_arr[ir,:]), gc.Vector(d_set_arr[ir,:])))

    for itri in range(0, nobj):
        p0 = gc.Point(msh.vertices[msh.faces[itri,0],:])
        p1 = gc.Point(msh.vertices[msh.faces[itri,1],:])
        p2 = gc.Point(msh.vertices[msh.faces[itri,2],:])
        triangle = gc.Triangle(p0, p1, p2)
        for ir in range (0, nrays):
            res_sca = triangle.intersect_v2(list_rays[ir])
            if res_sca[2] is not None:
                ds_sca = get_intersect_dataset(*res_sca)
                is_int_2d[itri,ir] = ds_sca['is_intersection'].values
                t_2d[itri,ir] = ds_sca['thit'].values
                p_2d[itri,ir,:] = ds_sca['phit'].values
                n_2d[itri,ir,:] = ds_sca['nhit'].values
                dpdu_2d[itri,ir,:] = ds_sca['dpdu'].values
                dpdv_2d[itri,ir,:] = ds_sca['dpdv'].values
                u_2d[itri,ir] = ds_sca['u'].values
                v_2d[itri,ir] = ds_sca['v'].values
            else:
                is_int_2d[itri,ir] = False
                t_2d[itri,ir] = None
                p_2d[itri,ir,:] = None
                n_2d[itri,ir,:] = None
                dpdu_2d[itri,ir,:] = None
                dpdv_2d[itri,ir,:] = None
                u_2d[itri,ir] = None
                v_2d[itri,ir] = None


    assert (np.array_equal(ds['thit'].values, t_2d, equal_nan=True))
    assert (np.array_equal(ds['phit'].values, p_2d, equal_nan=True))
    assert (np.array_equal(ds['nhit'].values, n_2d, equal_nan=True))
    assert (np.array_equal(ds['u'].values, u_2d, equal_nan=True))
    assert (np.array_equal(ds['v'].values, v_2d, equal_nan=True))
    assert (np.array_equal(ds['dpdu'].values, dpdu_2d, equal_nan=True))
    assert (np.array_equal(ds['dpdv'].values, dpdv_2d, equal_nan=True))


def test_triangle_1d_arr1():
    msh = gc.Sphere(1.).to_trianglemesh(reso_theta=5, reso_phi=5)
    o = gc.Point(-0.4, 0.4, 2.)
    d = gc.Vector(0., 0., -1.)
    r0 = gc.Ray(o, d)

    p0 = gc.Point(msh.vertices[msh.faces[:,0],:])
    p1 = gc.Point(msh.vertices[msh.faces[:,1],:])
    p2 = gc.Point(msh.vertices[msh.faces[:,2],:])
    triangles = gc.Triangle(p0, p1, p2)
    res = triangles.intersect_v2(r0, diag_calc=False)
    ds = get_intersect_dataset(*res)

    is_int_1d = np.full((msh.ntriangles), True, dtype=bool)
    t_1d = np.zeros((msh.ntriangles), dtype=np.float64)
    p_1d = np.zeros((msh.ntriangles,3), dtype=np.float64)
    n_1d = np.zeros_like(p_1d)
    dpdu_1d = np.zeros_like(p_1d)
    dpdv_1d = np.zeros_like(p_1d)
    u_1d = np.zeros_like(t_1d)
    v_1d = np.zeros_like(t_1d)

    for itri in range(0, msh.ntriangles):
        p0 = gc.Point(msh.vertices[msh.faces[itri,0],:])
        p1 = gc.Point(msh.vertices[msh.faces[itri,1],:])
        p2 = gc.Point(msh.vertices[msh.faces[itri,2],:])
        triangle = gc.Triangle(p0, p1, p2)
        res_sca = triangle.intersect_v2(r0)
        if res_sca[2] is not None:
            ds_sca = get_intersect_dataset(*res_sca)
            is_int_1d[itri] = ds_sca['is_intersection'].values
            t_1d[itri] = ds_sca['thit'].values
            p_1d[itri,:] = ds_sca['phit'].values
            n_1d[itri,:] = ds_sca['nhit'].values
            dpdu_1d[itri,:] = ds_sca['dpdu'].values
            dpdv_1d[itri,:] = ds_sca['dpdv'].values
            u_1d[itri] = ds_sca['u'].values
            v_1d[itri] = ds_sca['v'].values
        else:
            is_int_1d[itri] = False
            t_1d[itri] = None
            p_1d[itri,:] = None
            n_1d[itri,:] = None
            dpdu_1d[itri,:] = None
            dpdv_1d[itri,:] = None
            u_1d[itri] = None
            v_1d[itri] = None

    assert (np.array_equal(ds['thit'].values, t_1d, equal_nan=True))
    assert (np.array_equal(ds['phit'].values, p_1d, equal_nan=True))
    assert (np.array_equal(ds['nhit'].values, n_1d, equal_nan=True))
    assert (np.array_equal(ds['u'].values, u_1d, equal_nan=True))
    assert (np.array_equal(ds['v'].values, v_1d, equal_nan=True))
    assert (np.array_equal(ds['dpdu'].values, dpdu_1d, equal_nan=True))
    assert (np.array_equal(ds['dpdv'].values, dpdv_1d, equal_nan=True))

def test_triangle_1d_arr2():

    x_, y_, z_ = np.meshgrid(np.linspace(-0.4, 0.4, 5, np.float64),
                         np.linspace(-0.4, 0.4, 5, np.float64),
                         2., indexing='ij')
    o_set_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
    nrays = o_set_arr.shape[0]
    d_set_arr = np.zeros_like(o_set_arr)
    d_set_arr[:,0] = 0.
    d_set_arr[:,1] = 0.
    d_set_arr[:,2] = -1.
    o_set = gc.Point(o_set_arr)
    d_set = gc.Vector(d_set_arr)
    r_set = gc.Ray(o_set, d_set)

    p0 = gc.Point(np.array([0., 0., 1]))
    p1 = gc.Point(np.array([0.70710678, 0., 0.70710678]))
    p2 = gc.Point(np.array([0.21850801, 0.67249851, 0.70710678]))
    triangle = gc.Triangle(p0, p1, p2)

    res = triangle.intersect_v2(r_set, diag_calc=False)
    ds = get_intersect_dataset(*res)

    is_int_1d = np.full((nrays), True, dtype=bool)
    t_1d = np.zeros((nrays), dtype=np.float64)
    p_1d = np.zeros((nrays,3), dtype=np.float64)
    n_1d = np.zeros_like(p_1d)
    dpdu_1d = np.zeros_like(p_1d)
    dpdv_1d = np.zeros_like(p_1d)
    u_1d = np.zeros_like(t_1d)
    v_1d = np.zeros_like(t_1d)
    list_rays = []
    for ir in range (0, nrays):
        list_rays.append(gc.Ray(gc.Point(o_set_arr[ir,:]), gc.Vector(d_set_arr[ir,:])))

    for ir in range (0, nrays):
        res_sca = triangle.intersect_v2(list_rays[ir])
        if (res_sca[2] is not None):
            ds_sca = get_intersect_dataset(*res_sca)
            is_int_1d[ir] = ds_sca['is_intersection'].values
            t_1d[ir] = ds_sca['thit'].values
            p_1d[ir,:] = ds_sca['phit'].values
            n_1d[ir,:] = ds_sca['nhit'].values
            dpdu_1d[ir,:] = ds_sca['dpdu'].values
            dpdv_1d[ir,:] = ds_sca['dpdv'].values
            u_1d[ir] = ds_sca['u'].values
            v_1d[ir] = ds_sca['v'].values
        else:
            is_int_1d[ir] = False
            t_1d[ir] = None
            p_1d[ir,:] = None
            n_1d[ir,:] = None
            dpdu_1d[ir,:] = None
            dpdv_1d[ir,:] = None
            u_1d[ir] = None
            v_1d[ir] = None

    assert (np.array_equal(ds['thit'].values, t_1d, equal_nan=True))
    assert (np.array_equal(ds['phit'].values, p_1d, equal_nan=True))
    assert (np.array_equal(ds['nhit'].values, n_1d, equal_nan=True))
    assert (np.array_equal(ds['u'].values, u_1d, equal_nan=True))
    assert (np.array_equal(ds['v'].values, v_1d, equal_nan=True))
    assert (np.array_equal(ds['dpdu'].values, dpdu_1d, equal_nan=True))
    assert (np.array_equal(ds['dpdv'].values, dpdv_1d, equal_nan=True))


def test_triangle_1d_arr3():
    msh = gc.Sphere(1.).to_trianglemesh(reso_theta=4, reso_phi=4)
    x_, y_, z_ = np.meshgrid(np.linspace(0.0, 0.6, 4, np.float64),
                         np.linspace(0., 0.4, 4, np.float64),
                         2., indexing='ij')
    o_set_arr = np.vstack((x_.ravel(), y_.ravel(), z_.ravel())).T
    nrays = o_set_arr.shape[0]
    d_set_arr = np.zeros_like(o_set_arr)
    d_set_arr[:,0] = 0.
    d_set_arr[:,1] = 0.
    d_set_arr[:,2] = -1.
    o_set = gc.Point(o_set_arr)
    d_set = gc.Vector(d_set_arr)
    r_set = gc.Ray(o_set, d_set)

    p0 = gc.Point(msh.vertices[msh.faces[:,0],:])
    p1 = gc.Point(msh.vertices[msh.faces[:,1],:])
    p2 = gc.Point(msh.vertices[msh.faces[:,2],:])
    triangles = gc.Triangle(p0, p1, p2)
    res = triangles.intersect_v2(r_set, diag_calc=True)
    ds = get_intersect_dataset(*res)

    ndiag = nrays
    is_int_1d = np.full((ndiag), True, dtype=bool)
    t_1d = np.zeros((ndiag), dtype=np.float64)
    p_1d = np.zeros((ndiag,3), dtype=np.float64)
    n_1d = np.zeros_like(p_1d)
    dpdu_1d = np.zeros_like(p_1d)
    dpdv_1d = np.zeros_like(p_1d)
    u_1d = np.zeros_like(t_1d)
    v_1d = np.zeros_like(t_1d)
    list_rays = []
    for ir in range (0, nrays):
        list_rays.append(gc.Ray(gc.Point(o_set_arr[ir,:]), gc.Vector(d_set_arr[ir,:])))

    for idiag in range (0, ndiag):
        p0 = gc.Point(msh.vertices[msh.faces[idiag,0],:])
        p1 = gc.Point(msh.vertices[msh.faces[idiag,1],:])
        p2 = gc.Point(msh.vertices[msh.faces[idiag,2],:])
        triangle = gc.Triangle(p0, p1, p2)
        res_sca = triangle.intersect_v2(list_rays[idiag])
        if (res_sca[2] is not None):
            ds_sca = get_intersect_dataset(*res_sca)
            is_int_1d[idiag] = ds_sca['is_intersection'].values
            t_1d[idiag] = ds_sca['thit'].values
            p_1d[idiag,:] = ds_sca['phit'].values
            n_1d[idiag,:] = ds_sca['nhit'].values
            dpdu_1d[idiag,:] = ds_sca['dpdu'].values
            dpdv_1d[idiag,:] = ds_sca['dpdv'].values
            u_1d[idiag] = ds_sca['u'].values
            v_1d[idiag] = ds_sca['v'].values
        else:
            is_int_1d[idiag] = False
            t_1d[idiag] = None
            p_1d[idiag,:] = None
            n_1d[idiag,:] = None
            dpdu_1d[idiag,:] = None
            dpdv_1d[idiag,:] = None
            u_1d[idiag] = None
            v_1d[idiag] = None

    assert (np.array_equal(ds['thit'].values, t_1d, equal_nan=True))
    assert (np.array_equal(ds['phit'].values, p_1d, equal_nan=True))
    assert (np.array_equal(ds['nhit'].values, n_1d, equal_nan=True))
    assert (np.array_equal(ds['u'].values, u_1d, equal_nan=True))
    assert (np.array_equal(ds['v'].values, v_1d, equal_nan=True))
    assert (np.array_equal(ds['dpdu'].values, dpdu_1d, equal_nan=True))
    assert (np.array_equal(ds['dpdv'].values, dpdv_1d, equal_nan=True))