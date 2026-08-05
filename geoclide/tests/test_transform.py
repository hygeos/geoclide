#!/usr/bin/env python
"""Tests of the transform module."""

import math
from typing import cast

import numpy as np
import pytest

import geoclide as gc

V1 = [np.array([5.0, 0.0, 0.0]), np.array([5.0, 3.0, 1.0])]
V2 = [np.array([5.0, 1.0, 1.0]), np.array([5.0, 3.0, 2.0])]
ANGLES = [30.0, 45.0, 60.0]


@pytest.mark.parametrize("v_arr", V1)
def test_get_translate_tf(v_arr):
    t = gc.get_translate_tf(gc.Vector(v_arr))
    m = np.identity(4)
    m[0, -1] = v_arr[0]
    m[1, -1] = v_arr[1]
    m[2, -1] = v_arr[2]
    assert np.all(t.m == m)
    m_inv = np.identity(4)
    m_inv[0, -1] = -v_arr[0]
    m_inv[1, -1] = -v_arr[1]
    m_inv[2, -1] = -v_arr[2]
    assert np.all(t.m_inv == m_inv)


@pytest.mark.parametrize("v_arr", V2)
def test_get_scale_tf(v_arr):
    t = gc.get_scale_tf(gc.Vector(v_arr))
    m = np.identity(4)
    m[0, 0] = v_arr[0]
    m[1, 1] = v_arr[1]
    m[2, 2] = v_arr[2]
    assert np.all(t.m == m)
    m_inv = np.identity(4)
    m_inv[0, 0] = 1 * (1 / v_arr[0])
    m_inv[1, 1] = 1 * (1 / v_arr[1])
    m_inv[2, 2] = 1 * (1 / v_arr[2])
    assert np.all(t.m_inv == m_inv)


@pytest.mark.parametrize("angle", ANGLES)
def test_get_rotate_x_tf(angle):
    t = gc.get_rotate_x_tf(angle)
    sin_t = math.sin(angle * (math.pi / 180.0))
    cos_t = math.cos(angle * (math.pi / 180.0))
    m = np.identity(4)
    m[1, 1] = cos_t
    m[1, 2] = -1.0 * sin_t
    m[2, 1] = sin_t
    m[2, 2] = cos_t
    assert np.all(t.m == m)
    m_inv = np.transpose(m)
    assert np.all(t.m_inv == m_inv)


@pytest.mark.parametrize("angle", ANGLES)
def test_get_rotate_y_tf(angle):
    t = gc.get_rotate_y_tf(angle)
    sin_t = math.sin(angle * (math.pi / 180.0))
    cos_t = math.cos(angle * (math.pi / 180.0))
    m = np.identity(4)
    m[0, 0] = cos_t
    m[2, 0] = -1.0 * sin_t
    m[0, 2] = sin_t
    m[2, 2] = cos_t
    assert np.all(t.m == m)
    m_inv = np.transpose(m)
    assert np.all(t.m_inv == m_inv)


@pytest.mark.parametrize("angle", ANGLES)
def test_get_rotate_z_tf(angle):
    t = gc.get_rotate_z_tf(angle)
    sin_t = math.sin(angle * (math.pi / 180.0))
    cos_t = math.cos(angle * (math.pi / 180.0))
    m = np.identity(4)
    m[0, 0] = cos_t
    m[0, 1] = -1.0 * sin_t
    m[1, 0] = sin_t
    m[1, 1] = cos_t
    assert np.all(t.m == m)
    m_inv = np.transpose(m)
    assert np.all(t.m_inv == m_inv)


@pytest.mark.parametrize("angle", ANGLES)
def test_get_rotate_tf_1(angle):
    v_x = gc.Vector(1.0, 0.0, 0.0)
    tx = gc.get_rotate_x_tf(angle)
    t = gc.get_rotate_tf(angle, v_x)
    assert np.all(tx.m == t.m)
    assert np.all(tx.m_inv == t.m_inv)

    v_y = gc.Vector(0.0, 1.0, 0.0)
    ty = gc.get_rotate_y_tf(angle)
    t = gc.get_rotate_tf(angle, v_y)
    assert np.all(ty.m == t.m)
    assert np.all(ty.m_inv == t.m_inv)

    v_z = gc.Vector(0.0, 0.0, 1.0)
    tz = gc.get_rotate_z_tf(angle)
    t = gc.get_rotate_tf(angle, v_z)
    assert np.all(tz.m == t.m)
    assert np.all(tz.m_inv == t.m_inv)


@pytest.mark.parametrize("angle", ANGLES)
@pytest.mark.parametrize("v_arr", V1)
def test_get_rotate_tf_2(angle, v_arr):
    v = gc.normalize(
        gc.Vector(v_arr)
    )  # important to normalize for Rodrigues formula
    t = gc.get_rotate_tf(angle, v)
    # Use independant method -> Rodrigues rotation formula
    m = np.identity(4)
    cos_t = math.cos(angle * (math.pi / 180.0))
    sin_t = math.sqrt(1.0 - cos_t * cos_t)
    mat_a = np.identity(3)
    mat_b = np.zeros((3, 3))
    mat_b[0, 1] = -v.z
    mat_b[0, 2] = v.y
    mat_b[1, 0] = v.z
    mat_b[1, 2] = -v.x
    mat_b[2, 0] = -v.y
    mat_b[2, 1] = v.x
    mat_c = mat_b.dot(mat_b)
    m[0:3, 0:3] = mat_a + mat_b * sin_t + mat_c * (1 - cos_t)
    assert np.allclose(t.m, m, rtol=0.0, atol=1e-15)
    m_inv = np.transpose(m)
    assert np.allclose(t.m_inv, m_inv, rtol=0.0, atol=1e-15)


def test_transform():
    t1 = gc.Transform()
    assert np.all(t1.m == np.identity(4))
    assert np.all(t1.m_inv == np.identity(4))

    t2 = gc.get_translate_tf(gc.Vector(5.0, 5.0, 5.0))
    p1 = gc.Point(0.0, 0.0, 0.0)
    v1 = gc.normalize(gc.Vector(1.0, 1.0, 1.0))
    n1 = gc.Normal(0.0, 0.0, 1.0)
    b1 = gc.BBox(p1, p1 + v1)
    b1_bis = gc.BBox()
    b1_bis = b1_bis.union(t2(b1.p0))
    b1_bis = b1_bis.union(t2(b1.p1))
    b1_bis = b1_bis.union(t2(b1.p2))
    b1_bis = b1_bis.union(t2(b1.p3))
    b1_bis = b1_bis.union(t2(b1.p4))
    b1_bis = b1_bis.union(t2(b1.p5))
    b1_bis = b1_bis.union(t2(b1.p6))
    b1_bis = b1_bis.union(t2(b1.p7))
    assert t2(p1) == gc.Point(5.0, 5.0, 5.0)
    assert t2(v1) == v1
    assert t2(n1) == n1
    assert t2(b1).pmin == b1_bis.pmin
    assert t2(b1).pmax == b1_bis.pmax


def test_transform_wrong_matrix():
    """Only the (4,4) and (nt,4,4) shapes must be accepted."""
    wrong = [
        np.identity(3),
        np.zeros((4, 3)),
        np.zeros((2, 3, 4)),
        np.zeros((2, 4, 3)),
        np.zeros(4),
        np.zeros((2, 2, 4, 4)),
    ]
    for m in wrong:
        with pytest.raises(ValueError):
            gc.Transform(m)
        with pytest.raises(ValueError):
            gc.Transform(m_inv=m)
        with pytest.raises(ValueError):
            gc.Transform(np.identity(4), m)

    # the accepted shapes must not raise
    gc.Transform(np.identity(4))
    gc.Transform(np.zeros((3, 4, 4)) + np.identity(4))


def test_transform_1d_arr():
    tf = np.zeros((2, 4, 4), dtype=np.float64)
    tf[0, :, :] = np.identity(4)
    tf[1, :, :] = gc.get_translate_tf(gc.Vector(1.0, 2.0, 3.0)).m
    multi_tf = gc.Transform(tf)
    tf1 = gc.Transform(tf[0, :, :])
    tf2 = gc.Transform(tf[1, :, :])
    p_arr = np.array([0.0, 0.0, 0.0])
    v_arr = np.array([0.0, 0.0, 1.0])
    p1 = gc.Point(p_arr)
    v1 = gc.Vector(v_arr)
    r_set = gc.Ray(p1, v1)

    r_set_mtf = cast(np.ndarray, multi_tf(r_set))
    r_set_tf1 = tf1(r_set)
    r_set_tf2 = tf2(r_set)

    assert r_set_mtf[0].o == r_set_tf1.o
    assert r_set_mtf[1].o == r_set_tf2.o

    assert r_set_mtf[0].d == r_set_tf1.d
    assert r_set_mtf[1].d == r_set_tf2.d


def test_transform_2d_arr():
    tf = np.zeros((2, 4, 4), dtype=np.float64)
    tf[0, :, :] = np.identity(4)
    tf[1, :, :] = gc.get_translate_tf(gc.Vector(1.0, 2.0, 3.0)).m
    multi_tf = gc.Transform(tf)
    tf1 = gc.Transform(tf[0, :, :])
    tf2 = gc.Transform(tf[1, :, :])
    p_arr = np.array([[0.0, 0.0, 0.0], [2.0, 2.0, 2.0]])
    v_arr = np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
    p_set = gc.Point(p_arr)
    v_set = gc.Vector(v_arr)
    r_set = gc.Ray(p_set, v_set)

    r_set_mtf = cast(np.ndarray, multi_tf(r_set))
    r_set_tf1 = tf1(r_set)
    r_set_tf2 = tf2(r_set)

    assert np.all(r_set_mtf[0].o == r_set_tf1.o)
    assert np.all(r_set_mtf[1].o == r_set_tf2.o)

    assert np.all(r_set_mtf[0].d == r_set_tf1.d)
    assert np.all(r_set_mtf[1].d == r_set_tf2.d)


def test_transform_diag():
    tf = np.zeros((2, 4, 4), dtype=np.float64)
    tf[0, :, :] = np.identity(4)
    tf[1, :, :] = gc.get_translate_tf(gc.Vector(1.0, 2.0, 3.0)).m
    multi_tf = gc.Transform(tf)
    tf1 = gc.Transform(tf[0, :, :])
    tf2 = gc.Transform(tf[1, :, :])
    p_arr = np.array([[0.0, 0.0, 0.0], [2.0, 2.0, 2.0]])
    v_arr = np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
    p_set = gc.Point(p_arr)
    v_set = gc.Vector(v_arr)
    r_set = gc.Ray(p_set, v_set)

    p1 = gc.Point(p_arr[0, :])
    v1 = gc.Vector(v_arr[0, :])
    r1 = gc.Ray(p1, v1)
    p2 = gc.Point(p_arr[1, :])
    v2 = gc.Vector(v_arr[1, :])
    r2 = gc.Ray(p2, v2)

    r_set_mtf = cast(np.ndarray, multi_tf(r_set, diag_calc=True))
    r1_tf1 = tf1(r1)
    r2_tf2 = tf2(r2)

    assert r_set_mtf[0].o == r1_tf1.o
    assert r_set_mtf[1].o == r2_tf2.o

    assert r_set_mtf[0].d == r1_tf1.d
    assert r_set_mtf[1].d == r2_tf2.d


def test_get_translate_tf_arr():
    v_arr = np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
    v_set = gc.Vector(v_arr)
    multi_tf = gc.get_translate_tf(v_set)
    mtf_m = multi_tf.m
    mtf_m_inv = multi_tf.m_inv

    mf_m = np.zeros_like(mtf_m)
    mf_m_inv = np.zeros_like(mtf_m_inv)
    for i in range(0, 2):
        tfi = gc.get_translate_tf(gc.Vector(v_arr[i, :]))
        mf_m[i, :, :] = tfi.m
        mf_m_inv[i, :, :] = tfi.m_inv
    assert np.all(mtf_m == mf_m)
    assert np.all(mtf_m_inv == mf_m_inv)


def test_get_scale_tf_arr():
    v_arr = np.array([[1.0, 2.0, 3.0], [4.0, 0.5, 7.0]])
    v_set = gc.Vector(v_arr)
    multi_tf = gc.get_scale_tf(v_set)
    mtf_m = multi_tf.m
    mtf_m_inv = multi_tf.m_inv

    mf_m = np.zeros_like(mtf_m)
    mf_m_inv = np.zeros_like(mtf_m_inv)
    for i in range(0, 2):
        tfi = gc.get_scale_tf(gc.Vector(v_arr[i, :]))
        mf_m[i, :, :] = tfi.m
        mf_m_inv[i, :, :] = tfi.m_inv
    assert np.all(mtf_m == mf_m)
    assert np.all(mtf_m_inv == mf_m_inv)


def test_get_rotate_x_tf_arr():
    angles = np.array([45, 78, 115])
    multi_tf = gc.get_rotate_x_tf(angles)
    mtf_m = multi_tf.m
    mtf_m_inv = multi_tf.m_inv

    mf_m = np.zeros_like(mtf_m)
    mf_m_inv = np.zeros_like(mtf_m)
    for i in range(0, 3):
        tfi = gc.get_rotate_x_tf(angles[i])
        mf_m[i, :, :] = tfi.m
        mf_m_inv[i, :, :] = tfi.m_inv
    assert np.allclose(mtf_m, mf_m, 0.0, 1e-15)
    assert np.allclose(mtf_m_inv, mf_m_inv, 0.0, 1e-15)


def test_get_rotate_y_tf_arr():
    angles = np.array([45, 78, 115])
    multi_tf = gc.get_rotate_y_tf(angles)
    mtf_m = multi_tf.m
    mtf_m_inv = multi_tf.m_inv

    mf_m = np.zeros_like(mtf_m)
    mf_m_inv = np.zeros_like(mtf_m)
    for i in range(0, 3):
        tfi = gc.get_rotate_y_tf(angles[i])
        mf_m[i, :, :] = tfi.m
        mf_m_inv[i, :, :] = tfi.m_inv
    assert np.allclose(mtf_m, mf_m, 0.0, 1e-15)
    assert np.allclose(mtf_m_inv, mf_m_inv, 0.0, 1e-15)


def test_get_rotate_z_tf_arr():
    angles = np.array([45, 78, 115])
    multi_tf = gc.get_rotate_z_tf(angles)
    mtf_m = multi_tf.m
    mtf_m_inv = multi_tf.m_inv

    mf_m = np.zeros_like(mtf_m)
    mf_m_inv = np.zeros_like(mtf_m)
    for i in range(0, 3):
        tfi = gc.get_rotate_z_tf(angles[i])
        mf_m[i, :, :] = tfi.m
        mf_m_inv[i, :, :] = tfi.m_inv
    assert np.allclose(mtf_m, mf_m, 0.0, 1e-15)
    assert np.allclose(mtf_m_inv, mf_m_inv, 0.0, 1e-15)


def test_get_rotate_tf_arr1():
    v_arr = np.array([1.0, 2.0, 3.0])
    v = gc.Vector(v_arr)
    angles = np.array([45, 78])
    multi_tf = gc.get_rotate_tf(angles, v, diag_calc=True)
    mtf_m = multi_tf.m
    mtf_m_inv = multi_tf.m_inv

    mf_m = np.zeros_like(mtf_m)
    mf_m_inv = np.zeros_like(mtf_m)
    for i in range(0, 2):
        tfi = gc.get_rotate_tf(angles[i], v)
        mf_m[i, :, :] = tfi.m
        mf_m_inv[i, :, :] = tfi.m_inv
    assert np.allclose(mtf_m, mf_m, 0.0, 1e-15)
    assert np.allclose(mtf_m_inv, mf_m_inv, 0.0, 1e-15)


def test_get_rotate_tf_arr2():
    v_arr = np.array([[1.0, 2.0, 3.0], [4.0, 0.5, 7.0]])
    v_set = gc.Vector(v_arr)
    angle = 45
    multi_tf = gc.get_rotate_tf(angle, v_set, diag_calc=True)
    mtf_m = multi_tf.m
    mtf_m_inv = multi_tf.m_inv

    mf_m = np.zeros_like(mtf_m)
    mf_m_inv = np.zeros_like(mtf_m)
    for i in range(0, 2):
        tfi = gc.get_rotate_tf(angle, gc.Vector(v_arr[i, :]))
        mf_m[i, :, :] = tfi.m
        mf_m_inv[i, :, :] = tfi.m_inv
    assert np.allclose(mtf_m, mf_m, 0.0, 1e-15)
    assert np.allclose(mtf_m_inv, mf_m_inv, 0.0, 1e-15)


def test_get_rotate_tf_diag():
    v_arr = np.array([[1.0, 2.0, 3.0], [4.0, 0.5, 7.0]])
    v_set = gc.Vector(v_arr)
    angles = np.array([45, 78])
    multi_tf = gc.get_rotate_tf(angles, v_set, diag_calc=True)
    mtf_m = multi_tf.m
    mtf_m_inv = multi_tf.m_inv

    mf_m = np.zeros_like(mtf_m)
    mf_m_inv = np.zeros_like(mtf_m)
    for i in range(0, 2):
        tfi = gc.get_rotate_tf(angles[i], gc.Vector(v_arr[i, :]))
        mf_m[i, :, :] = tfi.m
        mf_m_inv[i, :, :] = tfi.m_inv
    assert np.allclose(mtf_m, mf_m, 0.0, 1e-15)
    assert np.allclose(mtf_m_inv, mf_m_inv, 0.0, 1e-15)


def test_flatten_tf_1d_arr1():
    angles = np.linspace(0.0, 180.0, 5)
    tr_arr = np.zeros((len(angles), 3), dtype=np.float64)
    tr_arr[:, 0] = 1
    tr_arr[:, 1] = 2
    tr_arr[:, 2] = 3
    multi_tf = gc.get_translate_tf(gc.Vector(tr_arr)) * gc.get_rotate_y_tf(
        angles
    )
    v_arr = np.array([1.0, 2.0, 3.0])
    p_arr = np.array([1.0, 1.0, 1.0])
    v = gc.Vector(v_arr)
    p = gc.Point(p_arr)
    n = gc.Normal(v_arr)

    v_arr_f = multi_tf(v, flatten=True).to_numpy()
    p_arr_f = multi_tf(p, flatten=True).to_numpy()
    n_arr_f = multi_tf(n, flatten=True).to_numpy()

    v_bis = cast(np.ndarray, multi_tf(v))
    p_bis = cast(np.ndarray, multi_tf(p))
    n_bis = cast(np.ndarray, multi_tf(n))
    nv = len(v_bis)
    v_arr_nf = np.zeros((nv, 3), dtype=np.float64)
    p_arr_nf = np.zeros_like(v_arr_nf)
    n_arr_nf = np.zeros_like(v_arr_nf)
    for i in range(0, nv):
        v_arr_nf[i, :] = v_bis[i].to_numpy()
        p_arr_nf[i, :] = p_bis[i].to_numpy()
        n_arr_nf[i, :] = n_bis[i].to_numpy()

    assert np.allclose(v_arr_nf, v_arr_f, 0.0, 1e-15)
    assert np.allclose(p_arr_nf, p_arr_f, 0.0, 1e-15)
    assert np.allclose(n_arr_nf, n_arr_f, 0.0, 1e-15)


def test_flatten_tf_1d_arr2():
    angle = 45.0
    tf = gc.get_translate_tf(gc.Vector(1, 2, 3)) * gc.get_rotate_y_tf(angle)
    v_arr = np.array([[1.0, 2.0, 3.0], [0.5, 5.0, 1.0]])
    p_arr = np.array([[1.0, 1.0, 1.0], [5.0, 1.0, 1.0]])
    v = gc.Vector(v_arr)
    p = gc.Point(p_arr)
    n = gc.Normal(v_arr)

    v_arr_f = tf(v, flatten=True).to_numpy()
    p_arr_f = tf(p, flatten=True).to_numpy()
    n_arr_f = tf(n, flatten=True).to_numpy()

    v_arr_nf = tf(v).to_numpy()
    p_arr_nf = tf(p).to_numpy()
    n_arr_nf = tf(n).to_numpy()

    assert np.allclose(v_arr_nf, v_arr_f, 0.0, 1e-15)
    assert np.allclose(p_arr_nf, p_arr_f, 0.0, 1e-15)
    assert np.allclose(n_arr_nf, n_arr_f, 0.0, 1e-15)


def test_flatten_tf_2d_arr():
    angles = np.linspace(0.0, 180.0, 5)
    nang = len(angles)
    tr_arr = np.zeros((nang, 3), dtype=np.float64)
    tr_arr[:, 0] = 1
    tr_arr[:, 1] = 2
    tr_arr[:, 2] = 3
    multi_tf = gc.get_translate_tf(gc.Vector(tr_arr)) * gc.get_rotate_y_tf(
        angles
    )
    v_arr = np.array([[1.0, 2.0, 3.0], [0.5, 5.0, 1.0]])
    p_arr = np.array([[1.0, 1.0, 1.0], [5.0, 1.0, 1.0]])
    v = gc.Vector(v_arr)
    p = gc.Point(p_arr)
    n = gc.Normal(v_arr)

    v_arr_f = multi_tf(v, flatten=True).to_numpy()
    p_arr_f = multi_tf(p, flatten=True).to_numpy()
    n_arr_f = multi_tf(n, flatten=True).to_numpy()

    v_bis = cast(np.ndarray, multi_tf(v))
    p_bis = cast(np.ndarray, multi_tf(p))
    n_bis = cast(np.ndarray, multi_tf(n))
    nv = len(cast(np.ndarray, v.x))
    v_arr_nf = np.zeros((nv * nang, 3), dtype=np.float64)
    p_arr_nf = np.zeros_like(v_arr_nf)
    n_arr_nf = np.zeros_like(v_arr_nf)
    for i in range(0, nang):
        ind0 = i * nv
        ind1 = ind0 + nv
        v_arr_nf[ind0:ind1, :] = v_bis[i].to_numpy()
        p_arr_nf[ind0:ind1, :] = p_bis[i].to_numpy()
        n_arr_nf[ind0:ind1, :] = n_bis[i].to_numpy()

    assert np.allclose(v_arr_nf, v_arr_f, 0.0, 1e-15)
    assert np.allclose(p_arr_nf, p_arr_f, 0.0, 1e-15)
    assert np.allclose(n_arr_nf, n_arr_f, 0.0, 1e-15)


def test_flatten_tf_diag():
    angles = np.array([45.0, 165.0])
    nang = len(angles)
    tr_arr = np.zeros((nang, 3), dtype=np.float64)
    tr_arr[:, 0] = 1
    tr_arr[:, 1] = 2
    tr_arr[:, 2] = 3
    multi_tf = gc.get_translate_tf(gc.Vector(tr_arr)) * gc.get_rotate_y_tf(
        angles
    )
    v_arr = np.array([[1.0, 2.0, 3.0], [0.5, 5.0, 1.0]])
    p_arr = np.array([[1.0, 1.0, 1.0], [5.0, 1.0, 1.0]])
    v = gc.Vector(v_arr)
    p = gc.Point(p_arr)
    n = gc.Normal(v_arr)

    v_arr_f = multi_tf(v, diag_calc=True, flatten=True).to_numpy()
    p_arr_f = multi_tf(p, diag_calc=True, flatten=True).to_numpy()
    n_arr_f = multi_tf(n, diag_calc=True, flatten=True).to_numpy()

    v_bis = cast(np.ndarray, multi_tf(v, diag_calc=True))
    p_bis = cast(np.ndarray, multi_tf(p, diag_calc=True))
    n_bis = cast(np.ndarray, multi_tf(n, diag_calc=True))
    nv = len(cast(np.ndarray, v.x))
    v_arr_nf = np.zeros((nv, 3), dtype=np.float64)
    p_arr_nf = np.zeros_like(v_arr_nf)
    n_arr_nf = np.zeros_like(v_arr_nf)
    for i in range(0, nang):
        v_arr_nf[i, :] = v_bis[i].to_numpy()
        p_arr_nf[i, :] = p_bis[i].to_numpy()
        n_arr_nf[i, :] = n_bis[i].to_numpy()

    assert np.allclose(v_arr_nf, v_arr_f, 0.0, 1e-15)
    assert np.allclose(p_arr_nf, p_arr_f, 0.0, 1e-15)
    assert np.allclose(n_arr_nf, n_arr_f, 0.0, 1e-15)
