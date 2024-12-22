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

V1 = [np.array([5.,0.,0.]), np.array([5.,3.,1.])]
V2 = [np.array([5.,1.,1.]), np.array([5.,3.,2.])]
ANGLES = [30., 45., 60.]


@pytest.mark.parametrize('v_arr', V1)
def test_get_translate_tf(v_arr):
    t = gc.get_translate_tf(gc.Vector(v_arr))
    m = np.identity(4)
    m[0,-1] = v_arr[0]
    m[1,-1] = v_arr[1]
    m[2,-1] = v_arr[2]
    assert (np.all(t.m == m))
    mInv = np.identity(4)
    mInv[0,-1] = -v_arr[0]
    mInv[1,-1] = -v_arr[1]
    mInv[2,-1] = -v_arr[2]
    assert (np.all(t.mInv == mInv))


@pytest.mark.parametrize('v_arr', V2)
def test_get_scale_tf(v_arr):
    t = gc.get_scale_tf(v_arr[0], v_arr[1], v_arr[2])
    m = np.identity(4)
    m[0,0] = v_arr[0]
    m[1,1] = v_arr[1]
    m[2,2] = v_arr[2]
    assert (np.all(t.m == m))
    mInv = np.identity(4)
    mInv[0,0] = 1 * (1/v_arr[0])
    mInv[1,1] = 1 * (1/v_arr[1])
    mInv[2,2] = 1 * (1/v_arr[2])
    assert (np.all(t.mInv == mInv))


@pytest.mark.parametrize('angle', ANGLES)
def test_get_rotateX_tf(angle):
    t = gc.get_rotateX_tf(angle)
    sin_t = math.sin(angle*(math.pi / 180.))
    cos_t = math.cos(angle*(math.pi / 180.))
    m = np.identity(4)
    m[1,1] = cos_t
    m[1,2] = -1.*sin_t
    m[2,1] = sin_t
    m[2,2] = cos_t
    assert (np.all(t.m == m))
    mInv = np.transpose(m)
    assert (np.all(t.mInv == mInv))


@pytest.mark.parametrize('angle', ANGLES)
def test_get_rotateY_tf(angle):
    t = gc.get_rotateY_tf(angle)
    sin_t = math.sin(angle*(math.pi / 180.))
    cos_t = math.cos(angle*(math.pi / 180.))
    m = np.identity(4)
    m[0,0] = cos_t
    m[2,0] = -1.*sin_t
    m[0,2] = sin_t
    m[2,2] = cos_t
    assert (np.all(t.m == m))
    mInv = np.transpose(m)
    assert (np.all(t.mInv == mInv))


@pytest.mark.parametrize('angle', ANGLES)
def test_get_rotateZ_tf(angle):
    t = gc.get_rotateZ_tf(angle)
    sin_t = math.sin(angle*(math.pi / 180.))
    cos_t = math.cos(angle*(math.pi / 180.))
    m = np.identity(4)
    m[0,0] = cos_t
    m[0,1] = -1.*sin_t
    m[1,0] = sin_t
    m[1,1] = cos_t
    assert (np.all(t.m == m))
    mInv = np.transpose(m)
    assert (np.all(t.mInv == mInv))


@pytest.mark.parametrize('angle', ANGLES)
def test_get_rotate_tf_1(angle):
    v_x = gc.Vector(1., 0., 0.)
    tx = gc.get_rotateX_tf(angle)
    t = gc.get_rotate_tf(angle, v_x)
    assert (np.all(tx.m == t.m))
    assert (np.all(tx.mInv == t.mInv))

    v_y = gc.Vector(0., 1., 0.)
    ty = gc.get_rotateY_tf(angle)
    t = gc.get_rotate_tf(angle, v_y)
    assert (np.all(ty.m == t.m))
    assert (np.all(ty.mInv == t.mInv))

    v_z = gc.Vector(0., 0., 1.)
    tz = gc.get_rotateZ_tf(angle)
    t = gc.get_rotate_tf(angle, v_z)
    assert (np.all(tz.m == t.m))
    assert (np.all(tz.mInv == t.mInv))


@pytest.mark.parametrize('angle', ANGLES)
@pytest.mark.parametrize('v_arr', V1)
def test_get_rotate_tf_2(angle, v_arr):
    v = gc.normalize(gc.Vector(v_arr)) # important to normalize for Rodrigues formula
    t = gc.get_rotate_tf(angle, v)
    # Use independant method -> Rodrigues rotation formula
    m = np.identity(4)
    cos_t = math.cos(angle*(math.pi / 180.))
    sin_t = math.sqrt(1.-cos_t*cos_t)
    matA = np.identity(3)
    matB = np.zeros((3,3))
    matB[0,1] = -v.z
    matB[0,2] =  v.y
    matB[1,0] =  v.z
    matB[1,2] = -v.x
    matB[2,0] = -v.y
    matB[2,1] =  v.x
    matC = matB.dot(matB)
    m[0:3,0:3] = (matA + matB*sin_t + matC*(1-cos_t))
    assert (np.allclose(t.m, m, rtol=0., atol=1e-15))
    mInv = np.transpose(m)
    assert (np.allclose(t.mInv, mInv, rtol=0., atol=1e-15))


def test_transform():
    t1 = gc.Transform()
    assert (np.all(t1.m == np.identity(4)))
    assert (np.all(t1.mInv == np.identity(4)))