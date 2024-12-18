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

P1 = [np.array([1.,2.,3.])]


def test_vector():
    v1 = gc.Vector(0.,0.,1.)
    v2 = gc.Vector(np.array([0.,0.,1.]))
    assert (v1 == v2)
    assert (np.all(v1.to_numpy() == np.array([0.,0.,1.])))
    assert (v1[0] == 0.)
    assert (v1[1] == 0.)
    assert (v1[2] == 1.)
    assert (v1[0] == v1.x)
    assert (v1[1] == v1.y)
    assert (v1[2] == v1.z)

@pytest.mark.parametrize('p_arr', P1)
def test_point(p_arr):
    p1 = gc.Point(p_arr[0], p_arr[1], p_arr[2])
    p2 = gc.Point(p_arr)
    assert (p1 == p2)
    assert (np.all(p1.to_numpy() == p_arr))
    assert (p1[0] == 1.)
    assert (p1[1] == 2.)
    assert (p1[2] == 3.)
    assert (p1[0] == p1.x)
    assert (p1[1] == p1.y)
    assert (p1[2] == p1.z)


def test_normal():
    n1 = gc.Normal(0.,0.,1.)
    n2 = gc.Normal(np.array([0.,0.,1.]))
    assert (n1 == n2)
    assert (np.all(n1.to_numpy() == np.array([0.,0.,1.])))
    assert (n1[0] == 0.)
    assert (n1[1] == 0.)
    assert (n1[2] == 1.)
    assert (n1[0] == n1.x)
    assert (n1[1] == n1.y)
    assert (n1[2] == n1.z)


def test_opeVector():
    v1 = gc.Vector(0., 0., 1.)
    v2 = gc.Vector(1., 0., 0.)
    assert (v1+v2 == gc.Vector(1., 0., 1.))
    assert (v1-v2 == gc.Vector(-1., 0., 1.))
    assert (v1*2 == gc.Vector(0., 0., 2.))
    assert (v1/2 == gc.Vector(0., 0., 0.5))
    v3 = v1+v2
    assert (v3.length() == math.sqrt(v3[0]**2+v3[1]**2+v3[2]**2))


def test_opePoint():
    p1 = gc.Point(1.,2.,3.)
    v1 = gc.Vector(0., 0., 1.)
    p2 = gc.Point(1.,1.,1.)
    assert (p1+v1 == gc.Point(1., 2., 4.))
    assert (p1-p2 == gc.Vector(0.0, 1.0, 2.0))
    assert (p1-v1 == gc.Point(1.0, 2.0, 2.0))
    assert (p1*2 == gc.Point(2.0, 4.0, 6.0))
    assert (p1/2 == gc.Point(0.5, 1.0, 1.5))


def test_opeNormal():
    n1 = gc.Normal(0.,0.,1.)
    n2 = gc.Normal(1., 0., 0.)
    assert (n1+n2 == gc.Normal(1.0, 0.0, 1.0))
    assert (n1-n2 == gc.Normal(-1.0, 0.0, 1.0))
    assert (n1*2 == gc.Normal(0.0, 0.0, 2.0))
    assert (n1/2 == gc.Normal(0.0, 0.0, 0.5))