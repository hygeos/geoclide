#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
import os
ROOTPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
import sys
sys.path.insert(0, ROOTPATH)
import geoclide as gc


def test_dot():
    assert (gc.dot(gc.Vector(0., 0., 1.), gc.Vector(0., 1., 0.)) == 0.)
    assert (gc.dot(gc.Vector(0., 0., 1.), gc.Vector(0., 0., 1.)) == 1.)
    v = gc.Vector(math.sqrt(2.)/2., 0., math.sqrt(2.)/2.)
    assert (gc.dot(gc.Vector(0., 0., 1.), v) == math.sqrt(2.)/2.)

    assert (gc.dot(gc.Normal(0., 0., 1.), gc.Vector(0., 1., 0.)) == 0.)
    assert (gc.dot(gc.Normal(0., 0., 1.), gc.Vector(0., 0., 1.)) == 1.)
    v = gc.Normal(math.sqrt(2.)/2., 0., math.sqrt(2.)/2.)
    assert (gc.dot(gc.Vector(0., 0., 1.), v) == math.sqrt(2.)/2.)

    assert (gc.dot(gc.Vector(0., 0., 1.), gc.Normal(0., 1., 0.)) == 0.)
    assert (gc.dot(gc.Vector(0., 0., 1.), gc.Normal(0., 0., 1.)) == 1.)
    v = gc.Normal(math.sqrt(2.)/2., 0., math.sqrt(2.)/2.)
    assert (gc.dot(gc.Vector(0., 0., 1.), v) == math.sqrt(2.)/2.)

    assert (gc.dot(gc.Normal(0., 0., 1.), gc.Normal(0., 1., 0.)) == 0.)
    assert (gc.dot(gc.Normal(0., 0., 1.), gc.Normal(0., 0., 1.)) == 1.)
    v = gc.Normal(math.sqrt(2.)/2., 0., math.sqrt(2.)/2.)
    assert (gc.dot(gc.Normal(0., 0., 1.), v) == math.sqrt(2.)/2.)


def test_cross():
    assert (gc.cross(gc.Vector(0., 0., 1.), gc.Vector(0., 1., 0.)) == gc.Vector(-1.0, 0.0, 0.0))
    assert (gc.cross(gc.Vector(0., 0., 1.), gc.Vector(0., 0., 1.)) == gc.Vector(0.0, 0.0, 0.0))
    v = gc.Vector(math.sqrt(2.)/2., 0., math.sqrt(2.)/2.)
    assert (gc.cross(gc.Vector(0., 0., 1.), v) ==  gc.Vector(0.0, math.sqrt(2.)/2., 0.0))

    assert (gc.cross(gc.Normal(0., 0., 1.), gc.Vector(0., 1., 0.)) == gc.Vector(-1.0, 0.0, 0.0))
    assert (gc.cross(gc.Normal(0., 0., 1.), gc.Vector(0., 0., 1.)) == gc.Vector(0.0, 0.0, 0.0))
    v = gc.Vector(math.sqrt(2.)/2., 0., math.sqrt(2.)/2.)
    assert (gc.cross(gc.Normal(0., 0., 1.), v) ==  gc.Vector(0.0, math.sqrt(2.)/2., 0.0))

    assert (gc.cross(gc.Vector(0., 0., 1.), gc.Normal(0., 1., 0.)) == gc.Vector(-1.0, 0.0, 0.0))
    assert (gc.cross(gc.Vector(0., 0., 1.), gc.Normal(0., 0., 1.)) == gc.Vector(0.0, 0.0, 0.0))
    v = gc.Normal(math.sqrt(2.)/2., 0., math.sqrt(2.)/2.)
    assert (gc.cross(gc.Vector(0., 0., 1.), v) ==  gc.Vector(0.0, math.sqrt(2.)/2., 0.0))


def test_normalize():
    length = math.sqrt(1. + 4. + 9.)
    a_1 = 1. / length
    a_2 = 2./ length
    a_3 = 3. / length
    assert (gc.normalize(gc.Vector(1.,2.,3.)) == gc.Vector(a_1, a_2, a_3))
    assert (gc.normalize(gc.Vector(-1.,-2.,-3.)) == gc.Vector(-a_1, -a_2, -a_3))

    assert (gc.normalize(gc.Normal(1.,2.,3.)) == gc.Normal(a_1, a_2, a_3))
    assert (gc.normalize(gc.Normal(-1.,-2.,-3.)) == gc.Normal(-a_1, -a_2, -a_3))


def test_coordinate_system():
    v1 = gc.Vector(0., 0., 1.)
    v2_m1, v3_m1 = gc.coordinate_system(v1, 'm1')
    assert (v2_m1 == gc.Vector(0.0, 1.0, -0.0))
    assert (v3_m1 == gc.Vector(-1.0, 0.0, 0.0))
    v2_m2, v3_m2 = gc.coordinate_system(v1, 'm2')
    assert (v2_m2 == gc.Vector(1.0, -0.0, -0.0))
    assert (v3_m2 == gc.Vector(-0.0, 1.0, -0.0))


def test_distance():
    p1 = gc.Point(0., 0., 0.)
    p2 = gc.Point(0., 0., 10.)
    assert (gc.distance(p1,p2) == 10.)
    p1 = gc.Point(1., 2., 1.9)
    p2 = gc.Point(5., 15., 3.)
    a_1 = math.sqrt(16 + 13**2 + 1.1**2)
    assert (gc.distance(p1,p2) == a_1)


def test_face_forward():
    n1 = gc.Normal(1., 0., 0.)
    v1 = gc.Vector(-1., 0., 0.)
    assert (gc.face_forward(v1, n1) == gc.Vector(1.0, -0.0, -0.0))
    n1 = gc.Normal(-1., 0., 0.)
    v1 = gc.Vector(-1., 0., 0.)
    assert (gc.face_forward(v1, n1) == gc.Vector(-1.0, 0.0, 0.0))


def test_max():
    v1 = gc.Vector(2.,3.,1.)
    assert (gc.max(v1) == 3.)


def test_min():
    v1 = gc.Vector(2.,3.,1.)
    assert (gc.min(v1) == 1.)


def test_argmax():
    v1 = gc.Vector(2.,3.,1.)
    assert (gc.argmax(v1) == 1)


def test_argmin():
    v1 = gc.Vector(2.,3.,1.)
    assert (gc.argmin(v1) == 2)


def test_permute():
    v1 = gc.Vector(2., 3., 1.)
    assert (gc.permute(v1, 1, 0, 2) == gc.Vector(3.0, 2.0, 1.0))
    assert (gc.permute(v1, np.array([1, 0, 2])) == gc.Vector(3.0, 2.0, 1.0))


def test_clamp():
    assert (gc.clamp(4., val_min=5., val_max=11.) == 5.)


def test_swap():
    a = 5.
    b = 12.
    assert (gc.swap(a, b) == (b, a))


def test_quadratic():
    a = 2.
    b = -5.
    c = 0.
    assert (gc.quadratic(a, b, c) == (True, 0.0, 2.5))

