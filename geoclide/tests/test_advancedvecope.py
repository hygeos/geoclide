#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import numpy as np
import os
ROOTPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
import sys
sys.path.insert(0, ROOTPATH)
import geoclide as gc

PHIS = np.linspace(0., 360., 6)
THETAS = np.linspace(-180., 180., 7)


@pytest.mark.parametrize('phi', PHIS)
@pytest.mark.parametrize('theta', THETAS)
def test_vec2ang(phi, theta):
    v1 = gc.ang2vec(theta=theta, phi=phi, vec_view='zenith')
    th, ph = gc.vec2ang(v1, vec_view='zenith')
    v2 = gc.ang2vec(theta=th, phi=ph, vec_view='zenith')
    assert (np.isclose(v1.x, v2.x, 0., 1e-14))
    assert (np.isclose(v1.y, v2.y, 0., 1e-14))
    assert (np.isclose(v1.z, v2.z, 0., 1e-14))

    v1 = gc.ang2vec(theta=theta, phi=phi, vec_view='nadir')
    th, ph = gc.vec2ang(v1, vec_view='nadir')
    v2 = gc.ang2vec(theta=th, phi=ph, vec_view='nadir')
    assert (np.isclose(v1.x, v2.x, 0., 1e-14))
    assert (np.isclose(v1.y, v2.y, 0., 1e-14))
    assert (np.isclose(v1.z, v2.z, 0., 1e-14))

