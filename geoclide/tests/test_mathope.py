#!/usr/bin/env python

import geoclide as gc


def test_clamp():
    assert gc.clamp(4.0, val_min=5.0, val_max=11.0) == 5.0


def test_quadratic():
    a = 2.0
    b = -5.0
    c = 0.0
    assert gc.quadratic(a, b, c) == (True, 0.0, 2.5)
