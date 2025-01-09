#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide.basic import Point, Ray, BBox
from geoclide.quadrics import Sphere
import xarray as xr
import numpy as np


def calc_intersection(shape, r1):
    """
    Performs intersection test between a shape and a ray and returns dataset

    Parameters
    ----------
    shape : BBox | Sphere
        The shape used for the intersection
    r1 : Ray
        The ray used for the iuntersection

    Returns
    ------
    out : xr.Dataset
        Look-up table with the intersection information
    
    Examples
    --------
    >>> import geoclide as gc
    >>> sphere = gc.Sphere(radius=1.) # sphere of radius 1
    >>> bbox = gc.BBox(p1=gc.Point(0., 0., 0.), p2=gc.Point(1.,1.,1.))
    >>> ray = gc.Ray(o=gc.Point(-2., 0., 0.8), d=gc.Vector(1.,0.,0.))
    >>> ds_sphere = gc.calc_intersection(sphere, ray)
    >>> ds_sphere
    <xarray.Dataset>
    Dimensions:          (xyz: 3, dim_0: 4, dim_1: 4)
    Coordinates:
    * xyz              (xyz) int64 0 1 2
    Dimensions without coordinates: dim_0, dim_1
    Data variables: (12/17)
        is_intersection  bool True
        ray_o            (xyz) float64 -2.0 0.0 0.8
        ray_d            (xyz) float64 1.0 0.0 0.0
        ray_mint         int64 0
        ray_maxt         float64 inf
        shape            <U6 'Sphere'
        ...               ...
        sphere_oTw_mInv  (dim_0, dim_1) float64 1.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 1.0
        sphere_wTo_m     (dim_0, dim_1) float64 1.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 1.0
        sphere_wTo_mInv  (dim_0, dim_1) float64 1.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 1.0
        thit             float64 1.4
        phit             (xyz) float64 -0.6 0.0 0.8
        nhit             (xyz) float64 -0.6 0.0 0.8
    >>> ds_box = gc.calc_intersection(bbox, ray)
    >>> ds_bbox
    <xarray.Dataset>
    Dimensions:          (xyz: 3)
    Coordinates:
    * xyz              (xyz) int64 0 1 2
    Data variables:
        is_intersection  bool True
        ray_o            (xyz) float64 -2.0 0.0 0.8
        ray_d            (xyz) float64 1.0 0.0 0.0
        ray_mint         int64 0
        ray_maxt         float64 inf
        shape            <U4 'BBox'
        bbox_pmin        (xyz) float64 0.0 0.0 0.0
        bbox_pmax        (xyz) float64 1.0 1.0 1.0
        thit             float64 2.0
        phit             (xyz) float64 0.0 0.0 0.8
    """
    if (not isinstance(r1, Ray)):
        raise ValueError('The parameter r1 must a Ray')

    if (isinstance(shape, BBox)):
        t0, t1, is_intersection = shape.intersect(r1)
        if is_intersection:
            if t0 > 0: thit = t0
            else: thit = t1
            phit = r1[thit]
            nhit = None # TODO compute the real normal
        else:
            thit = 0.
            phit = Point(0., 0., 0.)
            nhit = None
    elif(isinstance(shape, Sphere)):
        thit, dg, is_intersection = shape.intersect(r1)
        phit = dg.p
        nhit = dg.n
    else:
        raise ValueError('The only supported shape are: BBox and Sphere')
    
    ds = xr.Dataset(coords={'xyz':np.arange(3)})
    ds['is_intersection'] = is_intersection
    ds['is_intersection'].attrs = {'description':'tells if there is an intersection between the ray and the shape'}

    ds['ray_o'] = xr.DataArray(r1.o.to_numpy(), dims='xyz')
    ds['ray_o'].attrs = {'type': 'Point', 'description':'the x, y and z components of the ray point'}
    ds['ray_d'] = xr.DataArray(r1.d.to_numpy(), dims='xyz')
    ds['ray_d'].attrs = {'type': 'Vector', 'description':'the x, y and z components of the ray vector'}
    ds['ray_mint'] = r1.mint
    ds['ray_mint'].attrs = {'description':'the mint attribut of the ray'}
    ds['ray_maxt'] = r1.maxt
    ds['ray_maxt'].attrs = {'description':'the maxt attribut of the ray'}

    if (isinstance(shape, BBox)):
        ds['shape'] = 'BBox'
        ds['bbox_pmin'] = xr.DataArray(shape.pmin.to_numpy(), dims='xyz')
        ds['bbox_pmin'].attrs = {'type': 'Point', 'description':'the x, y and z components of the pmin BBox attribut'}
        ds['bbox_pmax'] = xr.DataArray(shape.pmax.to_numpy(), dims='xyz')
        ds['bbox_pmax'].attrs = {'type': 'Point', 'description':'the x, y and z components of the pmax BBox attribut'}
    if (isinstance(shape, Sphere)):
        ds['shape'] = 'Sphere'
        ds['sphere_radius'] = shape.radius
        ds['sphere_radius'].attrs = {'description':'the sphere radius attribut'}
        ds['sphere_z_min'] = shape.zmin
        ds['sphere_z_min'].attrs = {'description':'the sphere zmin attribut'}
        ds['sphere_z_max'] = shape.zmax
        ds['sphere_z_max'].attrs = {'description':'the sphere zmax attribut'}
        ds['sphere_phi_max'] = shape.phiMax
        ds['sphere_phi_max'].attrs = {'unit':'Radian', 'description':'the sphere phiMax attribut'}
        ds['sphere_oTw_m'] = xr.DataArray(shape.oTw.m)
        ds['sphere_oTw_m'].attrs = {'description':'the transformation matrix of the sphere oTw attribut'}
        ds['sphere_oTw_mInv'] = xr.DataArray(shape.oTw.mInv)
        ds['sphere_oTw_mInv'].attrs = {'description':'the inverse transformation matrix of the sphere oTw attribut'}
        ds['sphere_wTo_m'] = xr.DataArray(shape.wTo.m)
        ds['sphere_wTo_m'].attrs = {'description':'the transformation matrix of the sphere wTo attribut'}
        ds['sphere_wTo_mInv'] = xr.DataArray(shape.wTo.mInv)
        ds['sphere_wTo_mInv'].attrs = {'description':'the inverse transformation matrix of the sphere wTo attribut'}

    ds['thit'] = thit
    ds['thit'].attrs = {'description':'the t ray factor for the intersection point calculation'}
    ds['phit'] = xr.DataArray(phit.to_numpy(), dims='xyz')
    ds['phit'].attrs = {'type': 'Point', 'description':'the x, y and z components of the intersection point'}
    if (nhit is not None):
        ds['nhit'] = xr.DataArray(nhit.to_numpy(), dims='xyz')
        ds['nhit'].attrs = {'type': 'Normal', 'description':'the x, y and z components of the normal at the intersection point'}

    return ds