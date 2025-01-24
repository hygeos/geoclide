#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide.basic import Point, Ray, BBox
from geoclide.quadrics import Sphere
from geoclide.trianglemesh import Triangle, TriangleMesh
import xarray as xr
import numpy as np


def calc_intersection(shape, r1, method='v3'):
    """
    Performs intersection test between a shape and a ray and returns dataset

    Parameters
    ----------
    shape : BBox | Sphere | Triangle | TriangleMesh
        The shape used for the intersection
    r1 : Ray
        The ray used for the iuntersection
    method : str, optional
        Used in triangle intersection test. Only two choice : 'v2' and 'v3'.
        The 'v3' have more robustness tests, and 'v2' is faster.

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
        o                (xyz) float64 -2.0 0.0 0.8
        d                (xyz) float64 1.0 0.0 0.0
        mint             int64 0
        maxt             float64 inf
        shape            <U6 'Sphere'
        ...               ...
        oTw_mInv         (dim_0, dim_1) float64 1.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 1.0
        wTo_m            (dim_0, dim_1) float64 1.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 1.0
        wTo_mInv         (dim_0, dim_1) float64 1.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 1.0
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
        o                (xyz) float64 -2.0 0.0 0.8
        d                (xyz) float64 1.0 0.0 0.0
        mint             int64 0
        maxt             float64 inf
        shape            <U4 'BBox'
        pmin             (xyz) float64 0.0 0.0 0.0
        pmax             (xyz) float64 1.0 1.0 1.0
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
    elif(isinstance(shape, Sphere)      or
         isinstance(shape, Triangle)    or 
         isinstance(shape, TriangleMesh)):
        if (isinstance(shape, Triangle) or isinstance(shape, TriangleMesh)):
            thit, dg, is_intersection = shape.intersect(r1, method=method)
        else:
            thit, dg, is_intersection = shape.intersect(r1)
        phit = dg.p
        nhit = dg.n
    else:
        raise ValueError('The only supported shape are: BBox and Sphere')
    
    ds = xr.Dataset(coords={'xyz':np.arange(3)})
    ds['is_intersection'] = is_intersection
    ds['is_intersection'].attrs = {'description':'tells if there is an intersection between the ray and the shape'}

    ds['o'] = xr.DataArray(r1.o.to_numpy(), dims='xyz')
    ds['o'].attrs = {'type': 'Point', 'description':'the x, y and z components of the ray point'}
    ds['d'] = xr.DataArray(r1.d.to_numpy(), dims='xyz')
    ds['d'].attrs = {'type': 'Vector', 'description':'the x, y and z components of the ray vector'}
    ds['mint'] = r1.mint
    ds['mint'].attrs = {'description':'the mint attribut of the ray'}
    ds['maxt'] = r1.maxt
    ds['maxt'].attrs = {'description':'the maxt attribut of the ray'}

    if (isinstance(shape, BBox)):
        ds.attrs = {'shape': 'BBox'}
        ds['pmin'] = xr.DataArray(shape.pmin.to_numpy(), dims='xyz')
        ds['pmin'].attrs = {'type': 'Point', 'description':'the x, y and z components of the pmin BBox attribut'}
        ds['pmax'] = xr.DataArray(shape.pmax.to_numpy(), dims='xyz')
        ds['pmax'].attrs = {'type': 'Point', 'description':'the x, y and z components of the pmax BBox attribut'}
    if (isinstance(shape, Sphere)):
        ds.attrs = {'shape':  'Sphere'}
        ds['radius'] = shape.radius
        ds['radius'].attrs = {'description':'the sphere radius attribut'}
        ds['z_min'] = shape.zmin
        ds['z_min'].attrs = {'description':'the sphere zmin attribut'}
        ds['z_max'] = shape.zmax
        ds['z_max'].attrs = {'description':'the sphere zmax attribut'}
        ds['phi_max'] = shape.phiMax
        ds['phi_max'].attrs = {'unit':'Radian', 'description':'the sphere phiMax attribut'}
    if (isinstance(shape, Triangle)):
        ds.attrs = {'shape': 'Triangle'}
        ds['p0'] = xr.DataArray(shape.p0.to_numpy(), dims='xyz')
        ds['p0'].attrs = {'description': 'the triangle p0 attribut'}
        ds['p1'] = xr.DataArray(shape.p1.to_numpy(), dims='xyz')
        ds['p1'].attrs = {'description': 'the triangle p1 attribut'}
        ds['p2'] = xr.DataArray(shape.p2.to_numpy(), dims='xyz')
        ds['p2'].attrs = {'description': 'the triangle p2 attribut'}
    if (isinstance(shape, TriangleMesh)):
        ds.attrs = {'shape': 'TriangleMesh'}
        ds['v'] = xr.DataArray(np.array([pi.to_numpy() for pi in shape.vertices]), dims=['nvertices', 'xyz'])
        ds['v'].attrs = {'description': 'The vertices xyz coordinates.'}
        ds['vi'] = xr.DataArray(shape.vertices_index.reshape(shape.ntriangles, 3), dims=['ntriangles', 'p0p1p2'])
        ds['vi'].attrs = {'description': 'For each triangle, the index of vertices point p0, p1 and p2 (from variable v).'}
        ds.attrs.update({'ntriangles': shape.ntriangles,
                         'nvertices' : shape.nvertices})
    if (not isinstance(shape, BBox)):
        ds['wTo_m'] = xr.DataArray(shape.wTo.m)
        ds['wTo_m'].attrs = {'description':'the transformation matrix of the ' + str(ds.attrs['shape']).lower() + ' wTo attribut'}
        ds['wTo_mInv'] = xr.DataArray(shape.wTo.mInv)
        ds['wTo_mInv'].attrs = {'description':'the inverse transformation matrix of the ' + str(ds.attrs['shape']).lower() + ' wTo attribut'}
        ds['oTw_m'] = xr.DataArray(shape.oTw.m)
        ds['oTw_m'].attrs = {'description':'the transformation matrix of the ' + str(ds.attrs['shape']).lower() + ' oTw attribut'}
        ds['oTw_mInv'] = xr.DataArray(shape.oTw.mInv)
        ds['oTw_mInv'].attrs = {'description':'the inverse transformation matrix of the ' + str(ds.attrs['shape']).lower() +' oTw attribut'}



    ds['thit'] = thit
    ds['thit'].attrs = {'description':'the t ray factor for the intersection point calculation'}
    ds['phit'] = xr.DataArray(phit.to_numpy(), dims='xyz')
    ds['phit'].attrs = {'type': 'Point', 'description':'the x, y and z components of the intersection point'}
    if (nhit is not None):
        ds['nhit'] = xr.DataArray(nhit.to_numpy(), dims='xyz')
        ds['nhit'].attrs = {'type': 'Normal', 'description':'the x, y and z components of the normal at the intersection point'}

    return ds