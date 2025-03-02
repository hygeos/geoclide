#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide.basic import Vector, Point, Normal, Ray
from geoclide.vecope import normalize, cross, face_forward
import numpy as np
from geoclide.transform import Transform
import xarray as xr
from datetime import datetime
from geoclide.constante import VERSION


class Shape(object):
    '''
    The parent class of the different shape classes
    '''
    indShape = 0
    def __init__(self, ObjectToWorld, WorldToObject):
        if (not isinstance(ObjectToWorld, Transform) or not isinstance(WorldToObject, Transform)):
            raise ValueError('The parameters oTw and wTo must be both Transform')
        self.oTw = ObjectToWorld
        self.wTo = WorldToObject


def get_intersect_dataset(shape_name, r, t=None, is_intersection=False, u=None, v=None,
                          dpdu=None, dpdv=None, diag_calc=False):
    """
    Create dataset containing the intersection test information

    - The intersect method return (when parameter ds_output = False) of all quadrics, 
      Triangle and TriangleMesh classes gives directly the inputs of this function

    Parameters
    ----------
    shape_name : str
        The shape class name
    r : Ray
        The ray(s) used for the intersection test
    t : float | 1-D ndarray | 2-D ndarray, optional
        The t ray variable for its first intersection at the shape surface
    is_intersection : bool | 1-D ndarray | 2-D ndarray, optional
        If there is an intersection -> True, else False
    u : float | 1-D ndarray | 2-D ndarray, optional
        The u coordinate of the parametric representation
    v : float | 1-D ndarray | 2-D ndarray, optional
        The u coordinate of the parametric representation
    dpdu : 1-D ndarray | 2-D ndarray, optional
        The surface partial derivative of phit with respect to u
    dpdv : 1-D ndarray | 2-D ndarray, optional
        The surface partial derivative of phit with respect to v
    diag_cal : bool, optional
            This indicates whether a diagonal calculation has been performed

    Returns
    -------
    out : xr.Dataset
        Look-up table with the intersection information
    """
    if not isinstance(r, Ray): raise ValueError('The r parameter must be a Ray')
    if dpdu is None: dpdu = np.array([np.nan, np.nan, np.nan])
    if dpdv is None: dpdv = np.array([np.nan, np.nan, np.nan])

    is_r_arr = isinstance(r.o.x, np.ndarray)
    is_obj_arr = len(dpdu.shape) == 2
    ds = xr.Dataset(coords={'xyz':np.arange(3)})
    not_int = np.logical_not(is_intersection)
    
    if is_r_arr:
        nrays = len(r.o.x)
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims=['nrays', 'xyz'])
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims=['nrays', 'xyz'])
        mint = np.zeros(nrays, dtype=np.float64)
        maxt = np.zeros_like(mint)
        mint[:] = r.mint
        maxt[:] = r.maxt
        ds['mint'] = xr.DataArray(mint, dims=['nrays'])
        ds['maxt'] = xr.DataArray(maxt, dims=['nrays'])
    else:
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims=['xyz'])
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims=['xyz'])
        ds['mint'] = xr.DataArray(r.mint)
        ds['maxt'] = xr.DataArray(r.maxt)

    if (not isinstance(t, np.ndarray)):
        ds['is_intersection'] = xr.DataArray(is_intersection)
        ds['thit'] = t
        ds['u'] = u
        ds['v'] = v
        if t is None: p = np.array([np.nan, np.nan, np.nan])
        else:
            p = r[t]
            p = p.to_numpy()
        ds['phit'] = xr.DataArray(p, dims='xyz')
        if t is None: n = np.array([np.nan, np.nan, np.nan])
        else:
            n = face_forward(Normal(normalize(cross(Vector(dpdu), Vector(dpdv)))), -r.d)
            n = n.to_numpy()
        ds['nhit'] = xr.DataArray(n, dims='xyz')
        ds['dpdu'] = xr.DataArray(dpdu, dims='xyz')
        ds['dpdv'] = xr.DataArray(dpdv, dims='xyz')
    elif (is_obj_arr and not is_r_arr): # multiple obj and 1 ray
        ds.attrs.update({'nobj': len(t)})
        ds['is_intersection'] = xr.DataArray(is_intersection, dims=['nobj'])
        ds['thit'] = xr.DataArray(t, dims=['nobj'])
        u[not_int] = None
        v[not_int] = None
        ds['u'] = xr.DataArray(u, dims=['nobj'])
        ds['v'] = xr.DataArray(v, dims=['nobj'])
        p = r[t].to_numpy()
        ds['phit'] = xr.DataArray(p, dims=['nobj', 'xyz'])
        n = face_forward(Normal(normalize(cross(Vector(dpdu), Vector(dpdv)))), -r.d)
        n = n.to_numpy()
        n[not_int] = None
        ds['nhit'] = xr.DataArray(n, dims=['nobj', 'xyz'])
        dpdu[not_int] = None
        dpdv[not_int] = None
        ds['dpdu'] = xr.DataArray(dpdu, dims=['nobj', 'xyz'])
        ds['dpdv'] = xr.DataArray(dpdv, dims=['nobj', 'xyz'])
    elif (not is_obj_arr and is_r_arr): # 1 obj and multiple rays
        ds.attrs.update({'nrays': nrays})
        ds['is_intersection'] = xr.DataArray(is_intersection, dims=['nrays'])
        ds['thit'] = xr.DataArray(t, dims=['nrays'])
        u[not_int] = None
        v[not_int] = None
        ds['u'] = xr.DataArray(u, dims=['nrays'])
        ds['v'] = xr.DataArray(v, dims=['nrays'])
        p = r[t].to_numpy()
        ds['phit'] = xr.DataArray(p, dims=['nrays', 'xyz'])
        n = face_forward(Normal(normalize(cross(Vector(dpdu), Vector(dpdv)))), -r.d)
        n = n.to_numpy()
        n[not_int] = None
        ds['nhit'] = xr.DataArray(n, dims=['nrays', 'xyz'])
        dpdu = np.full((nrays, 3), dpdu, dtype=np.float64)
        dpdv = np.full((nrays, 3), dpdv, dtype=np.float64)
        dpdu[not_int] = None
        dpdv[not_int] = None
        ds['dpdu'] = xr.DataArray(dpdu, dims=['nrays', 'xyz'])
        ds['dpdv'] = xr.DataArray(dpdv, dims=['nrays', 'xyz'])
    elif (diag_calc or shape_name == 'TriangleMesh'): # multiple shape and rays but calc only the diagonal
        ndiag = len(t)
        if shape_name == 'TriangleMesh' and not diag_calc:
            dim_name = 'nrays'
            ds.attrs.update({'nrays': ndiag})
        else:
            dim_name = 'ndiag'
            ds.attrs.update({'nobj': ndiag, 'nrays': ndiag, 'ndiag':ndiag})
        ds['is_intersection'] = xr.DataArray(is_intersection, dims=[dim_name])
        ds['thit'] = xr.DataArray(t, dims=[dim_name])
        u[not_int] = None
        v[not_int] = None
        ds['u'] = xr.DataArray(u, dims=[dim_name])
        ds['v'] = xr.DataArray(v, dims=[dim_name])
        p = r[t].to_numpy()
        ds['phit'] = xr.DataArray(p, dims=[dim_name, 'xyz'])
        
        n = face_forward(Normal(normalize(cross(Vector(dpdu), Vector(dpdv)))), -r.d)
        n = n.to_numpy()
        n[not_int] = None
        ds['nhit'] = xr.DataArray(n, dims=[dim_name, 'xyz'])
        dpdu[not_int] = None
        dpdv[not_int] = None
        ds['dpdu'] = xr.DataArray(dpdu, dims=[dim_name, 'xyz'])
        ds['dpdv'] = xr.DataArray(dpdv, dims=[dim_name, 'xyz'])
    else: # multiple shape and rays, 2-D output
        nobj = t.shape[0]
        ds.attrs.update({'nobj': nobj, 'nrays': nrays})
        ds['is_intersection'] = xr.DataArray(is_intersection, dims=['nobj', 'nrays'] )
        ds['thit'] = xr.DataArray(t, dims=['nobj', 'nrays'])
        u[not_int] = None
        v[not_int] = None
        ds['u'] = xr.DataArray(u, dims=['nobj', 'nrays'])
        ds['v'] = xr.DataArray(v, dims=['nobj', 'nrays'])

        n = np.zeros((nobj, nrays,3), dtype=np.float64)
        p = np.zeros_like(n)

        if nobj >= nrays:
            ro = r.o.to_numpy()
            rd = r.d.to_numpy()
            n_bis = Normal(normalize(cross(Vector(dpdu), Vector(dpdv))))
            for ir in range (0, nrays):
                ri = Ray(Point(ro[ir,:]), Vector(rd[ir,:]), mint[ir], maxt[ir])
                n[:,ir,:] = face_forward(n_bis, -ri.d).to_numpy()
                p[:,ir,:] = ri[t[:,ir]].to_numpy()
        else:
            for iobj in range (0, nobj):
                n_bis = Normal(normalize(cross(Vector(dpdu[iobj,:]), Vector(dpdv[iobj,:]))))
                n[iobj,:,:] = face_forward(n_bis, -r.d).to_numpy()
                p[iobj,:,:] = r[t[iobj,:]].to_numpy()
        n[not_int,:] = None
        
        # n = face_forward(Normal(normalize(cross(Vector(dpdu), Vector(dpdv)))), -r.d)
        ds['phit'] = xr.DataArray(p, dims=['nobj', 'nrays', 'xyz'])
        ds['nhit'] = xr.DataArray(n, dims=['nobj', 'nrays', 'xyz'])
        
        dpdu_2d = np.repeat(dpdu[:,np.newaxis,:], nrays, axis=1)
        dpdv_2d = np.repeat(dpdv[:,np.newaxis,:], nrays, axis=1)
        dpdu_2d[not_int,:] = None
        dpdv_2d[not_int,:] = None
        ds['dpdu'] = xr.DataArray(dpdu_2d, dims=['nobj', 'nrays', 'xyz'])
        ds['dpdv'] = xr.DataArray(dpdv_2d, dims=['nobj', 'nrays', 'xyz'])

    ds['o'].attrs = {'type': 'Point', 'description':'the x, y and z components of the ray point'}
    ds['d'].attrs = {'type': 'Vector', 'description':'the x, y and z components of the ray vector'}
    ds['mint'].attrs = {'description':'the mint attribut of the ray'}
    ds['maxt'].attrs = {'description':'the maxt attribut of the ray'}
    ds['is_intersection'].attrs = {'description':'this variable tells if there is an intersection between the ray and the shape'}
    ds['thit'].attrs = {'description':'the t ray factor for the intersection point calculation'}
    ds['u'].attrs = {'description':'the u coordinate of the parametric representation'}
    ds['v'].attrs = {'description':'the v coordinate of the parametric representation'}
    ds['phit'].attrs = {'type': 'Point', 'description':'the x, y and z components of the intersection point'}
    ds['nhit'].attrs = {'type': 'Normal', 'description':'the x, y and z components of the normal at the intersection point'}
    ds['dpdu'].attrs = {'type': 'Vector', 'description':'the surface partial derivative of phit with respect to u'}
    ds['dpdv'].attrs = {'type': 'Vector', 'description':'the surface partial derivative of phit with respect to v'}
    ds.attrs = {'shape': shape_name}
    date = datetime.now().strftime("%Y-%m-%d")  
    ds.attrs.update({'date':date, 'version': VERSION})
    return ds


def get_bbox_intersect_dataset(bbox, r, t0=None, t1=None, is_intersection=False, diag_calc=False):
    """
    Create dataset containing the intersection test information

    - The intersect method return of BBox class gives the t0, t1 and is_intersection 
      inputs of this function

    Parameters
    ----------
    bbox : BBox
        The bounding box(es) used for the intersection test
    r : Ray
        The ray(s) used for the intersection test
    t0 : float | 1-D ndarray | 2-D ndarray
        The t ray variable for the first intersection
    t1 : float | 1-D ndarray | 2-D ndarray
        The t ray variable for the second intersection
    is_intersection : bool | 1-D ndarray | 2-D ndarray, optional
        If there is an intersection -> True, else False
    diag_cal : bool, optional
            This indicates whether a diagonal calculation has been performed

    Returns
    -------
    out : xr.Dataset
        Look-up table with the intersection information
    """
    is_r_arr = isinstance(r.o.x, np.ndarray)
    is_bbox_arr = isinstance(bbox.pmin.x, np.ndarray)

    ds = xr.Dataset(coords={'xyz':np.arange(3)})

    if is_r_arr:
        nrays = len(r.o.x)
        ro = r.o.to_numpy()
        rd = r.d.to_numpy()
        ds['o'] = xr.DataArray(ro, dims=['nrays', 'xyz'])
        ds['d'] = xr.DataArray(rd, dims=['nrays', 'xyz'])
        mint = np.zeros(nrays, dtype=np.float64)
        maxt = np.zeros_like(mint)
        mint[:] = r.mint
        maxt[:] = r.maxt
        ds['mint'] = xr.DataArray(mint, dims=['nrays'])
        ds['maxt'] = xr.DataArray(maxt, dims=['nrays'])
    else:
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims=['xyz'])
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims=['xyz'])
        ds['mint'] = xr.DataArray(r.mint)
        ds['maxt'] = xr.DataArray(r.maxt)

    if is_r_arr or is_bbox_arr:
        c1 = t0 > 0
        not_c1 = np.logical_not(c1)
        thit = np.full((t0.shape), np.nan, dtype=np.float64)
        c2 = np.logical_and(is_intersection, c1)
        c3 = np.logical_and(is_intersection, not_c1)
        if np.any(c2): thit[c2] = t0[c2]
        if np.any(c3): thit[c3] = t1[c3]

    if is_r_arr and is_bbox_arr and not diag_calc:
        nobj = len(bbox.p0.x)
        ds.attrs.update({'nobj': nobj, 'nrays': nrays})
        ds['is_intersection'] = xr.DataArray(is_intersection, dims=['nobj', 'nrays'])
        ds['thit'] = xr.DataArray(thit, dims=['nobj', 'nrays'])
        phit = np.zeros((nobj, nrays, 3), dtype=np.float64)
        for ir in range (0, nrays):
            ri = Ray(Point(ro[ir,:]), Vector(rd[ir,:]), mint[ir], maxt[ir])
            phit[:,ir,:] = ri[thit[:,ir]].to_numpy()
        ds['phit'] = xr.DataArray(phit, dims=['nobj', 'nrays', 'xyz'])
    elif is_r_arr or is_bbox_arr:
        if diag_calc :
            dim_name = 'ndiag'
            size = nrays
            ds.attrs.update({'nobj': size, 'nrays': size})
        elif is_r_arr:
            dim_name = 'nrays'
            size = nrays
        else:
            dim_name = 'nobj'
            size = len(bbox.p0.x)
        ds.attrs.update({dim_name: size})
        phit = r[thit]
        ds['is_intersection'] = xr.DataArray(is_intersection, dims=[dim_name])
        ds['thit'] = xr.DataArray(thit, dims=[dim_name])
        ds['phit'] = xr.DataArray(phit, dims=[dim_name, 'xyz'])
    else:
        if (t0 is None) : thit = None
        elif (t0 > 0) : thit = t0
        else : thit = t1
        phit = r[thit]
        ds['is_intersection'] = xr.DataArray(is_intersection)
        ds['thit'] = xr.DataArray(thit)
        ds['phit'] = xr.DataArray(phit, dims=['xyz'])
        
    ds['o'].attrs = {'type': 'Point', 'description':'the x, y and z components of the ray point'}
    ds['d'].attrs = {'type': 'Vector', 'description':'the x, y and z components of the ray vector'}
    ds['mint'].attrs = {'description':'the mint attribut of the ray'}
    ds['maxt'].attrs = {'description':'the maxt attribut of the ray'}
    ds['is_intersection'].attrs = {'description':'this variable tells if there is an intersection between the ray and the shape'}
    ds['thit'].attrs = {'description':'the t ray factor for the intersection point calculation'}
    ds['phit'].attrs = {'type': 'Point', 'description':'the x, y and z components of the intersection point'}
    ds.attrs = {'shape': bbox.__class__.__name__}
    date = datetime.now().strftime("%Y-%m-%d")  
    ds.attrs.update({'date':date, 'version': VERSION})
    return ds