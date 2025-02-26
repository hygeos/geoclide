#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide.basic import Vector, Point, Normal, Ray
from geoclide.vecope import normalize, cross, face_forward
import numpy as np
from geoclide.transform import Transform
import xarray as xr


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


class DifferentialGeometry(object):
    '''
    The general parametric description of most of the shapes

    - A point p is described by a function depending to variables u and v such that p=f(u,v)
    - From u and v we can get two directions (partial derivative of p) parallel to the shape surface

    Parameters
    ----------
    p : Point, optional
        The concerned position, at the surface of the given shape
    dpdu : Vector, optional
        The surface partial derivative of p with respect to u
    dpdv : vector, optional
        The surface partial derivative of p with respect to v
    u : float, optional
        The u coordinate of the parametric representation
    v : float, optional
        The v coordinate of the parametric representation
    ray_dir : float, optional
        The direction of the ray hitting the shape
    shape : Sphere | Triangle | ...
        The shape used
    '''
    def __init__(self, p=None, dpdu=None, dpdv=None, u=None, v=None, ray_dir=None, shape = None):
        if ( (p is None) and (dpdu is None) and (dpdv is None) and (u is None) and
             (v is None) and (ray_dir is None) and (shape is None) ):
            self.p = Point(0., 0., 0.)
            self.dpdu = Vector(0., 0., 0.)
            self.dpdv = Vector(0., 0., 0.)
            self.n = Normal(0., 0., 0.)
            self.u = 0.
            self.v = 0.
            self.ray_dir = Vector(0., 0., 0.)
            self.shape = None
        elif ( isinstance(p, Point)        and
               isinstance(dpdu, Vector)    and
               isinstance(dpdv, Vector)    and
               isinstance(ray_dir, Vector)):
            self.p = p
            self.dpdu = dpdu
            self.dpdv = dpdv
            self.ray_dir = ray_dir
            self.n = face_forward(Normal(normalize(cross(dpdu, dpdv))), -ray_dir)
            self.u = u
            self.v = v
            self.shape = shape
        else:
            raise ValueError('Problem with parameter(s)')


def get_intersect_dataset(shape_name, r, t=None, is_intersection=False, u=None, v=None,
                          dpdu=None, dpdv=None, diag_calc=False):
    """
    Create dataset containing the intersection test information

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
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims='xyz')
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims='xyz')
        ds['mint'] = r.mint
        ds['maxt'] = r.maxt
    elif (is_obj_arr and not is_r_arr): # multiple obj and 1 ray
        ds.attrs.update({'nobj': len(t)})
        ds['thit'] = xr.DataArray(t, dims=['nobj'])
        ds['u'] = xr.DataArray(u, dims=['nobj'])
        ds['v'] = xr.DataArray(v, dims=['nobj'])
        p = r[t]
        ds['phit'] = xr.DataArray(p, dims=['nobj', 'xyz'])
        n = face_forward(Normal(normalize(cross(Vector(dpdu), Vector(dpdv)))), -r.d)
        ds['nhit'] = xr.DataArray(n.to_numpy(), dims=['nobj', 'xyz'])
        ds['dpdu'] = xr.DataArray(dpdu, dims=['nobj', 'xyz'])
        ds['dpdv'] = xr.DataArray(dpdv, dims=['nobj', 'xyz'])
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims=['xyz'])
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims=['xyz'])
        ds['mint'] = r.mint
        ds['maxt'] = r.maxt 
    elif (not is_obj_arr and is_r_arr): # 1 obj and multiple rays
        ds.attrs.update({'nrays': len(r.o.x)})
        ds['thit'] = xr.DataArray(t, dims=['nrays'])
        ds['u'] = xr.DataArray(u, dims=['nrays'])
        ds['v'] = xr.DataArray(v, dims=['nrays'])
        p = r[t]
        ds['phit'] = xr.DataArray(p, dims=['nrays', 'xyz'])
        n = face_forward(Normal(normalize(cross(Vector(dpdu), Vector(dpdv)))), -r.d)
        ds['nhit'] = xr.DataArray(n.to_numpy(), dims=['nrays', 'xyz'])
        ds['dpdu'] = xr.DataArray(dpdu, dims=['xyz'])
        ds['dpdv'] = xr.DataArray(dpdv, dims=['xyz'])
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims=['nrays', 'xyz'])
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims=['nrays', 'xyz'])
        mint = np.zeros(ds.nrays, dtype=np.float64)
        maxt = np.zeros_like(mint)
        mint[:] = r.mint
        maxt[:] = r.maxt
        ds['mint'] = mint
        ds['maxt'] = maxt
    elif (diag_calc): # multiple shape and rays but calc only the diagonal
        size = len(t)
        ds.attrs.update({'nobj': size, 'nrays': size, 'ndiag':size})
        ds['thit'] = xr.DataArray(t, dims=['ndiag'])
        ds['u'] = xr.DataArray(u, dims=['ndiag'])
        ds['v'] = xr.DataArray(v, dims=['ndiag'])
        p = r[t]
        ds['phit'] = xr.DataArray(p, dims=['ndiag', 'xyz'])
        n = face_forward(Normal(normalize(cross(Vector(dpdu), Vector(dpdv)))), -r.d)
        ds['nhit'] = xr.DataArray(n.to_numpy(), dims=['ndiag', 'xyz'])
        ds['dpdu'] = xr.DataArray(dpdu, dims=['nobj', 'xyz'])
        ds['dpdv'] = xr.DataArray(dpdv, dims=['nobj', 'xyz'])
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims=['nrays', 'xyz'])
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims=['nrays', 'xyz'])
        mint = np.zeros(ds.diag, dtype=np.float64)
        maxt = np.zeros_like(mint)
        mint[:] = r.mint
        maxt[:] = r.maxt
        ds['mint'] = mint
        ds['maxt'] = maxt
    else: # multiple shape and rays, 2-D output
        nrays = t.shape[1]
        nobj = t.shape[0]
        ds.attrs.update({'nobj': nobj, 'nrays': nrays})
        ds['is_intersection'] = xr.DataArray(is_intersection, dims=['nobj', 'nrays'] )
        not_int = np.logical_not(is_intersection)
        ds['thit'] = xr.DataArray(t, dims=['nobj', 'nrays'])
        u[not_int] = None
        v[not_int] = None
        ds['u'] = xr.DataArray(u, dims=['nobj', 'nrays'])
        ds['v'] = xr.DataArray(v, dims=['nobj', 'nrays'])

        mint = np.zeros(nrays, dtype=np.float64)
        maxt = np.zeros_like(mint)
        mint[:] = r.mint
        maxt[:] = r.maxt
        ds['mint'] = mint
        ds['maxt'] = maxt
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims=['nrays', 'xyz'])
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims=['nrays', 'xyz'])
        

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

    ds['is_intersection'].attrs = {'description':'this variable tells if there is an intersection between the ray and the shape'}
    ds['thit'].attrs = {'description':'the t ray factor for the intersection point calculation'}
    ds['u'].attrs = {'description':'the u coordinate of the parametric representation'}
    ds['v'].attrs = {'description':'the v coordinate of the parametric representation'}
    ds['mint'].attrs = {'description':'the mint attribut of the ray'}
    ds['maxt'].attrs = {'description':'the maxt attribut of the ray'}
    ds['o'].attrs = {'type': 'Point', 'description':'the x, y and z components of the ray point'}
    ds['d'].attrs = {'type': 'Vector', 'description':'the x, y and z components of the ray vector'}
    ds['phit'].attrs = {'type': 'Point', 'description':'the x, y and z components of the intersection point'}
    ds['nhit'].attrs = {'type': 'Normal', 'description':'the x, y and z components of the normal at the intersection point'}
    ds['dpdu'].attrs = {'type': 'Vector', 'description':'the surface partial derivative of phit with respect to u'}
    ds['dpdv'].attrs = {'type': 'Vector', 'description':'the surface partial derivative of phit with respect to v'}
    ds.attrs = {'shape': shape_name}
    return ds

