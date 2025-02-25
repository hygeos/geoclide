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


def get_intersect_dataset(t, is_intersection, dpdu, dpdv, u, v, r, shape, diag_calc):
    """
    """

    is_r_arr = isinstance(r.o.x, np.ndarray)
    is_obj_arr = len(dpdu.shape) == 2
    ds = xr.Dataset(coords={'xyz':np.arange(3)})

    if (not isinstance(t, np.ndarray)):
        ds['thit'] = t
        ds['thit'].attrs = {'description':'the t ray factor for the intersection point calculation'}
        ds['u'] = u
        ds['u'].attrs = {'description':'The u coordinate of the parametric representation'}
        ds['v'] = v
        ds['v'].attrs = {'description':'The v coordinate of the parametric representation'}
        p = r[t]
        ds['phit'] = xr.DataArray(p, dims='xyz')
        ds['phit'].attrs = {'type': 'Point', 'description':'the x, y and z components of the intersection point'}
        n = face_forward(Normal(normalize(cross(Vector(dpdu), Vector(dpdv)))), -r.d)
        ds['nhit'] = xr.DataArray(n.to_numpy(), dims='xyz')
        ds['nhit'].attrs = {'type': 'Normal', 'description':'the x, y and z components of the normal at the intersection point'}
        ds['dpdu'] = xr.DataArray(dpdu, dims='xyz')
        ds['dpdu'].attrs = {'type': 'Vector', 'description':'The surface partial derivative of phit with respect to u'}
        ds['dpdv'] = xr.DataArray(dpdv, dims='xyz')
        ds['dpdv'].attrs = {'type': 'Vector', 'description':'The surface partial derivative of phit with respect to v'}
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims='xyz')
        ds['o'].attrs = {'type': 'Point', 'description':'the x, y and z components of the ray point'}
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims='xyz')
        ds['d'].attrs = {'type': 'Vector', 'description':'the x, y and z components of the ray vector'}
        ds['mint'] = r.mint
        ds['mint'].attrs = {'description':'the mint attribut of the ray'}
        ds['maxt'] = r.maxt
        ds['maxt'].attrs = {'description':'the maxt attribut of the ray'}
    elif (is_obj_arr and not is_r_arr): # multiple obj and 1 ray
        ds.attrs.update({'nobj': len(t)})
        ds['thit'] = xr.DataArray(t, dims=['nobj'])
        ds['thit'].attrs = {'description':'the t ray factor for the intersection point calculation'}
        ds['u'] = xr.DataArray(u, dims=['nobj'])
        ds['u'].attrs = {'description':'The u coordinate of the parametric representation'}
        ds['v'] = xr.DataArray(v, dims=['nobj'])
        ds['v'].attrs = {'description':'The v coordinate of the parametric representation'}
        p = r[t]
        ds['phit'] = xr.DataArray(p, dims=['nobj', 'xyz'])
        ds['phit'].attrs = {'type': 'Point', 'description':'the x, y and z components of the intersection point'}
        n = face_forward(Normal(normalize(cross(Vector(dpdu), Vector(dpdv)))), -r.d)
        ds['nhit'] = xr.DataArray(n.to_numpy(), dims=['nobj', 'xyz'])
        ds['nhit'].attrs = {'type': 'Normal', 'description':'the x, y and z components of the normal at the intersection point'}
        ds['dpdu'] = xr.DataArray(dpdu, dims=['nobj', 'xyz'])
        ds['dpdu'].attrs = {'type': 'Vector', 'description':'The surface partial derivative of phit with respect to u'}
        ds['dpdv'] = xr.DataArray(dpdv, dims=['nobj', 'xyz'])
        ds['dpdv'].attrs = {'type': 'Vector', 'description':'The surface partial derivative of phit with respect to v'}
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims=['xyz'])
        ds['o'].attrs = {'type': 'Point', 'description':'the x, y and z components of the ray point'}
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims=['xyz'])
        ds['d'].attrs = {'type': 'Vector', 'description':'the x, y and z components of the ray vector'}
        ds['mint'] = r.mint
        ds['mint'].attrs = {'description':'the mint attribut of the ray'}
        ds['maxt'] = r.maxt
        ds['maxt'].attrs = {'description':'the maxt attribut of the ray'}  
    elif (not is_obj_arr and is_r_arr): # 1 obj and multiple rays
        ds.attrs.update({'nrays': len(r.o.x)})
        ds['thit'] = xr.DataArray(t, dims=['nrays'])
        ds['thit'].attrs = {'description':'the t ray factor for the intersection point calculation'}
        ds['u'] = xr.DataArray(u, dims=['nrays'])
        ds['u'].attrs = {'description':'The u coordinate of the parametric representation'}
        ds['v'] = xr.DataArray(v, dims=['nrays'])
        ds['v'].attrs = {'description':'The v coordinate of the parametric representation'}
        p = r[t]
        ds['phit'] = xr.DataArray(p, dims=['nrays', 'xyz'])
        ds['phit'].attrs = {'type': 'Point', 'description':'the x, y and z components of the intersection point'}
        n = face_forward(Normal(normalize(cross(Vector(dpdu), Vector(dpdv)))), -r.d)
        ds['nhit'] = xr.DataArray(n.to_numpy(), dims=['nrays', 'xyz'])
        ds['nhit'].attrs = {'type': 'Normal', 'description':'the x, y and z components of the normal at the intersection point'}
        ds['dpdu'] = xr.DataArray(dpdu, dims=['xyz'])
        ds['dpdu'].attrs = {'type': 'Vector', 'description':'The surface partial derivative of phit with respect to u'}
        ds['dpdv'] = xr.DataArray(dpdv, dims=['xyz'])
        ds['dpdv'].attrs = {'type': 'Vector', 'description':'The surface partial derivative of phit with respect to v'}
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims=['nrays', 'xyz'])
        ds['o'].attrs = {'type': 'Point', 'description':'the x, y and z components of the ray point'}
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims=['nrays', 'xyz'])
        ds['d'].attrs = {'type': 'Vector', 'description':'the x, y and z components of the ray vector'}
        mint = np.zeros(ds.nrays, dtype=np.float64)
        maxt = np.zeros_like(mint)
        mint[:] = r.mint
        maxt[:] = r.maxt
        ds['mint'] = mint
        ds['mint'].attrs = {'description':'the mint attribut of the ray'}
        ds['maxt'] = maxt
        ds['maxt'].attrs = {'description':'the maxt attribut of the ray'}
    elif (diag_calc): # multiple shape and rays but calc only the diagonal
        size = len(t)
        ds.attrs.update({'nobj': size, 'nrays': size, 'ndiag':size})
        ds['thit'] = xr.DataArray(t, dims=['ndiag'])
        ds['thit'].attrs = {'description':'the t ray factor for the intersection point calculation'}
        ds['u'] = xr.DataArray(u, dims=['ndiag'])
        ds['u'].attrs = {'description':'The u coordinate of the parametric representation'}
        ds['v'] = xr.DataArray(v, dims=['ndiag'])
        ds['v'].attrs = {'description':'The v coordinate of the parametric representation'}
        p = r[t]
        ds['phit'] = xr.DataArray(p, dims=['ndiag', 'xyz'])
        ds['phit'].attrs = {'type': 'Point', 'description':'the x, y and z components of the intersection point'}
        n = face_forward(Normal(normalize(cross(Vector(dpdu), Vector(dpdv)))), -r.d)
        ds['nhit'] = xr.DataArray(n.to_numpy(), dims=['ndiag', 'xyz'])
        ds['nhit'].attrs = {'type': 'Normal', 'description':'the x, y and z components of the normal at the intersection point'}
        ds['dpdu'] = xr.DataArray(dpdu, dims=['nobj', 'xyz'])
        ds['dpdu'].attrs = {'type': 'Vector', 'description':'The surface partial derivative of phit with respect to u'}
        ds['dpdv'] = xr.DataArray(dpdv, dims=['nobj', 'xyz'])
        ds['dpdv'].attrs = {'type': 'Vector', 'description':'The surface partial derivative of phit with respect to v'}
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims=['nrays', 'xyz'])
        ds['o'].attrs = {'type': 'Point', 'description':'the x, y and z components of the ray point'}
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims=['nrays', 'xyz'])
        ds['d'].attrs = {'type': 'Vector', 'description':'the x, y and z components of the ray vector'}
        mint = np.zeros(ds.diag, dtype=np.float64)
        maxt = np.zeros_like(mint)
        mint[:] = r.mint
        maxt[:] = r.maxt
        ds['mint'] = mint
        ds['mint'].attrs = {'description':'the mint attribut of the ray'}
        ds['maxt'] = maxt
        ds['maxt'].attrs = {'description':'the maxt attribut of the ray'}
    else: # multiple shape and rays, 2-D output
        nrays = t.shape[1]
        nobj = t.shape[0]
        ds.attrs.update({'nobj': nobj, 'nrays': nrays})
        ds['is_intersection'] = xr.DataArray(is_intersection, dims=['nobj', 'nrays'] )
        ds['is_intersection'].attrs = {'description':'tells if there is an intersection between the ray and the shape'}
        not_int = np.logical_not(is_intersection)
        ds['thit'] = xr.DataArray(t, dims=['nobj', 'nrays'])
        ds['thit'].attrs = {'description':'the t ray factor for the intersection point calculation'}
        u[not_int] = None
        v[not_int] = None
        ds['u'] = xr.DataArray(u, dims=['nobj', 'nrays'])
        ds['u'].attrs = {'description':'The u coordinate of the parametric representation'}
        ds['v'] = xr.DataArray(v, dims=['nobj', 'nrays'])
        ds['v'].attrs = {'description':'The v coordinate of the parametric representation'}

        mint = np.zeros(nrays, dtype=np.float64)
        maxt = np.zeros_like(mint)
        mint[:] = r.mint
        maxt[:] = r.maxt
        ds['mint'] = mint
        ds['mint'].attrs = {'description':'the mint attribut of the ray'}
        ds['maxt'] = maxt
        ds['maxt'].attrs = {'description':'the maxt attribut of the ray'}
        ds['o'] = xr.DataArray(r.o.to_numpy(), dims=['nrays', 'xyz'])
        ds['o'].attrs = {'type': 'Point', 'description':'the x, y and z components of the ray point'}
        ds['d'] = xr.DataArray(r.d.to_numpy(), dims=['nrays', 'xyz'])
        ds['d'].attrs = {'type': 'Vector', 'description':'the x, y and z components of the ray vector'}

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
        ds['phit'].attrs = {'type': 'Point', 'description':'the x, y and z components of the intersection point'}
        ds['nhit'] = xr.DataArray(n, dims=['nobj', 'nrays', 'xyz'])
        ds['nhit'].attrs = {'type': 'Normal', 'description':'the x, y and z components of the normal at the intersection point'}
        dpdu_2d = np.repeat(dpdu[:,np.newaxis,:], nrays, axis=1)
        dpdv_2d = np.repeat(dpdv[:,np.newaxis,:], nrays, axis=1)
        dpdu_2d[not_int,:] = None
        dpdv_2d[not_int,:] = None
        ds['dpdu'] = xr.DataArray(dpdu_2d, dims=['nobj', 'nrays', 'xyz'])
        ds['dpdu'].attrs = {'type': 'Vector', 'description':'The surface partial derivative of phit with respect to u'}
        ds['dpdv'] = xr.DataArray(dpdv_2d, dims=['nobj', 'nrays', 'xyz'])
        ds['dpdv'].attrs = {'type': 'Vector', 'description':'The surface partial derivative of phit with respect to v'}
        

    ds.attrs = {'shape': shape.__class__.__name__}
    return ds

