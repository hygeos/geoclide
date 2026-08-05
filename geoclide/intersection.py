"""
Standardized ray-shape intersection tests.

This module provides the calc_intersection function, which
performs the intersection test between any geoclide shape (BBox,
Sphere, Spheroid, Disk, Triangle or TriangleMesh) and a ray or a
set of rays, and returns an xarray dataset gathering the
intersection results together with the shape attributes.
"""

from __future__ import annotations

from typing import cast

import xarray as xr

from geoclide.basic import BBox, Ray
from geoclide.quadrics import Disk, Sphere, Spheroid
from geoclide.shapes import Shape
from geoclide.trianglemesh import Triangle, TriangleMesh


def calc_intersection(
    shape: BBox | Sphere | Spheroid | Disk | Triangle | TriangleMesh,
    r: Ray,
    **kwargs,
) -> xr.Dataset:
    """
    Performs intersection test between a shape and a ray and returns
    dataset

    Parameters
    ----------
    shape : BBox, Sphere, Spheroid, Disk, Triangle or TriangleMesh
        The shape used for the intersection(s)
    r : Ray
        The ray(s) used for the intersection(s)
    **kwargs
        The keyword arguments are passed on to intersect method. The
        ds_output parameter is forced here to always be True.

    Returns
    -------
    Dataset
        Xarray dataset containing the intersection information.

        Key variables included:

        - **o**: The origin(s) of the ray(s) [xyz]
        - **d**: The direction(s) of the ray(s) [xyz]
        - **mint**: The mint attribute of the ray(s)
        - **maxt**: The maxt attribute of the ray(s)
        - **is_intersection**: If there is an intersection ->
          True, else False
        - **thit**: The t ray variable(s) of the intersection
          point(s)
        - **phit**: The intersection point(s) [xyz]
        - **nhit**: The surface normal(s) at the intersection
          point(s) [xyz] (not for a BBox)
        - **u**, **v**, **dpdu**, **dpdv**: The parametric
          coordinates and surface partial derivatives (not for a
          BBox)
        - the shape attributes (e.g. radius, z_min, z_max and
          phi_max for a sphere, or pmin and pmax for a bounding
          box)
        - **wto_m**, **wto_m_inv**, **otw_m**, **otw_m_inv**: The
          transformation matrices of the shape (not for a BBox)

    Examples
    --------
    >>> import geoclide as gc
    >>> sphere = gc.Sphere(radius=1.) # sphere of radius 1
    >>> bbox = gc.BBox(p1=gc.Point(0., 0., 0.), p2=gc.Point(1.,1.,1.))
    >>> ray = gc.Ray(o=gc.Point(-2., 0., 0.8), d=gc.Vector(1.,0.,0.))
    >>> ds_sphere = gc.calc_intersection(sphere, ray)
    >>> ds_sphere['thit'].values
    array(1.4)
    >>> ds_sphere['phit'].values
    array([-0.6,  0. ,  0.8])
    >>> ds_bbox = gc.calc_intersection(bbox, ray)
    >>> ds_bbox['thit'].values
    array(2.)
    >>> ds_bbox['phit'].values
    array([0. , 0. , 0.8])
    """
    if not isinstance(r, Ray):
        raise ValueError("The parameter r1 must a Ray")

    if "ds_output" in kwargs:
        kwargs.pop("ds_output", False)
    if (isinstance(shape, BBox)) or issubclass(shape.__class__, Shape):
        ds = cast(
            xr.Dataset, shape.intersect(r, ds_output=True, **kwargs)
        )
    else:
        raise ValueError(
            "The only supported shape are: BBox, Sphere, Spheroid, "
            "Disk, Triangle and TriangleMesh"
        )

    if isinstance(shape, BBox):
        ds["pmin"] = xr.DataArray(shape.pmin.to_numpy(), dims="xyz")
        ds["pmin"].attrs = {
            "type": "Point",
            "description": "the x, y and z components of the pmin BBox "
            "attribut",
        }
        ds["pmax"] = xr.DataArray(shape.pmax.to_numpy(), dims="xyz")
        ds["pmax"].attrs = {
            "type": "Point",
            "description": "the x, y and z components of the pmax BBox "
            "attribut",
        }
    if isinstance(shape, Sphere):
        ds["radius"] = shape.radius
        ds["radius"].attrs = {"description": "the sphere radius attribut"}
        ds["z_min"] = shape.zmin
        ds["z_min"].attrs = {"description": "the sphere zmin attribut"}
        ds["z_max"] = shape.zmax
        ds["z_max"].attrs = {"description": "the sphere zmax attribut"}
        ds["phi_max"] = shape.phi_max
        ds["phi_max"].attrs = {
            "unit": "Degree",
            "description": "the sphere phi_max attribut",
        }
    if isinstance(shape, Spheroid):
        ds["radius_xy"] = shape.alpha
        ds["radius_xy"].attrs = {
            "description": "the equatorial radius of the spheroid "
            "(alpha attribut)"
        }
        ds["radius_z"] = shape.gamma
        ds["radius_z"].attrs = {
            "description": "the distance between the spheroid center "
            "and pole (gamma attribut)"
        }
    if isinstance(shape, Disk):
        ds["radius"] = shape.radius
        ds["radius"].attrs = {"description": "the radius of the disk"}
        ds["inner_radius"] = shape.inner_radius
        ds["inner_radius"].attrs = {
            "description": "the inner radius of the disk "
            "(if > 0 -> annulus case)"
        }
        ds["phi_max"] = shape.phi_max
        ds["phi_max"].attrs = {
            "unit": "Degree",
            "description": "the disk phi_max attribut",
        }
        ds["z_height"] = shape.z_height
        ds["z_height"].attrs = {"description": "the disk z_height attribut"}
    if isinstance(shape, Triangle):
        ds["p0"] = xr.DataArray(shape.p0.to_numpy(), dims="xyz")
        ds["p0"].attrs = {"description": "the triangle p0 attribut"}
        ds["p1"] = xr.DataArray(shape.p1.to_numpy(), dims="xyz")
        ds["p1"].attrs = {"description": "the triangle p1 attribut"}
        ds["p2"] = xr.DataArray(shape.p2.to_numpy(), dims="xyz")
        ds["p2"].attrs = {"description": "the triangle p2 attribut"}
    if isinstance(shape, TriangleMesh):
        ds["vertices"] = xr.DataArray(
            shape.vertices, dims=["nvertices", "xyz"]
        )
        ds["vertices"].attrs = {"description": "The vertices xyz coordinates."}
        ds["faces"] = xr.DataArray(shape.faces, dims=["ntriangles", "p0p1p2"])
        ds["faces"].attrs = {
            "description": "For each triangle, the index of vertices "
            "point p0, p1 and p2 (from variable v)."
        }
        ds.attrs.update(
            {
                "ntriangles": shape.ntriangles,
                "nvertices": shape.nvertices,
            }
        )
    if not isinstance(shape, BBox):
        shape_name = str(ds.attrs["shape"]).lower()
        ds["wto_m"] = xr.DataArray(shape.wto.m)
        ds["wto_m"].attrs = {
            "description": "the transformation matrix of the "
            + shape_name
            + " wto attribut"
        }
        ds["wto_m_inv"] = xr.DataArray(shape.wto.m_inv)
        ds["wto_m_inv"].attrs = {
            "description": "the inverse transformation matrix of the "
            + shape_name
            + " wto attribut"
        }
        ds["otw_m"] = xr.DataArray(shape.otw.m)
        ds["otw_m"].attrs = {
            "description": "the transformation matrix of the "
            + shape_name
            + " otw attribut"
        }
        ds["otw_m_inv"] = xr.DataArray(shape.otw.m_inv)
        ds["otw_m_inv"].attrs = {
            "description": "the inverse transformation matrix of the "
            + shape_name
            + " otw attribut"
        }

    return ds
