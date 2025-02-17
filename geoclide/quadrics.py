#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide.shapes import Shape, DifferentialGeometry
from geoclide.mathope import clamp, quadratic
from geoclide.vecope import distance
from geoclide.basic import Ray, Vector, Point
from geoclide.transform import Transform, get_scale_tf
from geoclide.trianglemesh import create_sphere_trianglemesh, TriangleMesh
import math
import numpy as np
from geoclide.constante import TWO_PI


class Sphere(Shape):
    '''
    Creation of the class Sphere

    - without transformation the sphere is centered at the origin
    - z0, z1 and phi_max are needed parameters for the creation of any partial sphere

    Parameters
    ----------
    radius : float
        The radius of the sphere
    z_min : float, optional
        The minimum z value of the sphere where z0 is between [-radius, 0]
    z_max : float, optional
        The maximum z value of the sphere where z1 is between [0, radius]
    phi_max : float, optional
        The maximum phi value in degrees of the sphere, where phi is between 0 and 360°
    oTw : Transform, optional
        From object to world space or the transformation applied to the sphere
    wTo : Transform, optional
        From world to object space or the in inverse transformation applied to the sphere
    '''
    def __init__(self, radius, z_min=None, z_max=None, phi_max=360., oTw=None, wTo=None):
        if z_min is None: z_min = -radius
        if z_max is None: z_max = radius
        if z_max < z_min : raise ValueError ('zmax must be greater than zmin')
        if (phi_max < 0. or phi_max > 360.):
            raise ValueError ('The value of the parameter phi_max must in the range: [0, 360]')
        if wTo is None and oTw is None:
            wTo = Transform()
            oTw = Transform()
        elif (wTo is None and isinstance(oTw, Transform)): wTo = oTw.inverse()
        elif (isinstance(wTo, Transform) and oTw is None): oTw = wTo.inverse()
        if (not np.isscalar(radius) or
            not np.isscalar(z_min)  or
            not np.isscalar(z_max)  or
            not np.isscalar(phi_max) ):
            raise ValueError('The parameters radius, z_min, z_max and phi_max must be all scalars')
        Shape.__init__(self, ObjectToWorld = oTw, WorldToObject = wTo)
        self.radius = radius
        self.zmin = clamp(z_min, -self.radius, self.radius)
        self.zmax = clamp(z_max, -self.radius, self.radius)
        self.theta_min = math.acos(clamp(self.zmin/self.radius, -1, 1))
        self.theta_max = math.acos(clamp(self.zmax/self.radius, -1, 1))
        self.phi_max = phi_max

    def is_intersection_t(self, r1):
        """
        Test if a ray intersects the sphere / partial sphere

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test

        Returns
        -------
        thit : float
            The t ray variable for its first intersection at the shape surface
        is_intersection : bool
            If there is an intersection -> True, else False

        Examples
        --------
        >>> import geoclide as gc
        >>> sph1 = gc.Sphere(radius=1.) # sphere of radius 1
        >>> sph2 = gc.Sphere(radius=1., z_max=0.5) # partial sphere where portion above z=0.5 is removed
        >>> r = gc.Ray(o=gc.Point(-2., 0., 0.8), d=gc.Vector(1.,0.,0.))
        >>> sph1.is_intersection_t(r)
        (1.4000000000000004, True)
        >>> sph2.is_intersection(r) # here no intersection since the sphere part above z=0.5 is removed
        False
        """
        if not isinstance(r1, Ray): raise ValueError('The given parameter must be a Ray')
        ray = Ray(r1)
        ray.o = self.wTo[r1.o]
        ray.d = self.wTo[r1.d]

        # Compute quadratic sphere coefficients
        a = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z
        b = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z)
        c = ray.o.x*ray.o.x + ray.o.y*ray.o.y + ray.o.z*ray.o.z - \
            self.radius*self.radius

        # Solve quadratic equation
        exist, t0, t1 = quadratic(a, b, c)
        if (not exist): return None, False
        
        # Compute intersection distance along ray
        if (t0 > ray.maxt or t1 < ray.mint): return None, False
        thit = t0

        if (t0 < ray.mint):
            thit = t1
            if (thit > ray.maxt): return None, False

        # Compute sphere hit position and $\phi$
        phit = ray[thit]
        if (phit.x == 0 and phit.y == 0): phit.x = 1e-5 * self.radius
        phi = math.atan2(phit.y, phit.x)
        if (phi < 0): phi += TWO_PI

        # Test sphere intersection against clipping parameters
        phi_max_rad = math.radians(self.phi_max)
        if ((self.zmin > -self.radius and phit.z < self.zmin) or
            (self.zmax <  self.radius and phit.z > self.zmax) or
            (phi > phi_max_rad) ):
            if (thit == t1): return None, False
            if (t1 > ray.maxt): return None ,False
            thit = t1
            # Compute sphere hit position and $\phi$
            phit = ray[thit]
            if (phit.x == 0 and phit.y == 0): phit.x = 1e-5 * self.radius
            phi = math.atan2(phit.y, phit.x)
            if (phi < 0): phi += TWO_PI
            if ((self.zmin > -self.radius and phit.z < self.zmin) or
                (self.zmax <  self.radius and phit.z > self.zmax) or
                (phi > phi_max_rad) ):
                return None, False

        return thit, True
    
    def is_intersection(self, r1):
        """
        Test if a ray intersects the sphere / partial sphere

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test

        Returns
        -------
        out : bool
            If there is an intersection -> True, else False

        Examples
        --------
        >>> import geoclide as gc
        >>> sph1 = gc.Sphere(radius=1.) # sphere of radius 1
        >>> sph2 = gc.Sphere(radius=1., z_max=0.5) # partial sphere where portion above z=0.5 is removed
        >>> r = gc.Ray(o=gc.Point(-2., 0., 0.8), d=gc.Vector(1.,0.,0.))
        >>> sph1.is_intersection(r)
        True
        >>> sph2.is_intersection(r) # here no intersection since the sphere part above z=0.5 is removed
        False
        """
        _, is_intersection = self.is_intersection_t(r1)
        return is_intersection

    def intersect(self, r1):
        """
        Test if a ray intersects the sphere / partial sphere

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test

        Returns
        -------
        thit : float
            The t ray variable for its first intersection at the shape surface
        dg : DifferentialGeometry
            The parametric parameters at the intersection point
        is_intersection : bool
            If there is an intersection -> True, else False

        Examples
        --------
        >>> import geoclide as gc
        >>> sph1 = gc.Sphere(radius=1.) # sphere of radius 1
        >>> sph2 = gc.Sphere(radius=1., z_max=0.5) # partial sphere where portion above z=0.5 is removed
        >>> r = gc.Ray(o=gc.Point(-2., 0., 0.8), d=gc.Vector(1.,0.,0.))
        >>> sph1.intersect(r)
        (19.399999999999988,
        <geoclide.shapes.DifferentialGeometry at 0x7f589349cf40>,
        True)
        >>> thit, dg, is_int = sph1.intersect(r)
        >>> dg.p # the intersection point
        Point(-0.6000000000000121, 0.0, 0.8)
        >>> dg.n # The surface normal at the intersection point
        Normal(-0.6, 0.0, 0.8000000000000002)
        >>> sph2.intersect(r) # here no intersection since the sphere part above z=0.5 is removed
        (None, None, False)
        """
        if not isinstance(r1, Ray): raise ValueError('The given parameter must be a Ray')
        ray = Ray(r1)
        ray.o = self.wTo[r1.o]
        ray.d = self.wTo[r1.d]

        # Compute quadratic sphere coefficients
        a = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z
        b = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z)
        c = ray.o.x*ray.o.x + ray.o.y*ray.o.y + ray.o.z*ray.o.z - \
            self.radius*self.radius

        # Solve quadratic equation
        exist, t0, t1 = quadratic(a, b, c)
        if (not exist): return None, None, False
        
        # Compute intersection distance along ray
        if (t0 > ray.maxt or t1 < ray.mint): return None, None, False
        thit = t0

        if (t0 < ray.mint):
            thit = t1
            if (thit > ray.maxt): return None, None, False

        # Compute sphere hit position and $\phi$
        phit = ray[thit]
        phit *= self.radius / distance(phit, Point(0., 0., 0.))
        if (phit.x == 0 and phit.y == 0): phit.x = 1e-5 * self.radius
        phi = math.atan2(phit.y, phit.x)
        if (phi < 0): phi += TWO_PI

        # Test sphere intersection against clipping parameters
        phi_max_rad = math.radians(self.phi_max)
        if ((self.zmin > -self.radius and phit.z < self.zmin) or
            (self.zmax <  self.radius and phit.z > self.zmax) or
            (phi > phi_max_rad) ):
            if (thit == t1): return None, None, False
            if (t1 > ray.maxt): return None, None, False
            thit = t1
            # Compute sphere hit position and $\phi$
            phit = ray[thit]
            if (phit.x == 0 and phit.y == 0): phit.x = 1e-5 * self.radius
            phi = math.atan2(phit.y, phit.x)
            if (phi < 0): phi += TWO_PI
            if ((self.zmin > -self.radius and phit.z < self.zmin) or
                (self.zmax <  self.radius and phit.z > self.zmax) or
                (phi > phi_max_rad) ):
                return None, None, False

        # Find parametric representation of sphere hit
        u = phi / phi_max_rad
        theta = math.acos(clamp(phit.z / self.radius, -1, 1))
        v = (theta - self.theta_min) / (self.theta_max - self.theta_min)

        # Compute sphere $\dpdu$ and $\dpdv$
        zradius = math.sqrt(phit.x*phit.x + phit.y*phit.y)
        invzradius = 1 / zradius
        cosphi = phit.x * invzradius
        sinphi = phit.y * invzradius
        dpdu = Vector(-phi_max_rad * phit.y, phi_max_rad * phit.x, 0)
        dpdv = (self.theta_max-self.theta_min) * Vector(phit.z*cosphi, phit.z*sinphi, -self.radius*math.sin(theta)) 

        # Initialize _DifferentialGeometry_ from parametric information
        dg = DifferentialGeometry(self.oTw[phit], self.oTw[dpdu], self.oTw[dpdv],
                                  u, v, r1.d, self)

        return thit, dg, True

    def area(self):
        """
        compute the sphere / partial sphere area

        - !! the scale transform is not considered for the area calculation !!
        """
        return (math.radians(self.phi_max) * self.radius * (self.zmax-self.zmin)) # The sphere / partial sphere area
    
    def to_trianglemesh(self, reso_theta=None, reso_phi=None):
        """
        Convert the sphere to a triangle mesh

        Parameters
        ----------
        reso_theta : int, optional
            The number of lines around the polar theta angle, minimum accepted value is 3
        reso_phi : int, optional
            The number of lines around the azimuth phi angle, minimum accepted value is 3

        Returns
        -------
        mesh : TriangleMesh
            The sphere converted to a triangle mesh
        """
        theta_zmin = clamp(math.degrees(self.theta_min), 0., 360.)
        theta_zmax = clamp(math.degrees(self.theta_max), 0., 360.)
        theta_min = min(theta_zmin, theta_zmax)
        theta_max = max(theta_zmin, theta_zmax)
        return create_sphere_trianglemesh(self.radius, reso_theta, reso_phi, theta_min, theta_max,
                                          self.phi_max, self.oTw, self.wTo)

        

class Spheroid(Shape):
    '''
    Creation of the class Spheroid

    - without transformation the spheroid is centered at the origin
    - spheroid equation: x/(alpha**2) + y/(alpha**2) + z/(gamma**2) = 1,
      where alpha = radius_xy and gamma = radius_z
    - prolate -> radius_z > radius_xy
    - oblate -> radius_z < radius_xy

    Parameters
    ----------
    radius_xy : float
        The equatorial radius of the spheroid
    radius_z : float
        The pole radius of the spheroid (distance from center to pole along z axis)
    oTw : Transform, optional
        From object to world space or the transformation applied to the spheroid
    wTo : Transform, optional
        From world to object space or the in inverse transformation applied to the spheroid
    '''
    def __init__(self, radius_xy, radius_z, oTw=None, wTo=None):
        if wTo is None and oTw is None:
            wTo = Transform()
            oTw = Transform()
        elif (wTo is None and isinstance(oTw, Transform)): wTo = oTw.inverse()
        elif (isinstance(wTo, Transform) and oTw is None): oTw = wTo.inverse()
        if (not np.isscalar(radius_xy) or not np.isscalar(radius_z)):
            raise ValueError('The parameters alpha and gamma must be all scalars')
        Shape.__init__(self, ObjectToWorld = oTw, WorldToObject = wTo)
        self.alpha = radius_xy
        self.gamma = radius_z
        self.alpha2 = radius_xy*radius_xy
        self.gamma2 = radius_z*radius_z

    def is_intersection_t(self, r1):
        """
        Test if a ray intersects the spheroid

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test

        Returns
        -------
        thit : float
            The t ray variable for its first intersection at the shape surface
        is_intersection : bool
            If there is an intersection -> True, else False

        Examples
        --------
        >>> import geoclide as gc
        >>> oblate = gc.Spheroid(radius_xy=3., radius_z=1.5)
        >>> prolate = gc.Spheroid(radius_xy=1.5, radius_z=3.)
        >>> r1 = gc.Ray(o=gc.Point(2.5, 0., 10.), d=(gc.Vector(0., 0., -1.)))
        >>> r2 = gc.Ray(o=gc.Point(10., 0., 2.5), d=(gc.Vector(-1., 0., 0.)))
        >>> oblate.is_intersection_t(r1)
        (9.170843802411135, True)
        >>> oblate.is_intersection_t(r2)
        (None, False)
        >>> prolate.is_intersection_t(r1)
        (None, False)
        >>> prolate.is_intersection_t(r2)
        (9.170843802411135, True)
        """
        if not isinstance(r1, Ray): raise ValueError('The given parameter must be a Ray')
        ray = Ray(r1)
        ray.o = self.wTo[r1.o]
        ray.d = self.wTo[r1.d]

        # Compute quadratic sphere coefficients
        inv_alpha2 = 1./self.alpha2
        inv_beta2 = inv_alpha2 # ellipsoid special case where alpha=beta
        inv_gamma2 = 1./self.gamma2
        a = ray.d.x*ray.d.x*inv_alpha2 + ray.d.y*ray.d.y*inv_beta2 + ray.d.z*ray.d.z*inv_gamma2
        b = 2 * (ray.d.x*ray.o.x*inv_alpha2 + ray.d.y*ray.o.y*inv_beta2 + ray.d.z*ray.o.z*inv_gamma2)
        c = ray.o.x*ray.o.x*inv_alpha2 + ray.o.y*ray.o.y*inv_beta2 + ray.o.z*ray.o.z*inv_gamma2 - 1

        # Solve quadratic equation
        exist, t0, t1 = quadratic(a, b, c)
        if (not exist): return None, False
        
        # Compute intersection distance along ray
        if (t0 > ray.maxt or t1 < ray.mint): return None, False
        thit = t0

        if (t0 < ray.mint):
            thit = t1
            if (thit > ray.maxt): return None, False

        return thit, True
    
    def is_intersection(self, r1):
        """
        Test if a ray intersects the spheroid

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test

        Returns
        -------
        out : bool
            If there is an intersection -> True, else False

        Examples
        --------
        >>> import geoclide as gc
        >>> oblate = gc.Spheroid(radius_xy=3., radius_z=1.5)
        >>> prolate = gc.Spheroid(radius_xy=1.5, radius_z=3.)
        >>> r1 = gc.Ray(o=gc.Point(2.5, 0., 10.), d=(gc.Vector(0., 0., -1.)))
        >>> r2 = gc.Ray(o=gc.Point(10., 0., 2.5), d=(gc.Vector(-1., 0., 0.)))
        >>> oblate.is_intersection(r1)
        True
        >>> oblate.is_intersection(r2)
        False
        >>> prolate.is_intersection(r1)
        False
        >>> prolate.is_intersection(r2)
        True
        """
        _, is_intersection = self.is_intersection_t(r1)
        return is_intersection

    def intersect(self, r1):
        """
        Test if a ray intersects the spheroid

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test

        Returns
        -------
        thit : float
            The t ray variable for its first intersection at the shape surface
        dg : DifferentialGeometry
            The parametric parameters at the intersection point
        is_intersection : bool
            If there is an intersection -> True, else False

        Examples
        --------
        >>> import geoclide as gc
        >>> oblate = gc.Spheroid(radius_xy=3., radius_z=1.5)
        >>> prolate = gc.Spheroid(radius_xy=1.5, radius_z=3.)
        >>> r1 = gc.Ray(o=gc.Point(2.5, 0., 10.), d=(gc.Vector(0., 0., -1.)))
        >>> r2 = gc.Ray(o=gc.Point(10., 0., 2.5), d=(gc.Vector(-1., 0., 0.)))
        >>> oblate.intersect(r1)
        (9.170843802411135,
        <geoclide.shapes.DifferentialGeometry at 0x7f9d1e0c5810>,
        True)
        >>> oblate.intersect(r2)
        (None, None, False)
        >>> prolate.intersect(r1)
        (None, None, False)
        >>> prolate.intersect(r2)
        (9.170843802411135,
        <geoclide.shapes.DifferentialGeometry at 0x7f9d1e0c5450>,
        True)
        """
        if not isinstance(r1, Ray): raise ValueError('The given parameter must be a Ray')
        ray = Ray(r1)
        ray.o = self.wTo[r1.o]
        ray.d = self.wTo[r1.d]

        # Compute quadratic sphere coefficients
        inv_alpha2 = 1./self.alpha2
        inv_beta2 = inv_alpha2 # ellipsoid special case where alpha=beta
        inv_gamma2 = 1./self.gamma2
        a = ray.d.x*ray.d.x*inv_alpha2 + ray.d.y*ray.d.y*inv_beta2 + ray.d.z*ray.d.z*inv_gamma2
        b = 2 * (ray.d.x*ray.o.x*inv_alpha2 + ray.d.y*ray.o.y*inv_beta2 + ray.d.z*ray.o.z*inv_gamma2)
        c = ray.o.x*ray.o.x*inv_alpha2 + ray.o.y*ray.o.y*inv_beta2 + ray.o.z*ray.o.z*inv_gamma2 - 1 

        # Solve quadratic equation
        exist, t0, t1 = quadratic(a, b, c)
        if (not exist): return None, None, False
        
        # Compute intersection distance along ray
        if (t0 > ray.maxt or t1 < ray.mint): return None, None, False
        thit = t0

        if (t0 < ray.mint):
            thit = t1
            if (thit > ray.maxt): return None, None, False

        # Compute sphere hit position and $\phi$
        phit = ray[thit]
        if (phit.x == 0 and phit.y == 0): phit.x = 1e-5 * self.alpha
        phi = math.atan2(phit.y, phit.x) # because alpha=beta
        if (phi < 0): phi += TWO_PI

        # Find parametric representation of sphere hit
        u = phi / TWO_PI
        theta = math.acos(clamp(phit.z / self.gamma, -1, 1))
        v = 1 - (theta / math.pi)

        # Compute sphere dpdu and dpdv
        zradius = math.sqrt(phit.x*phit.x + phit.y*phit.y)
        invzradius = 1 / zradius
        cosphi = phit.x * invzradius
        sinphi = phit.y * invzradius
        fac = -math.pi*(self.alpha/self.gamma)*phit.z
        dpdu = Vector(-TWO_PI*phit.y, TWO_PI*phit.x, 0.)
        dpdv = Vector(fac*cosphi, fac*sinphi, math.pi*self.gamma*math.sin(theta))

        # Initialize _DifferentialGeometry_ from parametric information
        dg = DifferentialGeometry(self.oTw[phit], self.oTw[dpdu], self.oTw[dpdv],
                                  u, v, r1.d, self)

        return thit, dg, True
    
    def area(self):
        """
        compute the spheroid area

        - !! the scale transform is not considered for the area calculation !!
        """
        if (self.gamma < self.alpha): # oblate spheroid
            e = math.sqrt(1 - (self.gamma2/self.alpha2))
            area = TWO_PI*self.alpha2 + math.pi*(self.gamma2/e)*math.log((1+e)/(1-e))
        elif (self.gamma > self.alpha): # prolate
            e = math.sqrt(1 - (self.alpha2/self.gamma2))
            area = TWO_PI*self.alpha2*(1 + (self.gamma/(self.alpha*e))*math.asin(e))
        else: # sphere
            area = 4.*math.pi*self.alpha2
        return area
    
    def to_trianglemesh(self, reso_theta=None, reso_phi=None):
        """
        Convert the spheroid to a triangle mesh

        Parameters
        ----------
        reso_theta : int, optional
            The number of lines around the polar theta angle, minimum accepted value is 3
        reso_phi : int, optional
            The number of lines around the azimuth phi angle, minimum accepted value is 3

        Returns
        -------
        mesh : TriangleMesh
            The sphere converted to a triangle mesh
        """
        rescale_xyz = get_scale_tf(self.alpha, self.alpha, self.gamma)
        msh = create_sphere_trianglemesh(radius=1, reso_theta=reso_theta, reso_phi=reso_phi)
        vertices_t = np.zeros((msh.nvertices,3))
        for iver in range (0, msh.nvertices):
            vertices_t[iver,:] = rescale_xyz[Point(msh.vertices[iver,:])].to_numpy()
        return TriangleMesh(vertices_t, msh.faces, oTw=self.oTw, wTo=self.wTo)
    

class Disk(Shape):
    '''
    Creation of the class Disk

    Parameters
    ----------
    radius : float
        The disk radius
    inner_radius : float, optional
        This is the inner radius in the case of partial disk (or annulus)
    phi_max : float, optional
        The maximum phi value in degrees of the disk/annulus, where phi is between 0 and 360°
    z_height : float, optional
        the disk height along the z axis
    oTw : Transform, optional
        From object to world space or the transformation applied to the spheroid
    wTo : Transform, optional
        From world to object space or the in inverse transformation applied to the spheroid
    
    Notes
    -----
    Even if z_height is given, the origin for rotation transformation do not change.
    For exemple: z_height=5 then we apply a rotation of 90 degrees around the y axis, the disk
    we be rotated from (0.,0.,0.), meaning the disk we be moved from position (0.,0.,5.) to
    (5.,0.,0.).
    '''
    def __init__(self, radius, inner_radius=0., phi_max=360., z_height=0., oTw=None, wTo=None):
        if wTo is None and oTw is None:
            wTo = Transform()
            oTw = Transform()
        if (phi_max < 0. or phi_max > 360.):
            raise ValueError ('The value of the parameter phi_max must in the range: [0, 360]')
        elif (wTo is None and isinstance(oTw, Transform)): wTo = oTw.inverse()
        elif (isinstance(wTo, Transform) and oTw is None): oTw = wTo.inverse()
        if (inner_radius >= radius): raise NameError ('The parameter inner_radius must be < to radius')
        if (not np.isscalar(radius)): raise ValueError('The parameters radius must be a scalar')
        if (not np.isscalar(inner_radius)): raise ValueError('The parameters inner_radius must be a scalar')
        if (not np.isscalar(phi_max)): raise ValueError('The parameters phi_max must be a scalar')
        if (not np.isscalar(z_height)): raise ValueError('The parameters z_height must be a scalar')
        Shape.__init__(self, ObjectToWorld = oTw, WorldToObject = wTo)
        self.radius = radius
        self.inner_radius = inner_radius
        self.phi_max = phi_max
        self.z_height = z_height

    def is_intersection_t(self, r1):
        """
        Test if a ray intersects the disk

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test

        Returns
        -------
        thit : float
            The t ray variable for its first intersection at the shape surface
        is_intersection : bool
            If there is an intersection -> True, else False

        Examples
        --------
        >>> import geoclide as gc
        >>> r1 = gc.Ray(gc.Point(1.2,0.,10.), gc.Vector(0.,0.,-1.))
        >>> r2 = gc.Ray(gc.Point(0.2,0.,10.), gc.Vector(0.,0.,-1.))
        >>> r3 = gc.Ray(gc.Point(1.6,0.,10.), gc.Vector(0.,0.,-1.))
        >>> annulus = gc.Disk(radius=1.5, inner_radius=0.8)
        >>> annulus.is_intersection_t(r1) # hit point is between the inner radius and radius
        (10.0, True)
        >>> annulus.is_intersection_t(r2) # the ray passes through the annulus hole, no intersection
        (None, False)
        >>> annulus.is_intersection_t(r3) # the ray passes outside, no intersection
        (None, False)
        """
        if not isinstance(r1, Ray): raise ValueError('The given parameter must be a Ray')
        ray = Ray(r1)
        ray.o = self.wTo[r1.o]
        ray.d = self.wTo[r1.d]

        # no intersection in the case the ray is parallel to the disk's plane
        if (ray.d.z == 0): return None, False
        thit = (self.z_height - ray.o.z) / ray.d.z
        if (thit <= 0 or thit >= ray.maxt): return None, False

        # get the intersection point, and distance between disk center and the intersection
        phit = ray[thit]
        hit_radius2 = phit.x*phit.x + phit.y*phit.y

        # if the hit point is outside the disk then no intersection
        if (hit_radius2 > self.radius*self.radius): return None, False

        # annulus case
        # check that the hit point is not in the annulus hole
        if (self.inner_radius > 0.):
            if (hit_radius2 < self.inner_radius*self.inner_radius): return None, False
        
        # partial disk/annulus case
        # check phi value to see if the hit point is inside the partial disk/annulus
        if (self.phi_max < 360.):
            phi = math.atan2(phit.y, phit.x)
            if (phi < 0.): phi += TWO_PI
            if (phi > math.radians(self.phi_max)): return None, False
        
        return thit, True
    
    def is_intersection(self, r1):
        """
        Test if a ray intersects the disk

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test

        Returns
        -------
        is_intersection : bool
            If there is an intersection -> True, else False

        Examples
        --------
        >>> import geoclide as gc
        >>> r1 = gc.Ray(gc.Point(1.2,0.,10.), gc.Vector(0.,0.,-1.))
        >>> r2 = gc.Ray(gc.Point(0.2,0.,10.), gc.Vector(0.,0.,-1.))
        >>> r3 = gc.Ray(gc.Point(1.6,0.,10.), gc.Vector(0.,0.,-1.))
        >>> annulus = gc.Disk(radius=1.5, inner_radius=0.8)
        >>> annulus.is_intersection(r1) # hit point is between the inner radius and radius
        True
        >>> annulus.is_intersection(r2) # the ray passes through the annulus hole, no intersection
        False
        >>> annulus.is_intersection(r3) # the ray passes outside, no intersection
        False
        """
        _, is_intersection = self.is_intersection_t(r1)
        return is_intersection

    def intersect(self, r1):
        """
        Test if a ray intersects the disk

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test

        Returns
        -------
        thit : float
            The t ray variable for its first intersection at the shape surface
        dg : DifferentialGeometry
            The parametric parameters at the intersection point
        is_intersection : bool
            If there is an intersection -> True, else False

        Examples
        --------
        >>> import geoclide as gc
        >>> r1 = gc.Ray(gc.Point(1.2,0.,10.), gc.Vector(0.,0.,-1.))
        >>> r2 = gc.Ray(gc.Point(0.2,0.,10.), gc.Vector(0.,0.,-1.))
        >>> r3 = gc.Ray(gc.Point(1.6,0.,10.), gc.Vector(0.,0.,-1.))
        >>> annulus = gc.Disk(radius=1.5, inner_radius=0.8)
        >>> annulus.intersect(r1) # hit point is between the inner radius and radius
        >>> (10.0, <geoclide.shapes.DifferentialGeometry at 0x7f988d9cfce0>, True)
        >>> annulus.intersect(r2) # the ray passes through the annulus hole, no intersection
        (None, None, False)
        >>> annulus.intersect(r3) # the ray passes outside, no intersection
        (None, None, False)
        """
        if not isinstance(r1, Ray): raise ValueError('The given parameter must be a Ray')
        ray = Ray(r1)
        ray.o = self.wTo[r1.o]
        ray.d = self.wTo[r1.d]

        # no intersection in the case the ray is parallel to the disk's plane
        if (ray.d.z == 0): return None, None, False
        thit = (self.z_height - ray.o.z) / ray.d.z
        if (thit <= 0 or thit >= ray.maxt): return None, None, False

        # get the intersection point, and distance between disk center and the intersection
        phit = ray[thit]
        hit_radius2 = phit.x*phit.x + phit.y*phit.y

        # if the hit point is outside the disk then no intersection
        if (hit_radius2 > self.radius*self.radius): return None, None, False

        # annulus case
        # check that the hit point is not in the annulus hole
        if (hit_radius2 < self.inner_radius*self.inner_radius): return None, None, False
        
        # partial disk/annulus case
        # check phi value to see if the hit point is inside the partial disk/annulus
        phi = math.atan2(phit.y, phit.x)
        phi_max_rad = math.radians(self.phi_max)
        if (phi < 0.): phi += TWO_PI
        if (phi > phi_max_rad): return None, None, False

        # get the parameteric representation
        u = phi / phi_max_rad
        hit_radius = math.sqrt(hit_radius2)
        v = (self.radius-hit_radius) / (self.radius-self.inner_radius)
        dpdu = Vector(-phi_max_rad*phit.y, phi_max_rad*phit.x, 0.)
        dpdv = Vector(phit.x, phit.y, 0.) * ( (self.inner_radius-self.radius)/hit_radius )

        # Initialize _DifferentialGeometry_ from parametric information
        dg = DifferentialGeometry(self.oTw[phit], self.oTw[dpdu], self.oTw[dpdv],
                                  u, v, r1.d, self)
        
        return thit, dg, True
    
    def area(self):
        """
        compute the disk / annulus area

        - !! the scale transform is not considered for the area calculation !!
        """
        return 0.5*math.radians(self.phi_max)*(self.radius*self.radius - self.inner_radius*self.inner_radius)

