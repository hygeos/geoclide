#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geoclide.shapes import Shape, DifferentialGeometry
from geoclide.basic import Vector, Point, Ray
import geoclide.vecope as gv
import numpy as np
from geoclide.mathope import gamma_f64
from geoclide.transform import Transform


class Triangle(Shape):
    '''
    Creation of the class Triangle

    Parameters
    ----------
    p0 : Point
        The first point of the triangle
    p1 : Point
        The second point of the triangle
    p2 : Point
        The the third point of the triangle
    oTw : Transform, optional
        From object to world space or the transformation applied to the triangle
    wTo : Transform, optional
        From world to object space or the in inverse transformation applied to the triangle
    p0t : Point, optional
        If given circumvent the automatically computed p0t (p0 after applying transformation)
    p1t : Point, optional
        If given circumvent the automatically computed p1t (p1 after applying transformation)
    p2t : Point, optional
        If given circumvent the automatically computed p2t (p2 after applying transformation)
    '''
    def __init__(self, p0=None, p1=None, p2=None, oTw=None, wTo=None,
                 p0t=None, p1t=None, p2t=None):
        # Manage None cases
        if p0 is None : p0 = Point()
        if p1 is None : p1 = Point()
        if p2 is None : p2 = Point()
        if oTw is None and wTo is None:
            oTw = Transform()
            wTo = Transform()
            self.p0t = p0
            self.p1t = p1
            self.p2t = p2
        elif ( (oTw is None or isinstance(oTw, Transform)) and
               (wTo is None or isinstance(wTo, Transform)) ):
            if (oTw is None): oTw = wTo.inverse() # if oTw is None then wTo should be Transform
            if (wTo is None): wTo = oTw.inverse() # if wTo is None then oTw should be Transform
            if (p0t is None): self.p0t = oTw[p0]
            if (p1t is None): self.p1t = oTw[p1]
            if (p2t is None): self.p2t = oTw[p2]

        if (not isinstance(p0, Point) or not isinstance(p1, Point) or not isinstance(p2, Point)):
            raise ValueError('The parameters p0, p1 and p2 must be all Point')
        if ( (p0t is not None and not isinstance(p0t, Point)) or
             (p1t is not None and not isinstance(p1t, Point)) or
             (p2t is not None and not isinstance(p2t, Point)) ):
            raise ValueError('The parameters p0t, p1t and p2t must be all Point')
        Shape.__init__(self, ObjectToWorld = oTw, WorldToObject = wTo)
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        if (p0t is not None): self.p0t = p0t
        if (p1t is not None): self.p1t = p1t
        if (p2t is not None): self.p2t = p2t

    def is_intersection(self, r1, method='v3'):
        """
        Test if a Ray intersect with the triangle

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test
        method : str, optional
            Tow choice -> 'v2' (use mainly pbrt v2 intersection test method) or 'v3' (pbrt v3)

        Returns
        -------
        out : bool
            If there is an intersection -> True, else False
        """
        if method == 'v3':
            return self.is_intersection_v3(r1)
        elif method == 'v2':
            return self.is_intersection_v2(r1)
        else:
            raise ValueError("Only 'v2' and 'v3' are valid values for method parameter")   
    
    def intersect(self, r1, method='v3'):
        """
        Test if a Ray intersect with the triangle and return intersection information

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test
        method : str, optional
            Tow choice -> 'v2' (use mainly pbrt v2 intersection test method) or 'v3' (pbrt v3)
        
        Returns
        -------
        thit : float
            The t ray variable for its first intersection at the shape surface
        dg : DifferentialGeometry
            The parametric parameters at the intersection point
        is_intersection : bool
            If there is an intersection -> True, else False

        Notes
        -----
        By default the 'v3' method is used since there are more robustness tests.
        But the 'v2' method is at least twice faster than 'v3'.
        """
        if method == 'v3':
            return self.intersect_v3(r1)
        elif method == 'v2':
            return self.intersect_v2(r1)
        else:
            raise ValueError("Only 'v2' and 'v3' are valid values for method parameter")   
        
    def intersect_v2(self, r1):
        """
        Test if a Ray intersect with the triangle using mainly pbrt v2 method,
        and return intersection information
        
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
        """
        if not isinstance(r1, Ray): raise ValueError('The given parameter must be a Ray')
        ray = Ray(r1)
        p0 = self.p0t
        p1 = self.p1t
        p2 = self.p2t
        e1 = p1 - p0
        e2 = p2 - p0
        s1 = gv.cross(ray.d, e2)
        divisor = gv.dot(s1, e1)

        if (divisor == 0):
            return None, None, False
        invDivisor = 1./divisor

        # compute the first barycentric coordinate
        s = ray.o - p0
        b1 = gv.dot(s, s1) * invDivisor
        if (b1 < -0.00000001 or  b1 > 1.00000001):
            return None, None, False

        # compute the second barycentric coordinate
        s2 = gv.cross(s, e1)
        b2 = gv.dot(ray.d, s2) * invDivisor
        if (b2 < 0 or  b1+b2 > 1):
            return None, None, False

        # compute the time at the intersection point
        t = gv.dot(e2, s2) * invDivisor
        if (t < ray.mint or t > ray.maxt):
            return None, None, False

        # compute triangle partial derivatives
        uvs = np.array([[0., 0.], [1., 0.], [1., 1.]])

        # compute deltas for triangle partial derivatives
        du1 = uvs[0][0] - uvs[2][0]
        du2 = uvs[1][0] - uvs[2][0]
        dv1 = uvs[0][1] - uvs[2][1]
        dv2 = uvs[1][1] - uvs[2][1]
        dp1 = p0 - p2
        dp2 = p1 - p2
        determinant = du1 * dv2 - dv1 * du2

        if (determinant == 0):
            dpdu, dpdv = gv.coordinate_system(gv.normalize(gv.cross(e2, e1)))
        else:
            invdet = 1./determinant
            dpdu = ( dp1*dv2   - dp2*dv1) * invdet
            dpdv = (dp1*(-du2) + dp2*du1) * invdet
        
        # interpolate $(u,v)$ triangle parametric coordinates
        b0 = 1 - b1 - b2
        tu = b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0]
        tv = b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1]

        # fill the DifferentialGeometry and thit
        dg = DifferentialGeometry(ray[t], dpdu, dpdv, tu, tv, self)
        thit = t

        return thit, dg, True
    
    def intersect_v3(self, r1):
        """
        Test if a Ray intersect with the triangle using mainly pbrt v3 method,
        and return intersection information

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
        """
        if not isinstance(r1, Ray): raise ValueError('The given parameter must be a Ray')
        ray = Ray(r1)
        p0 = self.p0t
        p1 = self.p1t
        p2 = self.p2t

        # Get triangle vertices and translate them in based on ray origin
        p0t = p0 - ray.o
        p1t = p1 - ray.o
        p2t = p2 - ray.o

        kz = gv.vargmax(gv.vabs(ray.d))
        kx = kz + 1
        if(kx == 3): kx = 0
        ky = kx + 1
        if(ky == 3): ky = 0

        d = gv.permute(ray.d, kx, ky, kz)
        p0t = gv.permute(p0t, kx, ky, kz)
        p1t = gv.permute(p1t, kx, ky, kz)
        p2t = gv.permute(p2t, kx, ky, kz)
        
        sx = -d.x/d.z
        sy = -d.y/d.z
        sz = 1./d.z
        p0t.x += sx*p0t.z
        p0t.y += sy*p0t.z
        p1t.x += sx*p1t.z
        p1t.y += sy*p1t.z
        p2t.x += sx*p2t.z
        p2t.y += sy*p2t.z
        
        # Compute edge function coefficients
        e0 = (p1t.x * p2t.y) - (p1t.y * p2t.x)
        e1 = (p2t.x * p0t.y) - (p2t.y * p0t.x)
        e2 = (p0t.x * p1t.y) - (p0t.y * p1t.x)

        # Perform triangle edge and determinant tests
        if ((e0 < 0 or e1 < 0 or e2 < 0) and (e0 > 0 or e1 > 0 or e2 > 0)):
            return None, None, False
        det = e0 + e1 + e2
        if (det == 0): return None, None, False

        # Compute scaled hit distance to triangle and test against ray $t$ range
        p0t.z *=  sz
        p1t.z *=  sz
        p2t.z *=  sz

        tScaled = e0*p0t.z + e1*p1t.z + e2*p2t.z

        if ( (det < 0 and (tScaled >= 0 or tScaled < ray.maxt*det)) or
             (det > 0 and (tScaled <= 0 or tScaled > ray.maxt*det)) ):
            return None, None, False

        # Compute barycentric coordinates and t value for triangle intersection
        invDet = 1./det
        b0 = e0 * invDet
        b1 = e1 * invDet
        b2 = e2 * invDet
        t = tScaled * invDet
        
        # Ensure that computed triangle t is conservatively greater than zero
        maxZt = np.max(np.abs(np.array([p0t.z, p1t.z, p2t.z])))
        deltaZ = gamma_f64(3) * maxZt
        maxXt = np.max(np.abs(np.array([p0t.x, p1t.x, p2t.x])))
        maxYt = np.max(np.abs(np.array([p0t.y, p1t.y, p2t.y])))
        deltaX = gamma_f64(5) * (maxXt + maxZt)
        deltaY = gamma_f64(5) * (maxYt + maxZt)         
        deltaE = 2 * (gamma_f64(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt)
        maxE = np.max(np.abs(np.array([e0, e1, e2])))
        deltaT = 3 * (gamma_f64(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) * abs(invDet)
        if (t <= deltaT): return None, None, False

        # Compute triangle partial derivatives
        # Below the z components is not needed since we are in 2D with u in x and v un y
        dpdu = Vector()
        dpdv = Vector()
        uv0 = Point(0., 0., 0.)
        uv1 = Point(1., 0., 0.)
        uv2 = Point(1., 1., 0.)
        duv02 = uv0 - uv2
        duv12 = uv1 - uv2
        dp02 = p0 - p2
        dp12 = p1 - p2
        determinant = duv02.x*duv12.y - duv02.y*duv12.x
        degenerate = bool(abs(determinant) < 1e-8)

        if (not degenerate):
            invdet = 1./ determinant
            dpdu = (duv12.y*dp02 - duv02.y*dp12)*invdet
            dpdv = (-duv12.x*dp02 + duv02.x*dp12)*invdet

        if ( degenerate or gv.cross(dpdu, dpdv).length_squared() == 0):
            ng = gv.cross(p2-p0, p1-p0)
            if ( ng.length_squared() == 0 ):
                return None, None, False
            dpdu, dpdv = gv.coordinate_system(gv.normalize(ng))

        phit = b0*p0+b1*p1+b2*p2
        uvhit =b0*uv0 + b1*uv1 + b2*uv2
        thit = t
        dg = DifferentialGeometry(phit, dpdu, dpdv, uvhit.x, uvhit.y, self)
        
        return thit, dg, True
    
    def is_intersection_v2(self, r1):
        """
        Test if a Ray intersect with the triangle using mainly pbrt v2 method

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test
        
        Returns
        -------
        out : bool
            If there is an intersection -> True, else False
        """
        if not isinstance(r1, Ray): raise ValueError('The given parameter must be a Ray')
        ray = Ray(r1)
        p0 = self.p0t
        p1 = self.p1t
        p2 = self.p2t
        e1 = p1 - p0
        e2 = p2 - p0
        s1 = gv.cross(ray.d, e2)
        divisor = gv.dot(s1, e1)

        if (divisor == 0):
            return False
        invDivisor = 1./divisor

        # compute the first barycentric coordinate
        s = ray.o - p0
        b1 = gv.dot(s, s1) * invDivisor
        if (b1 < -0.00000001 or  b1 > 1.00000001):
            return False

        # compute the second barycentric coordinate
        s2 = gv.cross(s, e1)
        b2 = gv.dot(ray.d, s2) * invDivisor
        if (b2 < 0 or  b1+b2 > 1):
            return False

        # compute the time at the intersection point
        t = gv.dot(e2, s2) * invDivisor
        if (t < ray.mint or t > ray.maxt):
            return False

        return True
    
    def is_intersection_v3(self, r1):
        """
        Test if a Ray intersect with the triangle using mainly pbrt v3 method

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test
        
        Returns
        -------
        out : bool
            If there is an intersection -> True, else False
        """
        if not isinstance(r1, Ray): raise ValueError('The given parameter must be a Ray')
        ray = Ray(r1)
        p0 = self.p0t
        p1 = self.p1t
        p2 = self.p2t

        # Get triangle vertices and translate them in based on ray origin
        p0t = p0 - ray.o
        p1t = p1 - ray.o
        p2t = p2 - ray.o

        kz = gv.vargmax(gv.vabs(ray.d))
        kx = kz + 1
        if(kx == 3): kx = 0
        ky = kx + 1
        if(ky == 3): ky = 0

        d = gv.permute(ray.d, kx, ky, kz)
        p0t = gv.permute(p0t, kx, ky, kz)
        p1t = gv.permute(p1t, kx, ky, kz)
        p2t = gv.permute(p2t, kx, ky, kz)
        
        sx = -d.x/d.z
        sy = -d.y/d.z
        sz = 1./d.z
        p0t.x += sx*p0t.z
        p0t.y += sy*p0t.z
        p1t.x += sx*p1t.z
        p1t.y += sy*p1t.z
        p2t.x += sx*p2t.z
        p2t.y += sy*p2t.z
        
        # Compute edge function coefficients
        e0 = (p1t.x * p2t.y) - (p1t.y * p2t.x)
        e1 = (p2t.x * p0t.y) - (p2t.y * p0t.x)
        e2 = (p0t.x * p1t.y) - (p0t.y * p1t.x)

        # Perform triangle edge and determinant tests
        if ((e0 < 0 or e1 < 0 or e2 < 0) and (e0 > 0 or e1 > 0 or e2 > 0)):
            return False
        det = e0 + e1 + e2
        if (det == 0): return False

        # Compute scaled hit distance to triangle and test against ray $t$ range
        p0t.z *=  sz
        p1t.z *=  sz
        p2t.z *=  sz

        tScaled = e0*p0t.z + e1*p1t.z + e2*p2t.z

        if ( (det < 0 and (tScaled >= 0 or tScaled < ray.maxt*det)) or
             (det > 0 and (tScaled <= 0 or tScaled > ray.maxt*det)) ):
            return False

        # Compute barycentric coordinates and t value for triangle intersection
        invDet = 1./det
        t = tScaled * invDet
        
        # Ensure that computed triangle t is conservatively greater than zero
        maxZt = np.max(np.abs(np.array([p0t.z, p1t.z, p2t.z])))
        deltaZ = gamma_f64(3) * maxZt
        maxXt = np.max(np.abs(np.array([p0t.x, p1t.x, p2t.x])))
        maxYt = np.max(np.abs(np.array([p0t.y, p1t.y, p2t.y])))
        deltaX = gamma_f64(5) * (maxXt + maxZt)
        deltaY = gamma_f64(5) * (maxYt + maxZt)         
        deltaE = 2 * (gamma_f64(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt)
        maxE = np.max(np.abs(np.array([e0, e1, e2])))
        deltaT = 3 * (gamma_f64(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) * abs(invDet)
        if (t <= deltaT): return False

        # Compute triangle partial derivatives
        # Below the z components is not needed since we are in 2D with u in x and v un y
        dpdu = Vector()
        dpdv = Vector()
        uv0 = Point(0., 0., 0.)
        uv1 = Point(1., 0., 0.)
        uv2 = Point(1., 1., 0.)
        duv02 = uv0 - uv2
        duv12 = uv1 - uv2
        dp02 = p0 - p2
        dp12 = p1 - p2
        determinant = duv02.x*duv12.y - duv02.y*duv12.x
        degenerate = bool(abs(determinant) < 1e-8)

        if (not degenerate):
            invdet = 1./ determinant
            dpdu = (duv12.y*dp02 - duv02.y*dp12)*invdet
            dpdv = (-duv12.x*dp02 + duv02.x*dp12)*invdet

        if ( degenerate or gv.cross(dpdu, dpdv).length_squared() == 0):
            ng = gv.cross(p2-p0, p1-p0)
            if ( ng.length_squared() == 0 ):
                return False
        
        return True
     

class TriangleMesh(Shape):
    '''
    Creation of the class TriangleMesh

    Parameters
    ----------
    vi : np.ndarray
        The 1d ndarray of size (3*ntriangles) containing the vertices indices of triangles (see the parameter v). 
        The 3 first indices are the vertices (p0, p1 and p3) indices of the first triangle and so on. 
        It can be a 2d ndarray of shape (ntriangles, 3)
    v : np.ndarray
        The vertices xyz coordinates. It is a 1d array of size (nvertices) containing Point objects 
        where the first element is the coordinate of first vertex and so on. It can be a 2d float array 
        of size (nvertices, 3).
    oTw : Transform, optional
        From object to world space or the transformation applied to the triangle mesh
    wTo : Transform, optional
        From world to object space or the in inverse transformation applied to the triangle mesh
    '''
    def __init__(self, vi, v, oTw=None, wTo=None):
        if (  not isinstance(vi, np.ndarray)                                or
              not (len(vi.shape) == 1 or len(vi.shape) == 2)                or
              not (np.issubdtype(vi.dtype, int) or np.issubdtype(vi.dtype, np.integer))  ):
            raise ValueError('The parameter vi must be a 1d or 2d ndarray of intergers')
        if (  not ( isinstance(v, np.ndarray) )                                              or
              not ( (len(v.shape) == 1 and isinstance(v[0], Point)) or (len(v.shape) == 2) )  ):
            raise ValueError('The paramerter v must be a 1d ndarray of Point objects or a 2d ndarray')
        
        if (len(vi) == 2): vi = vi.flatten()
        self.vertices_index = vi
        self.nvertices = np.amax(vi) + int(1)
        if not isinstance(v[0], Point):
            v_bis = np.empty((self.nvertices), dtype=Point)
            for iv in range (0, self.nvertices):
                v_bis[iv] = Point(v[iv,:])
            self.vertices = v_bis
        else:
            self.vertices = v
        
        if oTw is None:
            self.vertices_t = self.vertices
        else:
            self.vertices_t = np.empty((self.nvertices), dtype=Point)
            for iv in range (0, self.nvertices):
                self.vertices_t[iv] = oTw[self.vertices[iv]]
        
        self.ntriangles = int(len(np.atleast_1d(self.vertices_index))/3)
        self.triangles = np.empty((self.ntriangles), dtype=Triangle)
            
        for itri in range(0, self.ntriangles):
            p0 = self.vertices[self.vertices_index[int(3*itri)]]
            p1 = self.vertices[self.vertices_index[int(3*itri) + int(1)]]
            p2 = self.vertices[self.vertices_index[int(3*itri) + int(2)]]
            p0t = self.vertices_t[self.vertices_index[int(3*itri)]]
            p1t = self.vertices_t[self.vertices_index[int(3*itri) + int(1)]]
            p2t = self.vertices_t[self.vertices_index[int(3*itri) + int(2)]]
            self.triangles[itri] = Triangle(p0, p1, p2, oTw, wTo, p0t, p1t, p2t)
        
        Shape.__init__(self, ObjectToWorld = self.triangles[0].oTw,
                       WorldToObject = self.triangles[0].wTo)
    
    def intersect(self, r1, method='v3'):
        """
        Test if a Ray intersect with the triangle mesh and return intersection information

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test
        method : str, optional
            Tow choice -> 'v2' (use mainly pbrt v2 triangle intersection test method) or 'v3' (pbrt v3)
        
        Returns
        -------
        thit : float
            The t ray variable for its first intersection at the shape surface
        dg : DifferentialGeometry
            The parametric parameters at the intersection point
        is_intersection : bool
            If there is an intersection -> True, else False
        """
        dg = None
        thit = float("inf")
        for itri in range(0, self.ntriangles):
            thit_bis, dg_bis, is_intersection_bis = self.triangles[itri].intersect(r1, method=method)
            if is_intersection_bis:
                if thit > thit_bis:
                    thit = thit_bis
                    dg = dg_bis
        if dg is None: return None, None, False
            
        return thit, dg, True
    
    def is_intersection(self, r1, method='v3'):
        """
        Test if a Ray intersect with the triangle mesh

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test
        method : str, optional
            Tow choice -> 'v2' (use mainly pbrt v2 triangle intersection test method) or 'v3' (pbrt v3)
        
        Returns
        -------
        out : bool
            If there is an intersection -> True, else False
        """
        for itri in range(0, self.ntriangles):
            if (self.triangles[itri].is_intersection(r1, method=method)):
                return True
        return False