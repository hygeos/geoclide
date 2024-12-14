#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math 
import numpy as np


class Vector(object):
    """
    Parameters
    ----------
    x : float | Point | Vector | Normal | np.ndarray
        If scalar -> x component of the vector.
        Else, circumvent the y and z parameters and take the components of the Point/Vector/Normal/np.ndarray.
    y : float
        The y component of the vector.
    z : float
        The z component of the vector.

    Exemples
    --------
    >>> v1 = Vector(0.,0.,1.)
    >>> v1
    >>> Vector(0,0,1)
    """
    def __init__(self, x = 0., y = 0., z = 0.):
        if ( np.isscalar(x) and
             np.isscalar(y) and
             np.isscalar(z) ):
            self.x = float(x)
            self.y = float(y)
            self.z = float(z)
        elif ( isinstance(x, Vector) or
               isinstance(x, Point) or
               isinstance(x, Normal) ):
            self.x = x.x
            self.y = x.y
            self.z = x.z
        elif ( isinstance(x, np.ndarray) and
               len(x) == 3 ):
            self.x = float(x[0])
            self.y = float(x[1])
            self.z = float(x[2])
        else: 
            raise ValueError("Wrong parameter value(s) for Vector")

    def __eq__(self, v2):
        if isinstance(v2, Vector):
            return (self.x==v2.x) and (self.y==v2.y) and (self.z==v2.z)
        else:
            raise NameError('Equality with a Vector must be only with another Vector')

    def __add__(self, v2):
        if isinstance(v2, Vector):
            return Vector(self.x+v2.x, self.y+v2.y, self.z+v2.z) 
        else:
            raise NameError('Addition with a Vector must be only with another Vector')

    def __sub__(self, v2):
        if isinstance(v2, Vector):
            return Vector(self.x-v2.x, self.y-v2.y, self.z-v2.z)
        else:
            raise NameError('Substraction with a Vector must be only with another Vector')

    def __truediv__(self, sca):
        if (np.isscalar(sca)):
            return Vector(self.x/sca, self.y/sca, self.z/sca) 
        else:
            raise NameError('A Vector can be divided only by a scalar')
    def __mul__(self, sca): 
        if (np.isscalar(sca)):
            return Vector(sca*self.x, sca*self.y, sca*self.z)
        else:
            raise NameError('A Vector can be multiplied only by a scalar')

    def __getitem__(self, ind):
        if ( not isinstance(ind, int) or
             not isinstance(ind, np.integer) ):
            IndexError("Only an integer is a valid index")
        if ind == 0:
            return self.x
        elif ind == 1:
            return self.y
        elif ind == 2 :
            return self.z
        else:
            IndexError(f"Index {ind} is out of range") 

    def __str__(self):
        return '(' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ')'
    
    def __repr__(self):
        return 'Vector(' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ')'

    def length(self):
        return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z) # L2 norm

    def to_numpy(self):
        return np.array([self.x, self.y, self.z], dtype=np.float64)
    

class Point(object):
    """
    Parameters
    ----------
    x : float | Point | Vector | Normal | np.ndarray
        If scalar -> x component of the point.
        Else, circumvent the y and z parameters and take the components of the Point/Vector/Normal/np.ndarray.
    y : float
        The y component of the point.
    z : float
        The z component of the point.

    Exemples
    --------
    >>> p1 = Point(0.,0.,1.)
    >>> p1
    >>> Point(0,0,1)
    """
    def __init__(self, x = 0., y = 0., z = 0.):
        if ( np.isscalar(x) and
             np.isscalar(y) and
             np.isscalar(z) ):
            self.x = float(x)
            self.y = float(y)
            self.z = float(z)
        elif ( isinstance(x, Vector) or
               isinstance(x, Point) or
               isinstance(x, Normal) ):
            self.x = x.x
            self.y = x.y
            self.z = x.z
        elif ( isinstance(x, np.ndarray) and
               len(x) == 3 ):
            self.x = float(x[0])
            self.y = float(x[1])
            self.z = float(x[2])
        else: 
            raise ValueError("Wrong parameter value(s) for Point")

    def __eq__(self, p2):
        if isinstance(p2, Point):
            return (self.x==p2.x) and (self.y==p2.y) and (self.z==p2.z)
        else:
            raise NameError('Equality with a Point must be only with another Point')

    def __add__(self, v):
        if isinstance(v, Vector):
            return Point(self.x+v.x, self.y+v.y, self.z+v.z)
        else:
            raise NameError('Addition with a Point must be only with a Vector')

    def __sub__(self, vp2):
        if isinstance(vp2, Vector):
            return Point(self.x-vp2.x, self.y-vp2.y, self.z-vp2.z)
        elif isinstance(vp2, Point):
            return Vector(self.x-vp2.x, self.y-vp2.y, self.z-vp2.z)
        else:
            raise NameError('Substraction with a Point must be with another Point or a Vector')

    def __truediv__(self, sca):
        if (np.isscalar(sca)):
            return Vector(self.x/sca, self.y/sca, self.z/sca) 
        else:
            raise NameError('A Point can be divided only by a scalar')

    def __mul__(self, sca): 
        if (np.isscalar(sca)):
            return Vector(sca*self.x, sca*self.y, sca*self.z)
        else:
            raise NameError('A Point can be multiplied only by a scalar')
        
    def __getitem__(self, ind):
        if ( not isinstance(ind, int) or
             not isinstance(ind, np.integer) ):
            IndexError("Only an integer is a valid index")
        if ind == 0:
            return self.x
        elif ind == 1:
            return self.y
        elif ind == 2 :
            return self.z
        else:
            IndexError(f"Index {ind} is out of range")

    def __str__(self):
        return '(' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ')'
    
    def __repr__(self):
        return 'Point(' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ')'

    def to_numpy(self):
        return np.array([self.x, self.y, self.z], dtype=np.float64)

class Normal(object):
    """
    Parameters
    ----------
    x : float | Point | Vector | Normal | np.ndarray
        If scalar -> x component of the normal.
        Else, circumvent the y and z parameters and take the components of the Point/Vector/Normal/np.ndarray.
    y : float
        The y component of the normal.
    z : float
        The z component of the normal.

    Exemples
    --------
    >>> n1 = Normal(0.,0.,1.)
    >>> n1
    >>> Normal(0,0,1)
    """
    def __init__(self, x = 0., y = 0., z = 0.):
        if ( np.isscalar(x) and
             np.isscalar(y) and
             np.isscalar(z) ):
            self.x = float(x)
            self.y = float(y)
            self.z = float(z)
        elif ( isinstance(x, Vector) or
               isinstance(x, Point) or
               isinstance(x, Normal) ):
            self.x = x.x
            self.y = x.y
            self.z = x.z
        elif ( isinstance(x, np.ndarray) and
               len(x) == 3 ):
            self.x = float(x[0])
            self.y = float(x[1])
            self.z = float(x[2])
        else: 
            raise ValueError("Wrong parameter value(s) for Normal")

    def __eq__(self, n2):
        if isinstance(n2, Normal):
            return (self.x==n2.x) and (self.y==n2.y) and (self.z==n2.z)
        else:
            raise NameError('Equality with a Normal must be only with another Normal')

    def __add__(self, n2):
        if isinstance(n2, Normal):
            return Normal(self.x+n2.x, self.y+n2.y, self.z+n2.z) 
        else:
            raise NameError('Addition with a Normal must be only with another Normal')

    def __sub__(self, n2):
        if isinstance(n2, Normal):
            return Vector(self.x-n2.x, self.y-n2.y, self.z-n2.z)
        elif n2 == 0:
            return (self.x==-1*self.x) and (self.y==-1*self.y) and (self.z==-1*self.z)
        else:
            raise NameError('Substraction with a Normal must be only with another Normal')

    def __truediv__(self, sca):
        if (np.isscalar(sca)):
            return Vector(self.x/sca, self.y/sca, self.z/sca) 
        else:
            raise NameError('A Point can be divided only by a scalar')

    def __mul__(self, sca): 
        if (np.isscalar(sca)):
            return Vector(sca*self.x, sca*self.y, sca*self.z)
        else:
            raise NameError('A Point can be multiplied only by a scalar')
        
    def __getitem__(self, ind):
        if ( not isinstance(ind, int) or
             not isinstance(ind, np.integer) ):
            IndexError("Only an integer is a valid index")
        if ind == 0:
            return self.x
        elif ind == 1:
            return self.y
        elif ind == 2 :
            return self.z
        else:
            IndexError(f"Index {ind} is out of range")

    def __str__(self):
        return '(' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ')'
    
    def __repr__(self):
        return 'Normal(' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ')'
    
    def length(self):
        return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z) # L2 norm

    def to_numpy(self):
        return np.array([self.x, self.y, self.z], dtype=np.float64)


class Ray(object):
    """
    Definition of ray:

    r(t) = o + t*d, where:
    - o is the origin point of the ray
    - d is the direction of the ray
    - t is a positive real scalar

    Parameters
    ----------
    o : Point | Ray
        Origin point of the ray.
        If the o parameter is a Ray -> circumvent all the parameters by the ray attributs
    d : Vector
        Direction of the ray
    mint : float
        The minimum t value
    maxt : float
        The maximum t value

    Examples
    --------
    >>> o = Point(0., 50., 2.)
    >>> d = Vector(0.,0.,1.)
    >>> r1 = Ray(o, d, mint=20, maxt=100)
    >>> r1
    r(t) = (0.0, 50.0, 2.0) + t*(0.0, 0.0, 1.0) with t ∈ [20,100[
    """
    __point_ini = Point()
    __vector_ini = Vector()
    def __init__(self, o = __point_ini, d = __vector_ini, mint = 0, maxt = float("inf")):
        if isinstance(o, Ray):
            self.o = o.o
            self.d = o.d
            self.mint = o.mint
            self.maxt = o.maxt
        else:
            if (not isinstance(o, Point)):
                raise ValueError("The parameter o must be a Point or a Ray")
            if (not isinstance(d, Vector)):
                raise ValueError("The parameter d must only be a Vector")
            if (not np.isscalar(mint) or not np.isscalar(maxt)):
                raise ValueError("The parameters mint and maxt must be both scalars")
            if (mint > maxt):
                raise ValueError("maxt must be greater than mint")
            self.o = o
            self.d = d
            self.mint = mint
            self.maxt = maxt

    def __getitem__(self, t):
        if (not np.isscalar(t)):
            raise NameError('The value must be a scalar')
        elif ( t < self.mint or t > self.maxt):
            raise NameError(f"The value {t} is out of bounds. It must be between {self.mint} and {self.maxt}")
        else:
            return (self.o + self.d*t)
        
    def __repr__(self):
        return f'r(t) = ({self.o.x}, {self.o.y}, {self.o.z}) + t*({self.d.x}, {self.d.y}, {self.d.z})' + \
               f' with t ∈ [{self.mint},{self.maxt}['
    
    def __str__(self):
        return f'({self.o.x}, {self.o.y}, {self.o.z}) + t*({self.d.x}, {self.d.y}, {self.d.z})' + \
               f' with t ∈ [{self.mint},{self.maxt}['

class BBox(object):
    '''
    Bounding Box

    Parameters
    ----------
    p1 : Point
        Frist point to use to create the BBox
    p2 : Point
        Second point to use to create the BBox

    Examples
    --------
    >>> p1 = Point(0., 0., 0.)
    >>> p2 = Point(1., 1., 1.)
    >>> b1 = BBox(p1, p2)
    >>> b1
    pmin=Point(0.0, 0.0, 0.0), pmax=Point(1.0, 1.0, 1.0)
    '''
    __point_ini_min = Point(float("-inf"), float("-inf"), float("-inf"))
    __point_ini_max = Point(float("inf"), float("inf"), float("inf"))
    def __init__(self, p1 = __point_ini_min, p2 = __point_ini_max):
        if ( isinstance(p1, Point)  and isinstance(p2, Point) ):
            self.pmin = Point(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z))
            self.pmax = Point(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z))
        else:
            raise ValueError('Bounding Box constructor accepts only points')
        
    def union(self, b):
        """
        Union with a Point or another BBox

        Parameters
        ----------
        b : Point | BBox
            The point or BBox to use for the union

        Returns
        -------
        b_union : BBox
            The new BBox after the union
        
        Examples
        --------
        >>> p1 = Point(0., 0., 0.)
        >>> p2 = Point(1., 1., 1.)
        >>> p3 = Point(1., 1., 3.)
        >>> b1 = BBox(p1, p2)
        >>> b1
        pmin=Point(0.0, 0.0, 0.0), pmax=Point(1.0, 1.0, 1.0)
        >>> b2 = b1.union(p3)
        >>> b2
        pmin=Point(0.0, 0.0, 0.0), pmax=Point(1.0, 1.0, 3.0)
        """
        b_union = BBox()
        if isinstance(b, Point):
            b_union.pmin.x = min(self.pmin.x, b.x)
            b_union.pmin.y = min(self.pmin.y, b.y)
            b_union.pmin.z = min(self.pmin.z, b.z)
            b_union.pmax.x = max(self.pmax.x, b.x)
            b_union.pmax.y = max(self.pmax.y, b.y)
            b_union.pmax.z = max(self.pmax.z, b.z)
        elif isinstance(b, BBox):
            b_union.pmin.x = min(self.pmin.x, b.pmin.x)
            b_union.pmin.y = min(self.pmin.y, b.pmin.y)
            b_union.pmin.z = min(self.pmin.z, b.pmin.z)
            b_union.pmax.x = max(self.pmax.x, b.pmax.x)
            b_union.pmax.y = max(self.pmax.y, b.pmax.y)
            b_union.pmax.z = max(self.pmax.z, b.pmax.z)
        else:
            raise NameError('The union must be with another BBox or Point')

        return b_union

    def __repr__(self):
        return f'pmin=Point({self.pmin.x}, {self.pmin.y}, {self.pmin.z}), pmax=Point({self.pmax.x}, {self.pmax.y}, {self.pmax.z})'
    
    def __str__(self):
        return f'pmin=({self.pmin.x}, {self.pmin.y}, {self.pmin.z}), pmax=({self.pmax.x}, {self.pmax.y}, {self.pmax.z})'
            
    def get_vertices(self):
        """
        Get the 8 vertices of a BBox as a list of points

        - p0=pmin, then next 3 points are in the XY plane at z=pmin.z the order being anti-clockwise
        - next 4 points are in the XY plane at z=pmax.z, starting with point p4 just above p0, so p6=pmax
        """
        return [
            Point(self.pmin.x,self.pmin.y,self.pmin.z),
            Point(self.pmax.x,self.pmin.y,self.pmin.z),
            Point(self.pmax.x,self.pmax.y,self.pmin.z),
            Point(self.pmin.x,self.pmax.y,self.pmin.z),
            Point(self.pmin.x,self.pmin.y,self.pmax.z),
            Point(self.pmax.x,self.pmin.y,self.pmax.z),
            Point(self.pmax.x,self.pmax.y,self.pmax.z),
            Point(self.pmin.x,self.pmax.y,self.pmax.z)
        ]

    def is_inside(self, P):
        """
        Test if Point P is included in BBox
        """
        return (P.x >= self.pmin.x) and (P.x <= self.pmax.x) and \
               (P.y >= self.pmin.y) and (P.y <= self.pmax.y) and \
               (P.z >= self.pmin.z) and (P.z <= self.pmax.z)


    def intersectP(self, r1) :
        """
        Test if a ray intersects the BBox

        There are 3 possibilities:
        - no intersection
        - only 1 intersection (case of ray located initially inside the BBox)
        - 2 intersections

        Parameters
        ----------
        r1 : Ray
            The ray to use for the intersection test

        Returns
        -------
        t0 : float
            The t ray variable for the first intersection.
            In case of only 1 intersection it represents nothing. 
        t1 : float
            The t ray variable for the second intersection.
            In case of only 1 intersection, t1 becomes the t ray variable for the first intersection.
        is_intersection : bool
            If there is at least 1 intersection -> True, else False.

        Examples
        --------
        >>> p1 = Point(0., 0., 0.)
        >>> p2 = Point(1., 1., 1.)
        >>> b1 = BBox(p1, p2)
        pmin=Point(0.0, 0.0, 0.0), pmax=Point(1.0, 1.0, 1.0)
        >>> p3 = Point(0.5, 0.5, 0.1)
        >>> v1 = Vector(0., 0., 1.)
        >>> r1 = Ray(p3, v1)
        >>> r1
        r(t) = (0.5, 0.5, 0.1) + t*(0.0, 0.0, 1.0) with t ∈ [0,inf[
        >>> t0, t1, is_intersection = b1.intersectP(r1)
        >>> t0, t1, is_intersection
        (0.0, 0.9, True)
        >>> r1[t1]
        Point(0.5, 0.5, 1.0)
        """
        t0 = 0.
        epsi = 1e-32 * 0.5
        t1 = np.inf
        gamma3 = (3*epsi)/(1 - 3*epsi)
        for i in range(3):
            if r1.d[i]!= 0 :
                invRayDir = 1. / r1.d[i]
            else: invRayDir=1e32
            tNear = (self.pmin[i] - r1.o[i]) * invRayDir
            tFar  = (self.pmax[i] - r1.o[i]) * invRayDir

            if (tNear > tFar):
                tmp  = tNear
                tNear= tFar
                tFar = tmp
            tFar *= 1 + 2*gamma3
            t0 = tNear if tNear > t0 else t0
            t1 = tFar  if  tFar < t1 else t1
            if (t0 > t1) : return 0, 0, False
        return t0, t1, True


# def CommonVertices(BBox1, BBox2):
#     ''' 
#     return a list of boolean checking if vertices are common between BBoxes 
#     '''
#     return np.array(list((map(lambda x: x in BBox2.get_vertices(), BBox1.get_vertices()))))

# def CommonFace(BBox1, BBox2, Fill_value=None):
#     ''' 
#     return the face index of th BBox1 which is common to BBox2 with
#     the convention of index from 0 to 5, for +X,-X,+Y,-Y,+Z,-Z faces
#     # (en.wikipedia.org/wiki/Cube_mapping)
#     return Fill if there is no common face
#     '''
#     ok = CommonVertices(BBox1, BBox2)
#     if ok.sum()==4:
#         n  = np.arange(8)[ok]
#         if   np.array_equal(n, np.array([1,2,5,6])):
#             return 0
#         elif np.array_equal(n, np.array([0,3,4,7])):
#             return 1
#         elif np.array_equal(n, np.array([2,3,6,7])):
#             return 2
#         elif np.array_equal(n, np.array([0,1,4,5])):
#             return 3
#         elif np.array_equal(n, np.array([4,5,6,7])):
#             return 4
#         elif np.array_equal(n, np.array([0,1,2,3])):
#             return 5
#         else: return Fill_value

#     else :
#         return Fill_value

