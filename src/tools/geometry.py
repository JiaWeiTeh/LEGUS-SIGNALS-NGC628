# -*- coding: utf-8 -*-
"""Geometry & coordinate helpers shared across the pipeline.

Consolidated in Patch 4 from duplicated copies (create_combined_table.py,
create_LEGUS_table.py, create_SIGNALS_table.py, the former my_functions.py, and
several plot scripts). `ang_dist` uses PyAstronomy — the version that built the
committed catalogue.
"""
import numpy as np
import PyAstronomy.pyasl
from shapely.geometry import Point, Polygon
from sys import maxsize


def r2d(rad):
    return rad*180/np.pi


def d2r(deg):
    return deg*np.pi/180


def d2a(deg):
    return deg*3600


def a2d(arc):
    return arc/3600


def a2r(arc):
    return d2r(a2d(arc))


def r2a(rad):
    return d2a(r2d(rad))


def pc2arc(radius, D = 9.9e6):
    # distance to NGC 628
    return r2a(radius/D)


def arc2pc(ang, D = 9.9e6):
    # distance to NGC 628
    return a2r(ang)*D


def dist(p1, p2):
    """Return distance between two points"""
    
    return np.linalg.norm(np.array(p1)-np.array(p2))


def ang_dist(p1,p2):
    """Return angular distance (arcsecond) between two points of RA and DEC"""
    """Equation from https://en.wikipedia.org/wiki/Angular_distance"""
    
    r1, d1 = p1
    r2, d2 = p2
    
    r1 = r2d(r1)
    d1 = r2d(d1)
    r2 = r2d(r2)
    d2 = r2d(d2)
    
    # in rad
    deg = PyAstronomy.pyasl.getAngDist(r1, d1, r2, d2)
    
    return d2a(deg)


def min_dist(polygon, point):
    """
    Given polygon and point, return the nearest distance to any of the lines.

    Parameters
    ----------
    polygon : list of datapoints representing the polygon.
        Polygon of interest.
        E.g., [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)].
    point : shapely.geometry.point.Point
        Point to test.
        E.g., [(0.5, 1.5)].

    """


    # a very simple check
    assert sum(polygon[0] == polygon[-1]) == 2, "Please input valid polygon"
    assert len(point) == 2, "Pleaase input valid point"
    #Modification based on: https://gis.stackexchange.com/questions/396/nearest-neighbor-between-point-layer-and-line-layer
    
    polygon = Polygon(polygon)
    point = Point(point)
    
    # pairs iterator:
    # http://stackoverflow.com/questions/1257413/1257446#1257446
    def pairs(lst):
        i = iter(lst)
        first = prev = i.__next__()
        for item in i:
            yield prev, item
            prev = item
        yield item, first
    
    # these methods rewritten from the C version of Paul Bourke's
    # geometry computations:
    # http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
    def magnitude(p1, p2):
        vect_x = p2.x - p1.x
        vect_y = p2.y - p1.y
        return np.sqrt(vect_x**2 + vect_y**2)
    
    def intersect_point_to_line(point, line_start, line_end):
        line_magnitude =  magnitude(line_end, line_start)
        u = ((point.x - line_start.x) * (line_end.x - line_start.x) +
             (point.y - line_start.y) * (line_end.y - line_start.y)) \
             / (line_magnitude ** 2)
    
        # closest point does not fall within the line segment, 
        # take the shorter distance to an endpoint
        if u < 0.00001 or u > 1:
            ix = magnitude(point, line_start)
            iy = magnitude(point, line_end)
            if ix > iy:
                return line_end
            else:
                return line_start
        else:
            ix = line_start.x + u * (line_end.x - line_start.x)
            iy = line_start.y + u * (line_end.y - line_start.y)
            return Point([ix, iy])
    
    nearest_point = None
    min_dist = maxsize
    
    for seg_start, seg_end in pairs(list(polygon.exterior.coords)[:-1]):
        line_start = Point(seg_start)
        line_end = Point(seg_end)
    
        intersection_point = intersect_point_to_line(point, line_start, line_end)
        cur_dist =  magnitude(point, intersection_point)
    
        if cur_dist < min_dist:
            min_dist = cur_dist
    
    return  min_dist


def point_inside_polygon(x,y,poly):
    n = len(poly)
    inside =False
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y
    return inside


def ang_dist_matrix(ra1, dec1, ra2, dec2):
    """Vectorised pairwise angular distance in arcsec.

    ra1/dec1: shape (N,) in DEGREES; ra2/dec2: shape (M,) in DEGREES.
    Returns an (N, M) matrix. Mirrors ang_dist() called as
    ang_dist((d2r(ra), d2r(dec)), ...) exactly (including the deg<->rad round-trip),
    so it matches the per-pair loop bit-for-bit while running on whole arrays.
    """
    a1 = r2d(d2r(np.asarray(ra1, float)))[:, None]
    b1 = r2d(d2r(np.asarray(dec1, float)))[:, None]
    a2 = r2d(d2r(np.asarray(ra2, float)))[None, :]
    b2 = r2d(d2r(np.asarray(dec2, float)))[None, :]
    return d2a(PyAstronomy.pyasl.getAngDist(a1, b1, a2, b2))
