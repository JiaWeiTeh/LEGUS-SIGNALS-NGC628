#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 15:30:41 2022

@author: Jia Wei Teh

This script contains functions which creates a combined catalogue of H2 region and star clusters,
with geometrical approach. 
    
"""

# libraries
import numpy as np
from shapely.geometry import Point, Polygon
from sys import maxsize
import PyAstronomy.pyasl
#--
import src.tools.create_LEGUS_table as create_LEGUS_table
import src.tools.create_SIGNALS_table as create_SIGNALS_table

# =============================================================================
# Create table that includes which cluster links to hii region and vice versa.
# =============================================================================
def create_combined_table():
    """
    Returns
    -------
    sc_cat : np.array
        star cluster catalogue.
    h2_cat : np.array
        hii region catalogue.
    """
    
    # For h2 catalogue
    # create copy of catalogue to modify
    h2_catalogue_unmod = create_SIGNALS_table.create_SIGNALS_table()
    rowdim, coldim = np.array(h2_catalogue_unmod).shape
    h2_cat = np.zeros(shape = (rowdim,coldim+2), dtype=np.ndarray)
    h2_cat[:,:-2] = h2_catalogue_unmod
        
    # For sc catalogue
    # create copy of catalogue to modify
    sc_catalogue_unmod = create_LEGUS_table.create_LEGUS_table()
    rowdim, coldim = np.array(sc_catalogue_unmod).shape
    sc_cat = np.zeros(shape = (rowdim,coldim+2), dtype=np.ndarray)
    sc_cat[:,:-2] = sc_catalogue_unmod

    # For H II regions, find out which star clusters are they overlapping with.
    for h in h2_cat:
        hID = h[0]
        p1 = (d2r(h[1]),d2r(h[2]))
        radius = pc2arc(h[12])
        sc_ID_in_this_h2, sc_d_in_this_h2 = [], []
        for s in sc_catalogue_unmod:
            scID = s[0]
            p2 = (d2r(s[3]),d2r(s[4]))
            dist = ang_dist(p1, p2)
            if (dist <= radius):
                sc_d_in_this_h2.append(dist)
                sc_ID_in_this_h2.append(scID)
        # record distance from each SC to H II region
        h[-2] = np.array(sc_d_in_this_h2)
        # record the list of ID of SC in this H II region
        h[-1] = np.array(sc_ID_in_this_h2)
    
    # For star clusters, find out which H II regions are they overlapping with.
    for s in sc_cat:
        scID = s[0]
        p2 = (d2r(s[3]),d2r(s[4]))
        h2_ID_in_this_sc, h2_d_in_this_sc = [], []
        for h in h2_cat:
            radius = pc2arc(h[12])
            hID = h[0]
            p1 = (d2r(h[1]),d2r(h[2]))
            dist = ang_dist(p1, p2)
            if (dist <= radius):
                h2_ID_in_this_sc.append(hID)
                h2_d_in_this_sc.append(dist)
        # record distance from each H II region to SC
        s[-2] = np.array(h2_d_in_this_sc)
        # record the list of ID of H II region in this SC
        s[-1] = np.array(h2_ID_in_this_sc)   
        
    # save table to avoid computing time in the future
    path2save = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/"
    # store as an explicit 2-element object array (robust to the two catalogues
    # having equal row counts, which np.array([...], dtype=object) cannot handle)
    combined = np.empty(2, dtype=object)
    combined[0], combined[1] = sc_cat, h2_cat
    np.save(path2save+"combined_catalogue.npy", combined)
   
    # return
    return sc_cat, h2_cat


# =============================================================================
# Helper functions
# =============================================================================

# Convert rad to deg

def r2d(rad):
    return rad*180/np.pi

# Convert deg to rad
def d2r(deg):
    return deg*np.pi/180

# Convert degree to arcsecond

def d2a(deg):
    return deg*3600

# Convert arcsec to deg

def a2d(arc):
    return arc/3600

# Convert arc to rad

def a2r(arc):
    return d2r(a2d(arc))

# Convert rad to arc

def r2a(rad):
    return d2a(r2d(rad))

# Convert pc to arcseconds

def pc2arc(radius, D = 9.9e6):
    # distance to NGC 628
    return r2a(radius/D)

# Convert arcseconds to pc

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

