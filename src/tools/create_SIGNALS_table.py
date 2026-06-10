#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 15:30:41 2022

@author: Jia Wei Teh

This script reads in the SIGNAL catalogue and creates a table for further
analysis. Note that this table contains all Class 1, 2, 3 and 4 HII regions. 
"""

# libraries
import numpy as np
from astropy.io import fits
#--
import src.tools.draw_FOV as draw_FOV


def create_SIGNALS_table():
    """
    This function reads in the SIGNAL catalogue and creates a table for further
    analysis.

    Returns
    -------
    h2_cat_ori : A list of arrays

    """
    
    def point_inside_polygon(x,y,poly):
        """
        Helper function which returns if points (x,y) is inside a polygon or not
        See: https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
    
        Parameters
        ----------
        x : float
            x location.
        y : float
            y location.
        poly : list
            Polygon list.
    
        Returns
        -------
        inside : boolean
            Is the given point within the polygon?
    
        """
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
    
    # path to catalogue
    path2h2 = r'/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/lib/SIGNALS/NGC628_catalog_WCS_corr.fits'
    
    # =============================================================================
    # Read file
    # =============================================================================
    # the fits
    fits_data = fits.open(path2h2)[0].data
    
    # create array to store catalogue. 
    h2_cat_1234 = []
    
    for row in fits_data:
        ra = float(row[1])
        dec = float(row[2])
        fov1, fov2 = draw_FOV.FOV()
        # check if HII region is in LEGUS FOV
        if not (point_inside_polygon(ra, dec, fov1) or 
            point_inside_polygon(ra, dec, fov2)):
            continue
        h2_cat_1234.append(np.array(row))
        
    return h2_cat_1234
    
    
    
    
    
    
    
    
    
    