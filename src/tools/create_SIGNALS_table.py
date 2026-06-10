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
# --- allow running this file directly: put repo root on sys.path ---
import os as _os, sys as _sys
_root = _os.path.dirname(_os.path.abspath(__file__))
while not _os.path.isdir(_os.path.join(_root, "src")) and _root != _os.path.dirname(_root):
    _root = _os.path.dirname(_root)
if _root not in _sys.path:
    _sys.path.insert(0, _root)
# ------------------------------------------------------------------
from src import paths
from src.tools.geometry import point_inside_polygon
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
    
    
    # path to catalogue
    path2h2 = paths.SIGNALS + "NGC628_catalog_WCS_corr.fits"
    
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
