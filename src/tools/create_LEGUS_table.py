#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 15:30:41 2022

@author: Jia Wei Teh

This script reads in the LEGUS catalogue and creates a table for further
analysis. Note that this table only contains Class 1, 2, 3 and 4 sources, 
and we consider only sources with nflt flags = 4 or 5.
"""

# libraries
import csv 
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
import numpy as np
#--
import src.tools.draw_FOV as draw_FOV


def create_LEGUS_table():
    """
    This function reads in the LEGUS catalogue and creates a table for further
    analysis.

    Returns
    -------
    sc_cat_1234 : A list of arrays

    """
    # path to catalogue
    path2sc = paths.LEGUS + "Catalog.csv"
    
    # create array to store catalogue. 
    sc_cat_1234 = []
    # counter to check if clusters are from centre/east pointing
    count = 1
    
    # =============================================================================
    # Read file
    # =============================================================================
    # FoV polygons are constant across rows -- build them once, not per row.
    fov1, fov2 = draw_FOV.FOV()
    with open(path2sc) as csv_file:
        reader = csv.reader(csv_file, delimiter=',')
        # skip 1st row
        names = next(reader)
        for row in reader:
            count += 1
            # remove replication and comments
            if row[0] == '#':
                continue
            # only consider nflt = 4 and 5
            if int(row[33]) < 4:
                continue
            # remove Class 0
            if int(row[34]) == 0:
                continue
            # centre pointings have index 1 to 3080
            # new name for east pointings, 1e -> 3081
            if count > 3080 :
                row[1] = str((int(row[1]) + 3080))
            # check if ra and dec is in the FOV.
            ra = float(row[4])
            dec = float(row[5])
            if not (point_inside_polygon(ra, dec, fov1) or
                point_inside_polygon(ra, dec, fov2)):
                continue
            # remove '#' column that indicates replication
            row = np.array(row[1:]).astype(float)
            # store into a list
            sc_cat_1234.append(row)
            
    # =============================================================================
    # Save table
    # =============================================================================
    # file
    path2csv = paths.LEGUS
    # svae
    np.savetxt(path2csv+'sc1234_valid_NGC628.csv', sc_cat_1234, delimiter=',')
    
    return sc_cat_1234
