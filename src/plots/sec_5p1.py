#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 11:46:51 2023

@author: Jia Wei Teh

This script looks at the regions that neither associated with 1,2,3,4 soruces.
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
import src.tools.draw_FOV as draw_FOV
import csv
#--
import src.tools.read_catalogue as read_catalogue
from src.tools.create_combined_table import d2r, a2d, pc2arc, arc2pc, ang_dist

# Read in catalogue
path2save = paths.DAT
_, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)

numberOfEmptyRegions = 0 
overlapWithClass0 = 0

# =============================================================================
# # create a catalogue that contains class 0
# =============================================================================
# First define a helper function



def create_LEGUS_table_wclass0():
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
    sc_cat_01234 = []
    # counter to check if clusters are from centre/east pointing
    count = 1
    
    # =============================================================================
    # Read file
    # =============================================================================
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
            # centre pointings have index 1 to 3080
            # new name for east pointings, 1e -> 3081
            if count > 3080 :
                row[1] = str((int(row[1]) + 3080))
            # check if ra and dec is in the FOV.
            ra = float(row[4])
            dec = float(row[5])
            fov1, fov2 = draw_FOV.FOV()
            if not (point_inside_polygon(ra, dec, fov1) or 
                point_inside_polygon(ra, dec, fov2)):
                continue
            # remove '#' column that indicates replication
            row = np.array(row[1:]).astype(float)
            # store into a list
            sc_cat_01234.append(row)
            
    # =============================================================================
    # Save table
    # =============================================================================
    # file
    path2csv = paths.LEGUS
    # svae
    np.savetxt(path2csv+'sc1234_valid_NGC628_wClass0.csv', sc_cat_01234, delimiter=',')
    
    return sc_cat_01234
        
sc_cat_01234 = create_LEGUS_table_wclass0()

for h2 in h2_catalogue:
    # do not include class 3/4 clusters
    h2class = read_catalogue.get_h2_param(h2, 'class')
    if h2class == 4 or h2class == 3 :
        continue
    if not read_catalogue.is_BPTStarforming(h2):
        continue
    assoc = read_catalogue.get_h2_param(h2, 'assoc', h2_catalogue)
    # if the hii region is empty;
    if not assoc.size: 
        numberOfEmptyRegions += 1
        # Now we check if they overlap with class 0
        hID = h2[0]
        p1 = (d2r(h2[1]),d2r(h2[2]))
        radius = pc2arc(h2[12])
        for s in sc_cat_01234:
            scID = s[0]
            scClass = s[33]
            if scClass == 0:
                p2 = (d2r(s[3]),d2r(s[4]))
                dist = ang_dist(p1, p2)
                if (dist <= radius):
                    overlapWithClass0 += 1
                    break
    
print(f"There are {numberOfEmptyRegions} HII regions with nothing in them, and {overlapWithClass0} of them has at least one Class 0.")
