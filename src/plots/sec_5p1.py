#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 11:46:51 2023

@author: Jia Wei Teh

This script looks at the regions that neither associated with 1,2,3,4 soruces.
"""


# libraries
import numpy as np
import src.tools.draw_FOV as draw_FOV
import csv
#--
import src.tools.read_catalogue as read_catalogue
from src.tools.create_combined_table import d2r, a2d, pc2arc, arc2pc, ang_dist

# Read in catalogue
path2save = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/"
_, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)

numberOfEmptyRegions = 0 
overlapWithClass0 = 0

# =============================================================================
# # create a catalogue that contains class 0
# =============================================================================
# First define a helper function

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


def create_LEGUS_table_wclass0():
    """
    This function reads in the LEGUS catalogue and creates a table for further
    analysis.

    Returns
    -------
    sc_cat_1234 : A list of arrays

    """
    # path to catalogue
    path2sc = r'/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/lib/LEGUS/Catalog.csv'
    
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
    path2csv = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/lib/LEGUS/"
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
        










