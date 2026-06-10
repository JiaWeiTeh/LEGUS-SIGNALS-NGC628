#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 23:45:23 2022

@author: Jia Wei Teh

This script saves class 0 star clusters into .reg file
"""

""" """

# libraries
from PyAstronomy.pyasl import coordsDegToSexa
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
import src.tools.draw_FOV as draw_FOV
import csv
#--
import src.tools.read_catalogue as read_catalogue
from src.tools.create_combined_table import d2r, a2d, pc2arc, arc2pc, ang_dist


# Read in catalogue
path2save = paths.DAT
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)
# path to region files
path2reg = paths.DAT


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

with open(path2reg+"scC0.reg","w+") as scC0:
                
    for sc in sc_cat_01234:
        scClass = read_catalogue.get_sc_param(sc, 'class', sctable = sc_cat_01234)
        if scClass != 0:
            continue
        scRA = read_catalogue.get_sc_param(sc, 'RA',  sctable = sc_cat_01234)
        scDEC = read_catalogue.get_sc_param(sc, 'DEC',  sctable = sc_cat_01234)
        # write into file
        ra, dec = coordsDegToSexa(
            float(scRA),
            float(scDEC),
            fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f')
            ).split()
        shape = "box"
        w = "0.5\""
        h = "0.5\""
        colour = "colour=red"
        width = "width=2"
        line = [shape, ra, dec, w, h, "#",colour,width,"\n"]
        line = " ".join(line)
        scC0.write(line)
        
# done
scC0.close()
