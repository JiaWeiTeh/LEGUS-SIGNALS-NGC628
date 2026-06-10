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
# --- allow running this file directly: put repo root on sys.path ---
import os as _os, sys as _sys
_root = _os.path.dirname(_os.path.abspath(__file__))
while not _os.path.isdir(_os.path.join(_root, "src")) and _root != _os.path.dirname(_root):
    _root = _os.path.dirname(_root)
if _root not in _sys.path:
    _sys.path.insert(0, _root)
# ------------------------------------------------------------------
from src import paths
from src.tools.geometry import (r2d, d2r, d2a, a2d, a2r, r2a,
                                   pc2arc, arc2pc, dist, ang_dist, ang_dist_matrix,
                                   min_dist, point_inside_polygon)
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

    # =========================================================================
    # Spatial match (vectorised). A star cluster belongs to an H II region if it
    # lies within that region's radius. Pairwise angular distances use the same
    # getAngDist as ang_dist (broadcast), so associations match the per-pair loop.
    # =========================================================================
    h2_RA  = np.array([h[1] for h in h2_cat], dtype=float)
    h2_DEC = np.array([h[2] for h in h2_cat], dtype=float)
    h2_ID  = np.array([h[0] for h in h2_cat], dtype=float)
    radius = np.array([pc2arc(h[12]) for h in h2_cat], dtype=float)   # arcsec, per region
    sc_RA  = np.array([s[3] for s in sc_catalogue_unmod], dtype=float)
    sc_DEC = np.array([s[4] for s in sc_catalogue_unmod], dtype=float)
    sc_ID  = np.array([s[0] for s in sc_catalogue_unmod], dtype=float)

    # (Nh2 x Nsc) angular-distance matrix (arcsec); within = distance <= region radius
    dmat = ang_dist_matrix(h2_RA, h2_DEC, sc_RA, sc_DEC)
    within = dmat <= radius[:, None]

    # For each H II region: the star clusters inside it (cluster order preserved)
    for i, h in enumerate(h2_cat):
        j = np.flatnonzero(within[i])
        h[-2] = dmat[i, j]
        h[-1] = sc_ID[j]
    # For each star cluster: the H II regions that contain it (region order preserved)
    for j, s in enumerate(sc_cat):
        i = np.flatnonzero(within[:, j])
        s[-2] = dmat[i, j]
        s[-1] = h2_ID[i]
        
    # save table to avoid computing time in the future
    path2save = paths.DAT
    # store as an explicit 2-element object array (robust to the two catalogues
    # having equal row counts, which np.array([...], dtype=object) cannot handle)
    combined = np.empty(2, dtype=object)
    combined[0], combined[1] = sc_cat, h2_cat
    np.save(path2save+"combined_catalogue.npy", combined)
   
    # return
    return sc_cat, h2_cat
