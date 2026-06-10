#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 13:37:51 2022

@author: Jia Wei Teh

This script extracts data from the output of cluster_slug
"""
# -*- coding: utf-8 -*-

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

# path to clusterslug libraries
path2lib = paths.SLUG_CLUSTER

# =============================================================================
# Create table
# =============================================================================
def create_clusterslug_table(stellarTrack = 'geneva',
                             includeClass4 = True,
                             ):
    """This function creates tables of PDF of age, mass and qh0 for clusters, 
    using cluster_slug. Input prameter includes the the stellar track and
    if we want to consider class 4 LEGUS sources."""
            
    # stellar track?
    if stellarTrack == 'geneva':
        track = '_modc020'
    elif stellarTrack == 'padova':
        track = '_modp020'
    else: 
        raise Exception('stellarTrack not found.')
    # do we include class 4?
    if includeClass4:
        classes = '_1234'
    else:
        classes = '_123'
    
    qtab = np.genfromtxt(path2lib + 'q0_table_phi073_apn1' + track + classes + ".tab", skip_header = 1)
    mtab = np.genfromtxt(path2lib + 'mass_table_phi073_apn1' + track + classes + ".tab", skip_header = 1)
    atab = np.genfromtxt(path2lib + 'age_table_phi073_apn1' + track + classes + ".tab", skip_header = 1)
    
    
    # save table to avoid computing time in the future
    path2save = paths.DAT
    # np.stack requires the three tables to share a shape (they always should) and
    # yields the identical (3, R, C) array, but fails loudly instead of silently
    # building a ragged object array if they ever diverge.
    np.save(path2save+"clusterslug_table_phi073_apn1_"+stellarTrack+classes+".npy", np.stack([qtab, mtab, atab]))
    
    return qtab, mtab, atab

#%%

# For AV
import numpy as np
import matplotlib.pyplot as plt
path2lib = paths.SLUG_CLUSTER
# stellar track?
track = '_modc020'
classes = '_1234'

avtab = np.genfromtxt(path2lib + 'av_table_phi073_apn1' + track + classes + ".tab", skip_header = 1)
    
plt.figure()
def create_overview(tabtype, ax):
    x = tabtype[0,:]
    pdfs = tabtype[1:,:]
    for i in range(500):
        ax.plot(x, -pdfs[i,:], alpha = .005)  
        
create_overview(avtab, plt.gca())






