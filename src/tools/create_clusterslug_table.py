#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 13:37:51 2022

@author: Jia Wei Teh

This script extracts data from the output of cluster_slug
"""
# -*- coding: utf-8 -*-

import numpy as np

# path to clusterslug libraries
path2lib = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/lib/SLUG2/cluster_slug/"

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
    path2save = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/"
    np.save(path2save+"clusterslug_table_phi073_apn1_"+stellarTrack+classes+".npy", np.array([qtab, mtab, atab]))
    
    return qtab, mtab, atab

#%%

# For AV
import numpy as np
import matplotlib.pyplot as plt
path2lib = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/lib/SLUG2/cluster_slug/"
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






