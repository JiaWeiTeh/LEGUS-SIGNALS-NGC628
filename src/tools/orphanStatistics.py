#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 18:41:52 2022

@author: Jia Wei Teh

This script prints out statistics about the orphan HII regions.
"""

# libraries
import numpy as np
import matplotlib.pyplot as plt
#--
import src.tools.plot_tools as plot_tools
import src.tools.read_catalogue as read_catalogue

# Read in catalogue
path2save = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/"
_, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)

# Loop through the catalogue and look for orphan HII regions. The regions are considered 
# orphans if they 1) have no star clusters, and 2) contain only Class 4 sources.
orphan_ID, populated_ID = [], []

for h2 in h2_catalogue:
    # do not include class 3/4 clusters
    h2class = read_catalogue.get_h2_param(h2, 'class')
    if h2class == 4 or h2class == 3 :
        continue
    if not read_catalogue.is_BPTStarforming(h2):
        continue
    # h2region data
    # We are only interested in regions associated with more than one 
    # of Class 1,2,3 sources. Regions with purely Class 4 are considered
    # orphan regions here.
    assoc = read_catalogue.get_h2_param(h2, 'assoc', h2_catalogue)
    # returns True if list is empty or if it is only Class 4
    isOrphan = all([read_catalogue.get_sc_param(sc, 'class') == 4 for sc in assoc])
    # record ID for further operations
    if isOrphan:
        orphan_ID.append(read_catalogue.get_h2_param(h2, 'ID'))
    else:
        populated_ID.append(read_catalogue.get_h2_param(h2, 'ID'))
        
# =============================================================================
# Plot1: Median radius distribution
# =============================================================================
fig, ax =  plt.subplots(1,1, figsize = (6,5), dpi = 300)
# bin range
bins = np.linspace(25, 180, 15)
# For orphan regions
plot_tools.plot_hist([read_catalogue.get_h2_param(h2, 'rad') for h2 in orphan_ID],
                     setticks  = [20, 2, 10, 5],
                     bins = bins,
                     histtype = 'step',
                     linewidth = 1,
                     linestyle = '--',
                     alpha = .8,
                     color = 'k',
                     label = 'Orphans, $R_{\\rm med} = %.2f \\rm pc$'%(np.median([read_catalogue.get_h2_param(h2, 'rad') for h2 in orphan_ID])),
                     )
# For populated regions
plot_tools.plot_hist([read_catalogue.get_h2_param(h2, 'rad') for h2 in populated_ID],
                     bins = bins,
                     linewidth = 1.5,
                     histtype = 'step',
                     color = 'b',
                     label = 'Populated, $R_{\\rm med} = %.2f \\rm pc$'%(np.median([read_catalogue.get_h2_param(h2, 'rad') for h2 in populated_ID])),
                     )
plt.legend(loc = 'upper right')




