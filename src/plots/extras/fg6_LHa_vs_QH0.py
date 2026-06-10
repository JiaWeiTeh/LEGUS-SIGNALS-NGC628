#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 20:48:30 2022

@author: Jia Wei Teh

This script plots the LHa of H II regions against QH0 of star clusters.
This one includes Class 1,2,3,4 Clusters. This is after the 1-dex constrain.
See Figure 6.
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
import matplotlib.pyplot as plt
import cmasher as cmr
from matplotlib.patches import Polygon
import numpy.ma as ma
#--
import src.tools.plot_tools as plot_tools
import src.tools.create_clusterslug_table as create_clusterslug_table

# Read in catalogue
path2save = paths.DAT
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)
try:
    qtab, mtab, atab = np.load(path2save+"clusterslug_table_phi073_apn1_geneva_1234.npy", allow_pickle = True)
except:
    create_clusterslug_table.create_clusterslug_table(stellarTrack = 'geneva')
clusterdata = np.load(path2save+"GenKroupc1234LHa_QH0_noCons.npy", allow_pickle = True)

# Unravel
f_esc_pdf, f_esc_percentiles, qh0_total_pdf, qh0_percentile, LHa_log, h2ID, median_mass = list(zip(*clusterdata))

# sigma values and error bar parameters
# 1-sigma
xmed = np.array(qh0_percentile)[:,2]
xmax = np.array(qh0_percentile)[:,3]
xmin = np.array(qh0_percentile)[:,1]
xerr_max = xmax - xmed
xerr_min = xmed - xmin
# 2-sigma
xmax_sig2 = np.array(qh0_percentile)[:,4]
xmin_sig2 = np.array(qh0_percentile)[:,0]
xerr_max_sig2 = xmax_sig2 - xmed
xerr_min_sig2 = xmed - xmin_sig2

# =============================================================================
# Place a filter on errorbar constrain, such that any data points with 
# left+right error greater than this value (in units of dex) are removed
# =============================================================================
lim = 0.5
# new value after error constrain
xmed_n = xmed[(xerr_max+xerr_min)<lim]
xerr_max_n = xerr_max[(xerr_max+xerr_min)<lim]
xerr_min_n = xerr_min[(xerr_max+xerr_min)<lim]
xerr_max_sig2_n = xerr_max_sig2[(xerr_max+xerr_min)<lim]
xerr_min_sig2_n = xerr_min_sig2[(xerr_max+xerr_min)<lim]
LHa_log_n = np.array(LHa_log)[(xerr_max+xerr_min)<lim]
h2ID_n = np.array(h2ID)[(xerr_max+xerr_min)<lim]
median_mass_n = np.array(median_mass)[(xerr_max+xerr_min)<lim]
qh0_total_pdf_n = np.array(qh0_total_pdf)[(xerr_max+xerr_min)<lim]
f_esc_pdf_n = np.array(f_esc_pdf)[(xerr_max+xerr_min)<lim]

# remove weird datapoints
polygon_data = [0, 25, 38, 42]
mask_array = np.zeros(len(xmed_n), dtype=int)
mask_array[polygon_data] = 1
xmed_n = ma.masked_array(xmed_n, mask = mask_array)
xerr_max_n = ma.masked_array(xerr_max_n, mask = mask_array)
xerr_min_n = ma.masked_array(xerr_min_n, mask = mask_array)
LHa_log_n = ma.masked_array(LHa_log_n, mask = mask_array)
median_mass_n = ma.masked_array(median_mass_n, mask = mask_array)



def draw_polygon(array, color, ax = None):
    """ This function draws polygon over forbidden zones"""
    # get current axis
    if ax is None:
        ax = plt.gca()
    # draw polygon
    p = Polygon(array, facecolor = color, alpha = 0.3)
    # add polygon
    ax.add_patch(p)
    return

# =============================================================================
# Plot 
# =============================================================================
fig, ax =  plt.subplots(1,1, figsize = (6,6), dpi = 300)

# Plot datapoints
plot_tools.plot_scatter(xmed_n, LHa_log_n, 
                  has_cmap = [(median_mass_n), cmr.chroma_r], 
                  clabel = [fig, "Log Cluster total mass ($M_\odot$)"],
                  xlabel = "Log $\\rm Q(H^0)$ (photon s$^{-1}$)",
                  ylabel = "Log $L_{\\rm H\\alpha}$ (erg s$^{-1}$)",
                   xlim = (47.5,52),
                  ylim = (36.4, 39.4),
                  edgecolor = 'k',
                  vmin = 2.7, 
                  vmax = 4.7,
                  zorder = 100,
                  marker = 'o',
                  s = 50,
                  setticks = [1,5,0.5,5],
                  label = 'NGC 628 (This work)'
                  )

# Plot errorbars
ax.errorbar(xmed_n, LHa_log_n, 
                xerr = [xerr_min_n, xerr_max_n],
                fmt = 'o',
                markersize =5,
                alpha = 0.7,
                zorder = 1,
                linewidth = 2,
                mec = 'k',
                capsize = 2,
                ecolor = 'darkgrey',
                )

plot_tools.plot_studies()

# miscellaneous
plt.text(48, 39, "``Forbidden Zone\"", fontstyle = "italic", size = 15)   
plt.text(48.3, 36.8, "$f_{\\rm esc}=0$", rotation = 58, size = 19) 
plt.legend(loc = 'lower right',
           facecolor = 'none',
           fontsize = 9,
           labelspacing = .5,
           frameon = True)

plot_tools.save("GenKroupc1234LHa_QH0_wCons_poster")









