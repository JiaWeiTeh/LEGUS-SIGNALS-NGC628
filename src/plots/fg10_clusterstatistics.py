#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 15:03:04 2022

@author: Jia Wei Teh

This script shows the relationship between inferred f_esc again properties
of star clusters and HII regions. See Figure 10.
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
from src.tools.stats import samplePDF
np.random.seed(0)  # reproducible Monte-Carlo (Patch 4: seeded)
import matplotlib.pyplot as plt
from tqdm import tqdm
#--
import src.tools.plot_tools as plot_tools
import src.tools.read_catalogue as read_catalogue

# Read in catalogue
path2save = paths.DAT
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)
clusterdata = np.load(path2save+"GenKroupc1234LHa_QH0_noCons.npy", allow_pickle = True)
qtab, mtab, atab = np.load(path2save+"clusterslug_table_phi073_apn1_geneva_1234.npy", allow_pickle = True)

# Unravel (see fg5_LHa_vs_QH0.py for more details)
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
LHa_log_n = np.array(LHa_log)[(xerr_max+xerr_min)<lim]
h2ID_n = np.array(h2ID)[(xerr_max+xerr_min)<lim]
qh0_total_pdf_n = np.array(qh0_total_pdf)[(xerr_max+xerr_min)<lim]
f_esc_pdf_n = np.array(f_esc_pdf)[(xerr_max+xerr_min)<lim]

# =============================================================================
# Additionally, we remove data points in the polygon regions.
# =============================================================================
# indeces datapoints that are in polygon after manual inspection
polygon_data = [0, 25, 38, 42]


# =============================================================================
# Now that we have double checked, we can update the list of datapoints
# which we will use to calculate f_esc. 
# =============================================================================
# Update list to remove points in polygons
h2ID_n = np.array([i for jj, i in enumerate(h2ID_n) if jj not in polygon_data])
f_esc_pdf_n = np.array([i for jj, i in enumerate(f_esc_pdf_n) if jj not in polygon_data])
qh0_total_pdf_n = np.array([i for jj, i in enumerate(qh0_total_pdf_n) if jj not in polygon_data])
LHa_log_n = np.array([i for jj, i in enumerate(LHa_log_n) if jj not in polygon_data])

# =============================================================================
# Define some function to help compute properties of clusters
# =============================================================================
def give_sortedMass(H2ID):
    """ Given hii region, return sorted median value of masses of clusters"""
    # use the H2 region to retrieve associated star clusters
    scIDlist = read_catalogue.get_h2_param(int(H2ID), "assoc")
    
    # create montecarlo samples
    montecarlo = np.zeros(shape = len(scIDlist), dtype = np.ndarray)
    for ii, ID in enumerate(scIDlist):
        range_log, mPDF = get_pdf(int(ID), mtab)
        montecarlo[ii] = list(np.percentile(samplePDF(range_log, mPDF), [50, 15.9, 84.1]))
    # sort 
    montecarlo = sorted(montecarlo)

    return np.array(montecarlo)

def give_sortedAge(H2ID):
    """ Given hii region, return sorted median value of masses of clusters"""
    # use the H2 region to retrieve associated star clusters
    scIDlist = read_catalogue.get_h2_param(int(H2ID), "assoc")
    
    # create montecarlo samples
    montecarlo = np.zeros(shape = len(scIDlist), dtype = np.ndarray)
    for ii, ID in enumerate(scIDlist):
        range_log, aPDF = get_pdf(int(ID), atab)
        montecarlo[ii] = list(np.percentile(samplePDF(range_log, aPDF), [50, 15.9, 84.1]))
    # sort 
    montecarlo = sorted(montecarlo)
    
    return np.array(montecarlo)
    
# =============================================================================
# Helper functions
# =============================================================================
# Helper function to create PDF
def get_pdf(ID, tabtype):
    """
    This function takes in cluster ID and returns the x and y values of 
    inferred PDF from clusterslug.

    Parameters
    ----------
    ID : int
        cluster ID.
    tabtype : quantity of interest
        qtab, mtab, or atab.

    Returns
    -------
    range_log : array
        x-values of PDF.
    pdf: 
        y-values of PDF.

    """
    # the range
    range_log = tabtype[0,:]
    # tab does not have ID. This will correlate the ID from catalogue into tab
    tab_id = sc_catalogue[:,0]
    # create pdf and range for given ID
    pdf = tabtype[1:,:][np.where(tab_id == ID)[0][0]]
    # return x and y values
    return range_log, pdf

# Helper function to sample PDF

# This function helps define errorbars in plots 
def errors(percentiles):
    """ return percentiles for plotting purposes"""
    mini2, mini, med, maxi, maxi2 = percentiles
    err_max = maxi - med
    err_min = med - mini
    err_max_sig2 = maxi2 - med
    err_min_sig2 = med - mini2
    return err_min_sig2, err_min, err_max, err_max_sig2

# =============================================================================
# We now are in position to loop through the HII regions and compute the statistics
# =============================================================================
sigma1 = []
ext = []
rGalac = []

for ii, (LHa, f_esc, ID) in enumerate(zip(LHa_log_n, f_esc_pdf_n, h2ID_n)):
    # get 2sigma, 1sigma and median
    percentiles = np.percentile(f_esc, (2.3, 15.9, 50, 84.1, 97.7))
    sigma1.append(percentiles)
    ext.append(read_catalogue.get_h2_param(ID, 'ext'))
    rGalac.append(read_catalogue.get_h2_param(ID, 'rGalac'))

reRun = False

if reRun:
    masslist = np.zeros(shape = len(h2ID_n), dtype = np.ndarray)
    agelist = np.zeros(shape = len(h2ID_n), dtype = np.ndarray)
    for ii, ID in enumerate(tqdm(h2ID_n)):
        masslist[ii] = give_sortedMass(ID)
        agelist[ii] = give_sortedAge(ID)
    np.save(path2save+"clusterMassMedians.npy", masslist)
    np.save(path2save+"clusterAgeMedians.npy", agelist)

else:
    masslist = np.load(path2save+"clusterMassMedians.npy", allow_pickle=True)
    agelist = np.load(path2save+"clusterAgeMedians.npy", allow_pickle=True)


# get lower values. If list length 4, return lowest two. if len = 5, return lower 2.
def get_lower(values):
    if len(values) == 1:
        return np.nan
    else:
        # values = values[values<7]
        return values[:int(np.floor(len(values)/2))]
def get_upper(values):
    if len(values) == 1:
        return np.nan
    else:
        # values = values[values<7]
        return values[-int(np.floor(len(values)/2)):]

# define plotting function locally to tweak stuffs
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
def plot_error(x, y, has_cmap = [], ax=None, clabel = [],
                 cticks = [],
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):
    """ Plots function. Accepts different axes """
    if ax is None:
        ax = plt.gca()
    if has_cmap:
        mappable = ax.scatter(x, y, c=has_cmap[0], cmap = has_cmap[1], **plt_kwargs)
        try:
            cbar = clabel[0].colorbar(mappable, pad = 1e-2)
            cbar.set_label(clabel[1], size = 12)
            cbar.ax.tick_params(labelsize=10)
        except IndexError:
            print("Add label for colorbar pls")
    else:
        ax.errorbar(x,y, **plt_kwargs)
    if cticks:
        cbar.ax.set_yticklabels(cticks)
    if xlabel is not None:
        ax.set_xlabel(xlabel, family = 'Times New Roman', fontsize = '20')
    if ylabel is not None:
        ax.set_ylabel(ylabel, family = 'Times New Roman', fontsize = '20')
    if title is not None:
        ax.set_title(title)
    if setticks:
        x, xinterval, y, yinterval = setticks
        ax.xaxis.set_major_locator(MultipleLocator(x))
        ax.xaxis.set_minor_locator(AutoMinorLocator(xinterval))
        ax.yaxis.set_major_locator(MultipleLocator(y))
        ax.yaxis.set_minor_locator(AutoMinorLocator(yinterval))
        ax.tick_params(axis='both', which = 'major', direction = 'in',length = 6, width = 1, labelsize = 18)
        ax.tick_params(axis='both', which = 'minor', direction = 'in',length = 4, width =1)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)  
    return(ax)

# Now, we plot
fig, ax =  plt.subplots(2,2, figsize = (10,10), dpi = 300)

def mini_helper(arr, ind):
    """This function returns the requried column. Only useful
    for this script because it is very much hardcoded."""
    new_array = np.ones(len(arr))
    for ii, [a] in enumerate(arr):
        new_array[ii] = a[ind]
    # return
    return new_array
        
# =============================================================================
# Age
# =============================================================================

# For age there is additional complication, and will be calculated as follows:
# Step 1: for each cluster, grab the median age
# Step 2: across all the median age, get the 25th percentile

def get_age_prop(mini_alist):
    """
    This function does what is mentioned above. 
    """
    
    # get the median value for each cluster
    medval = mini_alist[:,0]
    # get the 25th percentile for this list of median values
    return np.percentile(medval, 25)


# split into multiple and single clusters
# multiple
m_age_x = [ get_age_prop(a) for a in agelist]
m_age_y = [sig[2] for sig in sigma1]
# single
s_age_x = mini_helper([ a for a in agelist if len(a) == 1], 0)
s_age_y = [sig[2] for sig, a in zip(sigma1, agelist) if len(a) == 1]

plot_error(m_age_x, m_age_y,
        ylabel = "$f_{esc}$",             
        markersize = 10,
        xlabel = 'log $\\rm Age_{young} (yr)$',
          ylim = (-1, 1),
        fmt = 'yo',
         yerr = [
             [errors(sig)[1] for sig in sigma1]
             ,
             [errors(sig)[2] for sig in sigma1]
             ]
             ,
        setticks = [.1, 5, .5, 5],
        zorder = 1,
        linewidth = 1,
        mec = 'k',
        label = 'multiple star clusters',
        capsize = 2,
        ecolor = 'grey',
        ax = ax[0][0]
              )

plot_error(s_age_x, s_age_y,
        markersize = 10,
        linestyle = 'none',
        marker = 's',
        color = 'gray',
        yerr = [
            [errors(sig)[1] for sig, a in zip(sigma1, agelist) if len(a) == 1]
            ,
            [errors(sig)[2] for sig, a in zip(sigma1, agelist) if len(a) == 1]
            ]
            ,
        zorder = 1,
        linewidth = 1,
        label = 'single star cluster',
        mec = 'k',
        capsize = 2,
        ecolor = 'grey',
        ax = ax[0][0],
              )

# =============================================================================
# Mass
# =============================================================================
# split into multiple and single clusters
# multiple
m_mass_x = [ np.log10(sum(10**(m[:,0]))) for m in masslist]
m_mass_y = [sig[2] for sig in sigma1]
# single
s_mass_x = mini_helper([ m for m in masslist if len(m) == 1], 0)
s_mass_y = [sig[2] for sig, a in zip(sigma1, masslist) if len(a) == 1]


plot_error(m_mass_x, m_mass_y,
         markersize = 10,
         xlabel = 'log $\\rm Mass (M_\\odot)$',
         # setticks = [0.1, 5, .5, 5],
           ylim = (-1, 1),
         fmt = 'yo',
          yerr = [
              [errors(sig)[1] for sig in sigma1]
              ,
              [errors(sig)[2] for sig in sigma1]
              ]
              ,
        setticks = [.5, 5, .5, 5],
        zorder = 1,
        linewidth = 1,
        mec = 'k',
        label = 'multiple star clusters',
        capsize = 2,
        ecolor = 'grey',
        ax = ax[0][1],
              )


plot_error(s_mass_x, s_mass_y,
            markersize = 10,
            linestyle = 'none',
            marker = 's',
            color = 'gray',
            yerr = [
                [errors(sig)[1] for sig, a in zip(sigma1, agelist) if len(a) == 1]
                ,
                [errors(sig)[2] for sig, a in zip(sigma1, agelist) if len(a) == 1]
                ]
                ,
            zorder = 1,
            linewidth = 1,
            label = 'single star cluster',
            mec = 'k',
            capsize = 2,
            ecolor = 'grey',
            ax = ax[0][1],
              )

# hide axis
# ax[0][1].get_yaxis().set_ticklabels([])

# =============================================================================
# Extinction
# =============================================================================
m_ext_x = ext
m_ext_y = [sig[2] for sig in sigma1]
s_ext_x = [ ext[ii] for ii, a in enumerate(agelist) if len(a) == 1]
s_ext_y = [sig[2] for sig, a in zip(sigma1, agelist) if len(a) == 1]

plot_error(m_ext_x, m_ext_y,
             ylabel = "$f_{esc}$",             
             markersize = 10,
             xlabel = 'SIGNALS $\\rm E(B-V)$',
             # setticks = [0.1, 5, .5, 5],
               ylim = (-1, 1),
             fmt = 'yo',
              yerr = [
                  [errors(sig)[1] for sig in sigma1]
                  ,
                  [errors(sig)[2] for sig in sigma1]
                  ]
                  ,
            zorder = 1,
            linewidth = 1,
            setticks = [.1, 5, .5, 5],
            mec = 'k',
            label = 'multiple star clusters',
            capsize = 2,
            ecolor = 'grey',
            ax = ax[1][0],
              )

plot_error(s_ext_x, s_ext_y,
              markersize = 10,
              linestyle = 'none',
              marker = 's',
              color = 'gray',
              yerr = [
                  [errors(sig)[1] for sig, a in zip(sigma1, agelist) if len(a) == 1]
                  ,
                  [errors(sig)[2] for sig, a in zip(sigma1, agelist) if len(a) == 1]
                  ]
                  ,
            zorder = 1,
            linewidth = 1,
            label = 'single star cluster',
            mec = 'k',
            capsize = 2,
            ecolor = 'grey',
            ax = ax[1][0],
              )


# =============================================================================
# galactocentric radius
# =============================================================================
m_gr_x = rGalac
m_gr_y = [sig[2] for sig in sigma1]
s_gr_x = [ rGalac[ii] for ii, a in enumerate(agelist) if len(a) == 1]
s_gr_y = [sig[2] for sig, a in zip(sigma1, agelist) if len(a) == 1]


plot_error(m_gr_x, m_gr_y, 
             markersize = 10,
             xlabel = '$\\rm r_{galac} (kpc)$',
             # setticks = [0.1, 5, .5, 5],
               ylim = (-1, 1),
             fmt = 'yo',
              yerr = [
                  [errors(sig)[1] for sig in sigma1]
                  ,
                  [errors(sig)[2] for sig in sigma1]
                  ]
                  ,
            zorder = 1,
            linewidth = 1,
            setticks = [1, 5, .5, 5],
            mec = 'k',
            label = 'multiple star clusters',
            capsize = 2,
            ecolor = 'grey',
            ax = ax[1][1],
              )

plot_error(s_gr_x, s_gr_y, 
              markersize = 10,
              linestyle = 'none',
              marker = 's',
              color = 'gray',
              yerr = [
                  [errors(sig)[1] for sig, a in zip(sigma1, agelist) if len(a) == 1]
                  ,
                  [errors(sig)[2] for sig, a in zip(sigma1, agelist) if len(a) == 1]
                  ]
                  ,
            zorder = 1,
            linewidth = 1,
            label = 'single star cluster',
            mec = 'k',
            capsize = 2,
            ecolor = 'grey',
            ax = ax[1][1],
              )
 
ax[1][1].legend(loc = 'lower left', fontsize = 18, framealpha = 0.5)

plotlabels = ["a)", "b)", "c)", "d)"]
row = 0
col = 0
for char in plotlabels:
    ax[row][col].text(.05, 0.9, char,
            fontsize  = 20,
           transform=ax[row][col].transAxes
         )
    col += 1
    if col > 1:
        row += 1
        col = 0

plt.tight_layout(w_pad = -.3)
plot_tools.save("cluster_statistics")
