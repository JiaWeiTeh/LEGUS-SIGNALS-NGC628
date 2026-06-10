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
import matplotlib.pyplot as plt
import cmasher as cmr
#--
import src.tools.plot_tools as plot_tools
import src.tools.read_catalogue as read_catalogue

# Read in catalogue
path2save = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/"
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)
clusterdata = np.load(path2save+"GenKroupc1234LHa_QH0_noCons.npy", allow_pickle = True)

# For AV
path2lib = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/lib/SLUG2/cluster_slug/"
# stellar track?
track = '_modc020'
classes = '_1234'
avtab = np.genfromtxt(path2lib + 'av_table_phi073_apn1' + track + classes + ".tab", skip_header = 1)

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

def give_sigmaAV(H2ID):
    """ Given hii region, return sorted median value of masses of clusters"""
    # use the H2 region to retrieve associated star clusters
    scIDlist = read_catalogue.get_h2_param(int(H2ID), "assoc")
    
    # create montecarlo samples
    montecarlo = np.zeros(shape = (len(scIDlist), 3), dtype = np.ndarray)
    std = np.zeros(shape = (len(scIDlist), 1), dtype = np.ndarray)
    for ii, ID in enumerate(scIDlist):
        range_log, aPDF = get_pdf(int(ID), avtab)
        montecarlo[ii] = np.percentile(samplePDF(range_log, aPDF), [50, 15.9, 84.1])
    # return
    return np.median(montecarlo, axis = 0)
    
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
    range_log = tabtype[0,:][::-1]
    # tab does not have ID. This will correlate the ID from catalogue into tab
    tab_id = sc_catalogue[:,0]
    # create pdf and range for given ID
    pdf = tabtype[1:,:][np.where(tab_id == ID)[0][0]][::-1]  * -1
    # print(pdf)
    # return x and y values
    return range_log, pdf

# Helper function to sample PDF
def samplePDF(xvalues, pdf):
    """ This function creates pdf that is not spiky (smoothes out the pdf).
    Solves problems arising from the 'delta-function' monte carlo approach.
    Then, it picks 10**5 samples from the pdf."""
    
    def normalisePDF(a, b, m, c):
        """Given xs and ys value of a bin, find the normalisation factor"""
        area = 1/2 * (a[1]+b[1])*abs(a[0]-b[0])
        return 1/area
        
    def line_func(p1, p2):
        """ Given two points, return the line coefficients m, c in y = mx + c"""
        m = (p1[1]-p2[1])/(p1[0]-p2[0])
        c = p1[1] - m * p1[0]
        return m, c

    # The probability of landing into one of the bins
    probability = abs(xvalues[1]-xvalues[0]) * (pdf[1:]+pdf[:-1])/2
    # normalize to bypass error tolerence 
    probability = probability/np.sum(probability)
    # pick a choice. which bin (left)? Number of iterations?
    # left means the selected bin is identified by the left side of the bin.
    niter = int(1e5)
    selectbin = np.random.choice(xvalues[:-1], p = probability, size = niter)
    # pick random number between [0, 1)
    rand_uni = np.random.uniform(0,1,niter)
    # array for montecarlo pdf values
    montecarlos = np.zeros(shape = (len(selectbin)))
    # For each selected bin,    
    for ii, bins in enumerate(selectbin):
        binindex = np.where(xvalues == bins)[0][0]
        # get the two points representing this bin
        a = (xvalues[binindex], pdf[binindex])
        b = (xvalues[binindex+1],pdf[binindex+1])
        # get coeffs of the line equation, m and c in y = mx + c
        m, c = line_func(a, b)
        # normalisation factor
        norm = normalisePDF(a, b, m, c)
        # solve for the y-value in the bin that corresponds to the uniform distribution
        roots = np.roots([m/2,\
                  c,\
                  -( m / 2 * a[0]**2  + c * a[0]  + rand_uni[ii]/norm )])
        # the correct root is within the bin
        hasroot = False
        # search for the correct solution in given roots
        for r in roots:
            if r < b[0] and r > a[0]:
                hasroot = True
                montecarlos[ii] = r
        if not hasroot:  
            raise Exception("Root not found")
    # return
    return montecarlos

# This function helps define errorbars in plots 
def errors(percentiles):
    """ return percentiles for plotting purposes"""
    med, mini, maxi = percentiles
    err_max = maxi - med
    err_min = med - mini
    return err_min, err_max

# =============================================================================
# We now are in position to loop through the HII regions and compute the statistics
# =============================================================================
ext_SIGNALS = np.zeros(shape = len(h2ID_n), dtype = np.ndarray)
ext_SIGNALS_err = np.zeros(shape = len(h2ID_n), dtype = np.ndarray)
ext_LEGUS = np.zeros(shape = len(h2ID_n), dtype = np.ndarray)
ext_std_LEGUS = np.zeros(shape = len(h2ID_n), dtype = np.ndarray)

sigma1 = []
for f_esc in f_esc_pdf_n:
    median = np.percentile(f_esc, (50))
    sigma1.append(median)



for ii, ID in enumerate(h2ID_n):
    ext_SIGNALS[ii] = read_catalogue.get_h2_param(ID, 'ext')
    ext_SIGNALS_err[ii] = read_catalogue.get_h2_param(ID, 'ext_err')
        
reRun = False

if reRun:
    for ii, ID in enumerate(h2ID_n):
        # get 2sigma, 1sigma and median
        ext_LEGUS[ii] = give_sigmaAV(int(ID))
    np.save(path2save+"clusterAVMedians.npy", ext_LEGUS)
    # np.save(path2save+"clusterAVstd.npy", ext_std_LEGUS)

else:
    ext_LEGUS = np.load(path2save+"clusterAVMedians.npy", allow_pickle=True)
    # ext_std_LEGUS = np.load(path2save+"clusterAVstd.npy", allow_pickle=True)


# Function that does scatter plots
def plot_scatter(x, y, has_cmap = [], ax=None, clabel = [],
                 cticks = [],
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):
    
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
        ax.scatter(x,y, **plt_kwargs)
    if cticks:
        cbar.remove()
        cbar = clabel[0].colorbar(mappable, ticks = cticks, pad = 1e-2)
        cbar.set_label(clabel[1], fontsize = 13)
        cbar.ax.tick_params(labelsize=10)
    if xlabel is not None:
        ax.set_xlabel(xlabel, family = 'Times New Roman', fontsize = '12')
    if ylabel is not None:
        ax.set_ylabel(ylabel, family = 'Times New Roman', fontsize = '12')
    if title is not None:
        ax.set_title(title)
    if setticks:
        x, xinterval, y, yinterval = setticks
        ax.xaxis.set_major_locator(MultipleLocator(x))
        ax.xaxis.set_minor_locator(AutoMinorLocator(xinterval))
        ax.yaxis.set_major_locator(MultipleLocator(y))
        ax.yaxis.set_minor_locator(AutoMinorLocator(yinterval))
        ax.tick_params(axis='both', which = 'major', direction = 'in',length = 6, width = 1)
        ax.tick_params(axis='both', which = 'minor', direction = 'in',length = 4, width =1)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)  
    return(ax)

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
fig, ax =  plt.subplots(1, 1, figsize = (10,10), dpi = 300)

# hide axis
# ax[0][1].get_yaxis().set_ticklabels([])

# =============================================================================
# Extinction
# =============================================================================

def isSingle(ind):
    ID = h2ID_n[ind]
    # This function checks if the HII region is associated only with a single cluster
    return len(read_catalogue.get_h2_param(ID, 'assoc')) == 1

# m_ext_x = ext
# m_ext_y = [sig[2] for sig in sigma1]
s_ext_x = np.array([ legus[0] for ii, legus in enumerate(ext_LEGUS) if isSingle(ii)])
s_ext_y = np.array([ signals * 3.1 for ii, signals in enumerate(ext_SIGNALS) if isSingle(ii)])

s_ext_c = np.array([ s for ii, s in enumerate(sigma1) if isSingle(ii)])
m_ext_c = np.array([ s for ii, s in enumerate(sigma1) if not isSingle(ii)])

m_ext_x = np.array([ legus[0] for ii, legus in enumerate(ext_LEGUS) if not isSingle(ii)])
m_ext_y = np.array([ signals * 3.1 for ii, signals in enumerate(ext_SIGNALS) if not isSingle(ii)])


# A_V/ E(B-V)  = R_V ~ 3.1

plot_error(m_ext_x, m_ext_y,
             xlabel = "SLUG $A_{V, \star}$",             
             markersize = 8,
             ylabel = 'SIGNALS $A_{V, neb}$',
             fmt = 'yo',
              xerr = [
                  [errors(sig)[0] for ii, sig in enumerate(ext_LEGUS) if not isSingle(ii)]
                  ,
                  [errors(sig)[1] for ii, sig in enumerate(ext_LEGUS) if not isSingle(ii)]
                  ]
                  ,
             yerr = [ signals * 3.1 for ii, signals in enumerate(ext_SIGNALS_err) if not isSingle(ii)],
            zorder = 1,
            label  = 'multiple star cluster',
            linewidth = 1,
            setticks = [.5, 5, .5, 5],
            mec = 'k', 
            capsize = 2, 
            # c = m_ext_c,
            ecolor = 'grey',
            xlim = (0, 2),
            ylim = (0, 2)
              )

plot_scatter(m_ext_x, m_ext_y, 
             has_cmap = [m_ext_c, cmr.iceburn],
             s = 300,
             ec = 'k',
             zorder = 10,
             clabel = [fig, '$f_{esc}$'],
             vmin = -1,
             vmax = 1,
             )


plot_error(s_ext_x, s_ext_y,
             markersize = 8,
             marker = 's',
             linestyle = 'none',
              color = 'gray',
              xerr = [
                  [errors(sig)[0] for ii, sig in enumerate(ext_LEGUS) if isSingle(ii)]
                  ,
                  [errors(sig)[1] for ii, sig in enumerate(ext_LEGUS) if isSingle(ii)]
                  ]
                  ,
             yerr = [ signals * 3.1 for ii, signals in enumerate(ext_SIGNALS_err) if isSingle(ii)],
            zorder = 1,
            label  = 'single star cluster',
            linewidth = 1,
            mec = 'k',
            capsize = 2, 
            ecolor = 'grey',
              )

plot_scatter(s_ext_x, s_ext_y, 
             has_cmap = [s_ext_c, cmr.iceburn],
             s = 300,
             marker = 's',
             vmin = -1,
             vmax = 1,
             ec = 'k',
             zorder = 10,
             )

# plt.plot([0,2], [0,2], c = 'k', label = '1:1')
plt.plot([0, 2], [0, 2*2.27], c = 'k', linestyle = '--', label = '$A_{V, \rm neb}/A_{V, \star}=2.27$')
# plt.plot([0, 2], [0, 2*2], c = 'k', linestyle = '-', label = '$A_{V, \rm neb}/A_{V, \star}=2$')

# line fit simple
# a, b = np.polyfit(list(s_ext_x), list(s_ext_y), 1)
#add line of best fit to plot
# plt.plot(s_ext_x, a*s_ext_x+b, c = 'b', label = 'Simple linear fit')   




plt.legend()
plot_tools.save("cluster_EV_AV")

#%%


# include errroarbs

sigma1 = []
for f_esc in f_esc_pdf_n:
    median = np.percentile(f_esc, (50, 15.9, 84.1))
    sigma1.append(median)



# Here we plot the ratio against f_esc

m_ratio_x = m_ext_y/m_ext_x
s_ratio_x = s_ext_y/s_ext_x

# get uncertainty in a very wonky way by looking at the product from
# dividing largest value over smallest value.

m_std = ext_LEGUS - 


fig, ax =  plt.subplots(1, 1, figsize = (10,10), dpi = 300)

plot_error(m_ext_c, m_ratio_x,
             ylim = (0, 10),
             xlim = (-1, .8),
             xlabel = "$f_{esc}$",             
             markersize = 8,
             ylabel = '$\\frac{\\textsc{signals } A_{V, neb}}{\\textsc{slug } A_{V, \star}}$',
             fmt = 'yo',
              xerr = [
                  [errors(sig)[0] for ii, sig in enumerate(sigma1) if not isSingle(ii)]
                  ,
                  [errors(sig)[1] for ii, sig in enumerate(sigma1) if not isSingle(ii)]
                  ]
                  ,
            zorder = 1,
            label  = 'multiple star cluster',
            linewidth = 1,
            setticks = [.5, 5, 2, 4],
            mec = 'k', 
            capsize = 2, 
            ecolor = 'grey',
             )

plot_error(s_ext_c, s_ratio_x,
           markersize = 8,
             marker = 's',
             linestyle = 'none',
              color = 'gray',
              xerr = [
                  [errors(sig)[0] for ii, sig in enumerate(sigma1) if isSingle(ii)]
                  ,
                  [errors(sig)[1] for ii, sig in enumerate(sigma1) if isSingle(ii)]
                  ]
                  ,
            zorder = 1,
            label  = 'single star cluster',
            linewidth = 1,
            mec = 'k',
            capsize = 2, 
            ecolor = 'grey',
             )


plot_scatter(m_ext_c, m_ratio_x,
             s = 300,
             ec = 'k',
             zorder = 10,
             c = 'y',
             )

plot_scatter(s_ext_c, s_ratio_x,
             s = 300,
             marker = 's',
             c = 'gray',
             ec = 'k',
             zorder = 10,
             )



plt.axhline(2.27, 0, 10, linestyle = '--', zorder = 10, c = 'k', linewidth = 3,
            label = '$A_{V, \rm neb}/A_{V, \star}=2.27$'
            )

plt.legend()





































