#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 22:31:42 2022

@author: Jia Wei Teh

This script Plots f_esc vs different emissions lines obtained from SIGNALS survey.
"""
# libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats 
#--
import src.tools.plot_tools as plot_tools
import src.tools.read_catalogue as read_catalogue

# Read in catalogue
path2save = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/"
_, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)
clusterdata = np.load(path2save+"GenKroupc1234LHa_QH0_noCons.npy", allow_pickle = True)

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
# which we will use to calculate f_esc. Then, we can start looking at the
# emission line ratios from SIGNALS
# =============================================================================
# Update list to remove points in polygons
h2ID_n = np.array([i for jj, i in enumerate(h2ID_n) if jj not in polygon_data])
f_esc_pdf_n = np.array([i for jj, i in enumerate(f_esc_pdf_n) if jj not in polygon_data])
qh0_total_pdf_n = np.array([i for jj, i in enumerate(qh0_total_pdf_n) if jj not in polygon_data])
LHa_log_n = np.array([i for jj, i in enumerate(LHa_log_n) if jj not in polygon_data])

# Useful function to compute the errorbars. This will be fed into plot_error.
def errors(percentiles):
    """ return percentiles for plotting purposes"""
    mini2, mini, med, maxi, maxi2 = percentiles
    err_max = maxi - med
    err_min = med - mini
    err_max_sig2 = maxi2 - med
    err_min_sig2 = med - mini2
    return err_min_sig2, err_min, err_max, err_max_sig2

# Set up a list to record errorbar values
percentiles = []
for ii, (LHa, f_esc, ID) in enumerate(zip(LHa_log_n, f_esc_pdf_n, h2ID_n)):
    # get 2sigma, 1sigma and median
    percentiles.append(np.percentile(f_esc, (2.3, 15.9, 50, 84.1, 97.7)))
    
# Names of SIGNALS emission line ratios, their errors and SNR. For future referencing.
emissions = [
            'N2Ha', 'N2Ha_err', 'N2Ha_SNR',\
            'S2Ha', 'S2Ha_err', 'S2Ha_SNR',\
            'S2N2', 'S2N2_err', 'S2N2_SNR', \
            'O3Hb','O3Hb_err','O3Hb_SNR',\
            'O2Hb', 'O2Hb_err', 'O2Hb_SNR',\
            'O23Hb', 'O23Hb_err', 'O23Hb_SNR',\
            'O3O2', 'O3O2_err', 'O3O2_SNR', \
            'O3N2', 'O3N2_err', 'O3N2_SNR', \
            'O2N2', 'O2N2_err', 'O2N2_SNR', \
            'S2S2', 'S2S2_err', 'S2S2_SNR', \
                ]
# full line ratio with wavelength
real_names = ["log [NII]6583/Ha", 
              "log [SII]6716+6731/Ha",
              "log [SII]6716+6731/[NII]6583",
              "log [OIII]5007/Hb",
              "log [OII]3727/Hb",
              "log ([OII]3727+[OII]5007)/Hb",
              "log [OIII]5007/[OII]3727",
              "log [OIII]5007/[NII]6583",
              "log [OII]3727/[NII]6583",
              "log [SII]6716/[SII]6731",
              ]
# shortened line ratio, in latex readable format.
lazy_name = ["$\log{\\rm [NII]/H\\alpha}$", 
              "$\log{\\rm [SII]/H\\alpha}$",
              "$\log{\\rm [SII]/[NII]}$",
              "$\log{\\rm [OIII]/H\\beta}$",
              "$\log{\\rm [OII]/H\\beta}$",
              "$\log{\\rm [OIII]+[OII]/H\\beta}$",
              "$\log{\\rm [OIII]/[OII]}$",
              "$\log{\\rm [OIII]/[NII]}$",
              "$\log{\\rm [OII]/[NII]}$",
              "$\log{\\rm [SII]/[SII]}$",
              ]
# Create dictionary
mapping = {}
for idx in range(int(len(emissions)/3)):
    # mapping[emissions[int(3*idx)]] = real_names[idx]
    mapping[emissions[int(3*idx)]] = lazy_name[idx]

# =============================================================================
# Plot
# =============================================================================
# labels
charlabel = [chr(i) for i in range(ord('a'),ord('i')+1)]
charlabel = ['a', 'b', 'c', 'd', 'e', 'skipped', 'f', 'g', 'h', 'i']

# some helpful functions
def threshold(ID, snr):
    # returns True/False depending on if SNR is over threshold
    return read_catalogue.get_h2_param(ID, snr) >= 3

def secondmax(vals):
    # returns second highest value for aesthetic purposes
    vals = [v for v in vals if not np.isnan(v)]
    return sorted(vals)[-2]

def secondmin(vals):
    # return second lowest value for aesthetic purposes
    vals = [v for v in vals if not np.isnan(v)]
    return sorted(vals)[1]

# initialise row/column index for subplots
ax_idx_x = 0
ax_idx_y = 0
# plot
fig, axs =  plt.subplots(3,3, figsize = (7,5), dpi = 210)
# loop through emission line list
for ii, idx in enumerate(range(int(len(emissions)/3))):
    line = emissions[3*idx]
    err = emissions[3*idx+1]
    snr = emissions[3*idx+2]
    # skip one of the lines
    if line == "O23Hb":
        continue
    # get axis
    ax = axs[ax_idx_x,ax_idx_y]
    # get x, y values
    x = [np.percentile(f_esc, 50) for ID, f_esc in zip(h2ID_n, f_esc_pdf_n) if threshold(ID, snr)]
    y = [read_catalogue.get_h2_param(ID, line) for ID in h2ID_n if threshold(ID, snr)]
    # plot
    plot_tools.plot_error(
                    x,
                    y,
                    ax = ax,
                    label = mapping[line],
                    ylabel = mapping[line],
                    xlim = (-1, 1),
                    ylim = (secondmin(y)-.2, secondmax(y)+.2),
                    markersize = 5,
                    marker=  'o',
                    ls = 'none',
                    color  = 'blue',
                    yerr = [read_catalogue.get_h2_param(ID, err) for ID in h2ID_n if threshold(ID, snr)],
                    xerr = [
                    [errors(sig)[1] for ID, sig in zip(h2ID_n, percentiles) if threshold(ID, snr)]
                    ,
                    [errors(sig)[2] for ID, sig in zip(h2ID_n, percentiles) if threshold(ID, snr)]
                    ]
                    ,
                    linewidth = 1,
                    mec = 'k',
                    capsize = 2,
                    ecolor = 'lightgrey',
                  )
    # if SNR < 3 exists, 
    try:
        plot_tools.plot_error(
                    [np.percentile(f_esc, 50) for ID, f_esc in zip(h2ID_n, f_esc_pdf_n) if not threshold(ID, snr)],
                    [read_catalogue.get_h2_param(ID, line) for ID in h2ID_n if not threshold(ID, snr)],
                    ax = ax,
                    markersize = 5,
                    marker=  'o',
                    ls = 'none',
                    color  = 'cyan',
                    yerr = [read_catalogue.get_h2_param(ID, err) for ID in h2ID_n if not threshold(ID, snr)],
                    xerr = [
                    [errors(sig)[1] for ID, sig in zip(h2ID_n, percentiles) if not threshold(ID, snr)]
                    ,
                    [errors(sig)[2] for ID, sig in zip(h2ID_n, percentiles) if not threshold(ID, snr)]
                    ]
                    ,
                    linewidth = 1,
                    mec = 'k',
                    capsize = 2,
                    ecolor = 'lightgrey',
                  )
    except: pass

    # =============================================================================
    # Now, we calculate the spearmans correlation.
    # =============================================================================
    # calculate correlation
    x = np.array(x)
    y = np.array(y)
        
    # calculate spearmans coefficient rho. remove nans.
    newx = x[np.logical_not(np.isnan(y))]
    newy = y[~(np.isnan(y))]
    
    # OR only include rho for f_esc greater than 0?
    # newx = x[np.logical_and(np.logical_not(np.isnan(y)), x > 0)]
    # newy = y[np.logical_and(~(np.isnan(y)), x > 0)]
    
    # OR force every negative values to 0?
    # newx = x[np.logical_not(np.isnan(y))]
    # newy = y[~(np.isnan(y))]
    # newx[newx < 0] = 0 
    
    # Statistics
    rho, pval = scipy.stats.spearmanr(newx, newy)
    print("%.2f, %.4f"%(rho, pval))
    # set alphabet
    ax.text(0.05, 0.85, charlabel[ii]+")", transform=ax.transAxes)
    
    # operations to update rows/columns
    ax_idx_y += 1
    if ax_idx_y >= 3:
        ax_idx_y = 0
        ax_idx_x += 1
    if ax_idx_x == 2:
        axs[ax_idx_x,ax_idx_y].set_xlabel("$f_{\\rm esc}$")

    # Plot the mean value for negative fesc
    f_bins_neg = [-1, 0]
    xavg, yavg = [], []
    # Get bin values to calculate mean
    for bins in range(len(f_bins_neg)-1):
        idx = np.where(
                (x > f_bins_neg[bins]) & (x <= f_bins_neg[bins+1])
            )[0]
        yavg.append(np.nanmean(y[idx]))
        xavg.append(np.nanmean(x[idx]))
    # plot
    plot_tools.plot_scatter(xavg, yavg,
                      s = 100,
                      c = 'lightgrey',
                      edgecolor = 'k',
                      zorder = 10,
                      alpha = 1,
                      ax = ax)
    
    # Plot the mean value for positive fesc
    f_bins_pos = np.linspace(0, 1, 5)
    xavg, yavg = [], []
    # Get bin values to calculate mean
    for bins in range(len(f_bins_pos)-1):
        idx = np.where(
                (x > f_bins_pos[bins]) & (x <= f_bins_pos[bins+1])
            )[0]
        # print(x)
        # print(idx)
        yavg.append(np.nanmean(y[idx]))
        xavg.append(np.nanmean(x[idx]))
    # plot
    plot_tools.plot_scatter(xavg, yavg,
                      s = 100,
                      c = 'yellow',
                      edgecolor = 'k',
                      zorder = 10,
                      alpha = 1,
                      ax = ax)
    
    # control label size
    ax.yaxis.label.set_size(10)
    # control tick size
    ax.tick_params(axis='both', which='major', labelsize=10)
    # align ylabels
    ax.yaxis.set_label_coords(-.3, 0.5)
    
    
# housekeeping
plt.tight_layout(w_pad = 0)

plot_tools.save("lineratios")


