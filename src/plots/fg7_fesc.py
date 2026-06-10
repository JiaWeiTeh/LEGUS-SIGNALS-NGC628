#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 14:13:13 2022

@author: Jia Wei Teh

This script calculates the f_esc of HII regions. See Figure 9.
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
#--
import src.tools.plot_tools as plot_tools

# Read in catalogue
path2save = paths.DAT
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
median_mass_n = np.array(median_mass)[(xerr_max+xerr_min)<lim]
qh0_total_pdf_n = np.array(qh0_total_pdf)[(xerr_max+xerr_min)<lim]
f_esc_pdf_n = np.array(f_esc_pdf)[(xerr_max+xerr_min)<lim]

# =============================================================================
# Additionally, we remove data points in the polygon regions.
# =============================================================================
# Here we plot them to visualise.
fig, ax =  plt.subplots(1,1, figsize = (6,5), dpi = 300)

plot_tools.plot_scatter(xmed_n, LHa_log_n, 
                  c = 'w',
                  xlabel = "Log $\\rm Q(H^0)$ (photon s$^{-1}$)",
                  ylabel = "Log $L_{\\rm H\\alpha}$ (erg s$^{-1}$)",
                   xlim = (46,52),
                  ylim = (35.4, 39.4),
                  edgecolor = 'k',
                  vmin = 2.7, 
                  vmax = 4.7,
                  zorder = 100,
                  marker = 'o',
                  s = 50,
                  # setticks = [1,5,0.5,5],
                  label = 'NGC 628 (This work)'
                  )
# Plot datapoints that are in polygon after manual inspection
polygon_data = [0, 25, 38, 42]
for idx in polygon_data:
    plot_tools.plot_scatter(xmed_n[idx], LHa_log_n[idx], 
                  s = 200,
                  )
plot_tools.plot_studies()

# =============================================================================
# Now that we have double checked, we can update the list of datapoints
# which we will use to calculate f_Esc.
# =============================================================================
# Update list to remove points in polygons
h2ID_n = np.array([i for jj, i in enumerate(h2ID_n) if jj not in polygon_data])
f_esc_pdf_n = np.array([i for jj, i in enumerate(f_esc_pdf_n) if jj not in polygon_data])
qh0_total_pdf_n = np.array([i for jj, i in enumerate(qh0_total_pdf_n) if jj not in polygon_data])
LHa_log_n = np.array([i for jj, i in enumerate(LHa_log_n) if jj not in polygon_data])

# Seperate into different bins
lower, mid, higher = np.percentile(LHa_log_n, (16, 50, 84))
HaHigh = []
HaLow = []
PerHigh = []
PerLow = []
PerOverall = []

for ii, (LHa, f_esc, ID) in enumerate(sorted(zip(LHa_log_n, f_esc_pdf_n, h2ID_n))):
    # get 2sigma, 1sigma and median
    percentiles = np.percentile(f_esc, (2.3, 15.9, 50, 84.1, 97.7))
    if LHa < mid:
        PerLow.append(percentiles)
        HaLow.append(f_esc)
    else:
        PerHigh.append(percentiles)
        HaHigh.append(f_esc)
    PerOverall.append(percentiles)

def errors(percentiles):
    """ return percentiles for plotting purposes"""
    err = []
    for percen in percentiles:
        mini2, mini, med, maxi, maxi2 = percen
        err_max = maxi - med
        err_min = med - mini
        err_max_sig2 = maxi2 - med
        err_min_sig2 = med - mini2
        err.append([err_min_sig2, err_min, err_max, err_max_sig2])
    return err

# Sort by luminosity
sortedLHa = np.array(sorted(LHa_log_n))

# =============================================================================
# Plot
# =============================================================================

left, width = 0.1, 0.5
bottom, height = 0.1, 0.65
spacing = 0.005
plt.figure(figsize = (6,5), dpi = 200)
rect_scatter = [left, bottom, width, height]
ax_scatter = plt.axes(rect_scatter)

plot_tools.plot_error(sortedLHa[:len(PerLow)],
             np.array(PerLow)[:,2],
             setticks = [0.5, 5, 0.5, 5],
             xlabel = "Log $L_{\\rm H\\alpha}$ (erg s$^{-1}$)",
             ylabel = "$f_{\\rm esc}$",             
              ylim = (-1, 1),
             xlim = (37.3, 38.8),
             markersize = 10,
             fmt = 'co',
             yerr = [np.array(errors(PerLow))[:,1], np.array(errors(PerLow))[:,2]],
            zorder = 1,
            linewidth = 1,
            mec = 'k',
            capsize = 2,
            ecolor = 'lightgrey',
              )

plot_tools.plot_error(sortedLHa[len(PerLow):],
             np.array(PerHigh)[:,2],
             markersize = 10,
             fmt = 'bo',
             yerr = [np.array(errors(PerHigh))[:,1], np.array(errors(PerHigh))[:,2]],
            zorder = 1,
            linewidth = 1,
            mec = 'k',
            capsize = 2,
            ecolor = 'lightgrey',
              )


# get values for f_esc
f_esc_low, f_esc_med, f_esc_high = np.percentile(np.array(HaHigh).flatten(), (16, 50, 84))

print(f'Median values for high LHa bins are {lower}, and for fesc {f_esc_low}, {f_esc_med}, {f_esc_high}')

# high Ha
plot_tools.plot_error(higher,
        f_esc_med,
        yerr = [ [f_esc_med - f_esc_low], [f_esc_high - f_esc_med]],
         markersize = 20,
         marker=  'o',
         ls = 'none',
         color  = 'yellow',
        linewidth = 2,
        mec = 'k',
        capsize = 2,
        ecolor = 'k',
        alpha = .8,
                  )

# get values for f_esc
f_esc_low, f_esc_med, f_esc_high = np.percentile(np.array(HaLow).flatten(), (16, 50, 84))

print(f'Median values for low LHa bins are {lower}, and for fesc {f_esc_low}, {f_esc_med}, {f_esc_high}')

# High Ha
plot_tools.plot_error(lower,
        f_esc_med,
        yerr = [ [f_esc_med - f_esc_low], [f_esc_high - f_esc_med]],
         markersize = 20,
         marker=  'o',
         ls = 'none',
         color  = 'yellow',
        linewidth = 2,
        mec = 'k',
        capsize = 2,
        ecolor = 'k',
        alpha = .8,
                  )

plt.hlines(0, 37, 40, linestyles = '--', color = 'k', zorder = 0)

# =============================================================================
# Adds histogram
# =============================================================================
# parameter for axes
rect_histy = [left + width + spacing, bottom, 0.3, height]
# create axes for histogram
ax_histy = plt.axes(rect_histy)

plot_tools.plot_hist(np.array(f_esc_pdf_n).flatten(),
              bins = 100,
              range = [-3,1],
              histtype = 'step',
              linewidth = 3,
              alpha = .7,
              density = True,
              zorder = 10,
              # label = 'Overall H II regions',
              label = "Overall population",
              color = 'k',
              orientation = 'horizontal')


plot_tools.plot_hist(np.array(HaHigh).flatten(),
              bins = 100,
              range = [-3,1],
              histtype = 'step',
              density = True,
               # label = 'LHa high = [%.2f, %.2f]'%(high, max(LHa_log_n)),
              label = "High $L_{\\rm H\\alpha}$ bin",
              color = 'b',
              ylim = (-1, 1),
              linewidth = 2,
              xlim = (0, 2),
              orientation = 'horizontal',
              xlabel = "PDF",
              setticks = [1, 5, .5, 5]
              )
plot_tools.plot_hist(np.array(HaLow).flatten(),
              bins = 100,
              range = [-3,1],
              histtype = 'step',
              linewidth = 2,
               # label = 'LHa low = [%.2f, %.2f]'%(min(LHa_log_n), low),
              label = "Low $L_{\\rm H\\alpha}$ bin",
              color = 'c',
              density = True,
              orientation = 'horizontal')

# set ticks and labels
# plt.gca().yaxis.tick_right()
plt.gca().yaxis.set_visible(True)
plt.gca().set_yticklabels([])
plt.hlines(0, 0, 2, linestyles = '--', color = 'k', zorder = 0)
plt.legend(loc = 'lower right', 
            fontsize = 10,
            labelspacing = .5,
            frameon = True)

plot_tools.save("GenKroupc1234Fesc")




print("Here are the median and 1-sigma uncertainty values for Low-, High-, and Overall LHa bins")

vals = np.percentile(np.array(HaLow).flatten(), (16, 50, 84))
print('low')
print(plot_tools.latexReadable(*vals))
vals = np.percentile(np.array(HaHigh).flatten(), (16, 50, 84))
print('high')
print(plot_tools.latexReadable(*vals))
vals = np.percentile(np.array(f_esc_pdf_n).flatten(), (16, 50, 84))
print('overall')
print(plot_tools.latexReadable(*vals))
    
    
    

