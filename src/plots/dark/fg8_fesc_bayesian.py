#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 18:26:12 2022

@author: Jia Wei Teh

This script calculates the f_esc across population of HII regions using 
Bayesian technique. See Figure 8.
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
from src.tools.stats import medianPDF, prob2pdf
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from tqdm import tqdm
#--
import src.tools.plot_tools as plot_tools

plt.style.use("dark_background")

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
polygon_data = [0, 24, 38, 42]
for idx in polygon_data:
    plot_tools.plot_scatter(xmed_n[idx], LHa_log_n[idx], 
                  s = 200,
                  )
plot_tools.plot_studies()


# =============================================================================
# Now that we have double checked, we can update the list of datapoints
# which we will use to calculate f_esc.
# =============================================================================
# Update list to remove points in polygons
h2ID_n = np.array([i for jj, i in enumerate(h2ID_n) if jj not in polygon_data])
f_esc_pdf_n = np.array([i for jj, i in enumerate(f_esc_pdf_n) if jj not in polygon_data])
qh0_total_pdf_n = np.array([i for jj, i in enumerate(qh0_total_pdf_n) if jj not in polygon_data])
LHa_log_n = np.array([i for jj, i in enumerate(LHa_log_n) if jj not in polygon_data])

# To calculate f_esc, the following steps are taken
# Step1: What is the range of f_esc?
n_esc = 1000
f_esc_range = np.linspace(-1,1,n_esc)[:-1] # remove f_esc = 0 to avoid division by 0.
# Any np.log(value) lesser than this value equals 0.
very_small = -100

# The middle value seperating bins
mid = 10**np.percentile(LHa_log_n, 5)

# Initialise bins
overall = np.zeros(shape = (len(LHa_log_n), n_esc-1), dtype = np.ndarray)
# High LHa bins
HaHigh = np.zeros(shape = (len(LHa_log_n), (n_esc)-1), dtype = np.ndarray)
# Low LHa bins
HaLow = np.zeros(shape = (len(LHa_log_n), (n_esc)-1), dtype = np.ndarray)

# Step2: for each HII reigon, 
for ii, (LHa, QH0_pdf) in enumerate(tqdm(zip(10**(LHa_log_n), qh0_total_pdf_n))):
    # turn into distribution of QH0 into log values
    logPDF = np.log10(QH0_pdf)
    # points where we want to do kernal density plot
    QH0_range = np.linspace(min(logPDF), max(logPDF))
    # scipy KDE
    kde_sp = gaussian_kde(logPDF)
    QH0_pdf = kde_sp.pdf(QH0_range)
    # probability sum to 1
    QH0_pdf = QH0_pdf/sum(QH0_pdf)
    # QH0 implied by each f_esc, and its probability under the (KDE) QH0 distribution
    Qh0_condition = np.log10(7.31e11*LHa/(1-f_esc_range))
    probability = np.interp(Qh0_condition, QH0_range, QH0_pdf)
    # log-probability, floored at `very_small` (also where probability == 0)
    with np.errstate(divide='ignore'):
        logp = np.log(probability)
    temp = np.where((probability == 0) | (logp < very_small), very_small, logp)
    
    overall[ii] = temp
    # initiate with nan so we can np.nansum afterwards
    HaHigh[ii] = np.array([np.nan]*(n_esc-1))
    HaLow[ii] = np.array([np.nan]*(n_esc-1))
    
    if LHa < mid:
       HaLow[ii] = temp
    else:
       HaHigh[ii] = temp

# sum across columns to produce final probability
overall = np.sum(overall, axis = 0)
HaHigh = np.nansum(HaHigh, axis = 0)
HaLow = np.nansum(HaLow, axis = 0)


overall_pdf = prob2pdf(f_esc_range, np.exp(overall.astype(float)))
HaHigh_pdf = prob2pdf(f_esc_range, np.exp(HaHigh.astype(float)))
HaLow_pdf = prob2pdf(f_esc_range, np.exp(HaLow.astype(float)))


# =============================================================================
# plot
# =============================================================================
fig, ax =  plt.subplots(1,1, figsize = (6,5), dpi = 210)
# Overall
plot_tools.plot_plot(f_esc_range, overall_pdf, 
                xlabel = "$f_{\\rm esc}$",
                ylabel = "PDF",
                c = 'w',
                xlim = (-0.2, 1),
                ylim = (0, 9),
                linewidth = 2,
                zorder = 100,
                setticks = [0.2, 4, 2, 4],
                label = "Overall population",
                )
# High LHa bins
plot_tools.plot_plot(f_esc_range, HaHigh_pdf, 
                c = 'orange',
                alpha = 0.5,
                linewidth = 2,
              label = "High $L_{\\rm H\\alpha}$ bin",
                )
# Low LHa bins
plot_tools.plot_plot(f_esc_range, HaLow_pdf, 
               c = 'yellow',
               alpha = 1,        
              label = "Low $L_{\\rm H\\alpha}$ bin",
               linewidth = 2,
               )
# miscellaneous
plt.legend(loc = 'upper right', 
            fontsize = 11,
            labelspacing = .5,
            frameon = False)
plt.vlines(0, 0, 20, linestyles = '--', color = 'w')

# =============================================================================
# Print plot statistics in Latex-readable table format
# =============================================================================


print("Here are the median and 1-sigma uncertainty values for Low-, High-, and Overall LHa bins")
print(plot_tools.latexReadable(*medianPDF(f_esc_range, HaLow_pdf))+"&"+\
        plot_tools.latexReadable(*medianPDF(f_esc_range, HaHigh_pdf))+"&"+\
            plot_tools.latexReadable(*medianPDF(f_esc_range, overall_pdf))
      )

plot_tools.save("GenKroupc1234Fesc_Bayesian_dark")
