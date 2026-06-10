#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 21:39:12 2022

@author: Jia Wei Teh

This script contains function that plots individual distributions of QH0. 
See Appendix figure B1.
"""
# libraries
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
#--
import src.tools.plot_tools as plot_tools

# Read in data
path2save = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/"
clusterdata = np.load(path2save+"GenKroupc1234LHa_QH0_noCons.npy", allow_pickle = True)

# Unravel
_, _, qh0_total_pdf, _, _, _, _ = list(zip(*clusterdata))

# =============================================================================
# With data of 1e5 samples of QH0 value per HII region, we convert them
# into a continuous gaussian-like plot using kde estimates.
# =============================================================================
# function that draws a continuous pdf from histogram of 1e5 samples.
def create_pdf(prePDF):
    # turn into log values
    logPDF = np.log10(prePDF)
    # points where we want to do kernal density plot
    QH0_range = np.linspace(min(logPDF), max(logPDF))
    # scipy KDE
    kde_sp = gaussian_kde(logPDF)
    QH0_pdf = kde_sp.pdf(QH0_range)
    # probability sum to 1
    QH0_pdf = QH0_pdf/sum(QH0_pdf)
    # return
    return QH0_range, QH0_pdf

# =============================================================================
# Plot
# =============================================================================
fig, ax =  plt.subplots(1,3, figsize = (9, 3), dpi = 210)

# plot 1 (sharp, single-peak)
plot_tools.plot_plot(*create_pdf(qh0_total_pdf[142]),
               c = 'b',
               ax = ax[0],
               xlim = [38, 51],
               ylim = [-.01, 0.16],
               setticks=[5, 5, .05, 5],
               xlabel = "Log $\\rm Q(H^0)$ (photon s$^{-1}$)",
               ylabel = "PDF",
               )
# add straight lines at y=0 on both sides
plot_tools.plot_plot([38, min(create_pdf(qh0_total_pdf[142])[0])],
               [0,0],
               c = 'b',
               ax = ax[0],
    )
plot_tools.plot_plot([max(create_pdf(qh0_total_pdf[142])[0]), 51],
               [0,0],
               c = 'b',
               ax = ax[0],
    )
ax[0].text(.1, 0.85, "Single-peaked",
           transform=ax[0].transAxes
         )

# plot 2 (multi-peak)
plot_tools.plot_plot(*create_pdf(qh0_total_pdf[119]),
               c = 'b',
               ax = ax[1],
               xlim = [38, 51],
               ylim = [-.01, 0.1],
               setticks=[5, 5, .05, 5],
               xlabel = "Log $\\rm Q(H^0)$ (photon s$^{-1}$)"
               )
# add straight lines at y=0 on both sides
plot_tools.plot_plot([38, min(create_pdf(qh0_total_pdf[119])[0])],
               [0,0],
               c = 'b',
               ax = ax[1],
    )
plot_tools.plot_plot([max(create_pdf(qh0_total_pdf[119])[0]), 51],
               [0,0],
               c = 'b',
               ax = ax[1],
    )
ax[1].text(.1, 0.85, "Multi-peaked",
           transform=ax[1].transAxes
         )

# plot 3 (broad, single-peak)
plot_tools.plot_plot(*create_pdf(qh0_total_pdf[83]),
               c = 'b',
               ax = ax[2],
               xlim = [38, 51],
               ylim = [-.01, 0.13],
               setticks=[5, 5, .05, 5],
               xlabel = "Log $\\rm Q(H^0)$ (photon s$^{-1}$)"
               )
# add straight lines at y=0 on both sides
plot_tools.plot_plot([38, min(create_pdf(qh0_total_pdf[83])[0])],
               [0,0],
               c = 'b',
               ax = ax[2],
    )
plot_tools.plot_plot([max(create_pdf(qh0_total_pdf[83])[0]), 51],
               [0,0],
               c = 'b',
               ax = ax[2],
    )
ax[2].text(.1, 0.77, "Single-peaked\n(broad)",
           transform=ax[2].transAxes
         )

# aesthetic
plt.tight_layout()
# save
plot_tools.save("visualise_QH0")



