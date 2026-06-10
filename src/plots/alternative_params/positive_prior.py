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
from src.tools.stats import medianPDF, prob2pdf, samplePDF
np.random.seed(0)  # reproducible Monte-Carlo (Patch 4: seeded)
import matplotlib.pyplot as plt
import scipy.stats 
from scipy.stats import gaussian_kde
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


# Normal method

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
    # remove all negatives
    f_esc = np.array(f_esc[f_esc>=0])
    # get 2sigma, 1sigma and median
    percentiles = np.percentile(f_esc, (2.3, 15.9, 50, 84.1, 97.7))
    if LHa < mid:
        PerLow.append(percentiles)
        HaLow += f_esc.tolist()
    else:
        PerHigh.append(percentiles)
        HaHigh += f_esc.tolist()
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
f_esc_low, f_esc_med, f_esc_high = np.percentile(np.array(HaHigh
                                                          ), (16, 50, 84))

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
f_esc_low, f_esc_med, f_esc_high = np.percentile(np.array(HaLow), (16, 50, 84))

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

plot_tools.plot_hist(np.array(f_esc_pdf_n[f_esc_pdf_n>=0]).flatten(),
              bins = 100,
              range = [-3,1],
              histtype = 'step',
              linewidth = 2,
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
              alpha = .5,
              density = True,
               # label = 'LHa high = [%.2f, %.2f]'%(high, max(LHa_log_n)),
              label = "High $L_{\\rm H\\alpha}$ bin",
              color = 'b',
              ylim = (-1, 1),
              xlim = (0, 2),
              orientation = 'horizontal',
              xlabel = "PDF",
              setticks = [1, 5, .5, 5]
              )
plot_tools.plot_hist(np.array(HaLow).flatten(),
              bins = 100,
              range = [-3,1],
              histtype = 'step',
               # label = 'LHa low = [%.2f, %.2f]'%(min(LHa_log_n), low),
              label = "Low $L_{\\rm H\\alpha}$ bin",
              alpha = .5,
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

plot_tools.save("GenKroupc1234Fesc_positive")




print("Here are the median and 1-sigma uncertainty values for Low-, High-, and Overall LHa bins")

vals = np.percentile(np.array(HaLow), (16, 50, 84))
print('low')
print(plot_tools.latexReadable(*vals))
vals = np.percentile(np.array(HaHigh), (16, 50, 84))
print('high')
print(plot_tools.latexReadable(*vals))
flatpdf = np.array(f_esc_pdf_n).flatten()
positivepdf = flatpdf[flatpdf>=0]
vals = np.percentile(np.array(positivepdf), (16, 50, 84))
print('overall')
print(plot_tools.latexReadable(*vals))
    
    
    

# Bayesian method

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
for ii, (LHa, QH0_pdf, esc_pdf) in enumerate(tqdm(zip(10**(LHa_log_n), qh0_total_pdf_n, f_esc_pdf_n))):
    QH0_pdf = np.array(QH0_pdf[esc_pdf>=0])
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
                c = 'k',
                xlim = (-0.2, 1),
                ylim = (0, 12),
                linewidth = 3,
                alpha = .7,
                zorder = 100,
                setticks = [0.2, 4, 2, 4],
                label = "Overall population",
                )
# High LHa bins
plot_tools.plot_plot(f_esc_range, HaHigh_pdf, 
                c = 'b',
                linewidth = 2,
              label = "High $L_{\\rm H\\alpha}$ bin",
                )
# Low LHa bins
plot_tools.plot_plot(f_esc_range, HaLow_pdf, 
               c = 'c',
              label = "Low $L_{\\rm H\\alpha}$ bin",
               linewidth = 2,
               )
# miscellaneous
plt.legend(loc = 'upper right', 
            fontsize = 11,
            labelspacing = .5,
            frameon = False)
plt.vlines(0, 0, 20, linestyles = '--', color = 'k')

# =============================================================================
# Print plot statistics in Latex-readable table format
# =============================================================================


print("Here are the median and 1-sigma uncertainty values for Low-, High-, and Overall LHa bins")
print(plot_tools.latexReadable(*medianPDF(f_esc_range, HaLow_pdf))+"&"+\
        plot_tools.latexReadable(*medianPDF(f_esc_range, HaHigh_pdf))+"&"+\
            plot_tools.latexReadable(*medianPDF(f_esc_range, overall_pdf))
      )

plot_tools.save("GenKroupc1234Fesc_Bayesian_positive")

#%%


# Emission line

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
    # remove all negatives
    f_esc = np.array(f_esc[f_esc>=0])
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
    x = [np.percentile(f_esc[f_esc>=0], 50) for ID, f_esc in zip(h2ID_n, f_esc_pdf_n) if threshold(ID, snr)]
    y = [read_catalogue.get_h2_param(ID, line) for ID in h2ID_n if threshold(ID, snr)]
    # plot
    plot_tools.plot_error(
                    x,
                    y,
                    ax = ax,
                    # label = mapping[line],
                    ylabel = mapping[line],
                    xlim = (0, 1),
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
    f_bins_neg = [1, 0]
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
    
    # now set values
    ax.text(.53, .09, "$\\rho = %.2f$\n$p = %.4f$"%(rho, pval), 
            fontsize = 10,
            bbox=dict(facecolor='w', edgecolor='k', alpha = .8),
            zorder = 1000,
            transform=ax.transAxes)

    # control label size
    ax.yaxis.label.set_size(10)
    # control tick size
    ax.tick_params(axis='both', which='major', labelsize=10)
    # align ylabels
    ax.yaxis.set_label_coords(-.3, 0.5)
    
    
# housekeeping
plt.tight_layout(w_pad = 0)

plot_tools.save("lineratios_zero2one_positiveprior")




#%%

# cluster statistics

# =============================================================================
# Do you want to re-run the analysis?
# =============================================================================
        
reRun = False

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
    f_esc = np.array(f_esc[f_esc>=0])
    # get 2sigma, 1sigma and median
    percentiles = np.percentile(f_esc, (2.3, 15.9, 50, 84.1, 97.7))
    sigma1.append(percentiles)
    ext.append(read_catalogue.get_h2_param(ID, 'ext'))
    rGalac.append(read_catalogue.get_h2_param(ID, 'rGalac'))

if reRun:
    masslist = np.zeros(shape = len(h2ID_n), dtype = np.ndarray)
    agelist = np.zeros(shape = len(h2ID_n), dtype = np.ndarray)
    for ii, ID in enumerate(tqdm(h2ID_n)):
        masslist[ii] = give_sortedMass(ID)
        agelist[ii] = give_sortedAge(ID)
    np.save(path2save+"clusterMassMedians_positiveprior.npy", masslist)
    np.save(path2save+"clusterAgeMedians_positiveprior.npy", agelist)

else:
    masslist = np.load(path2save+"clusterMassMedians_positiveprior.npy", allow_pickle=True)
    agelist = np.load(path2save+"clusterAgeMedians_positiveprior.npy", allow_pickle=True)


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
def plot_error(x, y, c = [], ax=None, clabel = [],
                 cticks = [],
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):
    """ Plots function. Accepts different axes """
    if ax is None:
        ax = plt.gca()
    if c:
        mappable = ax.errorbar(x, y, c=c[0], cmap = c[1], **plt_kwargs)
        cbar = clabel[0].colorbar(mappable)
        cbar.set_label(clabel[1], fontsize = 10)
        cbar.ax.tick_params(labelsize=10)
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
          ylim = (0, 1),
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
           ylim = (0, 1),
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
               ylim = (0, 1),
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
               ylim = (0, 1),
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
 
ax[1][1].legend(loc = 'upper left', fontsize = 18, framealpha = 0.5)

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
plot_tools.save("cluster_statistics_postiveprior")
