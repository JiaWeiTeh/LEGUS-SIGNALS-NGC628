#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 23:36:35 2023

@author: Jia Wei Teh

This script includes operations from scratch that computes Table A1. 
This one plots a geneva track with class 1, 2, 3 clusters only.
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
from src.tools.read_catalogue import is_BPTStarforming
from src.tools.stats import medianPDF, prob2pdf, samplePDF
np.random.seed(0)  # reproducible Monte-Carlo (Patch 4: seeded)
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from tqdm import tqdm
#--
import src.tools.draw_FOV as draw_FOV
import src.tools.read_catalogue as read_catalogue
import src.tools.plot_tools as plot_tools
from src.tools.create_combined_table import a2d, pc2arc, point_inside_polygon, min_dist

# Read in catalogue
path2save = paths.DAT
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)
qtab, mtab, atab = np.load(path2save+"clusterslug_table_phi073_apn1_geneva_1234.npy", allow_pickle = True)


# =============================================================================
# Main functions
# =============================================================================
# function to evaluate f_esc
def add_QH0_mcmc(scIDlist, LHa):
    """ Given ID list of star clusters and corresponding H II region LHa, this function 
    returns the convovled (total) PDF and other useful informations.
    Method: For each cluster group, perform 1e5 Monte Carlo realisations. Each
    realisation, sample one value from PDF of every cluster. The final distribution
    will be the convolved PDF of the cluster group.
    
    Input in the form of [ID], LHa """
    # create montecarlo sampling
    montecarlo_across = np.zeros(shape = len(scIDlist), dtype = np.ndarray)
    for ii, ID in enumerate(scIDlist):
        range_log, qPDF = get_pdf(int(ID), qtab)
        montecarlo_across[ii] = 10**samplePDF(range_log, qPDF)
    # Sum across for summation of qh0
    qh0_total_pdf = np.sum(montecarlo_across, axis = 0)
    # percentile
    qh0_percentile = np.log10(np.percentile(qh0_total_pdf,  (2.3, 15.9, 50, 84.1, 97.7)))
    # QHa corresponding to Lha (Kenicutt)
    qHa = LHa *  7.31e11
    # escape fraction PDF
    f_esc_pdf = (np.array(qh0_total_pdf) - qHa)/np.array(qh0_total_pdf)
    # percentile
    f_esc_percentiles = np.percentile(f_esc_pdf, (15.9, 50, 84.1))
    # return the full 1e5 montecarlo samples and the percentile value
    return f_esc_pdf, f_esc_percentiles, qh0_total_pdf, qh0_percentile

# function to evaluate median value of combined mass
def add_mass_mcmc(scIDlist):
    """Similar to above, this script takes in the ID of clusters and returns the PDF
    after Monte Carlo algorithm to visualise."""
    # create montecarlo samples
    montecarlo = np.zeros(shape = len(scIDlist), dtype = np.ndarray)
    for ii, ID in enumerate(scIDlist):
        range_log, mPDF = get_pdf(int(ID), mtab)
        montecarlo[ii] = 10**samplePDF(range_log, mPDF)
    # Sum across for summation of mass
    summass = np.sum(montecarlo, axis = 0)
    # return the median value for colorbar
    return np.percentile(np.log10(summass), 50)

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

# A function that checks if an HII region is 100% in star-forming region

# =============================================================================
# Do you want to run the analysis or use pre-run results?
# Warninng: re-running takes ~50mins, plus the final result
# will differ slightly depending on how QH0 is being sampled. 
# =============================================================================
reRun = False

if reRun:
    # =============================================================================
    # With functions defined above, we can now begin evaluating convolved QH0 distributions.
    # For every HII region, we compare between the LHa and QH0 distributions to 
    # evaluate f_esc.
    # =============================================================================
    
    # Pre-computation. First, we do some checks.
    # For each HII region, make sure associated cluster is only associated with itself;
    # for those that do, put them into this list.
    one2onelist = []
    # Loop through
    for h in h2_catalogue:
        associatedSC = read_catalogue.get_h2_param(h, 'assoc')
        # if region has SC
        if len(associatedSC) > 0: 
            onlyHasOneRegion = all([len(read_catalogue.get_sc_param(s, 'assoc')) == 1\
                                for s in associatedSC])
            if onlyHasOneRegion: 
                one2onelist.append(int(read_catalogue.get_h2_param(h, 'ID')))
    
    # Field of view array
    fov1, fov2 = draw_FOV.FOV() 
    # Set up counter 
    removed_34 = 0
    removed_FOV = 0
    removed_121 = 0
    removed_LHa = 0
    
    # Main computation
    clusterdata = []
    # Now, loop through the catalogue 
    for ii, h in enumerate(tqdm(h2_catalogue)):
        # Select only Hii regions with associated clusters to decrease computation time
        scIDlist = read_catalogue.get_h2_param(h, 'assoc')
        # if the list is not empty,
        if scIDlist.size:
            # There are several process of which we use to remove unwanted data points.
            # Step1: remove class 3/4 Hii regions. 
            # 849
            h2Class = read_catalogue.get_h2_param(h, 'class')
            if h2Class in [3, 4]:
                removed_34 += 1
                continue
            
            # Step2: For each HII region, make sure associated cluster is only associated with itself.
            h2ID = read_catalogue.get_h2_param(h, 'ID') 
            if int(h2ID) not in one2onelist:
                removed_121 += 1
                continue
            
            # Step3: remove HII regions that are near the border (R<1.5R_norm)
            # 8
            h2Rad = read_catalogue.get_h2_param(h, 'rad')
            h2RA = read_catalogue.get_h2_param(h, 'RA')
            h2Dec = read_catalogue.get_h2_param(h, 'DEC')
            h2Rad_deg = a2d(pc2arc(h2Rad))
            # check
            if point_inside_polygon(h2RA, h2Dec, fov1):
                if min_dist(fov1, (h2RA, h2Dec))/h2Rad_deg < 1.5:
                    removed_FOV += 1
                    continue
            elif point_inside_polygon(h2RA, h2Dec, fov2):
                if min_dist(fov2, (h2RA, h2Dec))/h2Rad_deg < 1.5:
                    removed_FOV += 1
                    continue
        
            # Step4: remove bad LHa values (showing up as either 0 or np.nan) .
            # This is just a safety check.
            LHa = read_catalogue.get_h2_param(h, 'LHa')
            if LHa <= 0 or np.isnan(LHa):
                removed_LHa += 1
                continue
            
            # Step5: remove HII regions that are not 100% in the star-forming 
            # BPT regions
            if not is_BPTStarforming(h):
                continue
            
            # Step6: Only interested in HII regions that are associated with 
            # class 1, 2 and 3. Remove any class 4.
            for scID in scIDlist:
                if read_catalogue.get_sc_param(scID, 'class') == 4:
                    scIDlist = scIDlist.tolist()
                    scIDlist.remove(scID)
                    scIDlist = np.array(scIDlist)
                
            # Now we are sure that the list contains only class 1, 2, and 3.
            if scIDlist.size:  
               # log
               LHa_log = np.log10(LHa)        
               # Step5: evaluate f_esc
               f_esc_pdf, f_esc_percentiles, qh0_total_pdf, qh0_percentile = add_QH0_mcmc(scIDlist, LHa)
               # evaluate mass
               median_mass = add_mass_mcmc(scIDlist)
           
               clusterdata.append((f_esc_pdf, 
                                  f_esc_percentiles, 
                                  qh0_total_pdf,
                                  qh0_percentile, 
                                  LHa_log,
                                  h2ID,
                                  median_mass
                                  ))

    # Save data for future use to save computation time.
    np.save(path2save+"GenKroupc123LHa_QH0_noCons.npy", clusterdata)
   
# If user chooses not to rerun the analysis and use pre-run data,
elif not reRun:
    clusterdata = np.load(path2save+"GenKroupc123LHa_QH0_noCons.npy", allow_pickle = True)

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
plot_tools.plot_studies()

# =============================================================================
# Now that we have double checked, we can update the list of datapoints
# which we will use to calculate f_Esc.
# =============================================================================

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

plot_tools.save("GenKroupc123Fesc")


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
    
    

# =============================================================================
# Now that we have double checked, we can update the list of datapoints
# which we will use to calculate f_esc. This is for Bayesian
# =============================================================================

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
    # create temporary array to store cluster details
    temp = np.zeros(shape = (n_esc)-1)
    
    for jj, f_esc in enumerate(f_esc_range):
        # QH0 corresponding to a given LHa and f_esc
        Qh0_condition = np.log10(7.31e11*LHa/(1-f_esc))
        # individual probability
        probability = np.interp(Qh0_condition, QH0_range, QH0_pdf)
        if probability == 0 or np.log(probability) < very_small:
            temp[jj] = very_small
        else:
            temp[jj] = np.log(probability)
    
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
                ylim = (0, 9),
                linewidth = 2,
                zorder = 100,
                setticks = [0.2, 4, 2, 4],
                label = "Overall population",
                )
# High LHa bins
plot_tools.plot_plot(f_esc_range, HaHigh_pdf, 
                c = 'b',
                alpha = 0.5,
                linewidth = 1,
              label = "High $L_{\\rm H\\alpha}$ bin",
                )
# Low LHa bins
plot_tools.plot_plot(f_esc_range, HaLow_pdf, 
               c = 'c',
               alpha = 0.5,        
              label = "Low $L_{\\rm H\\alpha}$ bin",
               linewidth = 1,
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

plot_tools.save("GenKroupc123Fesc_Bayesian")
