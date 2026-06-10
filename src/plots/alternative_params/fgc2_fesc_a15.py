#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 14:28:51 2022

@author: Jia Wei Teh

This script plots the LHa of H II regions against QH0 of star clusters .
This one includes Class 1,2,3,4 Clusters. Includes clusters where r_dist < 1.5 r_h2. 

Because this is basically a re-run, this script will create the f_esc plot
from scratch. This means that the 
"""

# libraries
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import gaussian_kde
from matplotlib.patches import Polygon
import cmasher as cmr
#--
import src.tools.create_LEGUS_table as create_LEGUS_table
import src.tools.create_SIGNALS_table as create_SIGNALS_table
import src.tools.draw_FOV as draw_FOV
import src.tools.read_catalogue as read_catalogue
import src.tools.plot_tools as plot_tools
from src.tools.create_combined_table import a2d, d2r, ang_dist, pc2arc, point_inside_polygon, min_dist
from src.plots.fg8_fesc_bayesian import h2ID_n

# Step 1:
# =============================================================================
# Create table that includes which cluster links to hii region and vice versa.
# =============================================================================
# What is the scaling factor a?
# a = 1.5 -> R_effective = R_hii * a
a = 1.5
# Path to save data for future references.
path2save = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/"

# =============================================================================
# Do you want to run the analysis or use pre-run results?
# Warninng: re-running takes ~50mins, plus the final result
# will differ slightly depending on how QH0 is being sampled. 
# =============================================================================
reRun = False

def create_combined_table():
    """
    Returns
    -------
    sc_cat : np.array
        star cluster catalogue.
    h2_cat : np.array
        hii region catalogue.
    """
    
    # For h2 catalogue
    # create copy of catalogue to modify
    h2_catalogue_unmod = create_SIGNALS_table.create_SIGNALS_table()
    rowdim, coldim = np.array(h2_catalogue_unmod).shape
    h2_cat = np.zeros(shape = (rowdim,coldim+2), dtype=np.ndarray)
    h2_cat[:,:-2] = h2_catalogue_unmod
        
    # For sc catalogue
    # create copy of catalogue to modify
    sc_catalogue_unmod = create_LEGUS_table.create_LEGUS_table()
    rowdim, coldim = np.array(sc_catalogue_unmod).shape
    sc_cat = np.zeros(shape = (rowdim,coldim+2), dtype=np.ndarray)
    sc_cat[:,:-2] = sc_catalogue_unmod

    # For H II regions, find out which star clusters are they overlapping with.
    for h in h2_cat:
        hID = h[0]
        p1 = (d2r(h[1]),d2r(h[2]))
        radius = pc2arc(h[12])
        sc_ID_in_this_h2, sc_d_in_this_h2 = [], []
        for s in sc_catalogue_unmod:
            scID = s[0]
            p2 = (d2r(s[3]),d2r(s[4]))
            dist = ang_dist(p1, p2)
            if (dist <= a * radius):
                sc_d_in_this_h2.append(dist)
                sc_ID_in_this_h2.append(scID)
        # record distance from each SC to H II region
        h[-2] = np.array(sc_d_in_this_h2)
        # record the list of ID of SC in this H II region
        h[-1] = np.array(sc_ID_in_this_h2)
    
    # For star clusters, find out which H II regions are they overlapping with.
    for s in sc_cat:
        scID = s[0]
        p2 = (d2r(s[3]),d2r(s[4]))
        h2_ID_in_this_sc, h2_d_in_this_sc = [], []
        for h in h2_cat:
            radius = pc2arc(h[12])
            hID = h[0]
            p1 = (d2r(h[1]),d2r(h[2]))
            dist = ang_dist(p1, p2)
            if (dist <= a * radius):
                h2_ID_in_this_sc.append(hID)
                h2_d_in_this_sc.append(dist)
        # record distance from each H II region to SC
        s[-2] = np.array(h2_d_in_this_sc)
        # record the list of ID of H II region in this SC
        s[-1] = np.array(h2_ID_in_this_sc)   
        
    # save table to avoid computing time in the future
    # store as an explicit 2-element object array (robust to the two catalogues
    # having equal row counts, which np.array([...], dtype=object) cannot handle)
    combined = np.empty(2, dtype=object)
    combined[0], combined[1] = sc_cat, h2_cat
    np.save(path2save+"combined_catalogue_a15.npy", combined)
   
    # return
    return sc_cat, h2_cat


# Step 2: Compute f_esc
# Read in catalogue
if reRun:
    print("Computing...")
    sc_catalogue, h2_catalogue  = create_combined_table()
else:
    sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue_a15.npy", allow_pickle = True)
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
            # we only care about HII regions that made it through the truncated
            # sample.
            h2ID = read_catalogue.get_h2_param(h, 'ID')
            if h2ID not in h2ID_n:
                continue
        
            # Step4: remove bad LHa values (showing up as either 0 or np.nan) .
            # This is just a safety check.
            LHa = read_catalogue.get_h2_param(h, 'LHa')
            if LHa <= 0 or np.isnan(LHa):
                removed_LHa += 1
                continue
            
            # Step5: remove HII regions that are not 100% in the star-forming 
            # BPT regions
            if not read_catalogue.is_BPTStarforming(h):
                continue
            
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
    np.save(path2save+"GenKroupc1234LHa_QH0_noCons_a15.npy", clusterdata)
   
# If user chooses not to rerun the analysis and use pre-run data,
elif not reRun:
    clusterdata = np.load(path2save+"GenKroupc1234LHa_QH0_noCons_a15.npy", allow_pickle = True)

# Unravel
f_esc_pdf, f_esc_percentiles, qh0_total_pdf, qh0_percentile, LHa_log, h2ID, median_mass = list(zip(*clusterdata))
    
# =============================================================================
# Plot
# =============================================================================
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

# figure
fig, ax =  plt.subplots(1,1, figsize = (6,5), dpi = 300)

# Plot datapoints
plot_tools.plot_scatter(xmed, LHa_log, 
                  has_cmap = [median_mass, cmr.chroma_r], 
                  clabel = [fig, "Log Cluster total mass ($M_\odot$)"],
                  xlabel = "Log $\\rm Q(H^0)$ (photon s$^{-1}$)",
                  ylabel = "Log $L_{\\rm H\\alpha}$ (erg s$^{-1}$)",
                  xlim = (46,52),
                  ylim = (35.4, 39.4),
                  vmin = 2.7, 
                  vmax = 4.7,
                  edgecolor = 'k',
                  alpha = 1,
                  zorder = 100,
                  marker = 'o',
                  s = 50,
                  setticks = [1,5,0.5,5],
                  label = 'NGC 628 (This work)'
                  )

# Number of errorbars to plot selectively?
# This is to avoid crowding. Set large number (i.e., 500) to show all errorbars
n_err = 500

# Plot errorbars
ax.errorbar(xmed[:n_err], LHa_log[:n_err], 
                xerr = [xerr_min[:n_err], xerr_max[:n_err]],
                fmt = ' ',
                markersize = 5,
                alpha = 0.8,
                zorder = 1,
                linewidth = 2,
                mec = 'k',
                capsize = 2,
                ecolor = 'darkgrey',
                )

# miscellaneous
plt.text(47, 39, "Forbidden Zone", fontstyle = "italic", size = 15)   
plt.text(47.4, 36.1, "$f_{\\rm esc}=0$", rotation = 55, size = 17)
# Plot results from other studies
plot_tools.plot_studies()

legend = plt.legend(loc = 'lower right', 
           fontsize = 9,
           labelspacing = .5,
           framealpha = 1,
           frameon = True)

plot_tools.save("GenKroupc1234LHa_QH0_noCons_a15")

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
fig, ax =  plt.subplots(1,1, figsize = (6,5), dpi = 300)

# Plot datapoints
plot_tools.plot_scatter(xmed_n, LHa_log_n, 
                  has_cmap = [(median_mass_n), cmr.chroma_r], 
                  clabel = [fig, "Log Cluster total mass ($M_\odot$)"],
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
plt.text(47, 39, "Forbidden Zone", fontstyle = "italic", size = 15)   
plt.text(47.4, 36.1, "$f_{\\rm esc}=0$", rotation = 55, size = 17) 
plt.legend(loc = 'lower right', 
           fontsize = 9,
           labelspacing = .5,
           frameon = True)

plot_tools.save("GenKroupc1234LHa_QH0_wCons_a15")



# Seperate into different bins
mid = np.percentile(LHa_log_n, 50)
HaHigh = []
HaLow = []

for ii, (LHa, f_esc, ID) in enumerate(sorted(zip(LHa_log_n, f_esc_pdf_n, h2ID_n))):
    # get 2sigma, 1sigma and median
    percentiles = np.percentile(f_esc, (2.3, 15.9, 50, 84.1, 97.7))
    if LHa < mid:
        HaLow.append(f_esc)
    else:
        HaHigh.append(f_esc)

# =============================================================================
# grab some values for non-bayesian methods, for printing out in a table
# =============================================================================
tot = np.array(f_esc_pdf_n).flatten()
high = np.array(HaHigh).flatten()
low = np.array(HaLow).flatten()

print("Here are the median and 1-sigma uncertainty values for Low-, High-, and Overall LHa bins for non-bayesian")
print(plot_tools.latexReadable(*np.percentile(low, (15.9, 50, 84.1)))+"&"+\
        plot_tools.latexReadable(*np.percentile(high, (15.9, 50, 84.1)))+"&"+\
            plot_tools.latexReadable(*np.percentile(tot, (15.9, 50, 84.1)))
      )
    
# =============================================================================
# We can now calculate the escape fraction
# =============================================================================

# This is the bayesian method
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

def prob2pdf(xvalues, pdf):
    """Converts probability function into probability density function"""
    # sum probability
    sumprob = sum(abs(xvalues[1]-xvalues[0])*(pdf[1:]+pdf[:-1])/2)
    # normalisation constant
    norm = 1/sumprob
    # return
    return pdf * norm

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

def medianPDF(xvalues, pdf):
    """Given xvalues and a pdf, return median and 1-sigma uncertainty """
    pdfsum = np.cumsum(pdf)*(xvalues[1]-xvalues[0])
    percentiles = np.array([
        xvalues[np.argmax(np.greater(pdfsum, 0.159))],
        xvalues[np.argmax(np.greater(pdfsum, 0.5))],
        xvalues[np.argmax(np.greater(pdfsum, 0.841))]])
    return percentiles

print("Here are the median and 1-sigma uncertainty values for Low-, High-, and Overall LHa bins")
print(plot_tools.latexReadable(*medianPDF(f_esc_range, HaLow_pdf))+"&"+\
        plot_tools.latexReadable(*medianPDF(f_esc_range, HaHigh_pdf))+"&"+\
            plot_tools.latexReadable(*medianPDF(f_esc_range, overall_pdf))
      )

plot_tools.save("GenKroupc1234Fesc_Bayesian_a15")






