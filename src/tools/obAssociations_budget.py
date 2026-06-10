#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 19:14:16 2022

@author: Jia Wei Teh

This script calculates the full ionisation budget in the FoV.
"""

# libraries
import numpy as np
from tqdm import tqdm
#--
import src.tools.read_catalogue as read_catalogue

# Read in catalogue
path2save = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/"
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)
clusterdata = np.load(path2save+"GenKroupc1234LHa_QH0_noCons.npy", allow_pickle = True)
qtab, mtab, atab = np.load(path2save+"clusterslug_table_phi073_apn1_geneva_1234.npy", allow_pickle = True)

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
def samplePDF(xvalues, pdf, niter):
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
    niter = niter
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

# =============================================================================
# # For star clusters,
# =============================================================================

# Do you want to run the analysis or use pre-run results?
# Warninng: re-running takes ~2hrs, plus the final result
# will differ slightly depending on how QH0 is being sampled. 
reRun = False
withClass4 = True
onlyIncludeCoincide = True

if reRun:
    # how many iterations per cluster?
    niter = int(1e5)
    # create an array with length of a row
    montecarlo_across = np.zeros(shape = niter)
    
    if withClass4:
        if onlyIncludeCoincide:
            # include only those that coincide with regions
            for ii, sc in enumerate(tqdm(sc_catalogue)):
                # get ID
                scID = read_catalogue.get_sc_param(sc, 'ID')
                assoc = read_catalogue.get_sc_param(sc, 'assoc')
                if len(assoc) > 0:
                    range_log, qPDF = get_pdf(int(scID), qtab)
                    # create montecarlo sampling and sum with previous ones
                    montecarlo_across = np.sum([montecarlo_across, 10**samplePDF(range_log, qPDF, niter)], axis = 0)
            # Save data for future use to save computation time.
            np.save(path2save+"obAssociations_budget_wClass4_coincide.npy", montecarlo_across)
        else:
            for ii, sc in enumerate(tqdm(sc_catalogue)):
                # get ID
                scID = read_catalogue.get_sc_param(sc, 'ID')
                range_log, qPDF = get_pdf(int(scID), qtab)
                # create montecarlo sampling and sum with previous ones
                montecarlo_across = np.sum([montecarlo_across, 10**samplePDF(range_log, qPDF, niter)], axis = 0)
            # Save data for future use to save computation time.
            np.save(path2save+"obAssociations_budget_wClass4.npy", montecarlo_across)
   
    elif not withClass4:
        for ii, sc in enumerate(tqdm(sc_catalogue)):
            # get ID
            scID = read_catalogue.get_sc_param(sc, 'ID')
            isClass4 = int(read_catalogue.get_sc_param(sc, 'class')) == 4
            if not isClass4:
                range_log, qPDF = get_pdf(int(scID), qtab)
                # create montecarlo sampling and sum with previous ones
                montecarlo_across = np.sum([montecarlo_across, 10**samplePDF(range_log, qPDF, niter)], axis = 0)
        # Save data for future use to save computation time.
        np.save(path2save+"obAssociations_budget.npy", montecarlo_across)
        
        
# If user chooses not to rerun the analysis and use pre-run data,
elif not reRun:
    if withClass4:
        if onlyIncludeCoincide:
            montecarlo_across = np.load(path2save+"obAssociations_budget_wClass4_coincide.npy", allow_pickle = True)
        else:
            montecarlo_across = np.load(path2save+"obAssociations_budget_wClass4.npy", allow_pickle = True)

    else:
        _, montecarlo_across = np.load(path2save+"obAssociations_budget.npy", allow_pickle = True)
        
# Final cleanup
qh0_percentile = np.log10(np.percentile(montecarlo_across,  (15.9, 50, 84.1)))
print(qh0_percentile)

# =============================================================================
# # For HII regions,
# =============================================================================
# initialise LHa
LHa = 0
# Loop through every region
for h2 in h2_catalogue:
    # avoid double computing
    newLHa = read_catalogue.get_h2_param(h2, 'LHa')
    if not np.isnan(newLHa):
        LHa += newLHa
# Finally, convert to ionising photons, f_g = .73
qHa = LHa *  7.31e11 / .73
print(np.log10(qHa))


# =============================================================================
# # the escape fraction
# =============================================================================
f_esc = (10**qh0_percentile - qHa)/(10**qh0_percentile)
print(f'The escape frasction has median value of {f_esc[1]}')




