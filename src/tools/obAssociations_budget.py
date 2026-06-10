#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 19:14:16 2022

@author: Jia Wei Teh

This script calculates the full ionisation budget in the FoV.
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
from tqdm import tqdm
#--
import src.tools.read_catalogue as read_catalogue

# Read in catalogue
path2save = paths.DAT
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)
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

# =============================================================================
# # For star clusters,
# =============================================================================

# Do you want to run the analysis or use pre-run results?
# Note: re-running is fast (~20 s) now that samplePDF is vectorised -- it used to
# take ~2 hrs. With the np.random.seed(0) set above the result is reproducible
# run to run, and matches the stored pre-run file to ~0.002 dex in the QH0
# percentiles (the stored file predates the seeding, hence the small offset).
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
        # 1-D array, matching np.save above and the two sibling branches (the stray
        # `_,` unpack assumed a 2-element file and crashed on the saved 1-D array)
        montecarlo_across = np.load(path2save+"obAssociations_budget.npy", allow_pickle = True)
        
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
