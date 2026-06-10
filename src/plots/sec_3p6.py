#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 11:46:51 2023

@author: Jia Wei Teh

This script calculates the sum of all objects to check for OB associations.
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
from src.tools.stats import samplePDF
np.random.seed(0)  # reproducible Monte-Carlo (Patch 4: seeded)
import matplotlib.pyplot as plt
import cmasher as cmr
from tqdm import tqdm
#--
import src.tools.draw_FOV as draw_FOV
import src.tools.create_clusterslug_table as create_clusterslug_table
import src.tools.read_catalogue as read_catalogue
import src.tools.plot_tools as plot_tools
from src.tools.create_combined_table import a2d, pc2arc, point_inside_polygon, min_dist

# Read in catalogue
path2save = paths.DAT
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)
try:
    qtab, mtab, atab = np.load(path2save+"clusterslug_table_phi073_apn1_geneva_1234.npy", allow_pickle = True)
except:
    create_clusterslug_table.create_clusterslug_table(stellarTrack = 'geneva')
    
pdfVals = np.zeros(shape = int(1e5), dtype = np.ndarray)

# =============================================================================
# Main functions
# =============================================================================
# function to evaluate f_esc
def add_QH0_mcmc_all():
    
    global pdfVals
    
    """ Given ID list of star clusters and corresponding H II region LHa, this function 
    returns the convovled (total) PDF and other useful informations.
    Method: For each cluster group, perform 1e5 Monte Carlo realisations. Each
    realisation, sample one value from PDF of every cluster. The final distribution
    will be the convolved PDF of the cluster group.
    
    Input in the form of [ID], LHa """
    
    for sc in sc_catalogue:
        SCclass = read_catalogue.get_sc_param(sc, 'class')
        if read_catalogue.get_sc_param(sc, 'assoc').size:
            ID = read_catalogue.get_sc_param(sc, 'ID')
            range_log, qPDF = get_pdf(int(ID), qtab)
            montecarlo_across = 10**samplePDF(range_log, qPDF)
            pdfVals += montecarlo_across
    
    return pdfVals

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
# Warning: re-running takes ~50mins, plus the final result
# will differ slightly depending on how QH0 is being sampled. 
# =============================================================================
reRun = False

if reRun:
    # Loop through all star clusters
    qh0_total_pdf = add_QH0_mcmc_all()
    qh0_percentile = np.log10(np.percentile(qh0_total_pdf.astype(float),  (2.3, 15.9, 50, 84.1, 97.7)))
    # save
    np.savetxt(path2save+'qh0_total_pdf_all.npy', qh0_total_pdf, delimiter=',')
   
else:
    qh0_total_pdf = np.loadtxt(path2save+"qh0_total_pdf_all.npy")

# Field of view array
fov1, fov2 = draw_FOV.FOV() 
    
# Total LHa
total_LHa = 0

for h2 in h2_catalogue:
    h2ID = read_catalogue.get_h2_param(h2, 'ID')
    # There are several process of which we use to remove unwanted data points.
    # Step1: remove class 3/4 Hii regions. 
    h2Class = read_catalogue.get_h2_param(h2, 'class')
    # Step3: remove HII regions that are near the border (R<1.5R_norm)
    # 8
    h2Rad = read_catalogue.get_h2_param(h2, 'rad')
    h2RA = read_catalogue.get_h2_param(h2, 'RA')
    h2Dec = read_catalogue.get_h2_param(h2, 'DEC')
    h2Rad_deg = a2d(pc2arc(h2Rad))
    # check
    if point_inside_polygon(h2RA, h2Dec, fov1):
        if min_dist(fov1, (h2RA, h2Dec))/h2Rad_deg < 1.5:
            continue
    elif point_inside_polygon(h2RA, h2Dec, fov2):
        if min_dist(fov2, (h2RA, h2Dec))/h2Rad_deg < 1.5:
            continue
    # Step4: remove bad LHa values (showing up as either 0 or np.nan) .
    # This is just a safety check.
    LHa = read_catalogue.get_h2_param(h2, 'LHa')
    if LHa <= 0 or np.isnan(LHa):
        continue
    # checks BPT
    if is_BPTStarforming(h2ID):
        continue
    
    total_LHa += LHa

total_qHa = total_LHa *  7.31e11
f_esc = (np.array(qh0_total_pdf) - total_qHa)/np.array(qh0_total_pdf)
plt.hist(f_esc, bins = 100)

f_esc_percentiles = np.percentile(f_esc, (15.9, 50, 84.1))


        # # QHa corresponding to Lha (Kenicutt)
        # qHa = LHa *  7.31e11
        # # escape fraction PDF
        # f_esc_pdf = (np.array(qh0_total_pdf) - qHa)/np.array(qh0_total_pdf)
        # # percentile
        # f_esc_percentiles = np.percentile(f_esc_pdf, (15.9, 50, 84.1))
        # # return the full 1e5 montecarlo samples and the percentile value
        # return f_esc_pdf, f_esc_percentiles, qh0_total_pdf, qh0_percentile
