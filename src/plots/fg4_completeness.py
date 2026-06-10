#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 21:26:15 2022

@author: Jia Wei Teh

This script plots figure 4 in paper, showing the ratio of orphan regions
obtained using synthetic clusters from SLUG libraries. Contains also 
scripts for Section 3.4.
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
from src.tools.stats import QH02LHa
import matplotlib.pyplot as plt
from astropy.io import fits
#--
import src.tools.plot_tools as plot_tools
import src.tools.draw_FOV as draw_FOV

# Use centre or east pointing?
usePointing = ["centre", "east"]
usePointing = "centre"  # change this if needed

# Use padova or geneva tracks?
stellarTracks = {"padova": "_modc020", "geneva": "_"}
stellarTracks = stellarTracks["geneva"] # change this if needed

# Extract data from slug's cluster library
path2file = paths.SLUG_LIB
# read cluster properties
with fits.open(path2file+"1e5_SC_cluster"+stellarTracks+"prop.fits") as file:
    props_table = file[1].data
    props_header = file[1].header
# read cluster photometries
with fits.open(path2file+"1e5_SC_cluster"+stellarTracks+"phot.fits") as file:
    phots_table = file[1].data

# =============================================================================
# Cluster selection (vectorised): LEGUS criteria with M_V <= -6
# =============================================================================
# distance modulus to NGC 628
distance_modulus = 29.8
# per-cluster band magnitudes (Vega): UV=F275W, U=F336W, B=F435W, V=F555W, I=F814W
UVm = np.asarray(phots_table.field(35), float)
Um  = np.asarray(phots_table.field(36), float)
Bm  = np.asarray(phots_table.field(40), float)
Vm  = np.asarray(phots_table.field(41), float)
Im  = np.asarray(phots_table.field(42), float)
# 90% completeness limits (Adamo 2017), per pointing
if usePointing == "centre":
    UV = UVm <= (23.29 - distance_modulus)
    U  = Um  <= (23.91 - distance_modulus)
    B  = Bm  <= (24.93 - distance_modulus)
    V  = Vm  <= (25.05 - distance_modulus)
    I  = Im  <= (24.27 - distance_modulus)
elif usePointing == "east":
    UV = UVm <= (23.38 - distance_modulus)
    U  = Um  <= (23.48 - distance_modulus)
    B  = Bm  <= (25.26 - distance_modulus)
    V  = Vm  <= (25.22 - distance_modulus)
    I  = Im  <= (24.22 - distance_modulus)
nflt = UV.astype(int) + U.astype(int) + V.astype(int) + B.astype(int) + I.astype(int)
# criteria: V and (B or I); (UV or U) and V and nflt>=4; and M_V <= -6
eligible = (V & (B | I)) & ((UV | U) & V & (nflt >= 4)) & (Vm <= -6)
# indices of selected clusters (set -> O(1) membership in the categorisation loop)
class123 = set(np.flatnonzero(eligible).tolist())


# =============================================================================
# Cluster categorisation
# =============================================================================
# Bin clusters into detectable by LEGUS and not. Then further bin into above/below (up/down)
# of SIGNALS's LHa detection limit.
upLEGUS, downLEGUS, upNOT, downNOT = [], [], [], []
# SIGNALS detection limit?
limSIGNALS = 10**35.65
# ionising photon conversion

# Loop through clusters
for ii, (phot, prop) in enumerate(zip(phots_table, props_table)):
    QH0 = phot[4]
    mass = prop[6]
    age = prop[2]
    AV = prop[11]
    LHa = QH02LHa(QH0, AV)
    # only age <= 5Myr allowed
    if age > 10 * 10**6:
        continue
    # if the cluster is detected by LEGUS
    if ii in class123:
        # is it above LHa?
        if LHa >= limSIGNALS:
            upLEGUS.append([age, LHa])
        else:
            downLEGUS.append([age, LHa])
    else:
        if LHa >= limSIGNALS:
            upNOT.append([age, LHa])
        else:
            downNOT.append([age, LHa])
        
    # subset to avoid crowding
    if ii > 2e4:
        break
    
            
# make into array
upLEGUS = np.array(upLEGUS)
downLEGUS = np.array(downLEGUS)
upNOT = np.array(upNOT)
downNOT = np.array(downNOT)
        
# =============================================================================
# Cluster statistics
# =============================================================================
# orphan regions. detected regions with no star clusters
# ratio = (greys UP)/(total UP)
ratio = len(upNOT)/(len(upNOT)+len(upLEGUS))
print("orphan H II ratio: %.2f"%ratio)  

# orphan star clusters. detected star clusters with no region
# ratio = (blacks LOWER)/(total black)      
ratio = len(upLEGUS)/(len(downLEGUS)+len(upLEGUS))
print("orphan YSC ratio: %.2f"%ratio)      

# =============================================================================
# Plot
# =============================================================================
fig, ax = plt.subplots(1,1, figsize = (5,4), dpi = 300)
# scatter size
size = 100
# transparency
alpha = .9
# upNOT
# sample size
n = 200
# upNOT_red = np.array(random.sample(list(upNOT), n))
plot_tools.plot_scatter(
             (upNOT\
               [:, 0]\
                   /1e6
                  ), 
             (upNOT\
               [:, 1]\
                  ), 
                  ylabel = "$L_{\\rm H\\alpha}$ (erg s$^{-1}$)",
                  xlabel = "Age (Myr)",
                  s = size,
                  alpha = alpha,
                  c = 'grey',
                  label = 'Not classified by LEGUS',
                  marker = '.')
# downNOT
# sample size
n = 10000
# downNOT_red = np.array(random.sample(list(downNOT), n))
plot_tools.plot_scatter(
             (downNOT\
               [:, 0]\
                   /1e6
                  ), 
             (downNOT\
               [:, 1]\
                  ), 
                  alpha = alpha,
                  s = size,
                  c = 'grey',
                  marker = '.')
# upLEGUS
# sample size
n = 250
# upLEGUS_red = np.array(random.sample(list(upLEGUS), n))
plot_tools.plot_scatter(
             (upLEGUS\
               [:, 0]\
                   /1e6
                  ), 
             (upLEGUS\
               [:, 1]\
                  ), 
                  s = size,
                   edgecolor = 'k',
                  alpha = alpha,
                  c = 'b',
                  label = 'Classified by LEGUS',
                  marker = '.')
# downLEGUS
# sample size
n = 300
# downLEGUS_red = np.array(random.sample(list(downLEGUS), n))
plot_tools.plot_scatter(
             (downLEGUS\
               [:, 0]\
                   /1e6
                  ), 
             (downLEGUS\
               [:, 1]\
                  ), 
                  s = size,
                  c = 'b',
                  edgecolor = 'k',
                  alpha = alpha,
                  marker = '.')
# LHa detection limit
plt.hlines((limSIGNALS), 
            0, 1e7, 
            color = 'k',
            label = 'SIGNALS detection threshold ($L_{\\rm H\\alpha}$)')
# miscellaneous
plt.yscale('log')
plt.ylim(1e32,1e41)
plt.xlim(.1, 10)
plt.legend(loc = 'upper right', 
            fontsize = 10,
            labelspacing = .5,
            frameon = False)
plot_tools.save('SLUG_completeness')


#%%

# =============================================================================
# Bonus section for 3.4 (true clusters in CLASS 4 sources)
# =============================================================================
# libraries
# Qeustion 1: What is the number of orphan HII regions with class 4 in them?
import src.tools.read_catalogue as read_catalogue
import src.tools.plot_tools as plot_tools
from src.tools.create_combined_table import d2r, d2a, pc2arc, arc2pc, ang_dist

# Read in catalogue
path2save = paths.DAT
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)

# distance of orphan to the nearest star cluster
regionwithclass4list = []

for h2 in h2_catalogue:
    # do not include class 3/4 clusters
    h2class = read_catalogue.get_h2_param(h2, 'class', h2_catalogue)
    if h2class == 4 or h2class == 3 :
        continue
    if not read_catalogue.is_BPTStarforming(h2):
        continue
    scAssoc_list = read_catalogue.get_h2_param(h2, 'assoc', h2_catalogue)
    containsOnlyClass4 = True
    for sc in scAssoc_list:
        if read_catalogue.get_sc_param(sc, 'class', sc_catalogue) != 4:
            containsOnlyClass4 = False
    # now these are orphan regions
    if containsOnlyClass4:
        for sc in scAssoc_list:
            if read_catalogue.get_sc_param(sc, 'class', sc_catalogue) == 4:
                regionwithclass4list.append(read_catalogue.get_h2_param(h2, 'ID'))
                break

print(f"{len(regionwithclass4list)} orphan HII regions contain class 4 sources")

#%%

# Question 2: What is the total area of orphan HII region?

def givePolyArea(xvals, yvals):
    # return area of polygon given x and y arrays
    return 0.5*np.abs(np.dot(xvals,np.roll(yvals,1))-np.dot(yvals,np.roll(xvals,1)))

def giveCircleArea(r): 
    # area of circle given radius
    return np.pi * r**2

# Area of FoV in arcseconds^2
FoV1, FoV2 = draw_FOV.FOV()
AreaFoV1 = d2a(d2a(givePolyArea(*zip(*FoV1))))
AreaFoV2 = d2a(d2a(givePolyArea(*zip(*FoV2))))
h2regionsarea = []

# area of orphan HII regions
for h2 in h2_catalogue:
    # do not include class 3/4 clusters
    h2class = read_catalogue.get_h2_param(h2, 'class', h2_catalogue)
    if h2class == 4 or h2class == 3 :
        continue
    if not read_catalogue.is_BPTStarforming(h2):
        continue
    scAssoc_list = read_catalogue.get_h2_param(h2, 'assoc', h2_catalogue)
    containsOnlyClass4 = True
    for sc in scAssoc_list:
        if read_catalogue.get_sc_param(sc, 'class', sc_catalogue) != 4:
            containsOnlyClass4 = False
    # now these are orphan regions
    if containsOnlyClass4:
        for sc in scAssoc_list:
            if read_catalogue.get_sc_param(sc, 'class', sc_catalogue) == 4:
                regionwithclass4list.append(read_catalogue.get_h2_param(h2, 'ID'))
                break
        # radius
        h2rad = pc2arc(read_catalogue.get_h2_param(h2, 'rad', h2_catalogue))
        # area
        h2regionsarea.append(giveCircleArea(h2rad))
    
AreaOrphan = sum(h2regionsarea)
probOverlap = AreaOrphan/(AreaFoV1+AreaFoV2)

print(f"The fraction of orphan HII region area is {probOverlap}")

# What is the number of Class 4 objects?
NumberofClass4 = 0
for sc in sc_catalogue:
    if read_catalogue.get_sc_param(sc, 'class') == 4:
        NumberofClass4 += 1
        
print(f"The expected number of Class 4 is {NumberofClass4*probOverlap}.")

# The probablility can be calculated via https://www.gigacalculator.com/calculators/binomial-probability-calculator.php
# with n = 637, x = 36
