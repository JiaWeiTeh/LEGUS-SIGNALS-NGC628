#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 23:13:35 2022

@author: Jia Wei Teh

This script studies the effect of V band magnitude cut on the number 
of orphan regions in SLUG's synthesised clusters. See figure 9.
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
# Cluster selection and orphan ratio vs V-band magnitude cut (vectorised)
# =============================================================================
# magcut values
magCUT = np.linspace(-12, -3, 49)
# distance modulus to NGC 628
distance_modulus = 29.8

# Per-cluster band magnitudes (Vega): UV=F275W, U=F336W, B=F435W, V=F555W, I=F814W
UVm = np.asarray(phots_table.field(35), float)
Um  = np.asarray(phots_table.field(36), float)
Bm  = np.asarray(phots_table.field(40), float)
Vm  = np.asarray(phots_table.field(41), float)
Im  = np.asarray(phots_table.field(42), float)
QH0 = np.asarray(phots_table.field(4),  float)
age = np.asarray(props_table.field(2),  float)
AV  = np.asarray(props_table.field(11), float)

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

# number of detected filters per cluster
nflt = UV.astype(int) + U.astype(int) + V.astype(int) + B.astype(int) + I.astype(int)
# selection (cut-independent): V and (B or I), and (UV or U) and V and nflt>=4
eligible = (V & (B | I)) & ((UV | U) & V & (nflt >= 4))

# SIGNALS Halpha detection limit; detectable + young (<=10 Myr) clusters
limSIGNALS = 10**35.65
LHa = QH02LHa(QH0, AV)
up_detectable = (age <= 10 * 10**6) & (LHa >= limSIGNALS)

# For each V-band cut: orphan ratio = detectable clusters missed by the cut / all detectable
orphanratio = []
for cut in magCUT:
    class123 = eligible & (Vm <= cut)
    upLEGUS = (up_detectable &  class123).sum()
    upNOT   = (up_detectable & ~class123).sum()
    orphanratio.append(upNOT / (upNOT + upLEGUS))
        
# =============================================================================
# Plot
# =============================================================================
fig, ax = plt.subplots(1,1, figsize = (6,5), dpi = 300)

plot_tools.plot_plot(magCUT, np.array(orphanratio)*100,
               xlabel = "$M_V$ cut (mag)",
               ylabel = "Percentage of orphan H \\textsc{ii} regions \nin synthesised clusters (\%)",
               setticks = [1, 5, 10, 5],
               c = 'k',
               xlim = (-12, -3),
               )

plot_tools.plot_scatter(-6,
                  orphanratio[np.where(np.array(magCUT) == -6)[0][0]]*100,
                  s = 70,
                  zorder = 10,
                  label = 'LEGUS magnitude cut',
                  facecolor = 'yellow',
                  edgecolor = 'k'
                      )

plt.legend(loc = 'lower left', 
        fontsize = 14,
        labelspacing = .5,
        frameon = False)

plot_tools.save('vmag_cut')
