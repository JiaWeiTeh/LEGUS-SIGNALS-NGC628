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
import matplotlib.pyplot as plt
from tqdm import tqdm
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
path2file = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/lib/SLUG2/cluster_lib/"
# read cluster properties
with fits.open(path2file+"1e5_SC_cluster"+stellarTracks+"prop.fits") as file:
    props_table = file[1].data
    props_header = file[1].header
# read cluster photometries
with fits.open(path2file+"1e5_SC_cluster"+stellarTracks+"phot.fits") as file:
    phots_table = file[1].data

# =============================================================================
# Loop through different values of V-band magnitude cut
# =============================================================================
# magcut values
magCUT = np.linspace(-12, -3, 49) 
# Initialise array to store orphan ratio
orphanratio = []

for cut in tqdm(magCUT):
    # Initialise list to store Class 1,2,3 sources
    class123 = []
    # distance modulus to NGC 628
    distance_modulus = 29.8
    
    # Loop through clusters
    for ii, phot in enumerate(phots_table):
       
        # UV, U, B, V, I correspond to the following bands:
        # WFC3/F275W
        # WFC3/F336W
        # ACS/F435W
        # ACS/F555W
        # ACS/F814W
        # Values are 90% completeness, quoted from Adamo 2017.
        
        if usePointing == "centre":
            UV = phot[35] <= (23.29 - distance_modulus)
            U = phot[36] <= (23.91- distance_modulus)
            B =  phot[40] <= (24.93 - distance_modulus)
            V = phot[41] <= (25.05 - distance_modulus)
            I = phot[42] <= (24.27 - distance_modulus)
        elif usePointing == "east":
            UV = phot[35] <= (23.38 - distance_modulus)
            U = phot[36] <= (23.48- distance_modulus)
            B =  phot[40] <= (25.26 - distance_modulus)
            V = phot[41] <= (25.22 - distance_modulus)
            I = phot[42] <= (24.22 - distance_modulus)        
    
        # Check the number of filters, nflt
        nflt = (np.array([UV, U, V, B, I])).sum()
        
        # cluster selection criteria 1
        if not (V and (B or I)):
            continue
        
        # cluster selection criteria 2
        if ((UV or U) and V and nflt >= 4):
            # V-band cut for visual inspection 
            if phot[41] <= cut:
                # store index
                class123.append(ii)
                
    # =============================================================================
    # Cluster categorisation
    # =============================================================================
    # Bin clusters into detectable by LEGUS and not. Then further bin into above/below (up/down)
    # of SIGNALS's LHa detection limit.
    upLEGUS, downLEGUS, upNOT, downNOT = [], [], [], []
    # SIGNALS detection limit?
    limSIGNALS = 10**35.65
    # ionising photon conversion
    def QH02LHa(QH0, A_V):
        """
        Parameters
        ----------
        QH0 : float
            from SLUG. in photon/s
        A_V : float
            from SLUG. Stellar extinction.
    
        Returns
        -------
        LHa_cor : Extinction corrected Halpha Luminosity.
                Escape fraction considered.
        """
        # 0.27 escaped
        # convert from QH0 to LHa (Lha = QH0/7.31e11) (Kenicutt 1998)
        LHa = QH0*(1-0.27)/7.31e11
        # Nebula extinction from stellar extinction
        Neb_AV = A_V * 2.27
        # Extinction correction  A_V = -2.5log(L_ex/L_0)
        LHa_cor = 10**(Neb_AV/(-2.5)) * LHa
        
        return LHa_cor
    
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
                
    orphanratio.append(len(upNOT)/len(upNOT+upLEGUS))
        
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




















