#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 23:45:23 2022

@author: Jia Wei Teh

This script saves class 1, 2 HII regions into a .reg file. Only includes

"""

""" """

# libraries
from PyAstronomy.pyasl import coordsDegToSexa
import numpy as np
#--
import src.tools.read_catalogue as read_catalogue
from src.tools.create_combined_table import pc2arc

# Read in catalogue
path2save = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/"
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)
qtab, mtab, atab = np.load(path2save+"clusterslug_table_phi073_apn1_geneva_1234.npy", allow_pickle = True)
clusterdata = np.load(path2save+"GenKroupc1234LHa_QH0_noCons.npy", allow_pickle = True)
# Unravel (see fg5_LHa_vs_QH0.py for more details)
f_esc_pdf, f_esc_percentiles, qh0_total_pdf, qh0_percentile, LHa_log, h2ID, median_mass = list(zip(*clusterdata))
# path to region files
path2reg = r'/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/'

count = 0

with open(path2reg+"h2regionsEmpty.reg","w+") as h2regionsEmpty:
  
    # Loop through each HII region
    for h2 in h2_catalogue:
        # list of associated clusters
        h2class = read_catalogue.get_h2_param(h2, 'class')
        if h2class == 3 or h2class == 4:
            continue
        if not read_catalogue.is_BPTStarforming(h2):
            continue 
        scIDlist = read_catalogue.get_h2_param(h2, 'assoc')
        if not scIDlist.size:
            RA = read_catalogue.get_h2_param(h2, 'RA')
            DEC = read_catalogue.get_h2_param(h2, 'DEC')
            rad = pc2arc(read_catalogue.get_h2_param(h2, 'rad'))
                
            # write into file
            ra, dec = coordsDegToSexa(
                RA,
                DEC,
                fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f')
                ).split()
            shape = "circle"
            size = str(rad)+"\""
            colour = "colour=blue"
            width = "width=3"
            # dashlist = "dashlist=8 3"
            # dash = "dash=1"
            line = [shape, ra, dec, size, "#",colour, width,\
                    # dash, dashlist,\
                        "\n"]
            line = " ".join(line)
            h2regionsEmpty.write(line)
    
# done
h2regionsEmpty.close()
    
    
