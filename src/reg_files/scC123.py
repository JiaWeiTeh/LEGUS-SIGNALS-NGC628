#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 23:45:23 2022

@author: Jia Wei Teh

This script saves class 4 star clusters into .reg file
"""

""" """

# libraries
from PyAstronomy.pyasl import coordsDegToSexa
import numpy as np
#--
import src.tools.read_catalogue as read_catalogue

# Read in catalogue
path2save = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/"
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)
qtab, mtab, atab = np.load(path2save+"clusterslug_table_phi073_apn1_geneva_1234.npy", allow_pickle = True)
clusterdata = np.load(path2save+"GenKroupc1234LHa_QH0_noCons.npy", allow_pickle = True)
# Unravel (see fg5_LHa_vs_QH0.py for more details)
f_esc_pdf, f_esc_percentiles, qh0_total_pdf, qh0_percentile, LHa_log, h2ID, median_mass = list(zip(*clusterdata))
# path to region files
path2reg = r'/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/'

# list of associated clusters
assoc_sc_list = []

with open(path2reg+"scC123.reg","w+") as scC123:
                
    for sc in sc_catalogue:
        scClass = read_catalogue.get_sc_param(sc, 'class')
        if scClass == 4:
            continue
        scRA = read_catalogue.get_sc_param(sc, 'RA')
        scDEC = read_catalogue.get_sc_param(sc, 'DEC')
        # write into file
        ra, dec = coordsDegToSexa(
            float(scRA),
            float(scDEC),
            fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f')
            ).split()
        shape = "box"
        w = "0.5\""
        h = "0.5\""
        colour = "colour=yellow"
        width = "width=2"
        line = [shape, ra, dec, w, h, "#",colour,width,"\n"]
        line = " ".join(line)
        scC123.write(line)
        
# done
scC123.close()
    
    
