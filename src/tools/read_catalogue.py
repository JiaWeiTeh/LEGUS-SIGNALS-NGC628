#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 16:40:20 2022

@author: Jia Wei Teh

This script contains function which help reading catalogues.
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
#--
# Read in catalogue
path2save = paths.DAT
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)

def get_row(ID, table):
    """ Given ID, return the properties of the H2 region or SC """
    try:
        for row in table:
            if row[0] == ID:
                return row
        # np.where (code below) is cleaner but slower 
        # return table[np.where(np.array(table)[:,0] == ID)[0][0]]
    except Exception:
        return None

def get_sc_param(ID, colname, sctable = sc_catalogue):
    """Given ID/properties and header, return value.
    ID can either be the objects ID, or the entire (full)
    row of the object in the catalogue."""
    
    try:
        column_names = np.array(['ID','x','y', 'RA', 'DEC',\
                           'UV', 'dUV', 'U', 'dU',\
                           'B', 'dB', 'V', 'dV', 'I', 'dI', \
                           'CI', \
                           'age', 'maxage', 'minage',\
                           'mass', 'minmass', 'maxmass',\
                           'ext', 'maxext', 'minext',\
                           'chi2UV','chi2U','chi2B', 'chi2V', 'chi2I',\
                           'reduchi', 'Q', 'nfil', 'class', 'mclass',\
                           'dist','assoc'
                               ])  
        # if given the properties of ID
        if hasattr(ID, "__len__") and not isinstance(ID, str):
            return ID[np.where(column_names == colname)[0][0]]
        else:
            # get properties from ID number
            row = get_row(ID, sctable)
            return row[np.where(column_names == colname)[0][0]]
    
    except TypeError:
        raise Exception("ID does not exist")
        return
    except IndexError:
        raise Exception("Column name not cound. Please ensure header is either"+
              " of the following: 'ID','x','y', 'RA', 'DEC',"+
                           "'UV', 'dUV', 'U', 'dU',"+
                          " 'B', 'dB', 'V', 'dV', 'I', 'dI', "+
                          " 'CI', "+
                          " 'age', 'maxage', 'minage',"+
                          " 'mass', 'minmass', 'maxmass',"+
                          " 'ext', 'maxext', 'minext',"+
                          " 'chi2UV','chi2U','chi2B', 'chi2V', 'chi2I',"+
                          " 'reduchi', 'Q', 'nfil', 'class', 'mclass',"+
                          " 'dist','assoc'")
        

def get_h2_param(ID, colname, h2table = h2_catalogue):
    """Given ID/properties and header, return value.
    ID can either be the objects ID, or the entire (full)
    row of the object in the catalogue."""    
    
    try:
        column_names = np.array(['ID','RA','DEC', 'rGalac',\
                                 'LHa', 'HaDIG', 'class',\
                                 'I0', 'Amp', 'sig', 'alpha', 'R2',\
                                 'rad', 'ext', 'ext_err',\
                                 'N2Ha', 'N2Ha_err', 'N2Ha_SNR',\
                                 'S2Ha', 'S2Ha_err', 'S2Ha_SNR',\
                                 'S2N2', 'S2N2_err', 'S2N2_SNR', \
                                 'O3Hb','O3Hb_err','O3Hb_SNR',\
                                 'O2Hb', 'O2Hb_err', 'O2Hb_SNR',\
                                 'O23Hb', 'O23Hb_err', 'O23Hb_SNR',\
                                 'O3O2', 'O3O2_err', 'O3O2_SNR', \
                                 'O3N2', 'O3N2_err', 'O3N2_SNR', \
                                 'O2N2', 'O2N2_err', 'O2N2_SNR', \
                                 'S2S2', 'S2S2_err', 'S2S2_SNR', \
                                 'dist','assoc'
                                 ])  
        
        # if already given the properties of ID
        if hasattr(ID, "__len__") and not isinstance(ID, str):
            return ID[np.where(column_names == colname)[0][0]]
        else:
            # get properties from ID number
            row = get_row(ID, h2table)
            return row[np.where(column_names == colname)[0][0]]
    
    except TypeError:
        raise Exception("ID does not exist")
        return
    except IndexError:
        raise Exception("Column name not cound. Please ensure header is either"+
              " of the following: 'ID','RA','DEC', 'rGalac',"+
                                 "'LHa', 'HaDIG', 'class',"+
                                 "'I0', 'Amp', 'sig', 'alpha', 'R2',"+
                                 "'rad', 'ext', 'ext_err',"+
                                 "'N2Ha', 'N2Ha_err', 'N2Ha_SNR',"+
                                 "'S2Ha', 'S2Ha_err', 'S2Ha_SNR',"+
                                 "'S2N2', 'S2N2_err', 'S2N2_SNR',"+
                                 "'O3Hb','O3Hb_err','O3Hb_SNR',"+
                                 "'O2Hb', 'O2Hb_err', 'O2Hb_SNR',"+
                                 "'O23Hb', 'O23Hb_err', 'O23Hb_SNR',"+
                                 "'O3O2', 'O3O2_err', 'O3O2_SNR', "+
                                 "'O3N2', 'O3N2_err', 'O3N2_SNR', "+
                                 "'O2N2', 'O2N2_err', 'O2N2_SNR', "+
                                 "'S2S2', 'S2S2_err', 'S2S2_SNR', "+
                                 "'dist','assoc'")

# A function that checks if an HII region is 100% in star-forming region
def is_BPTStarforming(h2):
    # first, extract line ratios
    s2Ha = get_h2_param(h2, 'S2Ha')    
    o3Hb = get_h2_param(h2, 'O3Hb')    
    n2Ha = get_h2_param(h2, 'N2Ha')       
    # equation for BPT diagrams
    # Kauffmann+03
    def BPT_H2(a):
        return 0.61 / (a-0.05)+ 1.3
    # Kewley+01
    def BPT_COMP(b):
        return 0.61  / (b-0.47)+ 1.19
    # Main S2 line
    def BPT_Sii_AGN(AGN):
        return (0.72 / (AGN - 0.32) + 1.30)
    
    return all([o3Hb < BPT_H2(n2Ha), o3Hb < BPT_COMP(n2Ha), o3Hb < BPT_Sii_AGN(s2Ha)])