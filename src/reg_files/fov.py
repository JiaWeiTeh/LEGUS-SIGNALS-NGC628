#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 23:38:21 2022

@author: Jia Wei Teh

This script creates .reg files for LEGUS FOV for visualisation in Ds9
"""

# libraries
from PyAstronomy.pyasl import coordsDegToSexa
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
import src.tools.draw_FOV as draw_FOV

# path to store .reg files
path2reg = paths.DAT
# create .reg file
reg_file = open(path2reg+"LEGUS_newfov.reg","w+")

# outlines
outline1, outline2 = draw_FOV.FOV()

# operations to put lines into.reg file
# Outline 1
for i in range(len(outline1)-1):
    if i == (len(outline1)-1):
        ra1, dec1 = coordsDegToSexa(
            float(outline1[-1][0]),
            float(outline1[-1][1]),
            fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f')
            ).split()
        ra2, dec2 = coordsDegToSexa(
            float(outline1[0][0]),
            float(outline1[0][1]),
            fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f')
            ).split()
    else:
        ra1, dec1 = coordsDegToSexa(
        float(outline1[i][0]),
        float(outline1[i][1]),
        fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f')
        ).split()
    ra2, dec2 = coordsDegToSexa(
        float(outline1[i+1][0]),
        float(outline1[i+1][1]),
        fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f')
        ).split()
    
    shape = "line"
    colour = "colour=black"
    width = "width=2"
    line = [shape, ra1, dec1, ra2, dec2, "#",colour,
            width,"\n"]
    line = " ".join(line)
    reg_file.write(line)
    
# Outline 2
for i in range(len(outline2)-1):
    if i == (len(outline2)-1):
        ra1, dec1 = coordsDegToSexa(
            float(outline2[-1][0]),
            float(outline2[-1][1]),
            fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f')
            ).split()
        ra2, dec2 = coordsDegToSexa(
            float(outline2[0][0]),
            float(outline2[0][1]),
            fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f')
            ).split()
    else:
        ra1, dec1 = coordsDegToSexa(
        float(outline2[i][0]),
        float(outline2[i][1]),
        fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f')
        ).split()
    ra2, dec2 = coordsDegToSexa(
        float(outline2[i+1][0]),
        float(outline2[i+1][1]),
        fmt=('%02d:%02d:%06.3f ', '%s%02d:%02d:%06.3f')
        ).split()
    
    shape = "line"
    colour = "colour=black"
    width = "width=2"
    line = [shape, ra1, dec1, ra2, dec2, "#",colour,
            width,"\n"]
    line = " ".join(line)
    reg_file.write(line) 
# close file
reg_file.close()
