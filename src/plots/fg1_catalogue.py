#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 23:28:25 2022

@author: Jia Wei Teh

This script plots the FoV of the combined catalogue. See Figure 1.
"""
# libraries
import matplotlib.pyplot as plt
# --- allow running this file directly: put repo root on sys.path ---
import os as _os, sys as _sys
_root = _os.path.dirname(_os.path.abspath(__file__))
while not _os.path.isdir(_os.path.join(_root, "src")) and _root != _os.path.dirname(_root):
    _root = _os.path.dirname(_root)
if _root not in _sys.path:
    _sys.path.insert(0, _root)
# ------------------------------------------------------------------
from src import paths
import numpy as np
from matplotlib import transforms
import pylab
import scipy.ndimage
#--
import src.tools.draw_FOV as draw_FOV
import src.tools.read_catalogue as read_catalogue
import src.tools.plot_tools as plot_tools
from src.tools.create_combined_table import d2r, a2d, pc2arc, arc2pc, ang_dist

# Read in catalogue
path2save = paths.DAT
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)

# Initiate lists
scRA, scDEC = [], []
# orphan or associated?
h2RA_o, h2DEC_o, h2rad_o = [], [], []
h2RA_a, h2DEC_a, h2rad_a = [], [], []
# the top3 removed Hii regions
h2RA_t3, h2DEC_t3, h2rad_t3 = [], [], []

# =============================================================================
# Categorise regions.
# Region without cluster = no colour.
# Region with cluster = red colour.
# Show only class 1, 2 hii region.
# =============================================================================
# Plot H2 regions
for h2 in h2_catalogue:
    # do not include class 3/4 clusters
    h2class = read_catalogue.get_h2_param(h2, 'class', h2_catalogue)
    if h2class == 4 or h2class == 3 :
        continue
    if not read_catalogue.is_BPTStarforming(h2):
        continue
    # h2region data
    # We are only interested in regions associated with more than one 
    # of Class 1,2,3 sources. Regions with purely Class 4 are considered
    # orphan regions here.
    assoc = read_catalogue.get_h2_param(h2, 'assoc', h2_catalogue)
    # returns True if list is empty or if it is only Class 4
    containsOnlyClass4 = all([read_catalogue.get_sc_param(sc, 'class', sc_catalogue) == 4\
                          for sc in assoc])
    # categorise
    if containsOnlyClass4:
        h2RA_o.append(read_catalogue.get_h2_param(h2, 'RA', h2_catalogue))
        h2DEC_o.append(read_catalogue.get_h2_param(h2, 'DEC', h2_catalogue))
        h2rad_o.append(read_catalogue.get_h2_param(h2, 'rad', h2_catalogue))
    else:
        h2RA_a.append(read_catalogue.get_h2_param(h2, 'RA', h2_catalogue))
        h2DEC_a.append(read_catalogue.get_h2_param(h2, 'DEC', h2_catalogue))
        h2rad_a.append(read_catalogue.get_h2_param(h2, 'rad', h2_catalogue))    
# Plot star clusters
for sc in sc_catalogue:
    scclass = read_catalogue.get_sc_param(sc, 'class', sc_catalogue)
    if scclass != 4:
        scRA.append(read_catalogue.get_sc_param(sc, 'RA', sc_catalogue))
        scDEC.append(read_catalogue.get_sc_param(sc, 'DEC', sc_catalogue)) 

# =============================================================================
# Plot
# =============================================================================

fig, ax = plt.subplots(1,1, figsize = (10,8), dpi = 300)
# invert axis
ax.invert_xaxis()
# set limit of plot
xlims = (24.227, 24.135)
ylims = (15.74, 15.807)
# this is the real limit of the FITS image
real_center = (24.1884, 15.7742)
pixelsize = a2d(0.03962000000000001)
npixel2center = 6000
realx = (real_center[0] + npixel2center*pixelsize, real_center[0] - npixel2center*pixelsize)
realy = (real_center[1] - npixel2center*pixelsize, real_center[1] + npixel2center*pixelsize)

# TODO
# First, insert image of HST field.
# path2image = r"/Users/jwt/Documents/Thesis/Honours_Figures/LEGUS_fits.png"
# im = plt.imread(path2image)
# tr = scipy.ndimage.rotate(im, -6)
# implot = plt.imshow(tr, 
#                     extent = [*realx, *realy],
#                     # transform = tr,
#                     # cmap = 'gray',
#                     alpha = 0.2,
#                     )

plt.scatter(*real_center, s = 30)

# Then, plot h2 regions
def my_circle_scatter_radii(axes, x_array, y_array, radii_array, orphan = True,
                            highlightTop3 = False):
    
    """Draws circle of H II regions with physical size """
    if highlightTop3:
         for (x, y, r) in zip(x_array, y_array, radii_array):
                r = a2d(pc2arc(r))
                circle = pylab.Circle((x,y), radius=r, 
                                      fc = 'yellow',
                                      ec = 'yellow',
                                      linewidth = 0.3,
                                      alpha = 0.7
                                      )
                axes.add_patch(circle)
    else:
        # if orphan == True
        if orphan:
            for (x, y, r) in zip(x_array, y_array, radii_array):
                r = a2d(pc2arc(r))
                circle = pylab.Circle((x,y), radius=r, 
                                      fc = 'none',
                                      ec = 'k',
                                      linewidth = 0.3,
                                      alpha = 0.5,
                                      )
                axes.add_patch(circle)
        else:
            for (x, y, r) in zip(x_array, y_array, radii_array):
                r = a2d(pc2arc(r))
                circle = pylab.Circle((x,y), radius=r, 
                                      fc = 'red',
                                      ec = 'k',
                                      linewidth = 0.3,
                                      alpha = 0.5,
                                      )
                axes.add_patch(circle)        
    return True

my_circle_scatter_radii(ax, h2RA_o, h2DEC_o, h2rad_o, True)
my_circle_scatter_radii(ax, h2RA_a, h2DEC_a, h2rad_a, False)
# Plot circle?
# my_circle_scatter_radii(ax, h2RA_t3, h2DEC_t3, h2rad_t3, highlightTop3 = True)

# plot SCs
plot_tools.plot_scatter(scRA, scDEC, setticks = [0.02, 5, 0.02, 5],
                  s = 1,
                  facecolors = 'k',
                  edgecolors = 'k',
                  # saves= "old_fov",
                  label = 'Star clusters (1253)',
                  marker = '.')

# for legend purposes
plot_tools.plot_scatter(0, 0,
                  s = 10,
                   fc = 'red',
                    ec = 'k',
                    linewidth = 0.3,
                    alpha = 0.5,
                  label = 'Populated HII regions (169)',
                  marker = 'o')
plot_tools.plot_scatter(0, 0,
                  s = 10,
                   fc = 'none',
                    ec = 'k',
                    linewidth = 0.3,
                    alpha = 0.5,
                  label = 'Orphan HII regions (165)',
                  marker = 'o')

# miscellaneous
plt.ylabel("Dec (deg)", family = 'Times New Roman', fontsize = '13')
plt.xlabel("RA (deg)", family = 'Times New Roman', fontsize = '13')
# Plot field of view
draw_FOV.FOV(draw = True, newFOV = True)
plt.xlim(*xlims)
plt.ylim(*ylims)

lgnd = plt.legend(loc="upper left", scatterpoints=1, fontsize=15)
lgnd.legendHandles[0]._sizes = [50]
lgnd.legendHandles[1]._sizes = [80]
lgnd.legendHandles[2]._sizes = [80]


# plt.legend(loc = 'upper left', 
#             fontsize = 14,
#             labelspacing = .3,
#             frameon = True)


plt.hlines(15.8, 24.155, 24.13768, 'k')
plt.text(24.155, 15.802, "60\" (2.9 kpc)", fontsize = 15)


# =============================================================================
# Here we include the resolution check box
# =============================================================================

from matplotlib.patches import Rectangle

ax.add_patch(Rectangle((24.2218, 15.7835), .0024, .0027,
             edgecolor = 'k',
             facecolor = 'None',
             )
             )

# save
plot_tools.save("Catalogue")


# =============================================================================
# Some statistics
# =============================================================================

print(f'There are {len(h2RA_o)} orphan HII regions.')
print(f'There are {len(h2RA_a)} associated HII regions.')
print(f'There are {len(h2RA_o) + len(h2RA_o)} total HII regions.')

print(f'{len(h2RA_o)/(len(h2RA_o) + len(h2RA_a))}% are orphan HII regions')






