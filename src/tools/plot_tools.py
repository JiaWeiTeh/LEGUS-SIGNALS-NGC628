#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 21:27:32 2022

@author: Jia Wei Teh

This script contains some useful plotting functions.
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
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif', size=12)

# Function that does scatter plots
def plot_scatter(x, y, has_cmap = [], ax=None, clabel = [],
                 cticks = [],
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):
    
    if ax is None:
        ax = plt.gca()
    if has_cmap:
        mappable = ax.scatter(x, y, c=has_cmap[0], cmap = has_cmap[1], **plt_kwargs)
        try:
            cbar = clabel[0].colorbar(mappable, pad = 1e-2)
            cbar.set_label(clabel[1], size = 12)
            cbar.ax.tick_params(labelsize=10)
        except IndexError:
            print("Add label for colorbar pls")
    else:
        ax.scatter(x,y, **plt_kwargs)
    if cticks:
        cbar.remove()
        cbar = clabel[0].colorbar(mappable, ticks = cticks, pad = 1e-2)
        cbar.set_label(clabel[1], fontsize = 13)
        cbar.ax.tick_params(labelsize=10)
    if xlabel is not None:
        ax.set_xlabel(xlabel, family = 'Times New Roman', fontsize = '12')
    if ylabel is not None:
        ax.set_ylabel(ylabel, family = 'Times New Roman', fontsize = '12')
    if title is not None:
        ax.set_title(title)
    if setticks:
        x, xinterval, y, yinterval = setticks
        ax.xaxis.set_major_locator(MultipleLocator(x))
        ax.xaxis.set_minor_locator(AutoMinorLocator(xinterval))
        ax.yaxis.set_major_locator(MultipleLocator(y))
        ax.yaxis.set_minor_locator(AutoMinorLocator(yinterval))
        ax.tick_params(axis='both', which = 'major', direction = 'in',length = 6, width = 1)
        ax.tick_params(axis='both', which = 'minor', direction = 'in',length = 4, width =1)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)  
    if saves is not None:
        save(saves)
    return(ax)

# Function that does line plots
def plot_plot(x, y, has_cmap = [], ax=None, clabel = [],
                 cticks = [],
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):

    if ax is None:
        ax = plt.gca()
    if has_cmap:
        mappable = ax.plot(x, y, c=has_cmap[0], cmap = has_cmap[1], **plt_kwargs)
        try:
            cbar = clabel[0].colorbar(mappable)
            cbar.set_label(clabel[1], fontsize = 10)
            cbar.ax.tick_params(labelsize=10)
        except IndexError:
            print("Add label for colorbar pls")
    else:
        ax.plot(x,y, **plt_kwargs)
    if cticks:
        cbar.remove()
        cbar = clabel[0].colorbar(mappable, ticks = cticks)
        cbar.set_label(clabel[1], fontsize = 10)
        cbar.ax.tick_params(labelsize=10)
    if xlabel is not None:
        ax.set_xlabel(xlabel, family = 'Times New Roman', fontsize = '13')
    if ylabel is not None:
        ax.set_ylabel(ylabel, family = 'Times New Roman', fontsize = '13')
    if title is not None:
        ax.set_title(title)
    if setticks:
        x, xinterval, y, yinterval = setticks
        ax.xaxis.set_major_locator(MultipleLocator(x))
        ax.xaxis.set_minor_locator(AutoMinorLocator(xinterval))
        ax.yaxis.set_major_locator(MultipleLocator(y))
        ax.yaxis.set_minor_locator(AutoMinorLocator(yinterval))
        ax.tick_params(axis='both', which = 'major', direction = 'in',length = 6, width = 1)
        ax.tick_params(axis='both', which = 'minor', direction = 'in',length = 4, width =1)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)  
    if saves is not None:
        save(saves)
    return(ax)

# Function that does errorbar plots
def plot_error(x, y, c = [], ax=None, clabel = [],
                 cticks = [],
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):

    if ax is None:
        ax = plt.gca()
    if c:
        mappable = ax.errorbar(x, y, c=c[0], cmap = c[1], **plt_kwargs)
        cbar = clabel[0].colorbar(mappable)
        cbar.set_label(clabel[1], fontsize = 10)
        cbar.ax.tick_params(labelsize=10)
    else:
        ax.errorbar(x,y, **plt_kwargs)
    if cticks:
        cbar.ax.set_yticklabels(cticks)
    if xlabel is not None:
        ax.set_xlabel(xlabel, family = 'Times New Roman', fontsize = '13')
    if ylabel is not None:
        ax.set_ylabel(ylabel, family = 'Times New Roman', fontsize = '13')
    if title is not None:
        ax.set_title(title)
    if setticks:
        x, xinterval, y, yinterval = setticks
        ax.xaxis.set_major_locator(MultipleLocator(x))
        ax.xaxis.set_minor_locator(AutoMinorLocator(xinterval))
        ax.yaxis.set_major_locator(MultipleLocator(y))
        ax.yaxis.set_minor_locator(AutoMinorLocator(yinterval))
        ax.tick_params(axis='both', which = 'major', direction = 'in',length = 6, width = 1)
        ax.tick_params(axis='both', which = 'minor', direction = 'in',length = 4, width =1)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)  
    if saves is not None:
        save(saves)
    return(ax)

# Function that does histogram plots
def plot_hist(x, ax=None, 
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):

    if ax is None:
        ax = plt.gca()
    ax.hist(x, **plt_kwargs)
    if xlabel is not None:
        ax.set_xlabel(xlabel, family = 'Times New Roman', fontsize = '13')
    if ylabel is not None:
        ax.set_ylabel(ylabel, family = 'Times New Roman', fontsize = '13')
    if title is not None:
        ax.set_title(title, family = 'Times New Roman', fontsize = '13')
    if setticks:
        x, xinterval, y, yinterval = setticks
        ax.xaxis.set_major_locator(MultipleLocator(x))
        ax.xaxis.set_minor_locator(AutoMinorLocator(xinterval))
        ax.yaxis.set_major_locator(MultipleLocator(y))
        ax.yaxis.set_minor_locator(AutoMinorLocator(yinterval))
        ax.tick_params(axis='both', which = 'major', direction = 'in',length = 6, width = 1)
        ax.tick_params(axis='both', which = 'minor', direction = 'in',length = 4, width =1)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)  
    if saves is not None:
        save(saves)
    return(ax)

# function that does bar plots
def plot_bar(x, ax=None, 
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):

    if ax is None:
        ax = plt.gca()
    labels, counts = np.unique(x, return_counts=True)
    ax.bar(labels, counts, align='center')
    ax.set_xticks(labels)
    if xlabel is not None:
        ax.set_xlabel(xlabel, family = 'Times New Roman', fontsize = '13')
    if ylabel is not None:
        ax.set_ylabel(ylabel, family = 'Times New Roman', fontsize = '13')
    if title is not None:
        ax.set_title(title)
    if setticks:
        x, xinterval, y, yinterval = setticks
        ax.xaxis.set_major_locator(MultipleLocator(x))
        ax.xaxis.set_minor_locator(AutoMinorLocator(xinterval))
        ax.yaxis.set_major_locator(MultipleLocator(y))
        ax.yaxis.set_minor_locator(AutoMinorLocator(yinterval))
        ax.tick_params(axis='both', which = 'major', direction = 'in',length = 6, width = 1)
        ax.tick_params(axis='both', which = 'minor', direction = 'in',length = 4, width =1)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)  
    if saves is not None:
        save(saves)
    return(ax)

# save figure
def save(filename, ext = '.pdf', transparent = True):
    """Saves figure"""
    path = paths.FIGS
    savefile = path + filename + ext
    plt.savefig(savefile, bbox_inches='tight', transparent = transparent)
    print(f"Figure saved in {savefile}")

# plots results from previous f_esc studies.
def plot_studies(ax = None):
    """Plots results from other studies onto LHa vs QH0"""
    
    if ax is None:
        ax = plt.gca()
    
    # --------------Keniccutt+98----------------
    
    kenniy = np.linspace(20, 80, 200)
    kennix = kenniy + 11.8639
    
    plt.plot(kennix, kenniy, linestyle = '--', linewidth  = 2, alpha = 1, 
             zorder = 90,
              c='k', label = 'Case B, Kenicutt98')   
    
    # --------------LMC, McLeod+19----------------
    
    LMCy = [37.63, 37.27, 36.80,  36.49, 37.18, 37.17,
            37.85, 36.04,36.55,  35.91, 35.80]
    LMCx = [49.78, 49.52, 48.88,48.65, 49.53, 49.22,
            50.08, 48.61, 49.13,48.27, 48.06]
    
    plt.scatter(LMCx, LMCy, s=100,
                alpha = 0.80, 
                c = 'lightgrey',
                edgecolor = 'k',
                zorder = 80,
                marker = 'd',
                label = 'LMC (McLeod+19)')
    
    # --------------NGC 300----------------
    
    NGC300x = [38.04, 37.14, 38.28, 37.62, 37.33]
    
    NGC300y = [49.52, 49.15, 49.66, 49.11, 49.5]
    
    # make the scatter
    plt.scatter(NGC300y, NGC300x, s=80,
                alpha = 0.8,
                edgecolor = 'k',
                c = 'gainsboro',
                zorder = 80,
                label = 'NGC300 (McLeod+20)',
                marker = 's')
    
    # --------------NGC7793----------------
    
    NGC7793x = [51.16,
                50.58,
                50.52,
                50.34,
                49.70,
                49.18,
                49.19,
                51.11]
    
    NGC7793y = np.array([50.74,
                50.31,
                50.35,
                49.98,
                49.10, 
                49.14,
                48.75,
                50.53])
    
    NGC7793y = 10**(NGC7793y)
    NGC7793y = NGC7793y/7.31e11
    NGC7793y = np.log10(NGC7793y)
    
    # make the scatter
    plt.scatter(NGC7793x, NGC7793y, s=150,
                alpha = .8,
                edgecolor = 'k',
                zorder = 80,
                c = 'darkgrey',
                label = 'NGC7793 (Della Bruna+21)',
                marker = '*')
    return

# Function which takes in values and outputs in latex format for table in paper.
def latexReadable(minx, medx, maxx):
    """Given min, median and max value, returns latex readable string"""
    
    plus = np.round(maxx - medx, 2)
    minus = np.round(medx - minx, 2)
    medx = np.round(medx, 2)
    
    return "$"+ str(medx) + "^{+" + str(plus) + "}_{-" + str(minus) + "}$"

