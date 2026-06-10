# -*- coding: utf-8 -*-

"""This script returns the outline of the field-of-view"""

import matplotlib.pyplot as plt
import numpy as np

def FOV(newFOV = True, draw = False):
    """
    This function returns array of coordinates defining the FoV of the observation.
    Note: coordinates are determined by inspecting the FoV in Ds9.
    
    Parameters
    ----------
    newFOV : Boolean {True | False} 
        -- True: newly defined FOV. 
        -- False: return the coordinates for V-band (Grasha+15) 
    draw: Boolean {True | False}
        -- True: return a plot of FoV
        -- False: return an array of coordinates
    Returns
    -------
    See above.

    """
    # newly defined FoV (see paper)
    if newFOV:
        FOV1 = np.array([
            (24.134064, 15.78436),
            (24.156826, 15.797972),
            (24.1743368, 15.754506),
            (24.15680393, 15.74375047),
            (24.134064, 15.78436)
            ]
            )
        
        FOV2 = np.array([
            (24.163153, 15.784225),
            (24.174967, 15.75489),
            (24.189738, 15.763954),
            (24.200787, 15.744650),
            (24.198326, 15.761191), 
            (24.201340, 15.743681),
            (24.207205, 15.733435),
            (24.215758,  15.735718),
            (24.229897, 15.744400),
            (24.221962, 15.793204),
            (24.192917, 15.776170),
            (24.174423, 15.808493),
            (24.157491, 15.798369),
            (24.164417, 15.781386),
            (24.163153, 15.784225)
            ]
            )
    
    else:
    # old F555 FoV
        FOV1 = [
            (24.248023  ,15.755537),
            (24.208466  ,15.731233),
            (24.189738   ,15.763954),
            (24.156795   ,15.743745),
            (24.134064   ,15.784362),
            (24.174423   ,15.808504),
            (24.192894   ,15.776181),
            (24.226613   ,15.795955),
            (24.248023   ,15.755537)]
        FOV2 = [
            (24.189738   ,15.763954),
            (24.156795   ,15.743745),
            (24.134064   ,15.784362),
            (24.174423   ,15.808504),
            (24.192894   ,15.776181)]

    # return plot or array?
    if draw:
        plt.plot(*zip(*FOV1), 'k', lw=.5, alpha = 0.5)
        plt.plot(*zip(*FOV2), 'k', lw=.5, alpha = 0.5)
    else:
        return FOV1, FOV2
    
    