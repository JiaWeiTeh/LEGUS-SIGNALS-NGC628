#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 22:26:43 2022

@author: Jia Wei Teh

This script demonstrates the effect of oversegmentation via distance
of orphan to the nearest populated hiiregion. See figure 2. Also include
calculations for Section 3.1, 3.2.
"""
# libraries
import numpy as np
import matplotlib.pyplot as plt
#--
import src.tools.read_catalogue as read_catalogue
import src.tools.plot_tools as plot_tools
from src.tools.create_combined_table import d2r, pc2arc, arc2pc, ang_dist

# Read in catalogue
path2save = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/"
sc_catalogue, h2_catalogue  = np.load(path2save+"combined_catalogue.npy", allow_pickle = True)

# Create array to store distances to cluster for every HII region.
distanceDist = []

# =============================================================================
# Check distances from one orphan region to nearest populated region
# =============================================================================
for h2 in h2_catalogue:
    # check if falls into star forming region
    if not read_catalogue.is_BPTStarforming(h2):
        continue
    # For the first loop, we look only class 1, 2 orphan regions
    # region is orphan if it only contains Class 4. For this analysis
    # we are only interested in Class 1, 2, 3 sources.
    scAssoc_list = read_catalogue.get_h2_param(h2, 'assoc', h2_catalogue)
    containsOnlyClass4 = True
    # if empty, then containsOnlyClass4 = True regardless
    # if not, then containsOnlyClass4 = True only when there is class4, 
    # which essentially means orphan for what we care here.
    for sc in scAssoc_list:
        if read_catalogue.get_sc_param(sc, 'class', sc_catalogue) != 4:
            containsOnlyClass4 = False
    if containsOnlyClass4:
        h2Class = read_catalogue.get_h2_param(h2, 'class', h2_catalogue)
        if h2Class == 1 or h2Class == 2:
            p1 = (d2r(h2[1]),d2r(h2[2]))
            radius = pc2arc(h2[12])
            h2Dist = [] # temporary storage
            # For the second loop, we look at all populated regions regardless of 
            # morphology.
            for h2_2 in h2_catalogue:
                scAssoc_list_2 = read_catalogue.get_h2_param(h2_2, 'assoc', h2_catalogue)
                if len(scAssoc_list_2) != 0:
                    # For this analysis we are only interested in Class 1, 2, 3 LEGUS sources.
                    # Hence, do not include pure Class 4 associations in this section.
                    containsOnlyClass4 = True
                    for sc_2 in scAssoc_list_2:
                        if read_catalogue.get_sc_param(sc_2, 'class', sc_catalogue) != 4:
                            containsOnlyClass4 = False
                            break
                    if not containsOnlyClass4:
                        p2 = (d2r(h2_2[1]),d2r(h2_2[2]))
                        dist = ang_dist(p1, p2)
                        h2Dist.append(dist)
            # append
            distanceDist.append(sorted(h2Dist))
# sort by distance
nearest_distance = np.array(distanceDist)[:,0]
# into parsecs
nearest_distance = np.array([arc2pc(x) for x in nearest_distance])
# print
print(f'The number of orphan HII regions having a populated neighbour within 100 pc is {len(nearest_distance[nearest_distance<100])}')

# =============================================================================
# Plot
# =============================================================================
fig, ax1 = plt.subplots(1, 1, figsize = (6,5), dpi = 200)

plot_tools.plot_hist(nearest_distance, 
               ax = ax1,
               bins = 40,
               range = (0, 1000),
               xlim = (0, 1000),
               ylim = (-.5, 27),
               fill = False,
                setticks = [200, 4, 5, 5],
               histtype = 'step',
               xlabel = "Distance (pc)",
               ylabel = "Number",
               )

# location for the zoomed portion 
sub_axes = plt.axes([.5, .5, .30, .3]) 

# plot the zoomed portion
plot_tools.plot_hist(nearest_distance,  bins = np.linspace(0, 100, 26),
               ax = sub_axes,
         histtype='step',
         stacked=True,
          # weights = weights,
         fill=False, color = 'grey',
         title = "zoom-in ($<$ 100 pc)",
         range = [0,100],
         ylabel = "Number",
         xlabel = "Distance (pc)",
         )
# save
plot_tools.save('Oversegmentation')

#%%
# =============================================================================
# Bonus for astrometric uncertainty section (3.1)
# =============================================================================

# distance of orphan to the nearest star cluster
dist2neareststar = []

for h2 in h2_catalogue:
    # do not include class 3/4 clusters
    h2class = read_catalogue.get_h2_param(h2, 'class', h2_catalogue)
    if h2class == 4 or h2class == 3 :
        continue
    if not read_catalogue.is_BPTStarforming(h2):
        continue
    scAssoc_list = read_catalogue.get_h2_param(h2, 'assoc', h2_catalogue)
    containsOnlyClass4 = True
    # if empty, then containsOnlyClass4 = True regardless
    # if not, then containsOnlyClass4 = True only when there is class4, 
    # which essentially means orphan for what we care here.
    for sc in scAssoc_list:
        if read_catalogue.get_sc_param(sc, 'class', sc_catalogue) != 4:
            containsOnlyClass4 = False
    if containsOnlyClass4:
        p1 = (d2r(h2[1]),d2r(h2[2]))
        radius = pc2arc(h2[12])
        # temporary list
        temp_dist = []
        for s in sc_catalogue:
            if read_catalogue.get_sc_param(s, 'class') in [1, 2, 3]:
                p2 = (d2r(s[3]),d2r(s[4]))
                dist = ang_dist(p1, p2) - radius
                temp_dist.append(dist)
        dist2neareststar.append(sorted(temp_dist))

# nearest star?
dist2neareststar = np.array(dist2neareststar)[:,0]
# within 1"?
dist2neareststar_within1pc = dist2neareststar[(dist2neareststar<.1)]
# print
print(f'There are {len(dist2neareststar_within1pc)} HII regions with star clusters within .1" of their outer edge.')


#%%

# =============================================================================
# Bonus for single-star HII region section (3.2)
# =============================================================================

# radius list
radiuslist = []

for h2 in h2_catalogue:
    # do not include class 3/4 clusters
    h2class = read_catalogue.get_h2_param(h2, 'class', h2_catalogue)
    if h2class == 4 or h2class == 3 :
        continue
    if not read_catalogue.is_BPTStarforming(h2):
        continue
    scAssoc_list = read_catalogue.get_h2_param(h2, 'assoc', h2_catalogue)
    containsOnlyClass4 = True
    # if empty, then containsOnlyClass4 = True regardless
    # if not, then containsOnlyClass4 = True only when there is class4, 
    # which essentially means orphan for what we care here.
    for sc in scAssoc_list:
        if read_catalogue.get_sc_param(sc, 'class', sc_catalogue) != 4:
            containsOnlyClass4 = False
    if containsOnlyClass4:
        radiuslist.append(read_catalogue.get_h2_param(h2, 'rad'))

print(f"The median radius of orphan HII regions is {np.median(radiuslist)}")










