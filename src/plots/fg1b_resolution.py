# -*- coding: utf-8 -*-
"""This script provides visualisation of the resolution difference."""

import matplotlib.pyplot as plt
import cv2
import src.tools.plot_tools as plot_tools
import numpy as np

# =============================================================================
# Image 1 from HST
# =============================================================================
img2 = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/resolution_check2_2.png"
img2 = cv2.imread(img2)
y, x = 500, 550
size = int(400)
crop_img_2 = img2[x-size:x+size, y-size:y+size]
crop_img_2 = cv2.cvtColor(crop_img_2, cv2.COLOR_BGR2RGB)


# =============================================================================
# Image 2 from SIGNALS
# =============================================================================
img1 = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/resolution_check2_1.png"
img1 = cv2.imread(img1)
y, x = 550, 550
size = int(380)
crop_img_1 = img1[x-size:x+size, y-size:y+size]
crop_img_1 = cv2.cvtColor(crop_img_1, cv2.COLOR_BGR2RGB)


# =============================================================================
# Create subplots
# =============================================================================
f = plt.figure(figsize=(17, 8), dpi = 500)
ax1 = f.add_subplot(221)
ax2 = f.add_subplot(211)

# =============================================================================
# Now, for the HST plot
# =============================================================================

ax1.imshow(crop_img_1)

xgap = 9.999999999976694e-06
xtrue = 24.2232
gap = 200

ygap = 9.666666666661901e-06
ytrue = 15.7832
gap = 300

# just a holder
plot_tools.plot_scatter(100, 100, 
                        ax = ax1,
                        c = 'none',
                        s = 100,
                    # xlim = (0, 750),
                    # ylim = (750, 0),
                    setticks = [200, 3, 300, 3],
                    )

xlabels = np.array([0, 
            xtrue+gap*xgap, 
            xtrue, 
            xtrue-gap*xgap, 
            xtrue-2*gap*xgap, ])
xlabels = np.round(xlabels, 4)
ylabels = np.array([0, 
            ytrue+gap*ygap, 
            ytrue, 
            ytrue-gap*ygap, ])

ylabels = np.round(ylabels, 4)

ylabels = [ "$" + str(y) + '$' for y in ylabels]
xlabels = [ "$" + str(x) + '$' for x in xlabels]
ax1.set_yticklabels(list(ylabels))
ax1.set_xticklabels(xlabels)
ax1.tick_params(axis='both', which = 'major', direction = 'in',length = 6, width = 1,
               labelsize = 13)
ax1.tick_params(axis='both', which = 'minor', direction = 'in',length = 4, width =1,
               labelsize = 13)
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')    
ax1.set_xlabel('RA (deg)', fontsize = '15')

ax1.text(.09, .9, "LEGUS NUV (F275W) image", 
        fontsize = 13,
        bbox=dict(facecolor='w', edgecolor='k', alpha = .8),
        zorder = 1000,
        transform=ax1.transAxes)

# =============================================================================
# Now, for the SIGNALS plot
# =============================================================================

ax2.imshow(crop_img_2)

xgap = 9.199999999999875e-06
xtrue = 24.2232
gap = 200

ygap = 8.250000000002977e-06
ytrue = 15.7837
gap = 300   

# just a holder 
plot_tools.plot_scatter(100, 100, ax = ax2,
                    setticks = [200, 3, 300, 3],
                    c= 'none',
                    )

xlabels = np.array([0, 
            xtrue+gap*xgap, 
            xtrue, 
            xtrue-gap*xgap, 
            xtrue-2*gap*xgap, ])
xlabels = np.round(xlabels, 4)
ylabels = np.array([0, 
            ytrue+gap*ygap, 
            ytrue, 
            ytrue-gap*ygap, ])

ylabels = np.round(ylabels, 4)

ylabels = [ "$" + str(y) + '$' for y in ylabels]
xlabels = [ "$" + str(x) + '$' for x in xlabels]
ax2.set_yticklabels(list(ylabels))
ax2.set_xticklabels(xlabels)
ax2.tick_params(axis='both', which = 'major', direction = 'in',length = 6, width = 1,
               labelsize = 13)
ax2.tick_params(axis='both', which = 'minor', direction = 'in',length = 4, width =1,
               labelsize = 13)
ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')    
ax2.set_xlabel('RA (deg)', fontsize = '15')

ax2.text(.09, .9,  "SIGNALS H$\\alpha$ image", 
        fontsize = 13,
        bbox=dict(facecolor='w', edgecolor='k', alpha = .8),
        zorder = 1000,
        transform=ax2.transAxes)

    
plt.tight_layout()

plot_tools.save("resolution_check_comparison", ext = '.pdf')




#%%


# Version 2


# =============================================================================
# Image 1 from HST
# =============================================================================
img1 = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/resolution_check3_1.png"
img1 = cv2.imread(img1)
x, y = 600, 600
size = int(250)
img1 = img1[y-size:y+size, x-size:x+size]
crop_img_1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)


# Version 2


# =============================================================================
# Image 2 from SIGNALS
# =============================================================================
img2 = r"/Users/jwt/Documents/Code/LEGUS-SIGNALS-NGC628/src/dat/resolution_check3_2.png"
img2 = cv2.imread(img2)
x, y = 610, 610
size = int(230)
img2 = img2[y-size:y+size, x-size:x+size]
crop_img_2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)


# =============================================================================
# Create subplots
# =============================================================================
f = plt.figure(figsize=(17, 8), dpi = 500)
ax1 = f.add_subplot(221)
ax2 = f.add_subplot(211)


# =============================================================================
# First, image from HST
# =============================================================================

ax1.imshow(crop_img_1)

xgap = 0.0017000000000031434/300
xtrue = 24.224112
xcount = 100

ygap = 0.0016510000000007352/300
ytrue = 15.784068
ycount = 200



# # just a holder
# plot_tools.plot_scatter(100, 100, 
#                         ax = ax1,
#                         c = 'none',
#                         s = 100,
#                     setticks = [100, 2, 200, 2],
#                     )

# xlabels = np.array([0, 
#             xtrue+xcount*xgap, 
#             xtrue, 
#             xtrue-xcount*xgap, 
#             xtrue-2*xcount*xgap, 
#             xtrue-3*xcount*xgap, 
#             ])




plot_tools.plot_scatter(100, 100, 
                        ax = ax1,
                        c = 'none',
                        s = 100,
                    setticks = [200, 2, 200, 2],
                    )

xlabels = np.array([0, 
            xtrue+xcount*xgap, 
            xtrue-1*xcount*xgap, 
            xtrue-3*xcount*xgap, 
            ])



ylabels = np.array([0, 
            ytrue+2*ycount*ygap, 
            ytrue+ycount*ygap, 
            ytrue, 
            ])

ax1.set_yticklabels(["$%.4f$" % round(y, 4) for y in ylabels])
ax1.set_xticklabels(["$%.4f$" % round(x, 4) for x in xlabels])


ax1.tick_params(axis='both', which = 'major', direction = 'in',length = 6, width = 1,
                labelsize = 13)
ax1.tick_params(axis='both', which = 'minor', direction = 'in',length = 4, width =1,
                labelsize = 13)
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')    
ax1.set_xlabel('RA (deg)', fontsize = '15')
ax1.set_ylabel('Dec (deg)', fontsize = '15')

ax1.text(.09, .9, "LEGUS NUV (F275W) image", 
        fontsize = 13,
        bbox=dict(facecolor='w', edgecolor='k', alpha = .8),
        zorder = 1000,
        transform=ax1.transAxes)


# =============================================================================
# Then, image from SIGNALS
# =============================================================================

ax2.imshow(crop_img_2)

xgap = 0.0019682000000003086/300
xtrue = 24.2241647
xcount = 100

ygap = 0.001722199999999674/300
ytrue = 15.7839351
ycount = 200


# just a holder
plot_tools.plot_scatter(100, 100, 
                        ax = ax2,
                        c = 'none',
                        s = 100,
                    setticks = [200, 2, 200, 2],
                    )

xlabels = np.array([0, 
            xtrue+xcount*xgap, 
            xtrue-1*xcount*xgap, 
            xtrue-3*xcount*xgap, 
            ])


ylabels = np.array([0, 
            ytrue+2*ycount*ygap, 
            ytrue+ycount*ygap, 
            ytrue, 
            ])


ax2.set_yticklabels(["$%.4f$" % round(y, 4) for y in ylabels])
ax2.set_xticklabels(["$%.4f$" % round(x, 4) for x in xlabels])



ax2.tick_params(axis='both', which = 'major', direction = 'in',length = 6, width = 1,
               labelsize = 13)
ax2.tick_params(axis='both', which = 'minor', direction = 'in',length = 4, width =1,
               labelsize = 13)
ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')    
ax2.set_xlabel('RA (deg)', fontsize = '15')

ax2.text(.09, .9,  "SIGNALS H$\\alpha$ image", 
        fontsize = 13,
        bbox=dict(facecolor='w', edgecolor='k', alpha = .8),
        zorder = 1000,
        transform=ax2.transAxes)

    
plt.tight_layout()

plot_tools.save("resolution_check_comparison2", ext = '.pdf')






