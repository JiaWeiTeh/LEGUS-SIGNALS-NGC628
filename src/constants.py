# -*- coding: utf-8 -*-
"""Central physical constants for the LEGUS-SIGNALS-NGC628 analysis.

These are the literals that were previously hard-coded across the figure and tool
scripts, collected here so each is defined once and documented. Importing this module
has no side effects and no `src` dependencies.
"""

# Kennicutt (1998) Hα -> ionising-photon conversion:
#   Q(Hα) = L_Hα * QHA_PER_LHA      [photons s^-1 per erg s^-1]
QHA_PER_LHA = 7.31e11

# Nebular-to-stellar extinction ratio, A_V,neb / A_V,star.
AV_NEB_OVER_STAR = 2.27

# Adopted distance to NGC 628 (9.9 Mpc), expressed in parsecs.
DIST_NGC628_PC = 9.9e6
