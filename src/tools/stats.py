# -*- coding: utf-8 -*-
"""Statistics / PDF helpers shared across the pipeline.

Consolidated in Patch 4 from duplicated copies in the figure scripts and
obAssociations_budget. `samplePDF`'s `niter` is now a parameter (default 1e5) so the
single definition covers both the plot-script calls (2 args) and obAssociations (3 args).
"""
import numpy as np


def samplePDF(xvalues, pdf, niter=int(1e5)):
    """ This function creates pdf that is not spiky (smoothes out the pdf).
    Solves problems arising from the 'delta-function' monte carlo approach.
    Then, it picks 10**5 samples from the pdf.

    The per-draw root solve is vectorised; the random draws (choice + uniform) are
    unchanged, so with a fixed seed the result is reproducible and matches the former
    per-draw `np.roots` implementation to ~1e-12 (closed-form quadratic instead of
    companion-matrix eigenvalues)."""
    # The probability of landing into one of the (left-edge) bins
    probability = abs(xvalues[1]-xvalues[0]) * (pdf[1:]+pdf[:-1])/2
    # normalize to bypass error tolerence
    probability = probability/np.sum(probability)
    # random draws -- identical RNG consumption to the original per-draw loop
    selectbin = np.random.choice(xvalues[:-1], p = probability, size = niter)
    rand_uni = np.random.uniform(0, 1, niter)
    # endpoints of each drawn bin (selectbin entries are exact xvalues nodes)
    bi = np.searchsorted(xvalues, selectbin)
    a0 = xvalues[bi];      a1 = pdf[bi]
    b0 = xvalues[bi + 1];  b1 = pdf[bi + 1]
    # line through the bin (y = m x + c) and its normalisation factor
    m = (a1 - b1) / (a0 - b0)
    c = a1 - m * a0
    norm = 1.0 / (0.5 * (a1 + b1) * np.abs(a0 - b0))
    # solve (m/2) r^2 + c r - K = 0 for the root inside (a0, b0)
    K = m / 2 * a0**2 + c * a0 + rand_uni / norm
    montecarlos = np.empty(niter)
    flat = (m == 0)                       # degenerate flat bin -> linear solve
    montecarlos[flat] = K[flat] / c[flat]
    q = ~flat
    disc = np.sqrt(c[q]**2 + 2.0 * m[q] * K[q])
    r1 = (-c[q] + disc) / m[q]
    r2 = (-c[q] - disc) / m[q]
    montecarlos[q] = np.where((r1 > a0[q]) & (r1 < b0[q]), r1, r2)
    return montecarlos

def prob2pdf(xvalues, pdf):
    """Converts probability function into probability density function"""
    # sum probability
    sumprob = sum(abs(xvalues[1]-xvalues[0])*(pdf[1:]+pdf[:-1])/2)
    # normalisation constant
    norm = 1/sumprob
    # return
    return pdf * norm

def medianPDF(xvalues, pdf):
    """Given xvalues and a pdf, return median and 1-sigma uncertainty """
    pdfsum = np.cumsum(pdf)*(xvalues[1]-xvalues[0])
    percentiles = np.array([
        xvalues[np.argmax(np.greater(pdfsum, 0.159))],
        xvalues[np.argmax(np.greater(pdfsum, 0.5))],
        xvalues[np.argmax(np.greater(pdfsum, 0.841))]])
    return percentiles

def QH02LHa(QH0, A_V):
    """
    Parameters
    ----------
    QH0 : float
        from SLUG. in photon/s
    A_V : float
        from SLUG. Stellar extinction.

    Returns
    -------
    LHa_cor : Extinction corrected Halpha Luminosity.
            Escape fraction considered.
    """
    # 0.27 escaped
    # convert from QH0 to LHa (Lha = QH0/7.31e11) (Kenicutt 1998)
    LHa = QH0*(1-0.27)/7.31e11
    # Nebula extinction from stellar extinction
    Neb_AV = A_V * 2.27
    # Extinction correction  A_V = -2.5log(L_ex/L_0)
    LHa_cor = 10**(Neb_AV/(-2.5)) * LHa
    
    return LHa_cor
