# -*- coding: utf-8 -*-

""" This script stores useful functions"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
from shapely.geometry import Point, Polygon
from sys import maxsize
# Normalize
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif', size=12)

# def read(filename):
#     """ Reads file """
#     path = r'/Users/jwt/Documents/Paper/LyC/data/'
#     return genfromtxt(path+filename, delimiter=',')



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

def get_sc_param(ID, colname, sctable):
    """Given ID/properties and header, return value"""
    
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
        

def get_h2_param(ID, colname, h2table):
    """Given ID/properties and header, return value"""
    
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

""" Radian - Degree - Arcseconds - pc conversion """

# Convert rad to deg

def r2d(rad):
    return rad*180/np.pi

# Convert deg to rad
def d2r(deg):
    return deg*np.pi/180

# Convert degree to arcsecond

def d2a(deg):
    return deg*3600

# Convert arcsec to deg

def a2d(arc):
    return arc/3600

# Convert arc to rad

def a2r(arc):
    return d2r(a2d(arc))

# Convert rad to arc

def r2a(rad):
    return d2a(r2d(rad))

# Convert pc to arcseconds

def pc2arc(radius, D = 9.9e6):
    # distance to NGC 628
    return r2a(radius/D)

# Convert arcseconds to pc

def arc2pc(ang, D = 9.9e6):
    # distance to NGC 628
    return a2r(ang)*D

""" End of Radian - Degree - Arcseconds - pc conversion"""

def dist(p1, p2):
    """Return distance between two points"""
    
    return np.linalg.norm(np.array(p1)-np.array(p2)) 


def ang_dist(p1,p2):
    """Return angular distance (arcsecond) between two points of RA and DEC"""
    """Equation from https://en.wikipedia.org/wiki/Angular_distance"""
    
    r1, d1 = p1
    r2, d2 = p2
    
    # in rad
    A = np.sin(d1)*np.sin(d2)+np.cos(d1)*np.cos(d2)*np.cos(r1-r2)
    rad = np.arccos(A)
    
    return r2a(rad)

    
def plot_scatter(x, y, has_cmap = [], ax=None, clabel = [],
                 cticks = [],
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):
    """ Plots function. Accepts different axes """
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

    
def plot_plot(x, y, has_cmap = [], ax=None, clabel = [],
                 cticks = [],
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):
    """ Plots function. Accepts different axes """
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


    # TODO
def plot_error(x, y, c = [], ax=None, clabel = [],
                 cticks = [],
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):
    """ Plots function. Accepts different axes """
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



def plot_hist(x, ax=None, 
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):
    """ Plots histogram. Accepts different axes """
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

def plot_bar(x, ax=None, 
                 title = None, xlabel = None, ylabel = None,
                 setticks = [],
                 xlim = None, ylim = None,
                 saves = None,
                 **plt_kwargs):
    """ Plots bar. Accepts different axes """
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



def save(filename, ext = '.png', transparent = True):
    
    """Saves figure"""
    path = "/Users/jwt/Documents/Paper/LyC/figures/"
    savefile = path + filename + ext
    plt.savefig(savefile, bbox_inches='tight', transparent = transparent)
    print(f"Figure saved in {savefile}")


def Ha2Photons(LHa):
    
    """Given Ha luminosity return the number of photon required. According to 
    the conversion in Kennicutt+98."""
    
    return 7.31e11*LHa
    
def Photons2Ha(Photons):
    
    """ The other way"""
    
    return Photons/7.31e11

def BPT_N2(ax = None):
    
    """ Plots BPT diagram """
    
    if ax is None:
        ax = plt.gca()
    
    # Kauffmann+03
    def BPT_H2(a):
        return 0.61 / (a-0.05)+ 1.3
    
    xH2 = np.arange(-1.281,0,0.01)
    yH2 = [BPT_H2(i) for i in xH2]
        
    # Kewley+01
    def BPT_COMP(b):
        return 0.61  / (b-0.47)+ 1.19
    
    xComp = np.arange(-4, 0.4, 0.001)
    yComp = [BPT_COMP(i) for i in xComp]
    
    ax.plot(xH2,yH2,'b--', linewidth = 1)
    ax.plot(xComp,yComp,'k-', linewidth = 2)
    
    ax.text(-1.4,-2.5,'star-forming',family = 'Times New Roman', fontsize = '14')
    ax.text(-.01,-2.5,'comp',family = 'Times New Roman', fontsize = '14')
    ax.text(0.1,1.2,'AGN',family = 'Times New Roman', fontsize = '14')



def BPT_S2(ax = None):
    
    """ Plots BPT diagram """
    
    if ax is None:
        ax = plt.gca()
    
    # Main AGN line
    def BPT_Sii_AGN(AGN):
        return (0.72 / (AGN - 0.32) + 1.30)
    
    xH2 = np.arange(-10,0.3,0.01)
    yH2 = [BPT_Sii_AGN(i) for i in xH2]
        
    # Sy2/LINER BPT S2
    def BPT_Sii_SYLIN(SYLINER):
        return (1.89*SYLINER + 0.76)
    
    xSYLIN = np.arange(-0.315, 0.1654, 0.001)
    ySYLIN = [BPT_Sii_SYLIN(i) for i in xSYLIN]
    
    plt.plot(xH2, yH2, 'k-', linewidth = 2)
    plt.plot(xSYLIN, ySYLIN, 'b--')


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
    
    # --------------NGC 300, McLeod+20----------------
    
    NGC300x = [38.04, 37.14, 38.28, 37.62, 37.33]
    
    NGC300y = [49.52, 49.15, 49.66, 49.11, 49.5]
    
    # make the scatter
    plt.scatter(NGC300y, NGC300x, s=80,
                alpha = 0.8,
                edgecolor = 'k',
                c = 'gainsboro',
                zorder = 80,
                # facecolor = 'none',
                label = 'NGC300 (McLeod+20)',
                marker = 's')
    
    
    # --------------NGC7793, Della Bruna+21----------------
    
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
    NGC7793y = Photons2Ha(NGC7793y)
    NGC7793y = np.log10(NGC7793y)
    
    # make the scatter
    plt.scatter(NGC7793x, NGC7793y, s=150,
                alpha = .8,
                edgecolor = 'k',
                zorder = 80,
                c = 'darkgrey',
                            # facecolor = 'none',
                label = 'NGC7793 (Della Bruna+21)',
                marker = '*')

def line_func(p1, p2):
    """ Given two points, return the line coefficients m, c in y = mx + c"""
    
    m = (p1[1]-p2[1])/(p1[0]-p2[0])
    c = p1[1] - m * p1[0]
    
    return m, c

def samplePDF(xvalues, pdf):
    """ This function creates pdf that is not spiky. Solves problems arising
    from the 'delta-function' monte carlo approach."""
    
    def normalisePDF(a, b, m, c):
        """Given xs and ys value of a bin, find the normalisation factor"""
        area = 1/2 * (a[1]+b[1])*abs(a[0]-b[0])
        return 1/area
        
    # The probability of landing into one of the bins
    probability = abs(xvalues[1]-xvalues[0]) * (pdf[1:]+pdf[:-1])/2
    # normalize to bypass error tolerence 
    probability = probability/np.sum(probability)
    # pick a choice. which bin (left)? Number of iterations?
    # left means the selected bin is identified by the left side of the bin.
    niter = int(1e5)
    selectbin = np.random.choice(xvalues[:-1], p = probability, size = niter)
    # pick random number between [0, 1)
    rand_uni = np.random.uniform(0,1,niter)
    # array for montecarlo pdf values
    montecarlos = np.zeros(shape = (len(selectbin)))
    # For each selected bin,    
    for ii, bins in enumerate(selectbin):
        binindex = np.where(xvalues == bins)[0][0]
        # get the two points representing this bin
        a = (xvalues[binindex], pdf[binindex])
        b = (xvalues[binindex+1],pdf[binindex+1])
        # get coeffs of the line equation, m and c in y = mx + c
        m, c = line_func(a, b)
        # normalisation factor
        norm = normalisePDF(a, b, m, c)
        # solve for the y-value in the bin that corresponds to the uniform distribution
        roots = np.roots([m/2,\
                  c,\
                  -( m / 2 * a[0]**2  + c * a[0]  + rand_uni[ii]/norm )])
        # the correct root is within the bin
        hasroot = False
        # search for the correct solution in given roots
        for r in roots:
            if r < b[0] and r > a[0]:
                hasroot = True
                montecarlos[ii] = r
        if not hasroot:  
            raise Exception("Root not found")
                
    return montecarlos

def prob2pdf(xvalues, pdf):
    """Converts probability function into probability density function"""
    
    sumprob = sum(abs(xvalues[1]-xvalues[0])*(pdf[1:]+pdf[:-1])/2)
    
    norm = 1/sumprob
    
    return pdf * norm


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
    QH0_f = QH0*(1-0.27)
    # convert from QH0 to LHa (Lha = QH0/7.31e11) (Kenicutt 1998)
    LHa = Photons2Ha(QH0_f)
    # Nebula extinction from stellar extinction
    Neb_AV = A_V * 2.27
    # Extinction correction  A_V = -2.5log(L_ex/L_0)
    LHa_cor = 10**(Neb_AV/(-2.5)) * LHa
    
    return LHa_cor
    
def medianPDF(xvalues, pdf):
    
    """Given xvalues and a pdf, return median and 1-sigma uncertainty """
    
    pdfsum = np.cumsum(pdf)*(xvalues[1]-xvalues[0])
    percentiles = np.array([
        xvalues[np.argmax(np.greater(pdfsum, 0.159))],
        xvalues[np.argmax(np.greater(pdfsum, 0.5))],
        xvalues[np.argmax(np.greater(pdfsum, 0.841))]])
    
    return percentiles



def min_dist(polygon, point):
    """
    Given polygon and point, return the nearest distance to any of the lines.

    Parameters
    ----------
    polygon : list of datapoints representing the polygon.
        Polygon of interest.
        E.g., [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)].
    point : shapely.geometry.point.Point
        Point to test.
        E.g., [(0.5, 1.5)].

    """


    # a very simple check
    assert sum(polygon[0] == polygon[-1]) == 2, "Please input valid polygon"
    assert len(point) == 2, "Pleaase input valid point"
    #Modification based on: https://gis.stackexchange.com/questions/396/nearest-neighbor-between-point-layer-and-line-layer
    
    polygon = Polygon(polygon)
    point = Point(point)
    
    # pairs iterator:
    # http://stackoverflow.com/questions/1257413/1257446#1257446
    def pairs(lst):
        i = iter(lst)
        first = prev = i.__next__()
        for item in i:
            yield prev, item
            prev = item
        yield item, first
    
    # these methods rewritten from the C version of Paul Bourke's
    # geometry computations:
    # http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
    def magnitude(p1, p2):
        vect_x = p2.x - p1.x
        vect_y = p2.y - p1.y
        return np.sqrt(vect_x**2 + vect_y**2)
    
    def intersect_point_to_line(point, line_start, line_end):
        line_magnitude =  magnitude(line_end, line_start)
        u = ((point.x - line_start.x) * (line_end.x - line_start.x) +
             (point.y - line_start.y) * (line_end.y - line_start.y)) \
             / (line_magnitude ** 2)
    
        # closest point does not fall within the line segment, 
        # take the shorter distance to an endpoint
        if u < 0.00001 or u > 1:
            ix = magnitude(point, line_start)
            iy = magnitude(point, line_end)
            if ix > iy:
                return line_end
            else:
                return line_start
        else:
            ix = line_start.x + u * (line_end.x - line_start.x)
            iy = line_start.y + u * (line_end.y - line_start.y)
            return Point([ix, iy])
    
    nearest_point = None
    min_dist = maxsize
    
    for seg_start, seg_end in pairs(list(polygon.exterior.coords)[:-1]):
        line_start = Point(seg_start)
        line_end = Point(seg_end)
    
        intersection_point = intersect_point_to_line(point, line_start, line_end)
        cur_dist =  magnitude(point, intersection_point)
    
        if cur_dist < min_dist:
            min_dist = cur_dist
            nearest_point = intersection_point
    
    return  min_dist


def point_inside_polygon(x,y,poly):
    n = len(poly)
    inside =False
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y
    return inside
