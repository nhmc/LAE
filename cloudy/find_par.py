from __future__ import division

from math import log, sqrt, pi
from barak.utilities import adict
from barak.absorb import split_trans_name
from barak.io import parse_config, loadobj
from barak.interp import AkimaSpline, MapCoord_Interpolator
from cloudy.utils import read_observed
import numpy as np
import os
from glob import glob
from barak.plot import get_nrows_ncols, puttext


from matplotlib.ticker import AutoMinorLocator

import astropy.constants as c
import astropy.units as u
from astropy.table import Table
    
import pylab as plt

import sys

# dex 1 sigma error in UVB (and so nH)
Unorm_sig = 0.3

USE_HEXBIN = True



def make_cmap_red():
    from matplotlib.colors import LinearSegmentedColormap
    x = np.linspace(0,1,9)
    cm = plt.cm.Reds(x)
    r,g,b = cm[:,0], cm[:,1], cm[:,2]

    g[0] = 1
    b[0] = 1

    cdict = dict(red=zip(x, r, r), green=zip(x, g, g), blue=zip(x, b, b))
    
    return LinearSegmentedColormap('red_nhmc', cdict)

def make_cmap_blue():
    from matplotlib.colors import LinearSegmentedColormap
    x = np.linspace(0,1,15)
    cm = plt.cm.Blues(x)
    r,g,b = cm[:,0], cm[:,1], cm[:,2]

    g[0] = 1
    b[0] = 1
    r[1:10] = r[4:13]
    g[1:10] = g[4:13]
    b[1:10] = b[4:13]

    cdict = dict(red=zip(x, r, r), green=zip(x, g, g), blue=zip(x, b, b))
    
    return LinearSegmentedColormap('blue_nhmc', cdict)




def find_min_interval(x, alpha):
    """ Determine the minimum interval containing a given probability.

    x is an array of parameter values (such as from an MCMC trace).

    alpha (0 -> 1) is the desired probability encompassed by the
    interval.

    Inspired by the pymc function of the same name.
    """

    assert len(x) > 1
    x = np.sort(x)

    # Initialize interval
    min_int = None, None

    # Number of elements in trace
    n = len(x)

    # Start at far left
    end0 = int(n*alpha)
    start, end = 0, end0

    # Initialize minimum width to large value
    min_width = np.inf

    for i in xrange(n - end0):
        hi, lo = x[end+i], x[start+i]
        width = hi - lo
        if width < min_width:
            min_width = width
            min_int = lo, hi

    return min_int



def make_interpolators_uvbtilt(trans, simnames):
    """ Make interpolators including different UV slopes, given by the
    simulation names.

    simname naming scheme should be (uvb_k00, uvb_k01, uvb_k02, ...),

    uvb k values must be sorted in ascending order!
    """

    Models = []
    aUV = []
    for simname in simnames:
        # need to define prefix, SIMNAME
        gridname = os.path.join(simname, 'grid.cfg')
    
        #print 'Reading', gridname
        cfg = parse_config(gridname)
        aUV.append(cfg.uvb_tilt)
    
        name = os.path.join(simname, cfg.prefix + '_grid.sav.gz')
        #print 'Reading', name
        M = loadobj(name)
        M = adict(M)

        Uconst = (M.U + M.nH)[0]
        #print 'Uconst', Uconst, cfg.uvb_tilt
        assert np.allclose(Uconst, M.U + M.nH)
        Models.append(M)

    ##########################################################################
    # Interpolate cloudy grids onto a finer scale for plotting and
    # likelihood calculation
    ##########################################################################

    roman_map = {'I':0, 'II':1, 'III':2, 'IV':3, 'V':4, 'VI':5,
                 'VII':6, 'VIII':7, 'IX':8, 'X':9, '2':2}
    Ncloudy = {}
    Ncloudy_raw = {}
    #print 'Interpolating...'
    for tr in trans:
        shape = len(M.NHI), len(M.nH), len(M.Z), len(aUV)
        Nvals = np.zeros(shape)
        if tr == 'Tgas':
            for i,M in enumerate(Models):
                Nvals[:,:,:,i] = M['Tgas'][:,:,:,0]
        elif tr == 'NH':
            for i,M in enumerate(Models):
                logNHI = M.N['H'][:,:,:,0]
                logNHII = M.N['H'][:,:,:,1]
                logNHtot = np.log10(10**logNHI + 10**logNHII)
                Nvals[:,:,:,i] = logNHtot
        elif tr in ['CII*']:
            for i,M in enumerate(Models):
                Nvals[:,:,:,i] = M.Nex[tr][:,:,:]
        else:
            atom, stage = split_trans_name(tr)
            ind = roman_map[stage]
            for i,M in enumerate(Models):
                Nvals[:,:,:,i] = M.N[atom][:,:,:,ind]

        # use ndimage.map_coordinates (which is spline interpolation)
        coord = M.NHI, M.nH, M.Z, aUV
        try:
            Ncloudy[tr] = MapCoord_Interpolator(Nvals, coord)
        except:
            import pdb; pdb.set_trace()

        Ncloudy_raw[tr] = Nvals

    #print 'done'
    return Ncloudy, Ncloudy_raw, Models, np.array(aUV, np.float)


def triplot(names, vals, sigvals, fig, indirect={}, labels=None, fontsize=14):

    from barak.plot import hist_yedge, hist_xedge, puttext
    npar = len(names)
    bins = {}
    for n in names:
        x0, x1 = vals[n].min(), vals[n].max()
        dx = x1 - x0
        lo = x0 - 0.1*dx
        hi = x1 + 0.1*dx
        bins[n] = np.linspace(lo, hi, 20)

    axes = {}
    for i0,n0 in enumerate(names):
        for i1,n1 in enumerate(names):
            if i0 == i1:# or i1 < i0: # uncomment to keep just one triangle.
                continue
            ax = fig.add_subplot(npar,npar, i0 * npar + i1 + 1)
            ax.locator_params(tight=True, nbins=8)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            axes[(n0 + ' ' + n1)] = ax

            y,x = vals[n0], vals[n1]
            if USE_HEXBIN:
                ax.hexbin(x,y,cmap=CM, gridsize=40,linewidths=0.1)
            else:
                ax.plot(x,y,'r.', ms=0.5, mew=0)#, alpha=0.5)

            color = 'k' if n0 not in indirect else 'g'
            text = labels[n0] if labels is not None else n0
            puttext(0.05, 0.95, text, ax, color=color ,fontsize=fontsize, va='top')
            color = 'k' if n1 not in indirect else 'g'
            text = labels[n1] if labels is not None else n1
            puttext(0.95, 0.08, text, ax, color=color ,fontsize=fontsize, ha='right')
            # set limits
            y0, y1 = np.percentile(vals[n0], [5, 95])
            dy = y1 - y0
            ax.set_ylim(y0 - dy, y1 + dy)
            x0, x1 = np.percentile(vals[n1], [5, 95])
            dx = x1 - x0
            ax.set_xlim(x0 - dx, x1 + dx)

            c = 'k'
            if i0 == 0:
                ax.xaxis.set_tick_params(labeltop='on')
                ax.xaxis.set_tick_params(labelbottom='off')
                for t in ax.get_xticklabels():
                    t.set_rotation(60)
            elif i0 == npar-1 or (i0 == npar-2 and i1 == npar-1):
                hist_xedge(vals[n1], ax, color='forestgreen',
                           fill=dict(color='forestgreen',alpha=0.3),
                           bins=bins[n1], loc='bottom')
                ax.axvline(sigvals[n1][0], ymax=0.2, color=c, lw=0.5)
                ax.axvline(sigvals[n1][1], ymax=0.2, color=c, lw=0.5)
                cen = sum(sigvals[n1]) / 2.
                ax.axvline(cen, ymax=0.2, color=c, lw=1.5)
                for t in ax.get_xticklabels():
                    t.set_rotation(60)
            else:
                ax.set_xticklabels('')

            if not (i1 == 0 or (i0 == 0 and i1 == 1) or i1 == npar-1):
                ax.set_yticklabels('')

            if (i0 == 0 and i1 == 1) or i1 == 0:
                hist_yedge(vals[n0], ax, color='forestgreen',
                           fill=dict(color='forestgreen',alpha=0.3),
                           bins=bins[n0], loc='left')
                ax.axhline(sigvals[n0][0], xmax=0.2, color=c, lw=0.5)
                ax.axhline(sigvals[n0][1], xmax=0.2, color=c, lw=0.5)
                cen = sum(sigvals[n0]) / 2.
                ax.axhline(cen, xmax=0.2, color=c, lw=1.5)

            if i1 == npar - 1:
                ax.yaxis.set_tick_params(labelright='on')
                ax.yaxis.set_tick_params(labelleft='off')

            #ax.minorticks_on()

    return axes

if 1:
    print_header = False
    if len(sys.argv[1:]) > 0 and sys.argv[1] ==  '--header':
        print_header = True
if 1:
    ##################################################
    # Read configuration file, set global variables
    ##################################################
    testing = 0
    
    cfgname = 'model.cfg'
    # we only need the cfg file for the prefix of the cloudy runs and
    # the name of the file with the observed column densities.
    opt = parse_config(cfgname)

    simnames = sorted(glob(opt['simname']))
    #print opt['simname']
    #print simnames

    #CM = make_cmap_blue()  # plt.cm.binary
    #CM = make_cmap_red()  # plt.cm.binary
    CM = plt.cm.gist_heat_r # plt.cm.binary
    #CM = plt.cm.afmhot_r # plt.cm.binary

    #CM = plt.cm.bone_r # plt.cm.binary
    #CM = plt.cm.terrain_r # plt.cm.binary
    #CM = plt.cm.ocean_r # plt.cm.binary



    trans = 'Tgas', 'NH'


if 1:
    ################################################################
    # Read the cloudy grids and make the interpolators
    ################################################################ 
    Ncloudy, Ncloudy_raw, Models, aUV = make_interpolators_uvbtilt(
        trans, simnames)
    M = Models[0]
    #import pdb; pdb.set_trace()
    Uconst_vals = []
    for model in Models:
        Uconst_vals.append((model['U'] + model['nH'])[0])

    # note it's a function of aUV!
    Uconst = AkimaSpline(aUV, Uconst_vals)

    # Now find the parameter chains
    samples = loadobj('samples_mcmc.sav.gz')
    nwalkers, nsamples, npar = samples['chain'].shape
    parvals = samples['chain'].reshape(-1, npar)

    PAR = samples['par']
    assert PAR['names'][-1] == 'aUV'
    assert PAR['names'][-2] == 'Z'
    assert PAR['names'][-3] == 'nH'
    assert PAR['names'][-4] == 'NHI'

    aUV = parvals[:,-1]
    logZ = parvals[:,-2]
    lognH = parvals[:,-3]
    logNHI = parvals[:,-4]

    logU = Uconst(aUV) - lognH

    #import pdb; pdb.set_trace()

    # call the interpolators with these parameter values.
    logT = Ncloudy['Tgas'](parvals[:,-4:].T)
    logNtot = Ncloudy['NH'](parvals[:,-4:].T)
    # note this is log of D in kpc
    logD = logNtot - lognH - np.log10(c.kpc.to(u.cm).value)

    logP = logT + lognH

    #import pdb; pdb.set_trace()

    H_massfrac = 0.76  # (1 / mu)

    # Joe's mass calculation
    mass = 4./3. * pi * (3./4. * 10**logD * u.kpc)**3 * 10**lognH * \
           u.cm**-3 * u.M_p / H_massfrac

    # D = NH / nH

    logM = np.log10(mass.to(u.M_sun).value)


if 1:
    # print out the results and uncertainties
    vals = dict(U=logU, T=logT, N=logNtot, D=logD, P=logP, M=logM,
                nH=lognH, aUV=aUV, NHI=logNHI, Z=logZ)

    levels = 0.6827, 0.9545

    sigvals = {}
    for key in vals:
        sigvals[key] = find_min_interval(vals[key], levels[0])

    if print_header:
        print r'$\log(Z/Z_\odot)$&$\alpha_{UV}$ & $\log \nH$ & $\log U$& $\log \NHI$ & $\log \NH$& $\log T$ & $\log (P/k)$& $\log D$ & $\log M$ \\'
        print r'                 &              & (\cmmm)    &         & (\cmm)      &  (\cmm)   & (K)      &  (\cmmm K)  &  (kpc)   & (\msun)  \\'

    s = ''
    ALLPAR = 'Z aUV nH U NHI N T P D M'.split()
    for key in ALLPAR:
        sig = 0.5 * (sigvals[key][1] - sigvals[key][0])
        val = 0.5 * (sigvals[key][1] + sigvals[key][0])
        if key in {'nH', 'D', 'P'}:
            sig1 = np.hypot(sig, Unorm_sig) 
            s += '$%.2f\\pm%.2f(%.2f)$ &' % (val, sig1, sig)
        elif key == 'M':
            sig1 = np.hypot(sig, 2*Unorm_sig)
            s += '$%.2f\\pm%.2f(%.2f)$ &' % (val, sig1, sig)
        else:
            s += '$%.2f\\pm%.2f$ &' % (val, sig)

    print s[:-1] +  r'\\'

if 1:
    labels = dict(U='$U$', Z='$Z$', NHI='$N_\mathrm{HI}$', aUV=r'$\alpha_\mathrm{UV}$',
                  T='$T$', P='$P$', N='$N_\mathrm{H}$', D='$D$', M='$Mass$')

if 0:
    fig = plt.figure(figsize=(12,12))
    fig.subplots_adjust(left=0.05, bottom=0.05, top=0.94,right=0.94, wspace=1e-4,hspace=1e-4)
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)

    names = 'U Z NHI aUV T P N D M'.split() 

    #direct = 'U Z NHI aUV'.split()
    axes = triplot(names, vals, sigvals, fig, labels=labels)
    plt.savefig('par.png', dpi=200)

if 1:
    fig = plt.figure(figsize=(8,8))
    fig.subplots_adjust(left=0.095, bottom=0.105, top=0.94,right=0.94, wspace=1e-4,hspace=1e-4)
    plt.rc('xtick', labelsize=9.5)
    plt.rc('ytick', labelsize=9.5)

    names = 'U Z N aUV'.split()
    axes = triplot(names, vals, sigvals, fig, labels=labels, fontsize=16)

    axes['U Z'].set_ylabel('$\log_{10}U$') 
    axes['Z U'].set_ylabel('$\log_{10}[Z/Z_\odot]$') 
    axes['N U'].set_ylabel('$\log_{10}N_\mathrm{H}$') 
    axes['aUV U'].set_ylabel(r'$\log_{10}\alpha_\mathrm{UV}$') 
    axes['aUV U'].set_xlabel('$\log_{10}U$') 
    axes['aUV Z'].set_xlabel('$\log_{10}[Z/Z_\odot]$') 
    axes['aUV N'].set_xlabel('$\log_{10}N_\mathrm{H}$') 
    axes['N aUV'].set_xlabel(r'$\log_{10}\alpha_\mathrm{UV}$') 
    # special case:
    if os.path.abspath('.') == '/Users/ncrighton/Projects/MPIA_QSO_LBG/Cloudy/J0004_NHI_2/comp1/final':
        for k in ('N U', 'N Z', 'N aUV'):
            axes[k].set_ylim(17.3, 19.2)
        for k in ('U N', 'Z N', 'aUV N'):
            axes[k].set_xlim(17.3, 19.2)

            

    #plt.savefig('par2.pdf')
    plt.savefig('par2.png',dpi=250)
