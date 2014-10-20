"""
Make a likelihood function for use with emcee

Given input Z, nH, k_C, k_N, k_Al, aUV, NHI return ln
of the likelihood.

This module must define the following objects:

- a dictionary P with keys. The value of every key is a tuple with the
same length (the number of model parameters)

    name  : parameter names
    min   : minimum allowed parameter values
    max   : maximum allowed parameter values

- a model(par) function that generates the model of the data given
  an array of parameter values

- a ln_likelihood(par) function

- a get_initial_positions(nwalkers) function that generates an array
  of shape (nwalkers, npar) with parameter values for the initial
  walker positions.

- a plot_model(par) function that plots a model fit to the data given a
  set of parameter values, in the same order as they are listed in P.

- optionall a print_par(par) function.

AXIS ORDER for column densities:

NHI, nH, Z

Reverse this when using XY indexing (e.g. CloughTocher_Interpolator)
"""
from __future__ import division

from math import log, sqrt, pi
from barak.interp import CloughTocher2d_interpolator
from barak.utilities import adict
from barak.absorb import split_trans_name, get_ionization_energy
from barak.io import parse_config, loadobj
from barak.interp import MapCoord_Interpolator
from cloudy.utils import read_observed as read_Nvals
import numpy as np
import os
from glob import glob
from barak.plot import arrplot, get_nrows_ncols, puttext

import astropy.units as u

from scipy.ndimage import map_coordinates

use_ipot = False

log10_cm_per_Mpc = np.log10((1*u.Mpc).to(u.cm).value)

def make_interpolators_uvbtilt(trans, simnames):
    """ Make interpolators including different UV slopes, given by the
    simulation names.

    simname naming scheme should be (uvb_k00, uvb_k01, uvb_k02, ...),

    uvb k values must be sorted in ascending order!
    """

    Models = []
    aUV = []
    for simname in simnames:
        gridname = os.path.join(simname, 'grid.cfg')
    
        print 'Reading', gridname
        cfg = parse_config(gridname)
        aUV.append(cfg.uvb_tilt)
    
        name = os.path.join(simname, cfg.prefix + '_grid.sav.gz')
        print 'Reading', name
        M = loadobj(name)
        M = adict(M)

        Uconst = (M.U + M.nH)[0]
        print 'Uconst', Uconst, cfg.uvb_tilt
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
    print 'Interpolating...'
    for tr in trans + ['NH']:
        shape = len(M.NHI), len(M.nH), len(M.Z), len(aUV)
        Nvals = np.zeros(shape)
        if tr in ['CII*']:
            for i,M in enumerate(Models):
                Nvals[:,:,:,i] = M.Nex[tr][:,:,:]
        elif tr == 'NH':
            for i,M in enumerate(Models):
                logNHI = M.N['H'][:,:,:,0]
                logNHII = M.N['H'][:,:,:,1]
                logNHtot = np.log10(10**logNHI + 10**logNHII)
                Nvals[:,:,:,i] = logNHtot            
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

    print 'done'
    return Ncloudy, Ncloudy_raw, Models, np.array(aUV, np.float)

def model(par, for_plot=False):
    """
    This will interpolate between cloudy models and apply the relative
    abundance variation multipliers and return a bunch of column
    densities.

    par =  all the ks then NHI, nH, Z

    or

    par = all the ks then NHI, nH, Z, aUV

    Returns
    -------
    Nmodel : list
      The model log10 column density for each transition.
    """
    if for_plot:
        trvals = tr_plot
    else:
        trvals = trans
    try:
        coord = par[IND_PAR['NHI']], par[IND_PAR['nH']], par[IND_PAR['Z']], \
                par[IND_PAR['aUV']]
    except:
        import pdb; pdb.set_trace()

    Nmodel = []
    for tr in trvals:
        atom, stage = split_trans_name(tr) 
        N = Ncloudy[tr](coord)
        if atom in IND_KPAR:
            N += par[IND_KPAR[atom]]
        Nmodel.append(N)

    return np.array(Nmodel)


def ln_pdf_siglohi(x, x0, siglo, sighi):
    """ ln of the pdf of an observation centred on an observed value with
    different high and low sigmas.

    Assumes gaussian distribution, but with different sigmas on either
    side of the distribution.

    pdf =  k * exp(-1/2 * ((x - x0)/sighi)**2)  for  x > 0
           k * exp(-1/2 * ((x - x0)/siglo)**2)  for  x <= 0

    where k = 1 / (sqrt(pi/2) * (sighi + siglo))    

    x must be scalar
    """
    #ln_k = -np.log(np.sqrt(0.5 * pi) * (sighi + siglo))
    x = np.atleast_1d(x)
    out = np.empty(x.shape)
    c0 =  x > x0
    out[c0] = -0.5 * ((x[c0] - x0)/sighi)**2 #+ ln_k
    out[~c0] =  -0.5 * ((x[~c0] - x0)/siglo)**2 #+ ln_k
    if out.shape == (1,):
        return out[0]
    else:
        return out

def ln_pdf_uplim(x, x0, sighi):
    """ Find ln(probability) for a series of model values given
    an upper limit.

    Assumes probability drops off as a gaussian of sigma sighi
    """
    # if there is a lower limit too, this is the normalisation:
    
    #ln_k = -log(sqrt(0.5*pi) * sighi + (x0 - lolim))

    # but we'll say there is no lower limits, which means we can't
    # normalise.

    x = np.atleast_1d(x)
    out = np.empty(x.shape)
    c0 =  x > x0
    out[c0] = -0.5 * ((x[c0] - x0)/sighi)**2 #+ ln_k
    out[~c0] =  0 # + ln_k
    if out.shape == (1,):
        return out[0]
    else:
        return out

def ln_pdf_lolim(x, x0, siglo):
    """ Find ln of the pdf for an x value given an lower limit.

    Assumes probability drops off as a gaussian of sigma siglo
    """
    # if there is an upper limit too, this is the normalisation:
    
    #ln_k = -log(sqrt(0.5*pi) * siglo + (uplim - x0))

    # but we'll say there is no upper limit, which means we can't
    # normalise.

    x = np.atleast_1d(x)
    out = np.empty(x.shape)
    c0 =  x > x0
    out[c0] = 0 #+ ln_k
    out[~c0] =  -0.5 * ((x[~c0] - x0)/siglo)**2 #+ ln_k
    if out.shape == (1,):
        return out[0]
    else:
        return out

def ln_likelihood(par, per_obs=False):
    """ Uses obs, trans.

    if per_obs (default False), also return a separate likelihood for
    each observation.
    """
    # if we are outside the allowable parameter ranges, return 0
    # likelihood.
    #import pdb; pdb.set_trace()
    for i,p in enumerate(par):
        if not (P['min'][i] < p < P['max'][i]):
            return -np.inf
    # force the cloud thickness to be < 1 Mpc
    coord = par[IND_PAR['NHI']], par[IND_PAR['nH']], par[IND_PAR['Z']], \
            par[IND_PAR['aUV']]
    logNH = Ncloudy['NH'](coord) 
    if (logNH - par[IND_PAR['nH']]) > log10_cm_per_Mpc:
        return -np.inf
                                
        
    Nmodel = model(par)

    lnprobtot = np.zeros(np.asarray(par[0]).shape)

    for pname in priors:
        if pname.startswith('min ') or pname.startswith('max '):
            continue
        # only deals with two-sided gaussian priors at the moment
        pval, siglo, sighi = priors[pname]
        p = ln_pdf_siglohi(par[IND_PAR[pname]], pval, siglo, sighi)
        lnprobtot += p

    allprob = []
    for i,tr in enumerate(trans):
        Nobs, siglo, sighi = obs[tr]
        if siglo == 0:
            #print(tr, 'lower limit')
            p = ln_pdf_lolim(Nmodel[i], Nobs, SIG_LIMIT)
            lnprobtot += p
            if per_obs:
                allprob.append(p)
        elif sighi == 0:
            #print(tr, 'upper limit')
            p = ln_pdf_uplim(Nmodel[i], Nobs, SIG_LIMIT)
            lnprobtot += p
            if per_obs:
                allprob.append(p)
        else:
            #print(tr)
            siglo = max(siglo, MIN_SIG)
            sighi = max(sighi, MIN_SIG)
            p = ln_pdf_siglohi(Nmodel[i], Nobs, siglo, sighi)
            lnprobtot += p
            if per_obs:
                allprob.append(p)

    if per_obs:
        return lnprobtot, allprob
    else:
        return lnprobtot

def get_initial_positions(nwalkers):
    # Get initial parameter positions (guesses!) for each walker

    Npar = len(P['names'])

    # one possibility:
    
    # generate random values from a normal distribution with a 1
    # sigma width 5 times smaller than the limits for each parameter.
    #p0 = np.random.randn(nwalkers, Npar)
    #for i in range(Npar):
    #    p0[:, i] = P.true[i] + p0[:, i] * (P.max[i] - P.min[i]) / nsigma
    #    # clip so we are inside the parameter limits
    #    p0[:, i] = p0[:, i].clip(P.min[i], P.max[i])

    # another possibility:
    #
    # uniform distribution between parameter limits
    p0 = np.random.uniform(size=(nwalkers, Npar))
    p1 = np.random.normal(size=(nwalkers, Npar))
    for i in range(Npar):
        if P['names'][i] in priors:
            pval, siglo, sighi = priors[P['names'][i]]
            # gaussian
            p0[:, i] = pval + p1[:, i] * 0.5 * (siglo + sighi) 
        else:
            p0[:, i] = P['min'][i] + p0[:, i] * (P['max'][i] - P['min'][i])

    return p0

def plot_model(pars):
    """ Plot the observed values and errors, along with the predicted
    model values.
    """
    import matplotlib.pyplot as pl
    from barak.plot import draw_arrows, puttext

    pl.rc('xtick', labelsize=11)
    pl.rc('ytick', labelsize=11)

    fig = pl.figure(figsize=(6.0, 3.5))
    ax = fig.add_subplot(111)
    ms = 6
    ipot = [get_ionization_energy(t) for t in tr_plot]

    xvals = list(range(len(tr_plot)))
    for par in pars:
        Nmodel = model(par, for_plot=True)
        if use_ipot:
            if len(pars) == 1:
                ax.plot(ipot, Nmodel, 'r.-', lw=1, zorder=0)
            else:
                ax.plot(ipot, Nmodel, 'r-', lw=0.2, zorder=0)
        else:
            if len(pars) == 1:
                ax.plot(xvals, Nmodel, 'r.-', lw=1, zorder=0)
            else:
                ax.plot(xvals, Nmodel, 'r-', lw=0.2, zorder=0)

    y0,y1 = ax.get_ylim()
    dy = y1 - y0

    for i,tr in enumerate(tr_plot):
        if use_ipot:
            ind = ipot[i]
        else:
            ind = i
        if tr in obs:
            colour = 'k' if tr in trans else 'w'
            fmt = 'o' + colour
            val, siglo, sighi = obs[tr]
            if siglo == 0:
                draw_arrows(ind, val, direction='up', ax=ax,lw=1)
                ax.plot(ind, val, fmt,ms=ms)
            elif sighi == 0:
                draw_arrows(ind, val, direction='down', ax=ax,lw=1)
                ax.plot(ind, val, fmt,ms=ms)
            else:
                ax.plot([ind, ind], [val-siglo, val+sighi], 'k',lw=1)
                ax.plot(ind, val, fmt,ms=ms)
            ax.text(ind, val + 0.08*dy, tr,
                    fontsize=11, ha='center')
        else:
            puttext(ind, 0.05, tr, ax=ax, xcoord='data',
                    fontsize=11, ha='center')

    puttext(0.9,0.1, 'Model', ax, color='r', ha='right')

    #print np.asarray(pars).shape

    if use_ipot:
        ax.set_xlabel('Ionization potential (eV)')
        ax.set_xlim(ipot[0]-1, ipot[-1] + 1)
    else:
        ax.set_xlim(-0.5, xvals[-1] + 0.5)
        ax.set_xticks([])
        ax.set_xlabel(r'$\mathrm{Ionization\ potential} \longrightarrow$')

    ax.set_ylabel(r'$\log_{10}\ (N/\mathrm{cm}^{-2})$')
    fig.tight_layout()
    fig.savefig('fig/model.pdf')

    #return fig, ax
    return fig

def print_par(par):
    """ Print the maximum likelihood parameters and their
    uncertainties.
    """
    rec = []
    for i in range(len(P['names'])):
        p = P['ml'][i]
        pmed = P['median'][i]
        m1 = P['p1sig'][i]
        p0 = 0.5 * (m1[0] + m1[1])
        sig1 = 0.5 * (m1[1] - m1[0]) 
        m2 = P['p2sig'][i]
        j1 = P['p1sig_joint'][i]
        j2 = P['p2sig_joint'][i]
        rec.append( (P['names'][i], p0, sig1, m1[0], m1[1],
                     m2[0], m2[1], j1[0], j1[1],
                     j2[0], j2[1], pmed, p) )

    names = 'name,cen,sig,m1l,m1u,m2l,m2u,j1l,j1u,j2l,j2u,med,ml'
    rec = np.rec.fromrecords(rec, names=names)

    hd = """\
# name : parameter name
# cen  : central value (half way between the marginalised 1 sigma region)
# sig  : 1 sigma error around central value
# m1l  : 1 sigma lower level (marginalised over all other parameters)
# m1u  : 1 sigma upper level (marginalised)
# m2l  : 2 sigma lower level (marginalised) 
# m2u  : 2 sigma upper level (marginalised) 
# j1l  : 1 sigma lower level (joint with all other parameters) 
# j1u  : 1 sigma upper level (joint)
# j2l  : 2 sigma lower level (joint) 
# j2u  : 2 sigma upper level (joint) 
# ml   : maximum likelihood value
# med  : median value
"""
    from barak.io import writetxt
    writetxt('fig/pars.txt', rec, header=hd, fmt_float='.4g', overwrite=1)


if 1:
    ##################################################
    # Read configuration file, set global variables
    ##################################################
    cfgname = 'model.cfg'
    # we only need the cfg file for the prefix of the cloudy runs and
    # the name of the file with the observed column densities.
    opt = parse_config(cfgname)

    testing = 0
    
    MIN_SIG = float(opt['min_sig'])
    SIG_LIMIT = 0.05
    # H is an honorary alpha element here; it just means no offset is
    # added.
    ALPHA_ELEMENTS = 'Si O Mg S Ca Ne Ti H'.split()
    FEPEAK_ELEMENTS = 'Fe Cr Mn Co Ni'.split()


    simnames = sorted(glob(opt['simname']))
    assert len(simnames) > 0
    


if 1:
    ##################################################
    # Read the observed column densities and errors
    ##################################################
    obs = read_Nvals('observed_logN')
    print 'Observed transitions'
    print obs
    trans_obs = sorted(obs)
    # don't do DI or HI
    trans_obs.remove('HI')
    if 'DI' in trans_obs:
        trans_obs.remove('DI')

    priors = read_Nvals('priors')

    if 'NHI' not in priors:
        priors['NHI'] = obs['HI']
    print "Priors found:"
    print priors


    trans = list(trans_obs)
    fh = open('dont_use')
    dont_use =[]
    for row in fh:
        row = row.strip()
        if row == '' or row.startswith('#'):
            continue
        dont_use.append(row)
    fh.close()
    print "Don't use these transitions in fitting:"
    print dont_use
    for tr in dont_use:
        if tr in trans:
            trans.remove(tr)
    print 'Using these transitions'
    print trans


if 1:
    ################################################################
    # Read the cloudy grids and make the interpolators
    ################################################################ 

    tr_plot = opt['tr_plot'].split()
    #tr_plot = ('MgI CaII OI OII OIII OIV MgII FeII FeIII SiII AlII CII AlIII '
    #           'NI NII NIII NIV SiIII SII SIV SV '
    #           'SiIV CIII CIV NV OVI').split()

    ipot = [get_ionization_energy(t) for t in tr_plot]
    isort = np.argsort(ipot)
    tr_plot = [tr_plot[i] for i in isort]
    assert all(tr in tr_plot for tr in trans)
    Ncloudy, Ncloudy_raw, Models, aUV = make_interpolators_uvbtilt(
        tr_plot, simnames)
    M = Models[0]

if 0 and testing:
    # check they look ok

    nrows, ncols = get_nrows_ncols(len(trans) * 2)
    fig = pl.figure(figsize=(8.4, 8.4))
    Z = np.linspace(M.Z[0], M.Z[-1], 100)
    nH = np.linspace(M.nH[0], M.nH[-1], 101)
    nH1, Z1 = np.meshgrid(nH, Z, indexing='ij')
    for i,tr in enumerate(trans):
        ax0 = fig.add_subplot(nrows, ncols, 2*i + 1)
        ax1 = fig.add_subplot(nrows, ncols, 2*i + 2)
        arrplot(Ncloudy_raw[tr].T, x=M.nH, y=M.Z, ax=ax0, colorbar=0)
        z = Ncloudy[tr]((nH1, Z1))
        arrplot(z.T, x=nH, y=Z, ax=ax1, colorbar=0)
        ax0.set_title(tr)

    pl.show()


if 1:
    ######################################################
    # Work out which parameters we're estimating
    ######################################################

    vals = []
    # Only estimate the multipliers we can fit for, based on the observed
    # transitions
    atoms = set(split_trans_name(tr)[0] for tr in trans)
    #if any(atom in ALPHA_ELEMENTS for atom in atoms):
    #    vals.append( ('k_alpha', -1, 1) )
    if 'C' in atoms and 'k_C' not in dont_use:
        vals.append( ('k_C', priors['min k_C'], priors['max k_C']) )
    if 'Al' in atoms and 'k_Al' not in dont_use:
        vals.append( ('k_Al', priors['min k_Al'], priors['max k_Al']) )
    if 'N' in atoms and 'k_N' not in dont_use:
        vals.append( ('k_N' , priors['min k_N'], priors['max k_N']) )
    if any(atom in FEPEAK_ELEMENTS for atom in atoms) \
           and 'k_Fe' not in dont_use:
        vals.append( ('k_Fe', priors['min k_Fe'], priors['max k_Fe']) )

    # These limits will be set later on, based on the grids used to generate
    # the cloudy models.

    vmin = M['NHI'][0] if 'min NHI' not in priors else priors['min NHI']
    vmax = M['NHI'][-1] if 'max NHI' not in priors else priors['max NHI']
    vals.append( ('NHI', vmin, vmax) )
    vmin = M['nH'][0] if 'min nH' not in priors else priors['min nH']
    vmax = M['nH'][-1] if 'max nH' not in priors else priors['max nH']
    vals.append( ('nH', vmin, vmax) )
    vmin = M['Z'][0] if 'min Z' not in priors else priors['min Z']
    vmax = M['Z'][-1] if 'max Z' not in priors else priors['max Z']
    vals.append( ('Z', vmin, vmax) )
    vmin = aUV[0] if 'min aUV' not in priors else priors['min aUV']
    vmax = aUV[-1] if 'max aUV' not in priors else priors['max aUV']
    vals.append( ('aUV', vmin, vmax) )

    print 'min max priors:'
    print zip(*vals)

    P = {}
    P['names'], P['min'], P['max'] = zip(*vals)

    print P

    # dictionary that maps an input atom to the k parameter index,
    # used in model().
    IND_KPAR = {}

    # dictionary mapping parameter names to indices also used in
    # model()
    IND_PAR = {}

    # we take alpha elements to define the metallicity.

    #if 'k_alpha' in P['names']:
    #    i = P['names'].index('k_alpha')
    #    IND_KPAR.update(Si=i, O=i, Mg=i, S=i, Ca=i, Ne=i, Ar=i, Ti=i)
    # Fe peak elements
    if 'k_Fe' in P['names']:
        i = P['names'].index('k_Fe')
        IND_KPAR.update(Fe=i, Cr=i, Mn=i, Co=i, Ni=i)
        IND_PAR.update(k_Fe=i)
    if 'k_C' in P['names']:
        i = P['names'].index('k_C')
        IND_KPAR['C'] = i
        IND_PAR.update(k_C=i)
    if 'k_N' in P['names']:
        i = P['names'].index('k_N')
        IND_KPAR['N'] = i
        IND_PAR.update(k_N=i)
    if 'k_Al' in P['names']:
        i = P['names'].index('k_Al')
        IND_KPAR['Al'] = i
        IND_PAR.update(k_Al=i)

    IND_PAR.update({'Z'  : P['names'].index('Z'),
                    'nH' : P['names'].index('nH'),
                    'NHI' : P['names'].index('NHI')
                    })

    if 'aUV' in P['names']:
        IND_PAR['aUV'] = P['names'].index('aUV')

if 1 and testing:
    # test likelihood (no k or uvtilt)

    Z = np.linspace(M.Z[0], M.Z[-1], 100)
    nH = np.linspace(M.nH[0], M.nH[-1], 101)
    NHI = np.linspace(M.NHI[0], M.NHI[-1], 101)
    Z1, nH1 = np.meshgrid(Z, nH)
    lnprob, alllnprob = ln_likelihood((Z1,nH1,0,0,0,0,0), per_obs=1)
    
    nrows, ncols = get_nrows_ncols(len(trans))
    fig = pl.figure(figsize=(8.4, 8.4))
    for i,tr in enumerate(trans):
        ax = fig.add_subplot(nrows, ncols, i + 1)
        c0 = alllnprob[i] < 0
        if c0.sum():
            vmin = np.percentile(alllnprob[i][c0], 50)
        else:
            vmin = -0.1

        arrplot(alllnprob[i].T, x=nH, y=Z, ax=ax, colorbar=0,vmin=vmin)
        ax.set_title(tr)

if 1 and testing:
    # test total
    fig = pl.figure(figsize=(4.4, 4.4))
    ax = fig.add_subplot(111)
    c0 = alllnprob[i] < 0
    vmin = np.percentile(alllnprob[i][c0], 50)
    arrplot(lnprob.T, x=nH, y=Z, ax=ax, colorbar=0, vmin=vmin)
    pl.show()

if 1 and testing:
    # check the interpoaltion is working by testing a couple of points
    # from the input grid.
    M = Models[0]
    alpha = aUV[0]
    plot(M.NHI, M.N['C'][:,4,4,1], 'o-')
    # see which values we're using
    print M.nH[4], M.Z[4]
    coord = 0,0,0,14, M.nH[4], M.Z[4], alpha; model(coord)[2]
