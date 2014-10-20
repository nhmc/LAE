""" Perform calculations on Spectral Energy Distributions (SEDs).

Inspired by the SED module in astLib by Matt Hilton
(http://astlib.sourceforge.net)

- VEGA: The SED of Vega, used for calculation of magnitudes on the Vega system.
- AB: Flat spectrum SED, used for calculation of magnitudes on the AB system.
- SUN: The SED of the Sun.

"""
from __future__ import division, print_function, unicode_literals
try:
    unicode
except NameError:
    unicode = basestring = str
    xrange = range

from .io import readtabfits, loadtxt
from .constants import c, c_kms, Jy
from .utilities import get_data_path

import numpy as np
from numpy.random import randn

import matplotlib.pyplot as pl

import os, math
import warnings

DATAPATH = get_data_path()
PATH_PASSBAND = DATAPATH + '/passbands/'
PATH_EXTINCT = DATAPATH + '/atmos_extinction/'
PATH_TEMPLATE = DATAPATH + '/templates/'

def _listfiles(topdir):
    names = [n for n in os.listdir(topdir) if os.path.isdir(topdir + n)]
    files = dict([(n, []) for n in names])
    for name in sorted(names):
        for n in sorted(os.listdir(topdir + name)):
            if n != 'README' and not os.path.isdir(topdir + name + '/'+ n) and \
                   not n.startswith('effic') and \
                   not n.endswith('.py') and not n.endswith('.pdf'):
                files[name].append(n)
    return files

TEMPLATES = _listfiles(PATH_TEMPLATE)
PASSBANDS = _listfiles(PATH_PASSBAND)

def get_bands(instr=None, names=None, ccd=None):
    """ Get one or more passbands by giving the instrument and
    filename.

    If `names` is not given, then every passband for that instrument
    is returned.  Passband instruments and filenames are listed in the
    dictionary PASSBANDS. names can be a list, a single string, or a
    comma-separated string of values.

    Examples
    --------

    >>> sdss = get_bands('SDSS', 'u,g,r,i,z')  # get the SDSS passbands
    >>> U = get_bands('LBC', 'u')    # get the LBC U_spec filter
    """
    if instr is None:
        return _listfiles(PATH_PASSBAND)
    if isinstance(names, basestring):
        if ',' in names:
            names = [n.strip() for n in names.split(',')]
        else:
            return Passband(instr + '/' + names)
    elif names is None:
        names = PASSBANDS[instr]

    return [Passband(instr + '/' + n, ccd=ccd) for n in names]

def get_SEDs(kind=None, names=None):
    """ Get one or more SEDs based on their type and filename

    If `names` is not given, then every SED of that type is returned.
    SED types and filenames are listed in the dictionary TEMPLATES.

    Examples
    --------
    >>> pickles = get_SEDs('pickles')   # pickles stellar library SEDs
    >>> lbga = get_SEDs('LBG', 'lbg_abs.dat')  # LBG absorption spectrum
    """
    if kind is None:
        return _listfiles(PATH_TEMPLATE)
    if isinstance(names, basestring):
        if ',' in names:
            names = [n.strip() for n in names.split(',')]
        else:
            return SED(kind + '/' + names)
    elif names is None:
        names = TEMPLATES[kind]

    return [SED(kind + '/' + n) for n in names]


class Passband(object):
    """This class describes a filter transmission curve. Passband
    objects are created by loading data from from text files
    containing wavelength in angstroms in the first column, relative
    transmission efficiency in the second column (whitespace
    delimited). For example, to create a Passband object for the 2MASS
    J filter:

    passband = Passband('J_2MASS.res')

    where 'J_2MASS.res' is a file in the current working directory
    that describes the filter.

    The available passbands are in PASSBANDS.

    Attributes
    ----------
    wa : array of floats
      Wavelength in Angstroms
    tr : array of floats
      Normalised transmission, including atmospheric extinction and
      detector efficiency. May or may not include extinction from the
      optical path.
    effective_wa : float
      Effective wavelength of the passband.

    Methods
    -------
    plot
    """
    def __init__(self, filename, ccd=None):
        if not filename.startswith(PATH_PASSBAND):
            filepath = PATH_PASSBAND + filename
        else:
            filepath = filename
        if filepath.endswith('.fits'):
            try:
                import pyfits
            except ImportError:
                import astropy.io.fits as pyfits
            rec = pyfits.getdata(filepath, 1)
            self.wa, self.tr = rec.wa, rec.tr
        else:
            self.wa, self.tr = loadtxt(filepath, usecols=(0,1), unpack=True)
        # check wavelengths are sorted lowest -> highest
        isort = self.wa.argsort()
        self.wa = self.wa[isort]
        self.tr = self.tr[isort]

        # get the name of the filter/passband file and the name of the
        # directory in which it lives (the instrument).
        prefix, filtername = os.path.split(filename)
        _, instr = os.path.split(prefix)
        self.name = filename

        if instr == 'LBC' and ccd is None:
            if filtername.startswith('LBCB') or filtername in 'ug':
                ccd = 'blue'
            elif filtername.startswith('LBCR') or filtername in 'riz':
                ccd = 'red'
        elif instr == 'FORS' and ccd is None:
            warnings.warn('No cdd ("red" or "blue") given, assuming red.')
            ccd = 'red'

        self.atmos = self.effic = None
        if ccd is not None:
            # apply ccd/optics efficiency
            name = PATH_PASSBAND + instr + '/effic_%s.txt' % ccd
            wa, effic = loadtxt(name, usecols=(0,1), unpack=1)
            self.effic = np.interp(self.wa, wa, effic)
            self.tr *= self.effic

        extinctmap = dict(LBC='kpno_atmos.dat', FORS='paranal_atmos.dat',
                          HawkI='paranal_atmos.dat',
                          KPNO_Mosaic='kpno_atmos.dat',
                          CTIO_Mosaic='ctio_atmos.dat')

        if instr in extinctmap:
            # apply atmospheric extinction
            wa, emag = loadtxt(PATH_EXTINCT + extinctmap[instr], unpack=1)
            self.atmos = np.interp(self.wa, wa, 10**(-0.4 * emag))
            self.tr *= self.atmos

        # trim away areas where band transmission is negligibly small
        # (<0.01% of peak transmission).
        isort = self.tr.argsort()
        sortedtr = self.tr[isort]
        maxtr = sortedtr[-1]
        imax = isort[-1]
        ind = isort[sortedtr < 1e-4 * maxtr]
        if len(ind) > 0:
            i = 0
            c0 = ind < imax
            if c0.any():
                i = ind[c0].max()
            j = len(self.wa) - 1
            c0 = ind > imax
            if c0.any():
                j = ind[c0].min()
            i = min(abs(i-2), 0)
            j += 1
            self.wa = self.wa[i:j]
            self.tr = self.tr[i:j]
        if self.atmos is not None:
            self.atmos = self.atmos[i:j]
        if self.effic is not None:
            self.effic = self.effic[i:j]

        # normalise
        self.ntr = self.tr / np.trapz(self.tr, self.wa)

        # Calculate the effective wavelength for the passband. This is
        # the same as equation (3) of Carter et al. 2009.
        a = np.trapz(self.tr * self.wa)
        b = np.trapz(self.tr / self.wa)
        self.effective_wa = math.sqrt(a / b)

        # find the AB and Vega magnitudes in this band for calculating
        # magnitudes.
        self.flux = {}
        self.flux['Vega'] = VEGA.calc_flux(self)
        self.flux['AB'] = AB.calc_flux(self)

    def __repr__(self):
        return 'Passband "%s"' % self.name

    def plot(self, effic=False, atmos=False, ymax=None, **kwargs):
        """ Plots the passband. We plot the non-normalised
        transmission. This may or may not include ccd efficiency,
        losses from the atmosphere and telescope optics.
        """
        tr = self.tr
        if ymax is not None:
            tr = self.tr / self.tr.max() * ymax
        pl.plot(self.wa, tr, **kwargs)
        if self.effic is not None and effic:
            pl.plot(self.wa, self.effic,
                    label='applied ccd efficiency', **kwargs)
        if self.atmos is not None and atmos:
            pl.plot(self.wa, self.atmos,
                    label='applied atmospheric extinction', **kwargs)

        pl.xlabel("Wavelength ($\AA$)")
        pl.ylabel("Transmission")
        if atmos or effic:
            pl.legend()
        if pl.isinteractive():
            pl.show()

class SED(object):
    """A Spectral Energy Distribution (SED).

    Instantiate with either a filename or a list of wavelengths and fluxes.
    Wavelengths must be in Angstroms, fluxes in erg/s/cm^2/Ang.

    To convert from f_nu to f_lambda in erg/s/cm^2/Ang, substitute
    using::

     nu = c / lambda
     f_lambda = c / lambda^2 * f_nu

    Available SED template filenames are in TEMPLATES.
    """
    def __init__(self, filename=None, wa=[], fl=[], z=0., label=None):

        # filename overrides wave and flux keywords
        if filename is not None:
            if not filename.startswith(PATH_TEMPLATE):
                filepath = PATH_TEMPLATE + filename
            if filepath.endswith('.fits'):
                rec = readtabfits(filepath)
                wa, fl = rec.wa, rec.fl
            else:
                wa, fl = loadtxt(filepath, usecols=(0,1), unpack=1)
            if label is None:
                label = filename

        # We keep a copy of the wavelength, flux at z = 0
        self.z0wa = np.array(wa)
        self.z0fl = np.array(fl)
        self.z0fl_noextinct = np.array(fl)
        self.EBmV = None

        self.wa = np.array(wa)
        self.fl = np.array(fl)
        self.z = z
        self.label = label

        if abs(z) > 1e-6:
            self.redshift_to(z)


    def __repr__(self):
        return 'SED "%s"' % self.label

    def copy(self):
        """Copies the SED, returning a new SED object.
        """
        newSED = SED(wa=self.z0wa, fl=self.z0fl, z=self.z, label=self.label)
        return newSED

    def integrate(self, wmin=None, wmax=None):
        """ Calculates flux (erg/s/cm^2) in SED within given wavelength
        range."""
        if wmin is None:
            wmin = self.wa[0]
        if wmax is None:
            wmax = self.wa[-1]

        i,j = self.wa.searchsorted([wmin, wmax])
        fl = np.trapz(self.fl[i:j], self.wa[i:j])

        return fl

    def plot(self, log=False, ymax=None, **kwargs):
        fl = self.fl
        if ymax is not None:
            fl = self.fl / self.fl.max() * ymax

        if log:
            pl.loglog(self.wa, fl, **kwargs)
        else:
            pl.plot(self.wa, fl, **kwargs)
        pl.xlabel('Wavelength ($\AA$)')
        pl.ylabel('Flux (ergs s$^{-1}$cm$^{-2}$ $\AA^{-1}$)')
        #pl.legend()
        if pl.isinteractive():
            pl.show()

    def redshift_to(self, z, cosmo=None):
        """Redshifts the SED to redshift z. """
        # We have to conserve energy so the area under the redshifted
        # SED has to be equal to the area under the unredshifted SED,
        # otherwise magnitude calculations will be wrong when
        # comparing SEDs at different zs

        self.wa = np.array(self.z0wa)
        self.fl = np.array(self.z0fl)

        z0fluxtot = np.trapz(self.z0wa, self.z0fl)
        self.wa *= z + 1
        zfluxtot = np.trapz(self.wa, self.fl)
        self.fl *= z0fluxtot / zfluxtot
        self.z = z

    def normalise_to_mag(self, ABmag, band):
        """Normalises the SED to match the flux equivalent to the
        given AB magnitude in the given passband.
        """
        magflux = mag2flux(ABmag, band)
        sedflux = self.calc_flux(band)
        norm = magflux / sedflux
        self.fl *= norm
        self.z0fl *= norm

    def calc_flux(self, band):
        """Calculate the mean flux for a passband, weighted by the
        response and wavelength in the given passband.

        Returns the mean flux (erg/s/cm^2/Ang) inside the band.
        """
        if self.wa[0] > band.wa[0] or self.wa[-1] < band.wa[-1]:
            msg = "SED does not cover the whole bandpass, extrapolating"
            warnings.warn(msg)
            dw = np.median(np.diff(self.wa))
            sedwa = np.arange(band.wa[0], band.wa[-1]+dw, dw)
            sedfl = np.interp(sedwa, self.wa, self.fl)
        else:
            sedwa = self.wa
            sedfl = self.fl

        i,j = sedwa.searchsorted([band.wa[0], band.wa[-1]])
        fl = sedfl[i:j]
        wa = sedwa[i:j]

        dw_band = np.median(np.diff(band.wa))
        dw_sed = np.median(np.diff(wa))
        if dw_sed > dw_band and dw_band > 20:
            warnings.warn(
                'WARNING: SED wavelength sampling interval ~%.2f Ang, '
                'but bandpass sampling interval ~%.2f Ang' %
                (dw_sed, dw_band))
            # interpolate the SED to the passband wavelengths
            fl = np.interp(band.wa, wa, fl)
            band_tr = band.tr
            wa = band.wa
        else:
            # interpolate the band transmission to the SED
            # wavelength values.
            band_tr = np.interp(wa, band.wa, band.tr)

        # weight by response and wavelength, appropriate when we're
        # counting the number of photons within the band.
        flux = np.trapz(band_tr * fl * wa, wa) / np.trapz(band_tr * wa, wa)
        
        return flux

    def calc_mag(self, band, system="AB"):
        """Calculates magnitude in the given passband.

        Note that the distance modulus is not added.

        mag_sigma : float
          Add a gaussian random deviate to the magnitude, with sigma
          given by this value.

        `system` is either 'Vega' or 'AB'
        """
        f1 = self.calc_flux(band)

        if f1 > 0:
            mag = -2.5 * math.log10(f1/band.flux[system])
            if system == "Vega":
                # Add 0.026 because Vega has V=0.026 (e.g. Bohlin & Gilliland 2004)
                mag += 0.026
        else:
            mag = np.inf

        return mag

    def calc_colour(self, band1, band2, system="AB"):
        """Calculates the colour band1 - band2.

        system is either 'Vega' or 'AB'.

        mag_sigma : float
          Add a gaussian random deviate to each magnitude, with sigma
          given by this value.
        """
        mag1 = self.calc_mag(band1, system=system)
        mag2 = self.calc_mag(band2, system=system)

        return mag1 - mag2

    def apply_extinction(self, ext_type, EBmV):
        """ Return a new SED instance with the extinction curve
        applied at the same redshift as the template.

        Allowed extinction laws are:

           MW
           SMC
           LMC
           starburst

        See the extinction module for more information.
        """
        import extinction as ext
        ecurve = dict(MW= ext.MW_Cardelli89             ,
                      SMC=ext.SMC_Gordon03              ,
                      LMC=ext.LMC_Gordon03              ,
                      starburst=ext.starburst_Calzetti00
                      )

        sed = self.copy()
        sed.z0fl[:] = sed.z0fl_noextinct
        tau = ecurve[ext_type](self.z0wa, EBmV=EBmV).tau
        sed.z0fl *= np.exp(-tau)
        sed.EBmV = EBmV

        sed.redshift_to(sed.z)
        return sed

def mag2flux(ABmag, band):
    """ Converts given AB magnitude into flux in the given band, in
    erg/s/cm^2/Angstrom.

    Returns the flux in the given band.
    """
    # AB mag (See Oke, J.B. 1974, ApJS, 27, 21)
    # fnu in erg/s/cm^2/Hz
    fnu = 10**(-(ABmag + 48.6)/2.5)
    # convert to erg/s/cm^2/Ang
    flambda = fnu_to_flambda(band.effective_wa, fnu)

    return flambda

def flux2mag(flambda, band):
    """Converts flux in erg/s/cm^2/Angstrom into AB magnitudes.

    Returns the magnitude in the given band.
    """
    # convert to erg/s/cm^2/Hz
    fnu = flambda_to_fnu(band.effective_wa, flambda)
    mag = -2.5*math.log10(fnu) - 48.6
    return mag

def mag2Jy(ABmag):
    """Converts an AB magnitude into flux density in Jy (fnu).
    """
    flux_nu = 10**(-(ABmag + 48.6)/2.5) / Jy
    return flux_nu

def Jy2mag(fluxJy):
    """Converts flux density in Jy into AB magnitude (fnu).
    """
    ABmag = -2.5 * (np.log10(fluxJy * Jy)) - 48.6
    return ABmag

def fnu_to_flambda(wa, f_nu):
    """ Convert flux per unit frequency to a flux per unit wavelength.

    Parameters
    ----------
    wa : array_like
      Wavelength in Angstroms
    f_nu : array_like
      Flux at each wavelength in erg/s/cm^2/Hz

    Returns
    -------
    f_lambda : ndarray
      Flux at each wavelength in erg/s/cm^2/Ang
    """
    return c / (wa * 1e-8)**2 * f_nu * 1e-8

def flambda_to_fnu(wa, f_lambda):
    """ Convert flux per unit wavelength to a flux per unit frequency.

    Parameters
    ----------
    wa : array_like
      Wavelength in Angstroms
    f_lambda : array_like
      Flux at each wavelength in erg/s/cm^2/Ang

    Returns
    -------
    f_nu : ndarray
      Flux at each wavelength in erg/s/cm^2/Hz
    """
    return (wa * 1e-8)**2 * f_lambda * 1e8 / c


def qso_template(wa, z):
    """ Return a composite QSO spectrum at redshift z.

    This uses the SDSS composite at wa > 1680 and a smoothed version
    of the HST/COS EUV+FUV AGN composite spectrum shown in Figure 5
    from Shull, Stevans, and Danforth 2012 for wa < 1680.

    Parameters
    ----------
    wa : array_like, shape (N,)
      Wavelength array in Angstroms
    z : float
      Redshift.
    
    Returns
    -------
    f_lambda : ndarray, shape (N,)
      The QSO spectrum in F_lambda (the normalisation is arbitrary).
    """
    wa = np.array(wa, copy=False)
    wrest = wa / (1+z)
    i = wrest.searchsorted(1680)
    if i == len(wrest):
        return qso_template_uv(wa, z)
    elif i == 0:
        return qso_template_sdss(wa, z)
    else:
        fl = np.ones(len(wa), float)
        f = qso_template_uv(wa, z)
        fl[:i] = f[:i] / f[i]
        f = qso_template_sdss(wa, z)
        fl[i:] = f[i:] / f[i]

    return fl

def qso_template_sdss(wa, z):
    """ Return a composite visible QSO spectrum at redshift z.

    The SDSS composite spectrum as a function of F_lambda is returned
    at each wavelength wa. wa must be in Angstroms.

    Only good between 700 and 8000 Angstroms (rest frame).
    """
    T = readtabfits(DATAPATH + '/templates/qso/dr1QSOspec.fits')
    return np.interp(wa, T.wa*(1+z), T.fl)

def qso_template_uv(wa, z):
    """ Return a composite UV QSO spectrum at redshift z.

    Wavelengths must be in Angstroms.

    This is a smoothed version of the HST/COS EUV+FUV AGN composite
    spectrum shown in Figure 5 of Shull, Stevans, and Danforth 2012.

    Only good between 550 and 1730 Angstroms (rest frame)
    """
    T = readtabfits(DATAPATH + 'templates/qso/Shull_composite.fits')
    return np.interp(wa, T.wa*(1 + z), T.fl)


def make_constant_dv_wa_scale(wmin, wmax, dv):
    """ Make a wavelength scale with bin widths corresponding to a
    constant velocity width scale.

    Parameters
    ----------
    wmin, wmax : floats
      Start and end wavelength.
    dv : float
      Velocity pixel width.

    Returns
    -------
    wa : ndarray
      Wavelength scale. 

    See Also
    --------
    barak.spec.make_wa_scale
    """
    dlogw = np.log10(1 + dv/c_kms)
    # find the number of points needed.
    npts = int(np.log10(wmax / wmin) / dlogw)
    wa = wmin * 10**(np.arange(npts)*dlogw)
    return wa

def vel_from_wa(wa, wa0, redshift=0):
    """ Find a velocity scale from wavelengths, given a redshift and
    transition rest wavelength.

    Parameters
    ----------
    wa : array, shape (N,)
      Observed wavelength array.
    wa0 : float
      Transition rest wavelength, must be same units as `wa`.
    redshift : float
      Redshift. Default 0.

    Returns
    -------
    vel : arr of shape (N,)
      Velocities in km/s.
    """
    obswa = wa0 * (1 + redshift)
    return (wa / obswa - 1) * c_kms

# Data
VEGA = SED('reference/Vega_bohlin2006')
SUN = SED('reference/sun_stis')

# AB SED has constant flux density (f_nu) 3631 Jy, see
# http://www.sdss.org/dr5/algorithms/fluxcal.html

fnu = 3631 * Jy   # erg/s/cm^2/Hz
wa = np.logspace(1, 10, 1e5)   # Ang
AB = SED(wa=wa, fl=fnu_to_flambda(wa, fnu))

# don't clutter the namespace
del wa, fnu
