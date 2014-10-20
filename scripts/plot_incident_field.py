import sys
from glob import glob
import pylab as pl
plt = pl
import numpy as np

from cloudy.utils import read_starburst99


#names = sys.argv[1:]

usage = """usage: python2.7 plot_incident_field.py */cloudy*tbl
"""

import numpy as np

import astropy.units as u
import astropy.constants as c

from barak.absorb import get_ionization_energy
from barak.plot import puttext, axvlines

plt.rc('xtick.major',pad=8)
plt.rc('ytick.major',pad=7)

wa, logFlam = read_starburst99('/Users/ncrighton/Code/Python/cloudy/data/'
                         'starburst.spectrum1')

from barak.sed import flambda_to_fnu
logfnu = np.log10(flambda_to_fnu(wa, 10**logFlam))
wa1 = 10**np.linspace(np.log10(wa[0]), np.log10(wa[-1]), 10*len(wa))
logfnu1 = np.interp(wa1, wa, logfnu)
starburst = np.rec.fromarrays(
    [wa1, logfnu1], names='wa,logfnu')

trans = ('CII CIII CIV SiII SiIII SiIV HI AlII '
         'AlIII NV OVI MgII NII').split()  # MgI

trans = ('CII CIII CIV SiII SiIII SiIV HI AlII '
         'AlIII NV OVI MgII NII FeII FeIII MgI OI OII OIII '
         'CaII NI CI').split()  # MgI

IP = np.array(get_ionization_energy(trans))

import pylab as plt
from atpy import Table

prefix = "/Users/ncrighton/Projects/MPIA_QSO_LBG/Cloudy/J0004/comp1"

names = sorted(glob(prefix + "/uvb_k0[0246]/*tilted*tbl"))

fig1 = pl.figure(figsize=(5.1,4.7))
fig1.subplots_adjust(left=0.2, bottom=0.15, top=0.96,right=0.96)
fig2 = pl.figure(figsize=(5.1,4.7))
fig2.subplots_adjust(left=0.2, bottom=0.15, top=0.96,right=0.96)
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)


for ax in (ax1,ax2):
    N0 = -21
    for i,n in enumerate(names):
        T = Table(n)
        ind = T.energy.searchsorted(0.99)
        if i == len(names) - 1:
            ref = T
            ref.log10jnu = T.log10jnu - T.log10jnu[ind] + N0
        color = '0.5' if i != 2 else 'k'
        ax.semilogx(T.energy, T.log10jnu - T.log10jnu[ind] + N0, '-', color=color)

puttext(0.04, 0.07, '$\mathrm{Haardt\ & \ Madau}$\n$2012,\ z=2.5$',
        ax=ax1, ha='left')

# plot a QSO spectrum

X = 10**np.linspace(-1, 4, 1000)   # hv
# alpha = -1.7
# fnu = np.log10(X**alpha)
# ind = X.searchsorted(0.99)
# plt.semilogx(X, fnu - fnu[ind] + N0, '-' + 'c')
#alpha = -1.5
# Scott
alpha = -0.56
fnu = np.log10(X**alpha)
ind = X.searchsorted(0.99)
ax2.semilogx(X, fnu - fnu[ind] + N0, '--c', lw=2,
             label=r'$J_\nu\propto\ \nu^{%.3g}$' % alpha)


# plot the starburst spectrum
x = (starburst.wa * u.AA).to(u.rydberg, equivalencies=u.equivalencies.spectral())
x = x.value
X = x[::-1]
F = starburst.logfnu[::-1]
ind = X.searchsorted(1.02)

from barak.convolve import convolve_psf

ax2.semilogx(X, convolve_psf(F - F[ind] + N0,50), 'r', lw=2,
             label='$\mathrm{starburst}$')

# apply dust extinction

from barak.extinction import starburst_Calzetti00
ext = starburst_Calzetti00(starburst.wa, EBmV=0.2)
Fnu = np.log10(10**starburst.logfnu * np.exp(-ext.tau))
F = Fnu[::-1]
ax2.semilogx(X, convolve_psf(F - F[ind] + N0, 50), 'r', lw=0.5,
             label='$\mathrm{starburst},\ E(B-V)=0.2$')

ax2.legend(frameon=0, loc='lower left', fontsize=12)
#plt.grid(ls='solid', lw=0.5, color='0.7')
xvals = (IP * u.eV).to(u.Ry).value

y = np.interp(xvals, ref.energy, ref.log10jnu)

from barak.absorb import split_trans_name

ax1.vlines(xvals, y + 0.1, y + 0.3, lw=0.5, color='0.5')#, ls=':')
for i in xrange(len(xvals)):
    atom,stage = split_trans_name(trans[i])
    t = atom + ' ' + stage
    align = 'left'
    if trans[i] in 'SiIV AlIII'.split():
        align = 'center'
    ax1.text(xvals[i], y[i]+0.35, t, fontsize=9,
            rotation=90, ha=align, va='bottom')

c = 'k'
puttext(0.79, 0.59, r'$\alpha_\mathrm{UV}=1$', color=c,ax=ax1, rotation=-5)
puttext(0.78, 0.41, r'$\alpha_\mathrm{UV}=0$',  color=c,ax=ax1, rotation=-13)
puttext(0.76, 0.30, r'$\alpha_\mathrm{UV}=-1$',color=c,ax=ax1, rotation=-15)
puttext(0.76, 0.12, r'$\alpha_\mathrm{UV}=-2$',color=c,ax=ax1, rotation=-18)

for ax in (ax1,ax2):
    ax.set_xlabel('$\mathrm{Energy\ (Rydbergs)}$')
    ax.set_ylabel(r'$\log_{10}\ [J_{\nu}\ \mathrm{(erg/s/cm^2/Hz/ster)}]$')
    #ax.set_xlim(0.7, 14)
    ax.set_xlim(0.4, 14)
    ax.set_xticklabels(['', '', '$\mathdefault{1}$', '$\mathdefault{10}$'])
    #ax.set_ylim(-24.9, -20.1)
    ax.set_ylim(-24.9, -17.1)
    #plt.savefig('aUV_curve.png', dpi=200)

fig1.savefig('aUV_curve.pdf')
fig2.savefig('aUV_QSO_starburst.pdf')
#plt.show()
