import sys
from glob import glob
import pylab as pl
plt = pl
import numpy as np

#names = sys.argv[1:]

usage = """usage: python2.7 plot_incident_fieldpy */cloudy*tbl
"""

import numpy as np

import astropy.units as u
import astropy.constants as c

import barak.absorb
from barak.absorb import get_ionization_energy
from barak.plot import puttext, axvlines, make_log_xlabels

plt.rc('xtick.major',pad=8)
plt.rc('ytick.major',pad=7)

from barak.sed import flambda_to_fnu

get_ionization_energy('HI')

T = barak.absorb.ION_CACHE['table'].data
prob = 8

c0 = (T['P'] > prob) & (T['wrest'] > 600)
trans = T['species'][c0]
tr_wrest = T['wrest'][c0]
# trans = ('CII CIII CIV SiII SiIII SiIV HI AlII '
#          'AlIII NV OVI MgII NII').split()  # MgI

# trans = ('CI CII CIII CIV SiII SiIII SiIV HI AlII '
#          'AlIII NI NII NIII NIV NV MgII FeII FeIII MgI OI OII OIII OIV OV OVI '
#          'CaII MnII CrII NeVIII ZnII MgX').split()  # MgI

IP = np.array(get_ionization_energy(trans))
isort = IP.argsort()
IP = IP[isort]
trans = [trans[i] for i in isort]
tr_wrest = [tr_wrest[i] for i in isort]


import pylab as plt
from atpy import Table

prefix = "/Users/ncrighton/Projects/MPIA_QSO_LBG/Cloudy/J0004/comp1"

name, = glob(prefix + "/uvb_k04/*tilted*tbl")

fig1 = pl.figure(figsize=(7.1,4.7))
fig1.subplots_adjust(left=0.15, bottom=0.15, top=0.96,right=0.96)
ax = fig1.add_subplot(111)


N0 = -21
ref = Table(name)
ind = ref.energy.searchsorted(0.99)
ref.log10jnu = ref.log10jnu - ref.log10jnu[ind] + N0
ax.plot(np.log10(ref.energy), ref.log10jnu - ref.log10jnu[ind] + N0, '-', color='k')

xvals = np.log10((IP * u.eV).to(u.Ry).value)

y = np.interp(xvals, np.log10(ref.energy), ref.log10jnu)

from barak.absorb import split_trans_name


xdiff = np.diff(xvals)

direction = ['up']
previousdown = False
for i in xrange(len(xvals) - 1):
    if xdiff[i] < 0.02 and not previousdown:
        direction.append('down')
        previousdown = True
    else:
        direction.append('up')
        previousdown = False
direction = np.array(direction)


ind = np.arange(len(xvals))[direction =='up']
previouslong = False
length = np.array(['short'] * len(xvals))
for i in xrange(len(ind) - 1):
    if (xvals[ind[i]] - xvals[ind[i-1]]) < 0.04 and not previouslong:
        length[ind[i]] = 'long'
        previouslong = True
    else:
        previouslong = False


ind = np.arange(len(xvals))[direction =='down']
previouslong = False
for i in xrange(len(ind) - 1):
    i += 1
    if (xvals[ind[i]] - xvals[ind[i-1]]) < 0.04 and not previouslong:
        length[ind[i]] = 'long'
        previouslong = True
    else:
        previouslong = False

offset = 0.003
xoff = {t:offset for t in trans}
#xoff.update(SiIV=0, CIII=0.01)
for i in xrange(len(xvals)):
    tr = trans[i]
    atom,stage = split_trans_name(tr)
    t = atom + ' ' + stage
    halign = 'center'
    valign = 'bottom'
    dx = xoff[tr]
    dy0 = 0.1
    dy1 = 0.4
    if length[i] == 'long':
        dy1 = 1.2

    dy2 = dy1 + 0.05
    if direction[i] == 'down':
        dy0 = -dy0
        dy1 = -dy1
        dy2 = -dy2
        valign = 'top'

    color='b'
    if tr_wrest[i] > 1215.67:
        color = 'r'
    elif tr_wrest[i] > 912.:
        color = 'g'
    
    ax.text(xvals[i] + dx, y[i]+ dy2, t, fontsize=9,
            rotation=90, ha=halign, va=valign, color=color)
    ax.plot([xvals[i]]*2, [y[i] + dy0, y[i] + dy1], lw=0.5, color='0.5')

c = 'k'

ax.set_xlabel('$\mathrm{Energy\ (Rydbergs)}$')
ax.set_ylabel(r'$\log_{10}\ [J_{\nu}\ \mathrm{(erg/s/cm^2/Hz/ster)}]$')
ax.set_xlim(-0.6, 1.6)
make_log_xlabels(ax)
ax.set_ylim(-24.3, -16.8)

puttext(0.95, 0.85, '$\mathrm{Haardt\ & \ Madau}$\n$2012,\ z=2.5$',
        ax=ax, ha='right')
puttext(0.95, 0.75, '$\mathrm{Strongest\ }\lambda_0>\ 1216\ \AA$', ax=ax, ha='right',color='r',fontsize=10)
puttext(0.95, 0.7,  '$\mathrm{Strongest\ }\lambda_0>\ 912\ \AA$', ax=ax, ha='right',color='g' ,fontsize=10)
puttext(0.95, 0.65, '$\mathrm{Strongest\ }\lambda_0>\ 600\ \AA$', ax=ax, ha='right',color='b' ,fontsize=10)

puttext(0.95, 0.55, '$\mathrm{Verner\ }P\ >\ %g$' % prob, ax=ax, ha='right',fontsize=10)

fig1.savefig('incident_field_trans.pdf')
