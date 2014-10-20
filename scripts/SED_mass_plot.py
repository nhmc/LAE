import numpy as np
import pylab as plt
from astropy.io import fits
from barak.utilities import between

mag =dict(u=26.8, g=26.08, r=25.9)
rmg = mag['r'] - mag['g']
gmu = mag['g'] - mag['u']

fh = fits.open('../magphys/daCunha_UDF/UDF_mags_props_forNeil.fits')
d = fh[1].data

Rmag = 0.5*(d.Vmag + d.Imag)
Gmag = 0.5*(d.Bmag + d.Vmag)

# ~700 satisfy c0
c0 = between(d.z, 2.3, 2.7)
# 38 satisfy c0 and c1.
c1 = between(Rmag, mag['r'] -0.4, mag['r'] + 0.4)

Rmag0 = Rmag[c0&c1]
Gmag0 = Gmag[c0&c1]
# could also use Vmag cut instead?

#c2 = between(Rmag-d.Vmag, rmg -0.2, rmg + 0.4)
d0 = d[c0 & c1]


# magphys fitted vals  8.817E+000  9.362E+000  9.727E+000
M16, M50, M84 = 8.817,  9.362,  9.727
# 

#

if 1:
    from barak.plot import hist_xedge, make_log_xlabels, draw_arrows

    xvals = np.linspace(8,11)
    
    fig = plt.figure(figsize=(4.2,7))
    fig.subplots_adjust(top=0.97, left=0.2, bottom=0.12, hspace=1e-3)
    ax1 = plt.subplot(211)
    ax1.fill_between(xvals, 25.9-0.08, y2=25.9+0.08, lw=0, color='0.8')
    ax1.axhline(25.9, ls='--', color='k')
    ax1.plot(d0.mstar, Rmag0, 'ro')
    ax1.set_ylim(25.38, 26.39)
    ax1.set_xlim(8.3, 10.1)
    make_log_xlabels(ax1)
    ax1.set_xticklabels('')
    ax1.set_ylabel('$r$')
    ax1.text(8.6, 25.58, '$2.3 < z < 2.7$', fontsize=18)
    ax2 = plt.subplot(212)
    ax2.fill_between(xvals, rmg-0.16, y2=rmg+0.16, lw=0, color='0.8')
    ax2.axhline(rmg, ls='--', color='k')
    ax2.plot(d0.mstar, Rmag0-Gmag0, 'go')
    c2 = Rmag0-Gmag0 > 0.5
    #draw_arrows(d0.mstar[c2], 0.02, ax=ax2, ms=1, capsize=5, direction='up')
    ax2.plot([M16, M84],[-0.62]*2, color='0.5', lw=2)
    ax2.plot(M50, -0.62,'o', color='0.5', ms=10, mew=0)
    hist_xedge(d0.mstar, ax2,color='k',bins=10, height=0.35)
    ax2.set_ylim(-0.66, 0.095)
    ax2.set_xlim(8.3, 10.1)
    make_log_xlabels(ax2, yoff=-0.08)
    ax2.set_xlabel('$\mathrm{Stellar\ Mass\ (M_{\odot})}$')
    ax2.set_ylabel('$r-g$')
    #plt.show()

if 1:
    # for paper
    xvals = np.linspace(8,11)
    fig = plt.figure(figsize=(4.2, 4.))
    fig.subplots_adjust(top=0.95, left=0.2, bottom=0.18)
    ax = plt.subplot(111)
    #ax.fill_between(xvals, rmg-0.16, y2=rmg+0.16, lw=0, color='0.85')
    #ax.axhline(rmg, ls='-', color='0.5')
    #ax.plot(d0.mstar, Rmag0-Gmag0, 'bo')
    #c2 = Rmag0-Gmag0 > 0.5
    #draw_arrows(d0.mstar[c2], 0.02, ax=ax, ms=1, capsize=5, direction='up')
    #ax.plot([M16, M84],[-0.62]*2, color='r', lw=2)
    #ax.plot(M50, -0.62,'o', color='r', ms=10, mew=0)
    ax.text(9.4, 0.06, '$2.3 < z < 2.7$', fontsize=14)
    ax.text(9.4, 0.005, '$25.5 < r < 26.3$', fontsize=14)
    #hist_xedge(d0.mstar, ax,color='k',bins=10, height=0.35)
    hist_xedge(d0.mstar, ax,color='k',bins=10, height=0.8)
    ax.set_ylim(-0.66, 0.12)
    ax.set_xlim(8.3, 10.1)
    make_log_xlabels(ax, yoff=-0.08)
    ax.set_xlabel('$\mathrm{Stellar\ Mass\ (M_{\odot})}$')
    ax.set_ylabel('$r-g$')
    print np.percentile(d0.mstar, [16, 50, 84])

    #plt.savefig('../pics/SED_mass.png', dpi=200)
    #plt.savefig('../pics/SED_mass.pdf')
    #from subprocess import call
    #call('pdf2ps ../pics/SED_mass.pdf', shell=1)
    #call('ps2pdf ../pics/SED_mass.ps', shell=1)
    plt.show()
    
