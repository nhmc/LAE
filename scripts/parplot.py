import numpy as np
from StringIO import StringIO
from barak.io import readtxt
import pylab as plt

def readcol(col):
    temp = []
    for row in col:
        limit = 'none'
        erlo = erhi = 0
        r = row.strip('$')
        strong_prior = False
        if r.startswith('\\mathbf{'):
            r = r[9:-1]
            strong_prior = True
        if '(' in r:
            r = r.split('(')[0]
        if r.startswith('>'):
            limit = 'lower'
            val = r.strip('>')
        elif r.startswith('<'):
            limit = 'upper'
            val = r.strip('<')
        else:
            val,er = r.split('\\pm')
            erlo = erhi = float(er)

        temp.append((float(val), erlo, erhi, limit, strong_prior))
    return np.rec.fromrecords(temp, names='val,erlo,erhi,limit,strong_prior')


# You need to remove the backslashes at the end of the line!



tab=r"""
 1 & -112 &$-0.70\pm0.14$ &$-1.14\pm0.31$ &$-2.14\pm0.12$ &$\mathbf{14.88\pm0.01}$ &$18.18\pm0.16$              &$4.19\pm0.03$ &$-2.85\pm0.33(0.14)$&$1.33\pm0.32(0.12)$ &$-0.46\pm0.42(0.30)$ &$3.55\pm0.96(0.75)$
 2 & -84 &$-0.39\pm0.12$ &$-0.15\pm0.27$ &$-2.54\pm0.09$ &$\mathbf{15.23\pm0.01}$ &$18.02\pm0.14$              &$4.14\pm0.02$ &$-2.37\pm0.32(0.12)$&$1.78\pm0.32(0.10)$ &$-1.10\pm0.39(0.25)$ &$2.07\pm0.88(0.64)$
 3 & -30 &$-0.34\pm0.37$ &$\mathbf{-0.35\pm0.42}$ &$-1.83\pm0.19$ &$\mathbf{14.47\pm0.26}$ &$17.96\pm0.40$     &$4.24\pm0.10$ &$-3.10\pm0.38(0.23)$&$1.06\pm0.37(0.22)$ &$-0.45\pm0.65(0.57)$ &$3.40\pm1.64(1.53)$
 4 & 51 &$-0.35\pm0.11$ &$0.30\pm0.15$  &$-2.73\pm0.12$ &$\mathbf{16.43\pm0.09}$ &$18.97\pm0.17$              &$4.14\pm0.03$ &$-2.11\pm0.33(0.13)$&$2.03\pm0.32(0.12)$ &$-0.40\pm0.41(0.28)$ &$4.45\pm0.94(0.72)$
 5 & 72 &$-0.21\pm0.11$ &$0.37\pm0.14$  &$-2.70\pm0.11$ &$\mathbf{16.37\pm0.09}$ &$18.90\pm0.16$              &$4.11\pm0.03$ &$-2.15\pm0.32(0.12)$&$1.97\pm0.32(0.12)$ &$-0.42\pm0.39(0.26)$ &$4.34\pm0.89(0.66)$
 6 & 92 &$-0.26\pm0.11$ &$0.20\pm0.18$  &$-2.66\pm0.10$ &$\mathbf{16.02\pm0.10}$ &$18.62\pm0.13$              &$4.12\pm0.03$ &$-2.21\pm0.32(0.11)$&$1.91\pm0.32(0.12)$ &$-0.65\pm0.36(0.21)$ &$3.59\pm0.80(0.53)$
 7 & 235 &$-0.69\pm0.23$ &$\mathbf{0.01\pm0.38}$  &$-2.43\pm0.17$ &$\mathbf{14.70\pm0.01}$ &$17.64\pm0.25$     &$4.27\pm0.04$ &$-2.46\pm0.37(0.21)$&$1.80\pm0.35(0.17)$ &$-1.38\pm0.55(0.46)$ &$1.17\pm1.31(1.16)$
 8 & 324 &$-1.07\pm0.42$ &$\mathbf{-0.24\pm0.48}$ &$-1.63\pm0.44$ &$\mathbf{13.59\pm0.01}$ &$17.76\pm0.57$     &$4.49\pm0.12$ &$-3.34\pm0.55(0.46)$&$1.09\pm0.45(0.34)$ &$-0.22\pm1.08(1.04)$ &$3.60\pm2.71(2.64)$
"""

fh = StringIO(tab)

#T = readtxt(fh, sep='&', names=['comp', 'nH', 'Z', 'aUV', 'NHI',
#                                'U', 'T', 'NH', 'D', 'Pk'])

T = readtxt(fh, sep='&', names=['comp', 'vel', 'Z', 'aUV', 'U', 'NHI', 'NH',
                                'T', 'nH', 'Pk', 'D', 'M'])

from barak.utilities import adict
V = adict()
for n in ('nH', 'Z', 'aUV','NHI','NH', 'U','T','D','Pk', 'M'):
    V[n] = readcol(T[n])
    if n == 'T':
        v = V[n]['val']
        V[n]['erlo'] = (10**v - 10**(v - V[n]['erlo'])) / 1e3
        V[n]['erhi'] = (10**(v + V[n]['erhi']) - 10**v)/ 1e3
        V[n]['val'] = 10**v / 1e3

len(V.nH)

# make plots
dv = np.array(T['vel'])

#########################################
# make Z plot
#########################################


#pl.rc('text', usetex=1)

font1 = 9
font2 = 12 

from barak.plot import errplot, draw_arrows, make_log_ylabels, shade_to_line
fig = plt.figure(figsize=(4.4, 7))
left = 0.18
right = 0.84
fig.subplots_adjust(left=left, hspace=1e-5,top=0.95,bottom=0.08, right=right)


ax = fig.add_axes((left, 0.6, right-left, 0.35))
Z0, Z1 = -2.8 - 0.75, -2.8 + 0.75
r0,r1=-300,400
ax.fill([r0,r0,r1,r1], [Z0,Z1,Z1,Z0], lw=0, color='0.7')
ax.text(100, -2.4, 'IGM' ,ha='center',va='center')
ax.text(220, -0.35, 'ISM (cb58)',ha='center',va='center',fontsize=14)


Z0, Z1 = -0.11, -0.52
r0,r1=-300,400
ax.fill([r0,r0,r1,r1], [Z0,Z1,Z1,Z0], lw=0, color='orange', alpha=0.3)

v = V['Z']
errplot(dv, v['val'], (v['val']-v['erlo'],
                       v['val']+v['erhi']), ax=ax, lw=1, ms=10)
ax.set_ylabel('$Z/Z_\odot$')
ax.set_xlabel('Velocity (km/s)')
ax.xaxis.get_label().set_fontsize(font2+1)
for t in ax.get_xticklabels():
    t.set_fontsize(font1+2)
make_log_ylabels(ax)
#ax.set_xticklabels('')
ax.set_xlim(-140,360)
ax.set_ylim(-2.7,0.19)
for t in ax.get_yticklabels():
    t.set_fontsize(font1+2)
ax.yaxis.get_label().set_fontsize(font2+2)

#########################################
# make the remaining smaller plots
#########################################

ylim = {'aUV':(-1.7, 1.2),
        'nH':(-3.9, -1.8),
        'T':(9.5, 32.9),
        'D':(-2.3, 0.8),
        'NH':(16.8, 19.4),
        }

ylabels = dict(aUV=r'$\alpha_\mathrm{UV}$',
               T='$\mathrm{T/(10^3 K)}$',
               NH='$N_\mathrm{H}\ \mathrm{(cm^{-2})}}$',
               nH='$\mathrm{n_H\ (cm^{-3})}$',
               D='$\mathrm{D\ (kpc)}$',
               Pk='$\mathrm{P/k\ (cm^{-3}K)}$',
               )

names =  ('nH', 'aUV','T','NH','D','Pk')
for i in range(6):
    ax = fig.add_subplot(12,1,i+7)
    k = names[i]
    v = V[k]
    c0 = v['limit'] == 'none'
    v1 = v[c0]
    errplot(dv[c0], v1['val'], (v1['val']-v1['erlo'],
                                v1['val']+v1['erhi']), c='k',
            ax=ax, lw=0.5, ms=6)
    c0 =  v['limit'] == 'upper'
    if c0.any():
        v1 = v[c0]
        draw_arrows(dv[c0], v1['val'], ax=ax, direction='down', c='k',lw=0.5)
        errplot(dv[c0], v1['val'], 0, ax=ax, lw=0.5, ms=6, c='k')
    c0 =  v['limit'] == 'lower'
    if c0.any():
        v1 = v[c0]
        if k in ('D', 'NH'):
            # tweak so lower limit arrow still shows
            draw_arrows(dv[c0], v1['val']+0.55, ax=ax,direction='up',c='k',
                        ms=0.6,capsize=4,lw=0.5)
        else:
            draw_arrows(dv[c0], v1['val'], ax=ax,direction='up',c='k',
                        ms=0.6,capsize=4,lw=0.5)
            errplot(dv[c0], v1['val'], 0, ax=ax, lw=0.5, ms=6, c='k')
    c0 = v['strong_prior']
    v1 = v[c0]
    errplot(dv[c0], v1['val'], (v1['val']-v1['erlo'],v1['val']+v1['erhi']),
            ax=ax, lw=0.5, ms=6,c='r')


    ax.set_ylabel(names[i])
    if i < 5:
        ax.set_xticklabels('')
    else:
        ax.set_xlabel('Velocity (km/s)')
        ax.xaxis.get_label().set_fontsize(font2-1)
    

    if k in ('T', 'D'):
        ax.yaxis.tick_right()
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_label_position("right")
    #if k == 'T':
    #    ax.set_yticks([10, 12.5, 15, 17.5, 20])
    if k == 'aUV':
        ax.set_yticks([-1,0,1])
    if k in ylim:
        ax.set_ylim(*ylim[k])
    if k in ylabels:
        ax.set_ylabel(ylabels[k])
    if k in ('nH', 'NH', 'D', 'Pk'):
        make_log_ylabels(ax)
    #if k == 'nH':
        #ax1 = plt.twinx(ax)
        #y0,y1 = ax.get_ylim()
        #ax1.set_ylim(Uconst-y0, Uconst-y1)
        #ax1.set_ylabel('$U$')
        #make_log_ylabels(ax1)
        #for t in ax1.get_yticklabels():
        #    t.set_fontsize(font1)
        #ax1.yaxis.get_label().set_fontsize(font2)
    for t in ax.get_yticklabels():
        t.set_fontsize(font1)
    for t in ax.get_xticklabels():
        t.set_fontsize(font1)
    ax.yaxis.get_label().set_fontsize(font2)
    ax.set_xlim(-140,360)
    #ax.yaxis.grid()

plt.savefig('parplot.png',dpi=200)
plt.savefig('parplot.pdf')
