# Make a plot showing the total NHI

import barak.spec
from barak.absorb import calc_DLA_trans

sp = barak.spec.read('q0002m422_nhmc_cont.txt')

plt.rc('font', size=13)

fig = plt.figure(figsize=(4.4, 2.5))
fig.subplots_adjust(left=0.13, right=0.97, top=0.95, bottom=0.21)

ax = pl.gca()

sp.fl /= 1000
sp.er /= 1000
sp.co /= 1000

ax.axhline(0, color='0.7')


logN = 16.838, 16.938, 17.038
model = []
for N in logN:
    m,ticks = calc_DLA_trans(sp.wa, 2.46408, 7, logN=N, bHI=27)
    model.append(m)

plt.fill_between(sp.wa, sp.co * model[0], y2=sp.co * model[2], lw=0,
                 alpha=0.2, color='r', zorder=9)

c0 = between(sp.wa, 3161.52, 3163.82)

ax.plot(sp.wa[c0], sp.fl[c0], color='0.7', lw=0.5, drawstyle='steps-mid', zorder=10)
sp.fl[c0] = np.nan
ax.plot(sp.wa, sp.fl,color='0.4', lw=0.5, drawstyle='steps-mid', zorder=10)

ax.plot(sp.wa, sp.co, '--k', zorder=10, lw=0.5)

plt.plot(sp.wa, sp.co * model[1], 'r', lw=1, zorder=10)

ax.set_xlim(3158.1,3170.9)
ax.set_ylim(-.2,3.47)
ax.set_xlabel(r'Wavelength ($\AA$)')
ax.set_ylabel(r'Flux (Arbitrary)')

#ax.set_yticklabels('')

plt.savefig('ll.pdf')
plt.savefig('ll.png')

plt.show()
