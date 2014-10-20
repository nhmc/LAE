from model import model, tr_plot
from barak.absorb import get_ionization_energy
from barak.plot import make_log_xlabels

usepot = 1

plt.rc('font', size=14)
plt.rc('xtick.major', pad=4)
plt.rc('ytick.major', pad=6)

if usepot:
    ipot = [get_ionization_energy(t) for t in tr_plot]
    # convert to Rydbergs, which makes a nicer plot.
    X = np.log10(ipot) - np.log10(13.6)
else:
    X = np.arange(len(tr_plot))
    
#plot_model(coords[1])

NHI_ = 14.88
Z_ = -0.6
nH_ = -2.92
aUV_ = 0
ncoord = 40
aUVvals= np.linspace(-2, 1, ncoord)
Zvals = np.linspace(-2,0, ncoord)
nHvals = np.linspace(-4,-1.5, ncoord)
coords1 = [(NHI_, nH_, Z_, aUV) for aUV in aUVvals]
coords2 = [(NHI_, nH_, Z, aUV_) for Z in Zvals]
coords3 = [(NHI_, nH, Z_, aUV_) for nH in nHvals]

color = 'b'

iontext = []
text = []
lines = []
figs = []
for k in range(3):
    fig = plt.figure(figsize=(5.4, 3.2))
    figs.append(fig)
    ax = fig.add_subplot(111)
    ax.set_ylabel('$\log_{10}\ (\mathrm{N/cm}^{-2})$')
    ax.set_xlabel('$\mathrm{Ionization\ potential\ (Ryd)}$')
    ax.set_xlim(np.log10(0.9), np.log10(12))
    make_log_xlabels(ax, yoff=-0.07)
    ax.set_xticks([np.log10(1), np.log10(10)])
    ax.set_xticklabels(['1', '10'])
    fig.subplots_adjust(bottom=0.17, left=0.11, right=0.97, top=0.96)
    ax.set_ylim(4.5, 14.9)
    plt.autoscale(0)
    itexts = []
    m0 = model(coords1[k], for_plot=1)
    for i in range(len(m0)):
        dy = 0
        dx = 0
        if usepot:
            if tr_plot[i] == 'SiII':
                dy = 0.1
            elif tr_plot[i] == 'MgII':
                dx = -0.02
                dy = -0.1
            elif tr_plot[i] == 'SiIV':
                dx = -0.03
            elif tr_plot[i] == 'CIV':
                dy = -1.1
                dx = -0.02
            elif tr_plot[i] == 'OI':
                dx = -0.02
            elif tr_plot[i] == 'FeII':
                dx = +0.045
                dy = -0.4
            elif tr_plot[i] == 'FeIII':
                dx = +0.03
                dy = -1
            elif tr_plot[i] == 'AlIII':
                dx = -0.03
            elif tr_plot[i] == 'FeIII':
                dy = 0.5
        t = ax.text(X[i]+dx, m0[i] + 0.4 + dy, tr_plot[i], fontsize=10,
                    ha='center')
        itexts.append(t)

    t = puttext(0.95, 0.1, '', ax, ha='right', color=color, fontsize=20)
    text.append(t)
    iontext.append(itexts)
    l, = ax.plot([], [], '.-', color=color, ms=8)
    lines.append(l)


for i in range(len(coords1)):
    #m = model(coords1[i], for_plot=1)
    m = model(coords2[i], for_plot=1)
    #m = model(coords3[i], for_plot=1)
    lines[2].set_data(X, m)
    #text[2].set_text(r'$\alpha_\mathrm{UV}=%4.2f$' % coords1[i][-1])
    text[2].set_text(r'$\mathrm{[X/H]}=%.2g$' % coords2[i][-2])
    #text[2].set_text(r'$n_\mathrm{H}=%.2g$' % coords3[i][1])
    for j in range(len(m)):
        dy = 0
        dx = 0
        if usepot:
            if tr_plot[j] == 'SiII':
                dy = 0.1
            elif tr_plot[j] == 'MgII':
                dx = -0.02
                dy = -0.1
            elif tr_plot[j] == 'SiIV':
                dx = -0.03
            elif tr_plot[j] == 'CIV':
                dy = -1.1
                dx = -0.02
            elif tr_plot[j] == 'OI':
                dx = -0.02
            elif tr_plot[j] == 'FeII':
                dx = +0.045
                dy = -0.4
            elif tr_plot[j] == 'FeIII':
                dx = +0.03
                dy = -1
            elif tr_plot[j] == 'AlIII':
                dx = -0.03
            elif tr_plot[j] == 'FeIII':
                dy = 0.5
        iontext[2][j].set_y(m[j] + 0.4 + dy)
    fig.savefig('%03i.png' % i, dpi=150)


#plt.tight_layout()
#plt.show()
#fig.savefig('aUV.pdf')
