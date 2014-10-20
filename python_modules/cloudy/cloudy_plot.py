from matplotlib import colors, gridspec
from astro.plot import puttext, A4LANDSCAPE, axvlines
from astro.constants import Ryd_Ang, eV, Ryd, k
from astro.utilities import indexnear

#pl.cm.register_cmap(name='BRG', cmap=make_cmap(*cvals.brg))

# number of axes to plot
nplot=24

M = loadobj('qg_grid.sav.gz')
#if M.nH[-1] < 1e-15:
#    M.nH[-1] = 0

nNHI = len(M.NHI)
nnH = len(M.nH)
nZ = len(M.Z)

roman_map = dict(I=0, II=1, III=2, IV=3, V=4, VI=5, VII=6, VIII=7, IX=8, X=9)
roman = set('IVX')

xlabels = dict(Z='Z = log (X/H) - log (X/H)$_\odot$',
               nH='log n$_H$ (cm$^{-3}$)',
               NHI='log $N_{HI}$ (cm$^{-2}$)')

labels = dict(NHI='log N$_{HI}$=%.3g',
              nH='log n$_H$=%.3g', Z='log Z=%.3g') 

def split_trans(trans):
    i = 1
    while trans[i] not in roman:
        i+=1
    return trans[:i], trans[i:]

# Quantities that vary with each model are saved as a 3 or 4
# dimensional array. The first three array axes give values as a
# function of NHI, nH, Z.
#
# The final axis in each N arrays is the transition number. For
# example, to access SiII for all models, use:
#
#    models.N.Si[:,:,:,1] or models.N.Si[...,1]
#
# to access MgI for models run with NHI[1] and Z[1] values, but all nH
# values:
#
#    models.N.Mg[1,:,1,0]

def plot_mod(x, z, yvals, ylabel, ax, ind=0, cmap=pl.cm.rainbow,
             printlabel=True):
    """ Plot column-density-derived values yvals as a function of the
    x values (NHI, nH or Z), showing variation of quantity z by
    different coloured curves. ind is the index of the value used,
    which isn't varied.
    """ 
    # Want index order to be indtype, x, z. By default it's NHI, nH,
    # Z. Otherwise it has to change...
    if (x,z) == ('NHI','Z'):
        yvals = np.swapaxes(yvals, 0, 1) 
    elif (x,z) == ('Z','NHI'):
        yvals = np.swapaxes(yvals, 0, 1)
        yvals = np.swapaxes(yvals, 1, 2) 
    elif (x,z) == ('nH','NHI'):
        yvals = np.swapaxes(yvals, 0, 2) 
    elif (x,z) == ('NHI', 'nH'):
        yvals = np.swapaxes(yvals, 0, 2)
        yvals = np.swapaxes(yvals, 1, 2)
    elif (x,z) == ('Z','nH'):
        yvals = np.swapaxes(yvals, 1, 2) 
    
    norm = colors.normalize(M[z].min(), M[z].max())
    label_indices = set((0, len(M[z])//2, len(M[z])-1))
    for i in range(len(M[z])):
        # spring, summer, autumn, winter are all good
        c = cmap(norm(M[z][i]))
        label = None
        if i in label_indices:
            label = labels[z] % M[z][i]
        #ax.plot(M[x], yvals[ind,:,i], '-', lw=2.5, color='k') 
        ax.plot(M[x], yvals[ind,:,i], '-', lw=1.5, color=c, label=label)

    val, = list(set(['nH','NHI','Z']).difference([x,z]))
    if printlabel:
        ax.set_title(labels[val] % M[val][ind], fontsize='medium')
        ax.title.set_y(1.01)

    ax.set_xlabel(xlabels[x], fontsize='small')
    ax.set_ylabel(ylabel)
    ax.minorticks_on()
    ax.set_xlim(M[x][0]+1e-3, M[x][-1]-1e-3)


def cleanup(ax, ratio, nplot, gs):
    """ After all the plots are made, clean up the spacing between
    plots and labels, make legend and remove unnecessary x ticklabels
    and labels.
    """
    gs[0].tight_layout(fig,rect=[0, 0, 0.25, 1],pad=0.5)
    gs[0].update(hspace=1e-5, wspace=1e-5)
    if len(ax) > nplot/2:
        gs[1].tight_layout(fig,rect=[0.25, 0, 0.5, 1],pad=0.5)
        gs[1].update(hspace=1e-5, wspace=1e-5, bottom=gs[0].bottom,
                     top=gs[0].top)
    if len(ax) > nplot:
        gs[2].tight_layout(fig,rect=[0.5, 0, 0.75, 1],pad=0.5)
        gs[2].update(hspace=1e-5, wspace=1e-5, bottom=gs[0].bottom,
                     top=gs[0].top)
    if len(ax) > 3*nplot/2:
        gs[3].tight_layout(fig,rect=[0.75, 0, 1, 1],pad=0.5)
        gs[3].update(hspace=1e-5, wspace=1e-5, bottom=gs[0].bottom,
                     top=gs[0].top)

    for i in range(len(ax)):
        if i not in (10, 11, 22, 23, 34, 35, 46, 47) and i < (len(ax)-2):
            ax[i].set_xticklabels([])
            ax[i].set_xlabel('')

    # axes = [ax[-2]]
    # if len(ax) > nplot:
    #     axes.append(ax[nplot - 2])
    # for a in axes:
    #     a.set_xticklabels(['',u'\u22123.5','',u'\u22122.5','',u'\u22121.5','',
    #                        u'\u22120.5', ''])
        
    ax[1].legend(frameon=0, loc='best')

def make_gridfig(nplot):
    """ Make a bunch of plots in a A4 landscape figure in four
    columns.
    """
    fig = pl.figure(figsize=A4LANDSCAPE)
    # divide into 4 columns
    gs = [gridspec.GridSpec(nplot/4, 2) for i in range(4)]
    gs[0].update(left=0, right=0.25)
    gs[1].update(left=0.25, right=0.5)
    gs[2].update(left=0.5, right=0.75)
    gs[3].update(left=0.75, right=1)
    return fig, gs

def match_ylim_tweak_xlim(axes):
    """ Make all axes have the same y limits, and slightly tweak x
    limits.
    """
    ymin, ymax = 99, -99
    for a in axes:
        y0, y1 = a.get_ylim()
        if y0 < ymin:
            ymin = y0
        if y1 > ymax:
            ymax = y1
    for a in axes:
        a.set_ylim(ymin+1e-3, ymax-1e-3)
        x0,x1 = a.get_xlim()
        a.set_xlim(x0+1e-3, x1-1e-3)

    
if 1:
    ###################################################################
    # Make lots of plots of column density ratios as a function of nH
    # and NHI.
    ###################################################################

    #ratios = ("""SiIV/SiII SiIII/SiII SiIV/SiIII CIII/CII CIV/CII CIV/CIII
    #AlIII/AlII NV/NII OVI/OI OI/SiII CIV/SiIV  MgII/FeII FeII/SiII
    #OVI/SiII OVI/NV HI/HII SiII/HI CII/HI AlII/HI""").split()

    ratios = ("""SiIV/SiII SiIII/SiII SiIV/SiIII CIV/CII
    AlIII/AlII HI/HII SiII/HI CII/HI AlII/HI OI/HI""").split()

    # ratios not useful for constraining U - they don't change
    # monotonically with U, or change significantly with metallicity.
    # AlII/SiII AlII/NII CII/SiII MgII/SiII MgII/NII NII/SiII SiII/HI
     
    fig, gs = make_gridfig(nplot) 
    iax = 0
    ax = []
    dax = adict()
    for ratio in ratios:
        if iax == 2*nplot:
            cleanup(ax, ratio, nplot, gs)
            fig, gs = make_gridfig(nplot)
            iax = 0
            ax = []
        print ratio
        trans = ratio.split('/')
        atoms,nums = zip(*[split_trans(t) for t in trans])
     
        i0, i1 = (roman_map[n] for n in nums)
        yvals = M.N[atoms[0]][..., i0] - M.N[atoms[1]][..., i1]
        ylabel = r'log (N$_{%s}$ / N$_{%s}$)' % tuple(
            atom + n for atom,n in zip(atoms, nums))

        # if trans[1] == 'HI':
        #     # normalise metals
        #     yvals = yvals - M.Z
        #     ylabel = r'log (N$_{%s}$ / N$_{%s}$ / Z)' % tuple(
        #         atom + n for atom,n in zip(atoms, nums))
        yvals = yvals.clip(-10)
        ylabel = ''

        i = 3 
        if iax < nplot / 2:  i = 0
        elif iax < nplot:  i = 1
        elif iax < 3*nplot / 2:  i = 2
            
        ax.extend([pl.subplot(gs[i][iax % (nplot/2)]),
                   pl.subplot(gs[i][(iax+1) % (nplot/2)])])
        dax[ratio] = ax[-2]
        p = (True if iax in (0, nplot/2, nplot, 3.*nplot/2) else False)
        plot_mod('nH','Z', yvals, ylabel, ax[iax],
                 ind=indexnear(M.NHI, 16.5), printlabel=p)
        plot_mod('NHI','Z', yvals, ylabel, ax[iax+1],
                 ind=indexnear(M.nH, -2.5), printlabel=p)

        ax[iax+1].set_ylabel('')
        ax[iax+1].set_yticklabels([])
        puttext(0.9, 0.1, ratio, ax[iax+1], fontsize='large', ha='right')

        match_ylim_tweak_xlim(ax[iax:iax+2])
        iax += 2
    else:
        cleanup(ax, ratio, nplot, gs)

    # plot 

if 1:
    #################################
    # plot the observed ratios
    ####################################
    
    # each entry is the value, 1sig low, 1sig high

    obs = parse_config('observed')
    for k in obs:
        vals =  map(float, obs[k].split())
        if len(vals) == 2:
            obs[k] = vals[0], vals[1], vals[1]
        elif len(vals) == 3:
            obs[k] = tuple(vals)
        else:
            raise ValueError('Error parsing entry %s' % obs[k])
        
    def get_ratio(a, b):
        """ Measure minimum and maximum ratio for a/b"""
        return (a[0]-a[1]) - (b[0]+b[2]), a[0]-b[0], (a[0]+a[2]) - (b[0]-b[1])

    # obs_ratios = """ SiIV/SiII SiIII/SiII SiIV/SiIII CIV/CII CIV/SiIV
    # AlIII/AlII CII/SiII OI/SiII MgII/FeII FeII/SiII CII/HI SiII/HI AlII/HI
    # NV/CIV """.split()

    obs_ratios = """ SiIV/SiII SiIII/SiII SiIV/SiIII CIV/CII
    AlIII/AlII CII/HI SiII/HI AlII/HI OI/HI""".split()

    for ratio in obs_ratios:
        numer, denom = ratio.split('/')
        low, best, high = get_ratio(obs[numer], obs[denom])
        dax[ratio].fill_between([-5, 0.5], low, high, alpha=0.2, zorder=10)
        dax[ratio].plot([-5, 0.5], [best, best], 'k', zorder=11)

    for a in ax:
        a.axvline(-2.5, color='k', ls='--')

if 0:
    ###################################################################
    # Make a plot showing the column density as a function of nH
    ###################################################################
    pl.figure()
    IONS = ('MgII CII CIII CIV SiII SiIII SiIV FeII '
            'OI OVI HI AlI AlII AlIII NII NI').split()
    #trans = ('SiII SiIII SiIV CII CIII CIV AlII AlIII '
    #         'MgI MgII FeI CaII CaI FeII OI OVI NII NV NeVIII MgX HII').split()
    atoms,nums = zip(*[split_trans(t) for t in IONS])
    ax = pl.gca()
    colors = dict(Si='y', C='k', Al='c', O='r', N='g', Fe='orange',
                  Ne='pink',  Mg='b', H='m', Ca='purple')
    count = dict((k,0) for k in colors)
    ls = ['-', '--', '-.', ':']
    iNHI = nNHI // 2
    iZ = nZ // 2
    for atom, num in zip(atoms, nums):
        col = colors[atom]
        N = M.N[atom][iNHI, :, iZ, roman_map[num]]
        ax.plot(M.nH, N, lw=3, ls=ls[count[atom] % 3],color=col,label=atom+num)
        count[atom] += 1

    ax.set_xlabel(xlabels['nH'])
    ax.set_ylabel('log N (cm$^{-2}$)')
    ax.set_title('log Z = %.2f, log NHI = %.2f' % (M.Z[iZ], M.NHI[iNHI]),
                 fontsize='medium')
    ax.set_ylim(8, 23)
    ax.set_xlim(M.nH[0], M.nH[-1])
    
    pl.legend(frameon=0, ncol=2)


if 0:
    ####################################################################
    # plot the ionisation energies of ions over the incident
    # continuum.
    ####################################################################
    
    #from astro.pyvpfit import readatom
    #atom = readatom()
    #IONS = list(atom)

    ions = readtxt(astro.datapath+'linelists/Verner_table4.txt.gz',readnames=1)
    IONS = ('MgI MgII MgX CI CII CIII CIV SiII SiIII SiIV FeI FeII '
            'OI OII OIII OIV OV OVI OVII CaI CaII HI AlI AlII AlIII '
            'NII NI SI SII SIII NaI NeVIII').split()
    # IONS = ('MgII CII CIII CIV SiII SiIII SiIV FeII '
    #         'OI OVI HI AlI AlII AlIII NII NI').split()

    ions1 = ions[np.in1d(ions.name, IONS)]
    ions1.sort(order='ie')
    energy = ions1.ie
    fig = pl.figure(figsize=A4LANDSCAPE)
    fig.subplots_adjust(left=0.11, right=0.95)

    ax = pl.gca()
    ax.loglog(M.cont.ryd * Ryd / eV, M.cont.fnu, label='UVB z=2.2')
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel(r'$F_{\nu}$ (ergs/s/cm$^2$/Hz)')
    ax.set_xlim(3, 1e3)
    axvlines(energy, 0, 1, ax=ax)
    for i in range(len(ions1)):
        puttext(energy[i], 0.8 + 0.07*(i % 3), ions1.name[i], ax,
                xcoord='data', rotation=90, fontsize='small', ha='right')
    
    ax1 = ax.twiny()
    ax1.set_xscale('log')
    E0, E1 = ax.get_xlim()
    T0, T1 = 2*E0*eV/k , 2*E1*eV/k
    ax1.set_xlim(T0,T1)
    ax.set_ylim(1e-23, 3e-17)
    ax1.set_xlabel("required gas temperature = 2*Energy/k (K)")


#rc('legend', borderaxespad=1.5, fontsize='small')
#rc('font', size=16)

rc('legend', borderaxespad=0.5, fontsize=8)
rc('font', size=12)

pl.show()

# Conclusions: Ratios of transitions in the same species (Si, C) are
# very sensitive to the ionization parameter (thus density), but
# mostly insensitive to the HI column density over the range logN 14
# -> 18 and metallicity over the range logZ -2 -> 0.
