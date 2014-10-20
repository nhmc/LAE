from barak.constants import c_kms
from barak.io import parse_config
from barak.pyvpfit import readf26
from barak.convolve import convolve_constant_dv
from barak.absorb import readatom, find_tau, findtrans, split_trans_name
from barak.sed import vel_from_wa, make_constant_dv_wa_scale
from barak.utilities import indexnear, between
from barak.plot import puttext, get_fig_axes
import barak.spec
from cStringIO import StringIO
from subprocess import call
 
import numpy as np
import pylab as plt

import velplot.velplot
reload(velplot.velplot)
from velplot.velplot import plot_tick_vel, read_transitions

ATOMDAT = readatom()

lw_model = 0.7
lw_data = 1

vzero = -219

#plt.ioff()
plt.rc('font', size=11)
plt.rc('xtick', labelsize=8.5)
plt.rc('ytick', labelsize=8.5)
plt.rc('xtick.major', size=3)      # major tick size in points
plt.rc('xtick.minor', size=1.5)      # minor tick size in points
plt.rc('xtick.major', pad=3.5)     # distance to major tick label in points
plt.rc('ytick.major', size=3)      # major tick size in points
plt.rc('ytick.minor', size=1.5)      # minor tick size in points
plt.rc('ytick.major', pad=4.5)     # distance to major tick label in points

PREFIX = '/home/nhmc/projects/MPIA_QSO_LBG/data/'

config = StringIO("""\
trfilename = figmetals2_trans
#trfilename = transitions_h1
spfilename = /Users/ncrighton/Projects/MPIA_QSO_LBG/data/vpfit/J0004/q0002m422_nhmc_cont.txt
f26name = /Users/ncrighton/Projects/MPIA_QSO_LBG/data/vpfit/J0004/all.f26
Rfwhm = 6.6
vmax =  142
vmin = -458
wadiv = None
# this is the Lya redshift -219 km/s (mean Lya offset from Ha in Hashimoto)
#v0redshift = 2.46333
#redshift = 2.46333
# original Lya redshift
v0redshift = 2.46586
#redshift = 2.46586
modelfill = False
osc = False
residuals = False
""")

# VCOMP = [(-354, -316),
#          (-316, -283.5),
#          (-286, -210),
#          (-210,  -91),
#          (0,     47),
#          (92,    124),
#          ]

VCOMP = [
    (-354, -316.5),
    (-316.5, -283.5),
    (-283.5, -210),
    (-180,  -115),
    (0,     27),
    (92,    124),
         ]


unrelated = {'HI 1025': [(-500, -357), (60, 200)],
             'HI 921': [(-500, -385), (20, 200)],
             'HI 914': [(150+vzero, 200+vzero), (260+vzero, 400+vzero)],
             'MgII 2796': [(-280,-200),(60,200)],
             'MgII 2803': [(-300+vzero,-130+vzero),(250+vzero,305+vzero)],
             'FeIII 1122': [(-300+vzero,-100+vzero),(155+vzero,400+vzero)],
             'SiIII 1206': [(-300+vzero,-95+vzero),(-70+vzero,10+vzero),
                            (115+vzero,195+vzero), (200+vzero,220+vzero),
                            (240+vzero,400+vzero)],
             'CIII 977': [(-190+vzero,-155+vzero), (0+vzero,400+vzero),],
             'CII 1334': [(0,300)],
             'SiII 1526': [(60,300)],
             'NV 1238': [(90+vzero,200+vzero),(280+vzero, 403+vzero)],
             'OVI 1031': [(-500,-360),(-100, -20),(45,200)],
             'OVI 1037': [(-500,-290), (-100, +45), (70,200)],
             }

cfg = parse_config(config)
transitions = read_transitions(cfg.trfilename, ATOMDAT)

if 1:
    sp = barak.spec.read(cfg.spfilename)
    ndiv = 4.
    wa_dv = make_constant_dv_wa_scale(sp.wa[0], sp.wa[-1], cfg.Rfwhm / ndiv)

    expand_cont_adjustment = 5

    vp = readf26(cfg.f26name)
    lines = vp.lines[vp.lines.name != '<>']
    tau, ticks, alltau = find_tau(sp.wa, lines, ATOMDAT, per_trans=1)
    model = convolve_constant_dv(sp.wa, np.exp(-tau), wa_dv, ndiv)
    models = [convolve_constant_dv(sp.wa, np.exp(-t), wa_dv, ndiv) for t in
              alltau]
    adjust = [l for l in vp.lines if l['name'].strip() == '<>']
    print 'Nadjust', adjust
    if len(adjust) > 0:
        regions = vp.regions
        isort = regions.wmin.argsort()
        print 'Applying continuum adjustments for', len(adjust), 'regions'
        for val in adjust:
            wav0 = ATOMDAT['<>'].wa[0] * (1 + val['z'])
            i0 = regions.wmin[isort].searchsorted(wav0)
            i1 = regions.wmax[isort].searchsorted(wav0)
            #print i0, i1
            #import pdb; pdb.set_trace()
            assert i0 - 1 == i1
            linterm = val['b']
            level = val['logN']
        
            wa0 = regions.wmin[isort[i0 - 1]] - expand_cont_adjustment
            wa1 = regions.wmax[isort[i1]] + expand_cont_adjustment
            c0 = between(sp.wa, wa0, wa1)
        
            mult = level + linterm * (sp.wa[c0]/wav0 - 1)
            sp.co[c0] = sp.co[c0] * mult

if 1:
    fig,AX = get_fig_axes(10, 2, len(transitions), aspect=0.24, width=6.3)
    #fig,AX = get_fig_axes(7, 2, len(transitions), aspect=0.26, width=6.3)
    fig.subplots_adjust(hspace=1e-5, wspace=1e-5, bottom=0.06, right=0.99,
                        left=0.08, top=0.97)
    label = []
    for tr in transitions:
        ion,wa = tr[0].split()
        atom,stage = split_trans_name(ion)
        if tr[0] == 'HI 926':
            label.append('H I')
        elif tr[0] == 'DI 937':
            label.append('D I\nLy-$\epsilon$')
        else:
            label.append(atom + ' ' + stage + '\n' + wa)

    print vp.lines

    c0 = ((vp.lines.name == 'C IV') &
          ~np.in1d(vp.lines.zpar,('AR','AS','AQ','AU','AT', 'AV')))
    zvals = vp.lines[c0].z
    is_low = np.array([np.any(abs(z - zvals) < 1e-5) for z in lines.z])
    
    tvel = c_kms * (zvals - cfg.v0redshift) / (1 + cfg.v0redshift)
    print repr(tvel - vzero)

    zvals = vp.lines[vp.lines.name == 'O VI'].z
    tvel_o6 = c_kms * (zvals - cfg.v0redshift) / (1 + cfg.v0redshift)

    #import pdb; pdb.set_trace()
    
    # tvel = [
    #     -217.842, -179.138,  -149.848,  -128.839, -116.810,
    #     -99.2886, -85.8642,  -75.2292,  -59.2156,  -47.1599,
    #     -31.9920, -20.8340,    2.0049,   17.7830,   28.5051,
    #     40.796]

    for i,tr in enumerate(transitions):
        ax = AX['axes'][i]
        for x in VCOMP:
            ax.fill_between(np.array(x) - vzero, [1.5, 1.5], y2=[-0.5,-0.5],
                            color=(0.996, 0.9455, 0.848),
                            edgecolor='w', lw=1.5)

        if i in (0, 10) :
            ax.xaxis.tick_top()
            ax.xaxis.set_ticks_position('both') # THIS IS THE ONLY CHANGE
            #ax.text(-300, 0.2, 'Metallicity\n$0.1-0.4\ Z_{\odot}$', fontsize=10)
        #    ax.text(-55, 0.3, 'Metallicity\n$>0.2\ Z_{\odot}$', fontsize=10)
            ax.text(-340 - vzero, 0.3, '1', fontsize=8, color='brown')
            ax.text(-308 - vzero, 0.3, '2', fontsize=8, color='brown')
            ax.text(-254 - vzero, 0.3, '3', fontsize=8, color='brown')
            ax.text(-170 - vzero, 0.3, '4', fontsize=8, color='brown', ha='center')
            ax.text(-150 - vzero, 0.3, '5', fontsize=8, color='brown', ha='center')
            ax.text(-128 - vzero, 0.3, '6', fontsize=8, color='brown', ha='center')
            ax.text(5   - vzero, 0.3, '7', fontsize=8, color='brown')
            ax.text(100  - vzero, 0.3, '8', fontsize=8, color='brown')

        ax.axhline(0, color='gray', lw=0.5)
        ax.axhline(1, color='gray', lw=0.5, ls='dashed', zorder=3) 
        vel = vel_from_wa(sp.wa, tr[1]['wa'], cfg.v0redshift)
        c0 = between(vel, -1200, 1200)
        v,wa,fl,er,co = vel[c0], sp.wa[c0], sp.fl[c0], sp.er[c0], sp.co[c0]
        is_unrel = np.zeros(len(v), bool)
        print tr[0]
        if tr[0] in unrelated:
            for v0,v1 in unrelated[tr[0]]:
                i0,i1 = v.searchsorted([v0, v1])
                print 'unrelated!', v0,v1
                is_unrel[i0:i1] = True
                i0 -= 1
                i1 += 2
                ax.plot(v[i0:i1] - vzero, fl[i0:i1]/co[i0:i1], '-', lw=lw_data,
                        color='0.8', drawstyle='steps-mid')

        fl1 = np.array(fl)
        fl1[is_unrel] = np.nan
        #print stats(fl1)
        ax.plot(v - vzero, fl1/co, color='0.3', lw=lw_data,
                ls='steps-mid')


        ax.plot(v - vzero, er/co, color='g', lw=0.5)
        ax.plot(v - vzero, model[c0], 'r', lw=lw_model)
        for j,m in enumerate(models):
            if is_low[j] and lines['name'][j].strip() != 'O VI':
                ax.plot(v - vzero, m[c0], 'r', lw=0.3)
            else:
                l, = ax.plot(v - vzero, m[c0], 'r', lw=0.3)
                l.set_dashes([3,2,3,2])

        if not tr[0].startswith('OVI') and not tr[0].startswith('NV'):
            for j,v in enumerate(tvel):
                #if j in (0,11,13,14,15):
                #    T = plot_tick_vel(ax, v, 0, '', lw=1, col='b')
                if cfg.vmin < v < cfg.vmax:
                    T = plot_tick_vel(ax, v - vzero, 0.04,'',
                                      height=0.14, lw=1, col='b')
        else:
            for j,v in enumerate(tvel_o6):
                #if j in (0,11,13,14,15):
                #    T = plot_tick_vel(ax, v, 0, '', lw=1, col='b')
                if cfg.vmin < v < cfg.vmax:
                    T = plot_tick_vel(ax, v - vzero, 0.04,'',
                                      height=0.14, lw=1, col='r')

        puttext(0.03, 0.12, label[i], ax, ha='left', va='bottom')
        
        #ax.set_xticks([-50, -25, 0, 25, 50])
        ax.set_yticks([0, 0.5, 1])
        ax.minorticks_on()
        ax.set_xlim(cfg.vmin - vzero, cfg.vmax - vzero)
        ax.set_ylim(-0.1, 1.39)
        if i > 9:
            ax.set_yticklabels('')
        if i not in (0,9,10, 19):
            ax.set_xticklabels('')
        if i == 4:
            ax.set_ylabel('Transmission')
        
    fig.text(0.52, 0.01, 'Velocity offset (km s$^{-1}$)', ha='center')
    #fig.text(0.04, 0.4, 'Transmission', va='center', ha='center', rotation=90)
    #fig.suptitle('%s %s $R_\perp=%.1f$kpc $z_{gal}=%.3f$' % (
    #    p.qso, p.gname, p.rho, p.zgal))
    fig.canvas.print_figure('fig_metals.png', dpi=300)
    fig.canvas.print_figure('fig_metals.pdf')
    #fig.canvas.print_figure('talk_lores_' + name, dpi=200)
    #fig.canvas.print_figure('pics/' + 'talk_' + name, dpi=200)
    #fig.clf()
    #plt.show()
    if 0:
        call('pdf2ps fig_metals.pdf', shell=True)
        call('ps2pdf fig_metals.ps', shell=True)
        call('ps2eps -f fig_metals.ps', shell=True)
