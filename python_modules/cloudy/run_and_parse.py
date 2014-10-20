from .utils import calc_uvb, calc_local_jnu, read_starburst99, get_data_path,\
     tilt_spec, find_gamma

from barak.utilities import adict, between
from barak.constants import Ryd, Ryd_Ang, pi, hplanck
from barak.io import saveobj, parse_config, writetable

import numpy as np

from collections import OrderedDict
from subprocess import call
from glob import glob
import os, sys, time, multiprocessing, gzip
from cStringIO import StringIO

cfg_temp = """\
abundances=None
table=None
cuba_name=None
run_cloudy=True
uvb_tilt=None
distance_starburst_kpc=None
overwrite=False
logfnu912=None
fesc=None
grains=False
constantpressure=False
"""

cfg_defaults = parse_config(StringIO(cfg_temp))

def write_example_grid_config():
    fh = open('grid.cfg', 'w')
    fh.write("""\
# Path to the CUBA file giving the UV background. 
cuba_name = /home/nhmc/code/repo/QSOClustering-svn/misc/CUBA/Q1G01/bkgthick.out
# Prefix for final grid file
prefix = qg
# Redshift for the UV background
z = 2.2
# Minimum, maximum and step for neutral hydrogen column density (cm^-2)
logNHI =   14.0  19.0 1.0
# Minimum, maximum and step for metallicity
logZ   =   -2.0  0.0  1.0
# Minimum, maximum and step for hydrogen density (cm^-3)
lognH  =   -5.0  0.0  0.5
uvb_tilt = 0.0
# number of processors to use
nproc = 4
# overwrite any existing files?
overwrite = False

# uncomment below line to use a Cloudy table command to specify the
# incident continuum strength and shape
# table = ism

# distance from a typical starburst galaxy (will modify background
# spectrum).

# distance_starburst_kpc = 10

# uncomment line below to use abundances different to solar + no dust
# grains.

# abundances = ism

# uncomment below to stop running cloudy (just generate input or parse results)
# run_cloudy = False
""")
    fh.close()


def write_uvb(outname, energy, logjnu, overwrite):
    """ Create a Cloudy input file using the energy and logjnu
    arrays.

    Parameters
    ----------
    energy   energy in Rydbergs 
    Jnu      Intensity at each energy (erg/s/cm^2/Hz/ster)
    """
    if not os.path.lexists(outname) or overwrite:
        print 'Writing UVB cloudy input to', outname
        fh = open(outname, 'w')
        fh.write('interpolate (%.9f -30.0)\n' % 1e-8)
        for i in xrange(len(energy)):
            fh.write('continue (%s %10.6f)\n' % (energy[i], logjnu[i]))         
        fh.write('continue (%f -30.0)\n' % 7.4e6)
        fh.close()
    else:
        print outname, 'exists, skipping'


def write_input(outfilename, z, logNHI, metallicity, lognH,
                fluxfilename=None,
                table=None, logU=None, fnu912=None, title='model',
                abundances=None, iterate='to convergence', grains=False,
                CIE_temperature=None, contname=None, trimming='-10',
                constantpressure=False):
    """ Write an input file for Cloudy.

    By default solar abundances with no grains is used. The CMB
    continuum is added at redshift `z`, and an additional incident
    continuum is specified by `fluxfilename` or `table`. If the table
    command doesn't specify the normalisation of the continuum, then
    you need to also specify either logU or fnu912.

    Parameters
    ----------
    outfilename : str
      Name of the Cloudy input file to create.
    z : float
      Redshift
    logNHI  : float
      Log of neutral hydrogen column density (cm^-2)
    metallicity : float
      log10(multiplier) for Solar metal distribution.
    lognH : float
      Log of total hydrogen density (cm^-3)
    fluxfilename : str (optional)
      File giving the shape of the incident radiation field in Cloudy
      input format. See write_uvb(). Either this or `table` must be specified.
    table : str (optional)
      The tabulated Cloudy continuum shape to use. For example 'ism'
      or 'HM05 redshift 2.2'. Either this or fluxfilename must be
      specified.
    logU : float (optional)
      Log of the dimensionless ionisation parameter. Either this or
      fnu912 must be specified if fluxfilename is given.
    fnu912 : float (optional)
      Flux in erg/s/cm^2/Hz at 1 Rydberg (912 Angstroms).  Either this
      or logU must be specified if fluxfilename is given.
    CIE_temperature : float (optional)
      Temperature at which to fix the Cloudy run (this simulates
      collisionally ionised equilibrium, or CIE). The default is None,
      meaning the temperature is not fixed.
    title : str ('model')
      Title for the model written to the input file.
    contname : str (optional)
      If not None (default), save the indicent continuum in `contname`.
    trimming : float (optional, default -10)
      Ignore elements with log abundance smaller than this to
      increase running time.
    abundances : str (optional)
      If given, this abundance string given to the CLoudy
      abundance command.

    """
    
    # output string
    s = ['title %s' % title]
    s += ['radius 30']              # log inner radius (cm)
    s += ['print last']             # Don't print results from iterations
    s += ['CMB redshift %s' % z]
    
    if CIE_temperature is None:
        s += ['stop temperature = 10K linear']
    s += ['hden %s' % lognH]
    if abundances is not None:
        s += ['abundances %s' % abundances]
    # metallicity log10 scale factor.
    s += ['metals %s log' % metallicity]
    s += ['stop neutral column density %s' % logNHI]
    # Ignore elements with log abundance smaller than this and
    # increase running time. The default is -6/-10
    s += ['set trimming %s' % trimming]
    if CIE_temperature is not None:
        s += ['constant temperature, t=%sK' % CIE_temperature]
    
    if logU is not None:
        s += ['ionization parameter = ' % logU]
    elif fnu912 is not None:
        s += ['f(nu) = %8.3f' % fnu912]
    elif table is None:
        raise ValueError('You must specify either logU or fnu912')

    if fluxfilename is not None:
        fh = open(fluxfilename)
        for row in fh:
            r = row.strip()
            if r:
                s += [r]
        fh.close()
    elif table is not None:
        s += ['table %s' % table]
    else:
        raise ValueError('You must specify either fluxfilename or table')

    s += ['iterate ' + iterate]
    if contname is not None:
        s += ['save incident continuum "%s"' % contname]
    s += ['']

    fh = open(outfilename, 'w')
    fh.write('\n'.join(s))
    fh.close()
    
    return

def write_grid_input(cfg, fnu912=None, fluxfilename=None, table=None,
                     abundances=None, indir='input/'):
    """ Make a grid of Cloudy input files given a range of
    metallicities, neutral column densities and total column
    densities in a parameter file.

    Parameters
    ----------
    cfg is a dictionary with keys able to accessed by attribute. The
    following keys defined:
        z: redshift
        logNHI, lognH, logZ: arrays of values used to generate grid
        logT: Array of temperatures used.  If this is preset the lognH
            key is ignored.
    fnu912: Fnu at 912 Angstroms
    fluxfilename: filename defining Fnu as a function of energy
    table: Cloudy table command to specify incident continuum shape
    abundances: Cloudy abundances command, if absent solar abundances
                without grains are used.
    """
    if os.path.lexists(indir):
        if cfg.overwrite:
            print 'Removing existing directory', indir
            call('rm -rf ' + indir, shell=1)
        else:
            print indir, 'exists, skipping'
            return
    os.mkdir(indir)

    # get grid values from parameter file
    nmodels = len(cfg.logNHI) * len(cfg.lognH) * len(cfg.logZ)
    if cfg.run_cloudy:
        print 'Assuming 3 min each, %i models will take %.2f hours.' % (
            nmodels, 3./60.*nmodels/cfg.nproc)
        s = raw_input('About %.0fMB of disk space is needed. Continue? [y] '
                      % (50./1000.*nmodels))
        if s.strip().lower() in ('n', 'no'):
            sys.exit()

    contname = cfg.prefix + '_nuFnu.dat'
    for i,NHI in enumerate(cfg.logNHI):
        for j,nH in enumerate(cfg.lognH):
            for k,Z in enumerate(cfg.logZ):
                inname = '%02i_%02i_%02i.in' % (i, j, k)
                outname = inname[:-3] + '.out'
                write_input(indir + inname, cfg.z, NHI, Z, nH,
                            fluxfilename=fluxfilename,
                            table=table, abundances=abundances,
                            fnu912=fnu912, contname=contname,
                            trimming=cfg.trimming, iterate='2 times',
                            grains=cfg.grains)
                # only save incident continuum for the first model
                contname = None

def run_single_process(command):
    """ Call a cloudy command and print a message. Used by
    multiprocessing.pool.map() in run_grid().
    """
    try:
        call(command, shell=1)
        call('gzip ' + command.split()[-1], shell=1)
        sys.stdout.write('.')
        sys.stdout.flush()
    except KeyboardInterrupt:
        pass

def run_grid(indir='input/', outdir='output/', cloudy='cloudy.exe', nproc=4, overwrite=False):
    """ Run cloudy on a grid of models generated using make_grid_input().

    Output files are saved to `outdir`.
    """
    if os.path.lexists(outdir):
        if overwrite:
            print 'Removing existing directory', outdir
            call('rm -rf ' + outdir, shell=1)
        else:
            print outdir, 'exists, skipping'
            return
    os.mkdir(outdir)

    # make list of input arguments to run_single_process()
    args = []
    for inname in sorted(glob(indir + '??_??_??.in')):
        outname = outdir + inname.split('/')[-1][:-3] + '.out'
        args.append(('%s < %s > %s' % (cloudy, inname, outname)))

    pool = multiprocessing.Pool(processes=nproc)
    print 'Running Cloudy using %i processes' % nproc

    t1 = time.time()
    
    p = pool.map(run_single_process, args)

    print '\n%.2f min elapsed' % ((time.time() - t1)/ 60.)

def parse_abundances(i, rows):
    """ Parse a list of abundances in a cloudy file (gas or grain).

    Parameters
    ----------
    i : int
      starting row
    rows : list of str
      the cloudy output file as a list of rows.

    Returns
    -------
    i : int
      Row index of the line following the abundance list
    abun : record array
      Array giving the elements and log10(abundance) relative to H.
    """
    gas_abun = []
    temp = ''
    while rows[i].strip():
        temp += rows[i]
        i += 1
    temp = temp.split(':')
    elem = [temp[0].strip()]
    abun = []
    for t in temp[1:-1]:
        vals = t.split()
        abun.append(float(vals[0]))
        elem.append(vals[1])

    abun.append(float(temp[-1].strip()))

    return OrderedDict((e,a) for e,a in zip(elem, abun)), i


def parse_output(filename):
    """ parse a cloudy output file. Return a dictionary with the
    following info:

    filename:  Output filename 
    info:      Model grid info.
    redshift:  CMB Redshift.
    nH:        log10 total hydrogen density (cm^-3).
    N:         log10 column density (cm^-2) per atom and ion. Ordered
               lowest to highest ionisation level.
    Nex:       log10 column density (cm^-2) for excited states.
    U:         log10 ionisation parameter (dimensionless).
    Tgas:      log10 gas temperature (K) for HI and HII.
    NHI:       Cloudy stopped once this log10 NHI (cm^-2) was reached.
    Tstop:     Cloudy stopped once this log10 temperature (K) was reached.
    Z:         Metallicity.
    title:     Model title.
    gas_abun:  Gas phase abundances.
    dust_abun: Dust abundances.
    """
    atom_names = OrderedDict(Sodium     ='Na',
                             Chromium   ='Cr', 
                             Lithium    ='Li', 
                             Beryllium  ='Be', 
                             Nickel     ='Ni', 
                             Potassium  ='K', 
                             Argon      ='Ar', 
                             Boron      ='B', 
                             Scandium   ='Sc', 
                             Carbon     ='C', 
                             Sulphur    ='S', 
                             Chlorine   ='Cl', 
                             Phosphorus ='P', 
                             Zinc       ='Zn', 
                             Oxygen     ='O', 
                             Vanadium   ='V', 
                             Manganese  ='Mn', 
                             Silicon    ='Si', 
                             Calcium    ='Ca', 
                             Magnesium  ='Mg', 
                             Iron       ='Fe', 
                             Copper     ='Cu', 
                             Aluminium  ='Al', 
                             Helium     ='He', 
                             Neon       ='Ne', 
                             Cobalt     ='Co', 
                             Titanium   ='Ti', 
                             Nitrogen   ='N', 
                             Hydrogen   ='H', 
                             Fluorine   ='Fl')

    out = {}
    out['filename'] = filename.split('/')[-1]

    if filename.endswith('.gz'):
        fh = gzip.open(filename, 'rb')
    else:
        fh = open(filename)
    
    rows = fh.readlines()
    fh.close()

    # get input info
    i = 0
    while True:
        if i == len(rows):
            return out
        r = rows[i].split()
        if len(r) < 2:
            i += 1
            continue
        if ' '.join(r) == 'Gas Phase Chemical Composition':
            break
        elif r[1] == 'title':
            out['title'] = r[2].strip('"')
        elif r[1] == 'CMB':
            out['redshift'] = float(r[3])
        elif r[1] == 'hextra':
            out['heat'] = float(r[2])
        elif r[1] == 'hden':
            out['nH'] = float(r[2])
        elif r[1] == 'metals':
            if r[2:4] == ['and', 'grains']:
                out['Z'] = float(r[4])
            else:
                out['Z'] = float(r[2])
        elif r[1] == 'stop':
            if r[2] == 'neutral':
                out['NHIstop'] = float(r[5])
            elif r[2] == 'column':
                out['NHIstop'] = float(r[4])
            elif r[2] == 'temperature':
                out['Tstop'] = float(r[4].strip('K'))
        i += 1

    
    i += 1
    out['gas_abun'], i = parse_abundances(i, rows)
    i += 1

    if rows[i].strip() == 'Grain Chemical Composition':
        i += 1
        out['dust_abun'], i = parse_abundances(i, rows)
    else:
        out['dust_abun'] = {}

    while True:
        r = rows[i].split()
        if len(r) < 2:
            i += 1
            continue
        if r[1] == 'Log(U):':
            break
        i += 1
        if i == len(rows):
            return out
    
    out['U'] = float(r[2])

    while not (rows[i].startswith(' Hydrogen') and 
               rows[i][53:81] == 'Log10 Column density (cm^-2)'):
        i += 1
        if i == len(rows):
            return out

    N = OrderedDict()
    while rows[i].strip() != '' and not rows[i].startswith(' Exc state'):
        r = rows[i]
        if r[1].isalpha():
            j = 11
            name = atom_names[r[:j].strip()]
            N[name] = []
            while r[j:j+7].strip() not in ('', '(H2)'):
                N[name].append(float(r[j:j+7]))
                j += 7
        else:
            j = 0
            while r[j:j+7].strip():
                N[name].append(float(r[j:j+7]))
                j += 7
        i += 1
        if i == len(rows):
            return out

    out['N'] = N

    # read column densities for excited states
    assert rows[i].startswith(' Exc state')
    row = rows[i][len(' Exc state'):].strip()
    #print row, rows[i+1]
    vals = [row[j:j+14] for j in range(0, len(row), 14)]
    row = rows[i+1].strip()
    vals += [row[j:j+14] for j in range(0, len(row), 14)]
    #print vals

    out['Nex'] = dict((v[:4], float(v[4:])) for v in vals)
    #print repr(out['Nex'])
    while 'Log10 Mean Temperature (over radius)' not in rows[i]:
        i += 1
    # temperature of neutral and ionised H
    out['Tgas'] = float(rows[i][11:18]), float(rows[i][18:26])
    
    return out
    
def parse_grid(cfg, outdir='output/'):

    # read each model output
    names = sorted(glob(outdir + '??_??_??.out*'))
    assert names, outdir + '??_??_??.out*'
    print 'Reading output'
    models = []
    for n in names:
        models.append(parse_output(n))

    # Read the incident continuum
    energy, nuFnu = np.loadtxt(cfg.prefix + '_nuFnu.dat', unpack=1)
    nu = energy * Ryd / hplanck
    fnu = nuFnu / nu
    isort = energy.argsort()
    cont = np.rec.fromarrays([energy[isort], fnu[isort]], names='ryd,fnu')   
        
    # combine into large arrays

    grid = {}
    grid['cont'] = cont
    grid['NHI'] = cfg.logNHI
    grid['nH'] = cfg.lognH
    grid['Z'] = cfg.logZ
    grid['redshift'] = cfg.z

    # initialise keys that will contain one or more values for each
    # model.
    grid['N'] = OrderedDict()
    for atom in models[0]['N']:
        grid['N'][atom] = []

    grid['Nex'] = OrderedDict()
    for trans in models[0]['Nex']:
        grid['Nex'][trans] = []

    for key in ('gas_abun', 'dust_abun'):
        grid[key] = OrderedDict()
        for atom in models[0][key]:
            grid[key][atom] = []

    for k in ('U', 'Tgas', 'Tstop', 'filename'):
        grid[k] = []

    # initialisation finished, now copy the values from each model.

    for model in models:
        for key in ('U', 'Tgas', 'Tstop', 'filename'):
            try:
                grid[key].append(model[key])
            except:
                import pdb; pdb.set_trace()

        for key in ('N', 'Nex', 'gas_abun', 'dust_abun'):
            for atom in model[key]:
                try:
                    grid[key][atom].append(model[key][atom])
                except:
                    import pdb; pdb.set_trace()

    # finally convert lists to arrays and reshape to multiple dimensions.

    nNHI = len(cfg.logNHI)
    nnH = len(cfg.lognH)
    nZ = len(cfg.logZ)

    for key in ('U', 'Tstop', 'filename'):
        grid[key] = np.array(grid[key]).reshape(nNHI, nnH, nZ)

    for key in ('N', 'gas_abun', 'dust_abun'):
        for atom in grid[key]:
            grid[key][atom] = np.array(grid[key][atom]).reshape(
                nNHI, nnH, nZ, -1)

    key = 'Tgas'
    grid[key] = np.array(grid[key]).reshape(nNHI, nnH, nZ, -1)
    key = 'Nex'
    for trans in grid[key]:
        grid[key][trans] = np.array(grid[key][trans]).reshape(nNHI, nnH, nZ)

    # U values only vary with nH, so we don't need to keep a big grid
    # of them.
    U = np.array(grid['U'][0,:,0]).squeeze()
    if U.ndim == 0:
        U = float(U)
    del grid['U']
    grid['U'] = U 

    # abundance values only vary with Z, so we don't need to keep a big grid
    # of them either
    for key in ('gas_abun', 'dust_abun'):
        for atom in grid[key]:
            val = grid[key][atom][0,0,:].squeeze()
            if val.ndim == 0:
                val = float(val)
            grid[key][atom] = val

    if cfg.uvb_tilt isnot None:
        grid['uvb_tilt'] = cfg.uvb_tilt 

    grid['help'] = """\
cont:      Total incident continuum Fnu (ergs/cm^2/s/Hz), output by
           Cloudy's save incident continuum command, as a function
           of energy in Rydbergs.
filename:  Output filename. shape (M, N, P)
redshift:  CMB Redshift.
uvb_tilt:  Value of the alpha_UV tilt parameter.
nH:        log10 total hydrogen density (cm^-3), shape (N,)

N:         log10 column density (cm^-2) per atom and ion. Ordered
           lowest to highest ionisation level. Each ion is shape (M, N,
           P, Q), where Q is the number of ion species.

Nex:       log10 column density (cm^-2) for excited states.
U:         log10 ionisation parameter (dimensionless), shape (N,)
Tgas:      log10 gas temperature (K) for HI and HII.  shape (M, N, P, 2)
NHI:       Cloudy stopped once this log10 NHI (cm^-2) was reached, shape (M,)
Tstop:     Cloudy stopped once this log10 temperature (K) was reached.
           shape(M, N, P)
Z:         Metallicity, shape (P,)
gas_abun:  Gas-phase abundance, each element has shape (P,)
dust_abun: Abundances on dust, each element has shape (P,)

For Tgas, Tstop, filename and N values, the index order is ind_NHI,
ind_nH, ind_Z."""

    return grid

def savehdf5(filename, M, overwrite=cfg.overwrite):
    import h5py

    if os.path.exists(filename):
        print 'File exists, skipping.'
        return
    else:
        print 'Writing to {}'.format(filename)

    fh = h5py.File(filename, 'w')

    fh.attrs['help'] = M['help']
    fh.attrs['redshift'] = M['redshift']

    if 'uvb_tilt' in M:
        fh.attrs['uvb_tilt'] = M['uvb_tilt']

    for k in 'nH Z NHI U Tstop filename cont Tgas'.split():
        a = M[k]
        d = fh.create_dataset(k, a.shape, dtype=a.dtype, compression='gzip')
        d[:] = a

    for k in 'N Nex gas_abun dust_abun'.split():
        g = fh.create_group(k)
        for atom in M[k]:
            a = M[k][atom]
            d = g.create_dataset(atom, a.shape, dtype=a.dtype,
                                 compression='gzip')
            d[:] = a

    fh.close()


def read_config(name):
    """ read the configuration file, doing some extra processing
    beyond that done by parse_config().
    """
    cfg = parse_config(name, defaults=cfg_defaults)
    cfg.overwrite = bool(cfg.overwrite)
    cfg.nproc = int(cfg.nproc)
    cfg.z = float(cfg.z)
    for k in 'logNHI lognH logZ'.split():
        vmin, vmax, step = map(float, cfg[k].split()) 
        cfg[k] = np.arange(vmin, vmax + 0.5*step, step)

    return cfg

def main():
    if not os.path.lexists('grid.cfg'):
        print ('./grid.cfg file not found, writing an example grid.cfg to '
               'the current directory')
        write_example_grid_config()
        sys.exit()

    cfg = read_config('grid.cfg')

    print ''
    print 'Input values:'
    for k in sorted(cfg):
        print '  %s: %s' % (k, cfg[k])
    print ''

    if cfg.table is None:
        fluxname = cfg.prefix + '_temp_uvb.dat'
        if cfg.cuba_name is None:
            cfg.cuba_name = get_data_path() + 'UVB.out'

        uvb = calc_uvb(cfg.z, cfg.cuba_name, match_fg=False)

        writetable('cloudy_jnu_HM.tbl', [uvb['energy'], uvb['logjnu']],
                   overwrite=1,
                   units=['Rydbergs', 'erg/s/cm^2/Hz/ster'],
                   names=['energy', 'log10jnu'])

        if cfg.uvb_tilt is not None:
            if cfg.distance_starburst_kpc is not None:
                raise RuntimeError('Must only specify one of uvb_tilt and\
                distance_starburst_kpc!')

            # remember which bits had 1e-30
            clow = uvb['logjnu'] == -30
            # tilt the UV background between 1 and 10 Rydbergs
            logjnu = tilt_spec(cfg.uvb_tilt, uvb['energy'], uvb['logjnu'],
                               emin=1, emax=10)

            logjnu[clow] = -30

            print('Tilting UVB using parameter {}'.format(cfg.uvb_tilt))

            # now re-normalise to match the photoionization rate of the
            # default spectrum.

            gamma_default = find_gamma(uvb['energy'], 10**uvb['logjnu'])
            mult = gamma_default / find_gamma(uvb['energy'], 10**logjnu)
            print 'Scaling tilted Jnu by %.3g to match default gamma' % mult
            logjnu = logjnu + np.log10(mult)

            writetable('cloudy_jnu_tilted.tbl', [uvb['energy'], logjnu],
                       overwrite=1,
                       units=['Rydbergs', 'erg/s/cm^2/Hz/ster'],
                       names=['energy', 'log10jnu'])
            uvb['logjnu'] = logjnu

        elif cfg.distance_starburst_kpc is not None:
            wa, F = read_starburst99(get_data_path() + 'starburst.spectrum1')
            nu, logjnu = calc_local_jnu(wa, F, cfg.distance_starburst_kpc,
                                        cfg.fesc)
            energy = nu * hplanck / Ryd
            # use HM uvb energy limits
            cond = between(uvb['energy'], energy[0], energy[-1])
            logjnu1 = np.interp(uvb['energy'][cond], energy, logjnu)
            uvb['logjnu'][cond] = np.log10(10**uvb['logjnu'][cond] +
                                           10**logjnu1)
            writetable('cloudy_jnu_total.tbl', [uvb['energy'], uvb['logjnu']],
                       overwrite=1,
                       units=['Rydbergs', 'erg/s/cm^2/Hz/ster'],
                       names=['energy', 'log10jnu'])

        write_uvb(fluxname, uvb['energy'], uvb['logjnu'], cfg.overwrite)

        # Fnu at 1 Rydberg
        k = np.argmin(np.abs(uvb['energy'] - 1.))
        logfnu912 = np.log10(10**uvb['logjnu'][k] * 4 * pi)
    else:
        logfnu912 = cfg.logfnu912
        fluxname = None

    write_grid_input(cfg, fnu912=logfnu912, fluxfilename=fluxname, table=cfg.table,
                     abundances=cfg.abundances)

    if cfg.run_cloudy:
        run_grid(nproc=cfg.nproc, overwrite=cfg.overwrite)

    models = parse_grid(cfg)

    filename = cfg.prefix + '_grid.sav.gz'
    print 'Writing to', filename
    saveobj(filename, models, overwrite=cfg.overwrite)
    savehdf5(filename.replace('.sav.gz', '.hdf5'), models,
             overwrite=cfg.overwrite)


    
if __name__ == '__main__':
    main()
