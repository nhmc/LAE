
M = loadobj('qg_grid.sav.gz')

import h5py
fh = h5py.File('filename.hdf5', 'w')

fh.attrs['help'] = M['help']
fh.attrs['redshift'] = M['redshift']

if 'uvb_tilt' in M:
    fh.attrs['uvb_tilt'] = M['uvb_tilt']

for k in 'nH Z NHI U Tstop filename cont Tgas'.split():
    d = fh.create_dataset(k, M[k].shape, dtype=M[k].dtype, compression='gzip')
    d[:] = M[k]

for k in 'N Nex gas_abun dust_abun'.split():
   g = fh.create_group(k)
   for atom in M[k]:
       d = g.create_dataset(atom, M[k][atom].shape, dtype=M[k][atom].dtype,
                            compression='gzip')
       d[:] = M[k][atom]

fh.close()
