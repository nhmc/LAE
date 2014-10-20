from subprocess import call
import os
import sys

usage = """ python run_emcee.py dirname"""


dirname, = sys.argv[1:]

if os.path.exists(dirname):
    call('rm -rf %s' % dirname, shell=1)

call('mkdir %s' % dirname, shell=1)

call('run_mcmc', shell=1)
call('plot_mcmc samples_burn.sav.gz', shell=1)
call('plot_mcmc samples_mcmc.sav.gz', shell=1)
call('mv samples* %s' % dirname, shell=1)
call('cp priors %s' % dirname, shell=1)
call('cp dont_use %s' % dirname, shell=1)
call('cp model.py %s' % dirname, shell=1)
call('cp model.cfg %s' % dirname, shell=1)
call('cp observed_logN %s' % dirname, shell=1)

call('mv fig/* %s' % dirname, shell=1)
call('rmdir fig', shell=1)
