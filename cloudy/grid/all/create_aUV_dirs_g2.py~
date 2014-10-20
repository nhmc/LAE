from subprocess import call
import os
import numpy as np

#modify the redshift, prefix, and NHI, Z, nH, nproc values below.

REDSHIFT = 2.465
aUVvals = np.arange(-2, 1.01, 0.5)

template = """\
# Path to the CUBA file giving the UV background. 
#cuba_name = /Users/ncrighton/Code/Repo/QSOClustering-svn/misc/CUBA/Q1G01/bkgthick.out
# Prefix for final grid file
prefix = qg
# Redshift for the UV background
z = {}
# Minimum, maximum and step for neutral hydrogen column density (cm^-2)
logNHI =  13.0 17.5 0.5
# Minimum, maximum and step for metallicity
logZ   =  -3.0  0.25  0.5
# Minimum, maximum and step for hydrogen density (cm^-3)
lognH  =  -4.0 -1.0  0.4

trimming = -10

# number of processors to use
nproc = 16
# overwrite any existing files?
overwrite = True

uvb_tilt = {:.3f}

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
#run_cloudy = False
"""


qsub_template = '''\
#!/bin/bash
# First you set attributes for the queue. 
# name the job (optional)
#PBS -N cloudy
# rename the standard output file (optional)[jobname.ojobid is default]
#PBS -o jobout
# rename the error file (optional)[jobname.ejobid is default]
#PBS -e joberr
# select the queue you want (gstar or sstar)
#PBS -q sstar
# set your resource requests ... 
# here we are asking for 16 cpu cores on 1 node
#PBS -l nodes=1:ppn=16
# also ~2gb RAM per cpu core
#PBS -l pmem=2048mb
# and a walltime of 1 day (format is days:hours:mins:secs)
#PBS -l walltime=01:00:00:00

date
source ~/.bashrc

# Now you can print some information about the nodes and GPUs that 
# have been assigned (optional)
echo "Deploying job to CPUs ..."
cat "$PBS_NODEFILE"

# you may then want to change to the directory where your code sits
cd {}
# which will now be the working directory ...
echo "Working directory:"
echo "$PBS_O_WORKDIR"
echo "y" > gridin
run_cloudy_grid < gridin > log
rm gridin
date
'''

prefix = os.path.abspath('./')

for i,aUV in enumerate(aUVvals):
    dirname = 'uvb_k{:02d}'.format(i) 
    call('mkdir ' + dirname, shell=1) 
    fh = open(dirname + '/grid.cfg', 'wb')
    s = template.format(REDSHIFT, aUV)
    fh.write(s)
    fh.close()

    fh = open(dirname + '/qsub.sh', 'wb')
    s = qsub_template.format(prefix + '/' + dirname)
    fh.write(s)
    fh.close()
