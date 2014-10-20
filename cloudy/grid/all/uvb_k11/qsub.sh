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
cd /lustre/projects/p028_swin/ncrighton/cloudy/J0004_NHI_2/all/uvb_k11
# which will now be the working directory ...
echo "Working directory:"
echo "$PBS_O_WORKDIR"
echo "y" > gridin
run_cloudy_grid < gridin > log
rm gridin
date
