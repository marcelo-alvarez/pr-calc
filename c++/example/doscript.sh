#!/bin/bash -l
  
#SBATCH -L SCRATCH
#SBATCH -C haswell
#SBATCH -J pr-calc-test
#SBATCH -N 40
#SBATCH -t 01:30:00
#SBATCH -p regular
#SBATCH -o pr-calc.log

cd $SLURM_SUBMIT_DIR

#Turn off implicit threadaing in Python, R
export OMP_NUM_THREADS=1

module load python
module load gsl
#module load openmpi

box=4e3
ngrid=2048
nprocs=64
#256
zmax=15
run=$box\_$ngrid

HOSTS=$(scontrol show hostnames $SBATCH_NODELIST | tr '\n' ,)

#Run the runscript
./doscript.sh $box $ngrid $nprocs $zmax $run
