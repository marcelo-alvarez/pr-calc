#!/bin/bash

# test successive runs on Cori at NERSC
# assuming colossus is installed
# and user is on an interactive node

module load gsl python
source activate py3

# set root location
cd ..
export PRCALC_ROOT=$PWD
cd example

# run test
./runscript.sh 1
