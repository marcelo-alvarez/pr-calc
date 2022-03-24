#!/bin/bash
  
# test successive runs on Cori at NERSC
# assuming colossus is installed
# and user is on an interactive node

# set root location
cd ..
export PRCALC_ROOT=$PWD
cd example

# run script 
./runscript.sh 32
