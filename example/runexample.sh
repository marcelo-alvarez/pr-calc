#!/bin/bash

#
# This script requires that the root directory of pr-calc is set in
# the environment variable PRCALC_ROOT
#

# load modules bash scripting helper functions 
module load python gsl
EXAMPLE_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PRCALC_ROOT=$EXAMPLE_DIR/..
source $PRCALC_ROOT/scripts/helpers.sh

# parse command line 
nprocs=$1 # number of processes passed as first argument 
if [ ! $nprocs > 0 ] ; then
    echo "nprocs = $nprocs not > 0, exiting"
    exit
fi

# set run name 
rname=$boxsize\_$Nboxres\_$lambda\_$zeta\_$mmin

# set parameters 
mmin=2e9    # minimum source halo mass in Msun
zeta=30     # ionizing efficiency 
lambda=140  # mean free path in Mpc/h
seed=18937  # set the random seed for the ICs
nR=100      # number of R bins in zreion table 
ndeltaR=101 # number of delta_R bins in zreion table (not used currently)
muKmin=-10  # minimum of kSZ map image colorbar in muK
muKmax=10   # maximum of kSZ map image colorbar in muK

# directories 
rundir=$PRCALC_ROOT/example                    # main run directory 
icsdir=$PRCALC_ROOT/example/ICs                # location of ICs and P(k)
cpbindir=$PRCALC_ROOT/c++/bin                  # location of compiled binaries
paramdir=$PRCALC_ROOT/example/parameterfiles   # location of parameter files
pytabdir=$PRCALC_ROOT/py/tables                # location of Pk.py, fcoll_zreion.py
pyscrdir=$PRCALC_ROOT/scripts/                 # location of map2pk.py, showmap.py
pydirdir=$PRCALC_ROOT/py                       # location of setParams.py

# python scripts 
pysetparams=$pydirdir/setParams.py             # parameter file generator
pytable=$pytabdir/Pk.py                        # power spectrum generator
pyfcoll=$pytabdir/fcoll_zreion.py              # fcoll and zreion table generator
pymap2pk=$pyscrdir/map2pk.py                   # map power spectrum calculator
pyshowmap=$pyscrdir/map2pdf.py                 # map displayer

# compiled c++ binaries 
icsbinary=$cpbindir/ics                        # ICs compiled binary
d2zbinary=$cpbindir/delta2zreion               # delta2zreion compiled binary
fsmbinary=$cpbindir/flatskymap                 # flat sky map compiled binary
fsmbinary=$cpbindir/flatskymap                 # flat sky map compiled binary 

# parameter file names 
iniparams=$paramdir/globalParams.ini           # uber ini parameter file  
colparams=$paramdir/param.col                  # Colossus parameter file
fsmparams=$paramdir/param.fsm                  # flat sky map parameter file 
icsparams=$paramdir/param.ics                  # ICs parameter file
d2zparams=$paramdir/param.d2z                  # delta2zreion parameter file 
fsmparams=$paramdir/param.fsm                  # flat sky map parameter file 

# data file names 
pkfile=pkfile.txt                              # P(k) file is $PWD/ICs/$pkfile
deltafile=delta                                # ICs data file is $PWD/ICs/$deltafile
kszmapbin=$rundir/maps/$rname.kszmap           # kSZ map 
kszcls=$rundir/$rname.ksz_cl.txt               # kSZ angular power spectrum
kszmapimg=$rundir/maps/$rname                  # kSZ map pdf  

printf "--- Setting cosmo params and writing parameter files"
read -r boxsize Nboxres <<< $(python $pysetparams $iniparams $paramdir)

printf "--- Generating linear power spectrum at z=0"
python $pytable $colparams $icsdir/$pkfile 

printf "--- Generating initial conditions"
srun -n $nprocs $icsbinary $icsparams -p $pkfile -o $deltafile -b $boxsize -n $Nboxres -v -s $seed

printf "--- Generating zreion tables"
python $pyfcoll $colparams $mmin $zeta $nR $ndeltaR

printf "--- Running delta2zreion"
srun -n $nprocs $d2zbinary $d2zparams -f $rname -m $mmin -z $zeta -r $lambda -v

printf "--- Testing flat sky map"
srun -n $nprocs $fsmbinary $fsmparams -i $rname -o $rname -v

printf "--- Output angular power spectrum to file"
python $pymap2pk $kszmapbin $kszcls 

printf "--- Show ksz map"
python $pyshowmap $kszmapbin $kszmapimg $muKmin $muKmax ksz 

