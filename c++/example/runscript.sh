#!/bin/bash

bin=../bin
nprocs=$1

mmin=3e9
zeta=100
echo $mmin
echo $zeta

if [ ! $nprocs > 0 ] ; then
    echo "nprocs = $nprocs not > 0, exiting"
    exit
fi

source ../../scripts/banner.sh
seed=18937
pkfile=wmap5_0_m.pk

mybanner "RFAST TEST WITH NPROCS = $nprocs"

mybanner "Testing initial conditions"
srun -n $nprocs $bin/ics parameterfiles/param.ics -p $pkfile -o delta -b 4e3 -n 2048 -v -s $seed

mybanner "Testing delta2zreion"
python /global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/py/tables/fcoll_zreion.py $mmin $zeta
srun -n $nprocs $bin/delta2zreion parameterfiles/param.d2z -f test -m $mmin -z $zeta -r 300 -v

mybanner "Testing all sky map"
#srun -n $nprocs $bin/allskymap parameterfiles/param.asm -i test -o test -v

mybanner "Testing flat sky map"
srun -n $nprocs $bin/flatskymap parameterfiles/param.fsm -i test -o test -v

mybanner "Output angular power spectrum to file"
/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/scripts/map2pk.py maps/test.kszmap test.ksz_cl.txt

# show the ksz map and output it to a pdf file
mybanner "Show ksz map"
/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/scripts/showmap.py maps/test.kszmap test.kszmap.pdf

