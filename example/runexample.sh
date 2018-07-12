#!/bin/bash

bin=../bin
nprocs=$1

if [ ! $nprocs > 0 ] ; then
    echo "nprocs = $nprocs not > 0, exiting"
    exit
fi

source ../scripts/banner.sh
seed=18937
pkfile=0.044_0.226_0.73_0.7_0.96_0.8.pk
#pkfile=0.048_0.262_0.69_0.68_0.97_0.82.pk
mybanner "RFAST TEST WITH NPROCS = $nprocs"

mybanner "Testing initial conditions"
mpirun -np $nprocs $bin/ics parameterfiles/param.icsw -p $pkfile -o delta -b 8e3 -n 512 -v -s seed

mybanner "Testing delta2zreion"
mpirun -np $nprocs $bin/delta2zreion parameterfiles/param.d2zw -f testw -m 1e9 -z 500 -r 32 -v

mybanner "Testing all sky map"
mpirun -np $nprocs $bin/allskymap parameterfiles/param.asmw -i testw -o testw -v

mybanner "Testing flat sky map"
#mpirun -np $nprocs $bin/flatskymap parameterfiles/param.map -i testw -o testw -v

mybanner "Testing Limber"
mpirun -np $nprocs $bin/limber parameterfiles/param.lmbw -i testw -o testw -v 

