#!/bin/bash

bin=../bin
nprocs=$1

if [ ! $nprocs > 0 ] ; then
    echo "nprocs = $nprocs not > 0, exiting"
    exit
fi

source ../scripts/banner.sh
seed=18937
pkfile=wmap5_0_m.pk

mybanner "RFAST TEST WITH NPROCS = $nprocs"

mybanner "Testing initial conditions"
mpirun -np $nprocs $bin/ics parameterfiles/param.ics -p $pkfile -o delta -b 8e3 -n 512 -v -s seed

mybanner "Testing delta2zreion"
mpirun -np $nprocs $bin/delta2zreion parameterfiles/param.d2z -f test -m 1e9 -z 500 -r 32 -v

mybanner "Testing all sky map"
mpirun -np $nprocs $bin/allskymap parameterfiles/param.asm -i test -o test -v

mybanner "Testing flat sky map"
mpirun -np $nprocs $bin/flatskymap parameterfiles/param.fsm -i test -o test -v


