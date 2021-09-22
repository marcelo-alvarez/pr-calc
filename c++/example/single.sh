#!/bin/bash

bin=../bin

source ../../scripts/banner.sh

box=$1
ngrid=$2
nprocs=$3
z=$4
base=$5

#zf=printf "08.5f" $z

#cp parameterfiles/param.fsm_gen parameterfiles/param.fsm_$zf
#replace.pl FINALREDSHIFT_REPLACE $zf parameterfiles/param.fsm_$zf
#replace.pl NGRID_REPLACE $ngrid parameterfiles/param.fsm_$zf
#replace.pl BOXSIZE_REPLACE $box parameterfiles/param.fsm_$zf

mybanner "Testing all sky map"
#srun -n $nprocs $bin/allskymap parameterfiles/param.asm -i $base -o $base

mybanner "Testing flat sky map"
srun -n $nprocs $bin/flatskymap parameterfiles/param.fsm -i $base -o $base -v
~                                                                              
