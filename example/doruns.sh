#!/bin/bash

bin=../bin

source ../../scripts/banner.sh

box=$1
ngrid=$2
nprocs=$3
zmin=$4
run=$5

seed=25645
pkfile=wmap5_0_m.pk

mybanner "RFAST TEST WITH NPROCS = $nprocs"

#make linear density
mybanner "Testing initial conditions"
srun -n $nprocs $bin/ics parameterfiles/param.ics -p $pkfile -o delta -b $box -n $ngrid -v -s $seed

#cp parameterfiles/param.d2z_gen parameterfiles/param.d2z
#replace.pl BOXSIZE_REPLACE $box parameterfiles/param.d2z
#replace.pl NGRID_REPLACE $ngrid parameterfiles/param.d2z

for mmin in 2e9 3e9 4e9; do
        for zeta in 50 75 100; do
		echo writing out zreion tables
		python /global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/py/tables/fcoll_zreion.py $mmin $zeta
		for lambda in 200 300 400; do

                        base=$run\_$lambda\_$zeta\_$mmin

                        echo Processing $lambda$zeta$mmin ...

                        mybanner "Testing delta2zreion"
                        srun -n $nprocs $bin/delta2zreion parameterfiles/param.d2z -z $zeta -r $lambda -m $mmin -f $base -v

                        ./single.sh $box $ngrid $nprocs $zmin $base
#                       mybanner "Testing all sky map"
#                       srun -n $nprocs $bin/allskymap parameterfiles/param.asm -i $base -o $base -v

                       mybanner "Testing flat sky map"
                       srun -n $nprocs $bin/flatskymap parameterfiles/param.fsm -i $base -o $base -v
                        #read the ksz map and output angular power spectrum to a pdf file
                        mybanner "Output angular power spectrum to file"
                        /global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/scripts/map2pk.py maps/test.kszmap $base.ksz_cl.txt

                        # show the ksz map and output it to a pdf file
                        mybanner "Show ksz map"
                        /global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/scripts/showmap.py maps/test.kszmap $base.kszmap.pdf

                        echo Work done. Returning to upper directory

                        mkdir -p Combined_Output
                        mv $base.ksz_cl.txt output/$base.history output/$base.stdout Combined_Output
#                       rm $base.*

                done
        done
done

