#!/bin/bash

logspace=$1

for f in `/bin/ls 8e3*.cl`
do 
    base=`basename $f .cl`
    echo $logspace | bincl $f $base.bcl
done

