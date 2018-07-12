#!/bin/bash

for f in `/bin/ls 8e3*.dtbmap`
do
    base=`basename $f .dtbmap`
    ~/src/taumaps/map2cc.py $base.kszmap $base.dtbmap $base\_ksz-dtb.cl
    ~/src/taumaps/map2cc.py $base.kszmap $base.taumap $base\_ksz-tau.cl
    ~/src/taumaps/map2cc.py $base.dtbmap $base.taumap $base\_dtb-tau.cl
done
