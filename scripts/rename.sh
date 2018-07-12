#!/bin/bash

for f in `/bin/ls *xu.dtbmap`
do
    base=`basename $f xu.dtbmap`

    field=ksz
    if [ -e $base\xu.$field\map ] ; then mv $base\xu.$field\map $base\ux.$field\map; fi
    if [ -e $base\yu.$field\map ] ; then mv $base\yu.$field\map $base\uy.$field\map; fi
    if [ -e $base\zu.$field\map ] ; then mv $base\zu.$field\map $base\uz.$field\map; fi

    field=dtb
    if [ -e $base\xu.$field\map ] ; then mv $base\xu.$field\map $base\ux.$field\map; fi
    if [ -e $base\yu.$field\map ] ; then mv $base\yu.$field\map $base\uy.$field\map; fi
    if [ -e $base\zu.$field\map ] ; then mv $base\zu.$field\map $base\uz.$field\map; fi

    field=tau
    if [ -e $base\xu.$field\map ] ; then mv $base\xu.$field\map $base\ux.$field\map; fi
    if [ -e $base\yu.$field\map ] ; then mv $base\yu.$field\map $base\uy.$field\map; fi
    if [ -e $base\zu.$field\map ] ; then mv $base\zu.$field\map $base\uz.$field\map; fi

done
