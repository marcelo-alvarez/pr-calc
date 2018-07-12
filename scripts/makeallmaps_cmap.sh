#!/bin/bash

fname=4e3_2048
size=2048
bytes=`echo $size | awk '{print $1*$1*4}'`
echo $bytes

min=-20
max=30
for f in `/bin/ls $fname*.kszmap`
do
    base=`basename $f .kszmap`
    tail -c $bytes $f > foo
    colormap foo $base.ppm $min $max 5 $size $size
    convert $base.ppm $base\_ksz.jpg
done

min=0
max=4e4
for f in `/bin/ls $fname*.dtbmap`
do
    base=`basename $f .dtbmap`
    tail -c $bytes $f > foo
    colormap foo $base.ppm $min $max 1 $size $size
    convert $base.ppm $base\_dtb.jpg
done

min=0.08
max=0.11
for f in `/bin/ls $fname*.taumap`
do
    base=`basename $f .taumap`
    tail -c $bytes $f > foo
    colormap foo $base.ppm $min $max 3 $size $size
    convert $base.ppm $base\_tau.jpg
done

rm foo *.ppm