#!/bin/bash

command=map2pdf.py

#### kSZ ####
field=ksz
min=-20; max=20
for f in `/bin/ls *0.06*.$field\map`
do base=`basename $f .$field\map`; $command $base.$field\map $base\_$field.pdf $min $max $field; done
min=-20; max=20
for f in `/bin/ls *0.09*.$field\map`
do base=`basename $f .$field\map`; $command $base.$field\map $base\_$field.pdf $min $max $field; done

#### 21-cm ####
field=dtb
kmin=0; max=4e4
for f in `/bin/ls *0.06*.$field\map`
do base=`basename $f .$field\map`; $command $base.$field\map $base\_$field.pdf $min $max $field; done

kmin=0; max=4e4
for f in `/bin/ls *0.09*.$field\map`
do base=`basename $f .$field\map`; $command $base.$field\map $base\_$field.pdf $min $max $field; done

#### tau ####
field=tau
min=0.04; max=0.08
for f in `/bin/ls *0.06*.$field\map`
do base=`basename $f .$field\map`; $command $base.$field\map $base\_$field.pdf $min $max $field; done

min=0.07; max=0.11
for f in `/bin/ls *0.09*.$field\map`
do base=`basename $f .$field\map`; $command $base.$field\map $base\_$field.pdf $min $max $field; done
