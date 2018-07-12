#!/bin/bash

for f in `/bin/ls 8e3*z_ksz-dtb.cl`
do

  base=`basename $f z_ksz-dtb.cl`

  code=_ksz-dtb
  outfile=$base$code.cl
  addthem.awk $base\x$code.cl $base\y$code.cl $base\z$code.cl > $outfile

  code=_ksz-tau
  outfile=$base$code.cl
  addthem.awk $base\x$code.cl $base\y$code.cl $base\z$code.cl > $outfile

  code=_dtb-tau
  outfile=$base$code.cl
  addthem.awk $base\x$code.cl $base\y$code.cl $base\z$code.cl > $outfile

done

