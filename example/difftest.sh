#!/bin/bash

# test successive runs on Cori at NERSC
# assuming colossus is installed
# and user is on an interactive node

# set root location
cd ..
export PRCALC_ROOT=$PWD
cd example

# run script 
./runscript.sh 1

mv ICs/delta ICs/delta0
mv output/test.zreion test.zreion0
mv maps/test.kszmap maps/test.kszmap0
mv test.ksz_cl.txt test.ksz_cl.txt0

# run script 
./runscript.sh 1

echo "testing IC differences:      "; diff ICs/delta ICs/delta0
echo "testing zreion differences:  "; diff output/test.zreion test.zreion0
echo "testing ksz map differences: "; diff maps/test.kszmap maps/test.kszmap0
echo "testing ksz cls differences: "; diff test.ksz_cl.txt test.ksz_cl.txt0
