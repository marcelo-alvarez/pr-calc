#!/bin/bash

./clean
cd ..
tar -zcvf pr-calc.tgz pr-calc -X pr-calc/exclude
scp pr-calc.tgz malvarez@gw.cita.utoronto.ca:

