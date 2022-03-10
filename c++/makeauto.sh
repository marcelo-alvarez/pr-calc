#!/bin/bash

source clean_automake.sh

touch NEWS README AUTHORS ChangeLog COPYING

aclocal 
automake --add-missing
autoconf

rm -rf COPYING NEWS README AUTHORS ChangeLog
