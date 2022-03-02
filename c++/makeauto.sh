#!/bin/bash

rm -f stamp-h1 aclocal.m4 Makefile.in autoscan.log config.status configure Makefile config.log config.h

aclocal 
automake 
autoconf

rm -rf aclocal.m4 autom4te.cache 
