Introduction
------------

Codes for generating reionization histories and synthetic sky maps as described in 
  Alvarez, M. A. and Abel, T., 2012, ApJ, 747, 126
and
  Alvarez, M. A., 2017, ApJ, 824, 118
  
Installation
------------

In order to compile you need to have the following installed:
  - MPI 
  - GNU Scientific Library (gsl)
  - FFTW2 (fftw-2.1.5)
  - FFTW3 (fftw-3.3.3)
  - CFITSIO 
  - Colossus
If you don't have FFTW2, FFTW3, or CFITSIO, these will be automatically
downloaded and installed in your home directory at installation.

Download::
	git clone https://github.com/marcelo-alvarez/pr-calc $installpath/pr-calc
	
Modify pr-calc/c++/Settings.<yourmachinename> to reflect the local paths to these 
libraries and include files, then::
	cd $installpath/pr-calc/c++
    ./install_pr-calc <yourmachinename>

Executables will be located in pr-calc/bin.

Example
-------

There is an example script to run the codes in the package, using srun, in the example directory::
	cd $installpath/pr-calc/example
	./runexample.sh <nprocs>
where <nprocs> is the number of MPI processes to run with. This will run several tasks:
1. Create initial linear density field (ics)
2. Generate reionization field from density field (delta2zreion)
3. Generate an all-sky map of 21-cm, tau_es, and dT_kSZ from density 
   and reionization field (allskymap)
4. Generate a rectangular map of a portion of the sky for the same
   quantities in the flat sky approximation (flatskymap)







