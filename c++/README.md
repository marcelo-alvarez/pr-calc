1. Compilation

   In order to compile you need to have the following installed:
      - MPI 
      - GNU Scientific Library (gsl)
      - FFTW2 (fftw-2.1.5)
      - FFTW3 (fftw-3.3.3)
      - CFITSIO 
   If you don't have FFTW2, FFTW3, or CFITSIO, these will be automatically
   downloaded and installed in your home directory

   Modify Settings.<yourmachinename> to reflect the local paths to these
   libraries and include files, then do
	./install_cal-pr <yourmachinename>

   The executables will be located in bin/ 

2. Usage for the main code, delta2zreion

   A) Command line options

      usage: delta2zreion <parameterfile> [options]

      OPTIONS:

	-h show this message
     	-f <File base name> [default = 'output'; output files are prepended
      	   with this string in a directory called output in the working
      	   directory (e.g. output/output.history, output/output.zreion, etc)]
   	-r <Rmfp> [default = 256; this is the mean free path in units of Mpc/h]
   	-m <Mmin> [default = 1e8; this is the minimum halo mass
      	   capable of hosting ionizing sources, in Msun]
   	-z <zeta> [default =  10; this is the ionizing efficiency of halos
      	   i.e. zeta = M_HII / Mhalo * rho_matter / rho_hydrogen]
   	-v verbose mode -- prints info to terminal [default = OFF]

   B) Parameter file

      The parameter file contains settings which don't change with the initial
      input file. See the file 'example/parameterfiles/param.zreion' for more 
      info.

3. Example run

   There is an example script to run the codes in the package in  example/ 
      cd example
      ./runexample.sh <nprocs>
   where <nprocs> is the number of MPI processes to run with. This should
   do 5 different tasks:
      1. Create initial linear density field (ics)
      2. Generate reionization field from density field (delta2zreion)
      3. Generate an all-sky map of 21-cm, tau_es, and T_kSZ from density 
         and reionization field (allskymap)
      4. Generate a rectangular map of a portion of the sky for the same
         quantities in the flat sky approximation (flatskymap)
      5. Use the Limber approximation to generate the power spectrum of
         the kSZ (limber)






