#include <math.h>
#include "proto_d2z.h"
#include "allvars_d2z.h"
#include <unistd.h>

void usage(){

  if(myid==0){
    printf("\n usage: d2z <parameterfile> [options]\n");
    
    printf("\n OPTIONS:\n");
    printf("\n   -h show this message");
    printf("\n   -f <File base name> [default = 'output'; all output files will be prepended");
    printf("\n      with this string and placed in a directory called output in the working");
    printf("\n      directory (e.g. output/output.history, output/output.zreion, etc)]");
    printf("\n   -r <Rmfp> [default = 256; this is the mean free path in units of Mpc/h]");
    printf("\n   -m <Mmin> [default = 1e8; this is the minimum halo mass");
    printf("\n      capable of hosting ionizing sources, in Msun]");
    printf("\n   -z <zeta> [default =  10; this is the efficiency of halos in ionizing");
    printf("\n      their surroundings -- i.e. zeta = M_HII / Mhalo * rho_matter / rho_hydrogen]");
    printf("\n   -l switch on lean mode which reduces memory footprint by ~4 bytes per cell but doubles");
    printf("\n      number of FFTs that must be performed [default = OFF]");
    printf("\n   -v switch on verbose mode to print stdout to terminal instead of .stdout file [default = OFF]"); 
    printf("\n   -e switch on verbose mode to print stderr to terminal instead of .stderr file [default = OFF]"); 
    printf("\n\n");
  }

  MPI_Finalize();
  exit(0);

}

void CommandLine(int argc, char *argv[])
{

  int c;

  opterr=0;
  
  if(argc<2) usage();

  sprintf(clParameters.Paramfile,"%s",argv[1]);
  sprintf(clParameters.Basename,"output/output");

  if(myid==0) system("[ -d output ] || mkdir output");
  
  clParameters.Rmax  = 256;
  clParameters.zeta  = 10;
  clParameters.M_min = 1e8;
  clParameters.verbose = 0;
  clParameters.everbose = 0;
  clParameters.lean = 0;

  while ((c = getopt (argc-1, &argv[1], "hvelf:r:m:z:")) != -1)
    switch (c)
      {
      case 'h':
	usage();
	break;
      case 'f':
	sprintf(clParameters.Basename,"output/%s",optarg);
	break;
      case 'r':
	clParameters.Rmax  = atof(optarg);
	break;
      case 'm':
	clParameters.M_min = atof(optarg);
	break;
      case 'z':
	clParameters.zeta  = atof(optarg);
	break;
      case 'v':
	clParameters.verbose = 1;
	break;
      case 'e':
	clParameters.everbose = 1;
	break;
      case 'l':
	clParameters.lean = 1;
	break;
      case '?':
	if (optopt == 'r' || optopt == 'f' || optopt == 'm' || optopt == 'z'){
	  if(myid==0) fprintf (stderr, "\n Option -%c requires an argument.\n", optopt);
	  usage();
	}
	else if (isprint (optopt)){
	  if(myid==0) fprintf (stderr, "\n Unknown option `-%c'.\n", optopt);
	  usage();
	}
	else{
	  if(myid==0) fprintf (stderr,
		   "Unknown option character `\\x%x'.\n",
		   optopt);
	  usage();
	}
	return;
      default:
	usage();
      }

  if(clParameters.verbose == 0){
    // Redirect stdout to output file
    char fname[256];
    sprintf(fname,"%s.stdout",clParameters.Basename);
    freopen(fname,"w",stdout);
  }

  if(clParameters.everbose == 0){
    // Redirect stderr to output file
    char fname[256];
    sprintf(fname,"%s.stderr",clParameters.Basename);
    freopen(fname,"w",stderr);
  }

  // Don't buffer stdout and stderr
  setvbuf(stdout, NULL, _IONBF, 0);
  setvbuf(stderr, NULL, _IONBF, 0);

  if(myid==0){
    printf("\n Command line: %s %s",argv[0],argv[argc-1]);
    for(int i=1;i<argc-1;i++) printf(" %s",argv[i]); printf("\n");  
  }

  return;

}

