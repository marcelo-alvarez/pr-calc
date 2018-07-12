#include <math.h>
#include "proto_ics.h"
#include "allvars_ics.h"
#include <unistd.h>

void usage(){

  if(myid==0){
    printf("\n usage: ics <parameterfile> [options]\n");
    
    printf("\n OPTIONS:\n");
    printf("\n   -h show this message");
    printf("\n   -v switch on verbose mode [default = OFF]"); 
    printf("\n   -p <power spectrum filename> [default = 'power.dat']");
    printf("\n      this will be found in the directory ICs (e.g. ICs/power.dat)");
    printf("\n   -o <output basename> [default = 'delta']");
    printf("\n      this will be placed in the directory ICs (e.g. ICs/delta)");
    printf("\n   -b <boxsize in Mpc/h> [default = 1e2]");
    printf("\n   -s <random seed> [default = 13579]");
    printf("\n   -n <N> [default = 512]");
    printf("\n\n");
  }

  MPI_Finalize();
  exit(0);

}

void CommandLine(int argc, char **argv)
{

  int c;

  opterr=0;
  
  if(argc<2) usage();

  sprintf(clParameters.ParameterFile,"%s",argv[1]);

  sprintf(clParameters.PowerSpectrumFile,"ICs/power.dat");
  sprintf(clParameters.Base,"ICs/delta");

  if(myid==0) system("[ -d ICs ] || mkdir ICs");

  clParameters.BoxSize = 1e2;
  clParameters.Seed    = 10;
  clParameters.N       = 512;
  clParameters.verbose = 0;

  while ((c = getopt (argc-1, &argv[1], "hvp:o:b:s:n:")) != -1){
    switch (c)
      {
      case 'h':
	usage();
	break;
      case 'v':
	clParameters.verbose = 1;
	break;
      case 'p':
	sprintf(clParameters.PowerSpectrumFile,"ICs/%s",optarg);
	break;
      case 'o':
	sprintf(clParameters.Base,"ICs/%s",optarg);
	break;
      case 'b':
	clParameters.BoxSize = atof(optarg);
	break;
      case 's':
	clParameters.Seed = atoi(optarg);
	break;
      case 'n':
	clParameters.N = atoi(optarg);
	break;
      case '?':
	if (optopt == 'p' || optopt == 'o' || optopt == 'b' || optopt == 's' || optopt == 'n'){
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
  }
  if(clParameters.verbose == 0){
    // Redirect stdout to output file
    char fname[256];
    sprintf(fname,"%s.stdout",clParameters.Base);
    freopen(fname,"w",stdout);
  }

  sprintf(DeltaFile,"%s",clParameters.Base);

  // Don't buffer stdout and stderr
  setvbuf(stdout, NULL, _IONBF, 0);
  setvbuf(stderr, NULL, _IONBF, 0);

  if(myid==0){
    printf("\n Command line:");
    for(int i=0;i<argc;i++) printf(" %s",argv[i]); printf("\n");  
  }

  return;

}

