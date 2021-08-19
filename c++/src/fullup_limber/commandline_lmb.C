#include <math.h>
#include "proto_lmb.h"
#include "allvars_lmb.h"
#include <unistd.h>

void usage(){

  if(myid==0){
    printf("\n usage: lmb <parameterfile> [options]\n");
    
    printf("\n OPTIONS:\n");
    printf("\n   -h show this message");
    printf("\n   -i <Input base name> [default = 'output'; all input files assumed to be prepended");
    printf("\n      with this string and in a directory called output in the working");
    printf("\n      directory (e.g. output/output.zreion, etc)]");
    printf("\n   -o <Output base name> [default = 'out'; all output files will be prepended");
    printf("\n      with this string and in a directory called curlpk in the working");
    printf("\n      directory (e.g. pk/out.ps, ps/out.info, etc)]");
    printf("\n   -v switch on verbose mode to print info to terminal instead of log file [default = OFF]"); 
    printf("\n   -u switch on delta_x = delta = 0 [default = OFF]");     
    printf("\n   -x switch on delta_x = 0         [default = OFF]");     
    printf("\n   -d switch on delta   = 0         [default = OFF]");     
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
  sprintf(clParameters.BaseIn,"output/output");
  sprintf(clParameters.BaseOut,"ps/out");

  if(myid==0) system("[ -d ps ] || mkdir ps");
  
  clParameters.verbose  = 0;
  clParameters.uniformu = 0;
  clParameters.uniformx = 0;
  clParameters.uniformd = 0;

  while ((c = getopt (argc, argv, "hvuxdz:i:o:")) != -1)
    switch (c)
      {
      case 'h':
	usage();
	break;
      case 'i':
	sprintf(clParameters.BaseIn,"output/%s",optarg);
	break;
      case 'o':
	sprintf(clParameters.BaseOut,"ps/%s",optarg);
	break;
      case 'v':
	clParameters.verbose = 1;
	break;
      case 'u':
	clParameters.uniformu = 1;
	break;
      case 'x':
	clParameters.uniformx = 1;
	break;
      case 'd':
	clParameters.uniformd = 1;
	break;
      case '?':
	if (optopt == 'i'){
	  if(myid==0) fprintf (stderr, "\n Option -%c requires an argument.\n", optopt);
	  usage();
	}
	if (optopt == 'o'){
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
    sprintf(fname,"%s.stdout",clParameters.BaseOut);
    freopen(fname,"w",stdout);
  }

  // Don't buffer stdout and stderr
  setvbuf(stdout, NULL, _IONBF, 0);
  setvbuf(stderr, NULL, _IONBF, 0);

  if(myid==0){
    printf("\n Command line:");
    for(int i=0;i<argc;i++) printf(" %s",argv[i]); printf("\n");  
  }

  return;

}

