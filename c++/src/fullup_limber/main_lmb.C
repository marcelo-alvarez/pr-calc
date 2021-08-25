#include "allvars_lmb.h"
#include "proto_lmb.h"
#include <math.h>

/////////////////////////////////////////////////////////////////////////
// THIS PROGRAM MAKES AN OPTICAL DEPTH MAP FROM A REIONIZATION REDSHIFT
// FIELD, ASSUMING THAT THE MATTER DENSITY IS UNIFORM.
//
//                                                AUTHOR: MARCELO ALVAREZ  
//                                             LAST EDIT:        01.25.13
/////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // Parse command line
  CommandLine(argc, &argv[0]);

  // Read parameter file
  ReadParameterFile();
  int N=Parameters.N;
  Parameters.BoxSize/=Parameters.h;
  nl=Parameters.nl;

  // Growth factor at initial time
  DInit=growth(Parameters.zInit,Parameters.Omegam,Parameters.Omegal,
	       Parameters.w);

  // Growth factor at z>>1
  float zhigh=100;
  float Dhigh=growth(zhigh,Parameters.Omegam,Parameters.Omegal,
	       Parameters.w);
  float Df = (1+zhigh)*Dhigh/DInit;

  // Set cl pre-factor
  double Mpc=3.08567758e24;
  double ccgs=2.998e10;
  double thomson=6.652458734e-25;
  double h0=Parameters.h*100.*1e5/Mpc;
  double omegam=Parameters.Omegam;
  double tcmb=2.726e6;
  double Yp=0.25;
  double mp=1.67262178e-24;
  double gcgs=6.674e-8;
  double rhoc=3.*h0*h0/8/M_PI/gcgs;
  double rhob=Parameters.Omegab*rhoc;

  double nh0 = (1-Yp)*rhob/mp;
  double nhe0 = Yp*rhob/4/mp;
  double ne0 = nh0 + nhe0;

  // clprefac_prp is in units of muK^2/Mpc 
  clprefac_prp = pow(thomson,2)*pow(ne0,2)*pow(tcmb,2)/
    4/M_PI/h0/sqrt(omegam)*ccgs*3.086e24; 

  // clprefac_par is in units of muK^2/Mpc^3 
  clprefac_par = pow(thomson,2)*pow(ne0,2)*pow(tcmb,2)/
    2/M_PI*h0*sqrt(omegam)/ccgs*pow(3.086e24,3);

  // clprefac_dop is in units of muK^2/Mpc^5
  clprefac_dop = pow(thomson,2)*pow(ne0,2)*pow(tcmb,2)/
    2/M_PI*pow(h0,3)*pow(omegam,1.5)/pow(ccgs,3)*pow(Mpc,5)*pow(Df,2);

  // clprefac_bnd is in units of muK^2/Mpc^4
  clprefac_bnd = pow(thomson,2)*pow(ne0,2)*pow(tcmb,2)/
    4/M_PI*pow(h0,2)*omegam/pow(ccgs,2)*pow(Mpc,4)*pow(Df,2);

  if(myid==0) printf("\n clprefac_prp, clprefac_par, clprefac_dop = %le %le %le\n",clprefac_prp,clprefac_par,clprefac_dop);

  // Set redshift to radius table
  Redshift2RadiusTable = new double[NZTABLE];
  SetRedshift2RadiusTable(Parameters.h, Parameters.Omegam, Parameters.Omegal, Redshift2RadiusTable);

  // Make FFTW plans
  plan  = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,N, N, N,
				  FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  iplan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,N, N, N,
				  FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
  
  rfftwnd_mpi_local_sizes(plan, &local_nx, &local_x_start,
			  &local_ny_after_transpose,
			  &local_y_start_after_transpose,
			  &total_local_size);

  // Allocate arrays
  AllocateArrays();

  // Read input files
  char ReionFile[256];
  sprintf(ReionFile,"%s.zreion",clParameters.BaseIn);

  ReadReionFile(ReionFile);
  ReadDeltaFile(Parameters.DeltaFile);

  GetMeanSigma(delta,1,N,Nlocal);
  GetMeanSigma(zreion,1,N,Nlocal);

  // Obtain power spectrum for density contrast at input time
  PowerSpectrum(delta,ps_delta,1);

  // Set wavenumber to pdelta table
  Wavenumber2PdeltaTable = new double[Parameters.nk];
  Wavenumbers4PdeltaTable = new double[Parameters.nk];
  SetWavenumber2P1DTable(Parameters.nk, Wavenumber2PdeltaTable, Wavenumbers4PdeltaTable, ps_delta, ook2);

  // Set history table
  Redshift2HistoryTable = new double[NZTABLE];
  SetRedshift2HistoryTable(Parameters.N, Nlocal, zreion, Redshift2HistoryTable);

  // Set Dg table
  Redshift2DgTable = new double[NZTABLE];
  SetRedshift2DgTable(Redshift2HistoryTable,Redshift2DgTable);

  // Set tau table
  Redshift2TauTable = new double[NZTABLE];
  SetRedshift2TauTable(Parameters.h, Parameters.Omegab, Parameters.Omegam, Redshift2TauTable, Redshift2HistoryTable);

  // Find minimum and maximum l-values, based upon zmin, zmax, kmin, and kmax
  float zmin=Parameters.zmin;
  float zmax=Parameters.zmax;
  float rmin=Redshift2Float(zmin, Redshift2RadiusTable);
  float rmax=Redshift2Float(zmax, Redshift2RadiusTable);
  float kmin=Parameters.kmin;
  float kmax=Parameters.kmax;

  lmin=rmax*kmin*2;
  lmax=rmin*kmax*0.5;

  // Loop over redshifts
  int nz=Parameters.nz;
  dz=(zmax-zmin)/nz;
  for(int iz=0;iz<nz;iz++){
    
    float Redshift = zmin + (iz+0.5)*dz;

    if(myid==0) printf("\n\n Doing redshift %d of %d z=%f --------- \n",iz+1,nz,Redshift);
  
    // Get scale factor and hubble parameter (in speed of light/Mpc) at current time
    aCurr=1/(1+Redshift);
    HCurr=Parameters.h*sqrt(Parameters.Omegam*pow((1+Redshift),3)+1-Parameters.Omegam);
    HCurr*=100/2.998e5;
    A = tcmb*ne0*thomson*pow(1+Redshift,2)*3.086e24; // prefactor in units of muK/Mpc

    // Current growth factor
    DCurr=growth(Redshift,Parameters.Omegam,Parameters.Omegal,Parameters.w);		 
    
    // Use linear theory to go from density to velocity
    Delta2Velocity();

    PowerSpectrum(prp_x,ps_prpx,1);
    PowerSpectrum(prp_y,ps_prpy,1);
    PowerSpectrum(prp_z,ps_prpz,1);

    PowerSpectrum(delta,ps_delta,1);

    // Generate 'specific momentum' field
    GenerateMomentum();

    // Seperate momentum into parallel and perpendicular components
    SeparateMomentum();

    // Differentiate parallel momentum
    DifferentiateMomentum();

    // Obtain power spectrum for each component
    PowerSpectrum(prp_x,ps_prpx,0);
    PowerSpectrum(prp_y,ps_prpy,0);
    PowerSpectrum(prp_z,ps_prpz,0);

    PowerSpectrum(delta,ps_delta,1);

    FILE *fps=fopen("ps.dat","w");
    for(int i=0;i<Parameters.nk;i++){
      if(inbin[i]>0) {
	//float kcurr = sqrt(1/ook2[i]);
	float kcurr = kmean[i];
	float kdoh = kcurr*sqrt(
				(ps_prpx[i]+ps_prpy[i]+ps_prpz[i])*
				pow(kcurr,3)/2/M_PI/M_PI
				)/HCurr;	  	  
	float delta2curr = pow(kcurr,3)*ps_delta[i]*pow(DCurr/DInit,2)/2/M_PI/M_PI;
	fprintf(fps,"%e %e %e\n",kcurr,kdoh,delta2curr);
      }
    }
    fclose(fps);
    
    if(iz>0) PowerSpectrum(par,ps_par,0);

    WritePS();

    // Set wavenumber to perpendicular ps table
    Wavenumber2PprpTable  = new double[Parameters.nk];
    Wavenumbers4PprpTable = new double[Parameters.nk];
    SetWavenumber2P3DTable(Parameters.nk, Wavenumber2PprpTable, Wavenumbers4PprpTable, 
                            ps_prpx, ps_prpy, ps_prpz, ook2);

    // Set wavenumber to parallel ps table
    Wavenumber2PparTable  = new double[Parameters.nk];
    Wavenumbers4PparTable = new double[Parameters.nk];
    SetWavenumber2P1DTable(Parameters.nk, Wavenumber2PparTable, Wavenumbers4PparTable, 
                            ps_par, ook2);

    // Accumulate l^2Cl/(2pi)
    if(iz==0) {
      AccumulateCl(0);
    } else{
      AccumulateCl(1);
    }

    // Delete wavenumber to prp and par table arrays
    delete[] Wavenumber2PprpTable;
    delete[] Wavenumbers4PprpTable;

    delete[] Wavenumber2PparTable;
    delete[] Wavenumbers4PparTable;

  }

  WriteCl();

  // Finalize and return

  if(myid==0) printf("\n\n");
  fclose(stdout);
  MPI_Finalize();    
  exit(0);

}

