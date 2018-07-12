#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto_map.h"
#include "allvars_map.h"

void MakeMaps()
{

  int i,j,k,index;

  float xm,ym,zm;    // position in map coordinates
  float x,y,z;       // position in cube coordinates

  float theta,phi;   // angle in cube coordinates
  float deg2rad=2.*M_PI/360.;

  float observer_offset = 1.5e4; // position of observer is (offset,offset,offset)

  // Local copies of parameters

  float BoxSize=Parameters.BoxSize;

  int N=Parameters.N;
  int NPixels=Parameters.NPixels;  

  float theta0=Parameters.theta*deg2rad;
  float phi0=Parameters.phi*deg2rad;
  float fov=Parameters.fov*deg2rad;

  float theta1=M_PI/2.-theta0;
  float dtheta=fov/NPixels;

  float z0=Parameters.InitialRedshift;
  float z1=Parameters.FinalRedshift;
  int Nz=Parameters.NRedshifts;

  float dz = (z1-z0)/Nz;
  float dzl = z0/Nz;

  int Nxmin = Nlocal * myid;

  float Omegam = Parameters.Omegam;
  float Omegab = Parameters.Omegab;
  float Omegal = 1-Omegam;
  float      h = Parameters.h;
  float  zInit = Parameters.zInit;

  float Tcmb = 2.726e6; // Tcmb in muK
  float hoc = h * 1e2 / 3e5; // H0/c in units of 1/Mpc
  float Tb0   = 7.3e3*(Omegab*h*h/0.02)*sqrt((0.15/Omegam/h/h)); // Tb in muK
  float dzdchi0 = pow(Omegam,0.5) * hoc;
  
  // local copies of maps

  float *taumapl = new float[NPixels*NPixels];
  float *kszmapl = new float[NPixels*NPixels];
  float *dtbmapl = new float[NPixels*NPixels];

  // Set Redshift to Radius Table
  
  Redshift2RadiusTable = new double[NZTABLE];
  SetRedshift2RadiusTable(h, Omegam, Omegal, Redshift2RadiusTable);

  // make redshift spacings such that dz = dz/dr * cell size

  float CellSize = BoxSize / N ;
  
  // first determine number of spacings in redshift
  int Nza = 0;
  float zcur = z0;
  float ra = Redshift2Float(zcur, Redshift2RadiusTable);	
  while(zcur < z1){
    float dza = hoc * sqrt(Omegam * pow(1+zcur,3) + 1 - Omegam) * CellSize;
    zcur+=dza;
    Nza++;
  }

  if(myid==0) printf("\n Nza = %d r0 = %e rf = %e\n",Nza,ra,Redshift2Float(zcur, Redshift2RadiusTable));

  float *redshiftsl       = new float[Nz]();
  float *growthfactorsl   = new float[Nz]();
  float *ionizedfractionl = new float[Nz]();
  float *meanhistory      = new float[Nz]();

  float *redshifts        = new float[Nza]();
  float *growthfactors    = new float[Nza]();
  float *ionizedfraction  = new float[Nza]();
  float *dkszdtau         = new float[Nza]();
  float *dtbcur           = new float[Nza]();

  redshifts[0]=z0;
  ionizedfraction[0]=1;
  zcur=z0;
  for(k=1;k<Nza;k++){
    float dza = hoc * sqrt(Omegam * pow(1+zcur,3) + 1 - Omegam) * CellSize;
    zcur+=dza;
    redshifts[k] = zcur;
    growthfactors[k]=growth(redshifts[k],Parameters.Omegam,Parameters.Omegal, Parameters.w);
    meanhistory[k] = Redshift2Float(zcur,Redshift2HistoryTable);
  }
  
  for(k=0;k<Nz;k++){
    redshiftsl[k] = (k-0.5)*dzl;
    growthfactorsl[k]=growth(redshiftsl[k],Parameters.Omegam,Parameters.Omegal,
    	       Parameters.w);
    ionizedfractionl[k] = 0;
    if(redshiftsl[k] < 5.5) ionizedfractionl[k] = 1;
  }

  float taul = tau(redshiftsl,ionizedfractionl,Nz,Omegam,Omegab,h);

  // All processes loop over all lines of sight, but only the segments of lines of 
  // sight that overlap the domain belonging to the process are accumulated.
  // At the end, the maps from all the processes are added to obtain the final
  // map

  if(clParameters.verbose==1 && myid==0) printf("\n");
  for(i=0; i<NPixels; i++){

    if(i%10==0 && clParameters.verbose==1 && myid==0) printf("\n Projecting %d of %d rows",i,NPixels);

    for(j=0; j<NPixels; j++){

      theta = M_PI/2. + (i-0.5)*dtheta-fov/2.;
      phi   = (j-0.5)*dtheta-fov/2.;

      // The pixel initially oriented with respect to the x-axis

      xm = sin(theta)*cos(phi);
      ym = sin(theta)*sin(phi);
      zm = cos(theta);

      // Rotate in the x-z plane counterclockwise by theta1=pi/2-theta0

      float xp = xm*cos(theta1) - zm*sin(theta1);
      float yp = ym;
      float zp = xm*sin(theta1) + zm*cos(theta1);

      // Now rotate in the x-y plane counterclockwise by phi0
      // so that x,y,z represents the unit vector of the pixel

      x = xp*cos(phi0) - yp*sin(phi0);
      y = xp*sin(phi0) + yp*cos(phi0);
      z = zp;

      // Loop over all redshifts along the line of sight for projected maps
      
      for(k=0; k<Nza; k++){
	ionizedfraction[k]=0;
	dkszdtau[k]=0;
	dtbcur[k]=0;
      }
      for(k=0; k<Nza; k++){
	
	float redshift = redshifts[k];
	float DCurr = growthfactors[k];
	float r = Redshift2Float(redshift, Redshift2RadiusTable);	
	
	// Note the negative sign here
	// dksz/dtau = - Tcmb * v / c = kszfac * vi
	// Tcmb * v / c = Tcmb * H(z)/(1+z)*D(z)/D(zi) * vi
	// kszfac = - Tcmb * v / c / vi = - Tcmb * H(z)/(1+z)*D(z)/D(zi)

	float hofz = hoc * sqrt(Omegam * pow(1+redshift,3) + 1 - Omegam);
	float kszfac = - Tcmb * hofz / (1+redshift) * DCurr / DInit;
	               
	// The point (xr,yr,zr) represents the comoving position along the line
	// of sight with respect to the observer, obtained by multiplying the
	// unit vector, (x,y,z), with the radius, r.

	float xr = x*r + observer_offset ;
	float yr = y*r + observer_offset ;
	float zr = z*r + observer_offset ;

	if(xr<0) xr += (int(abs(xr/BoxSize))+1)*BoxSize;
	if(yr<0) yr += (int(abs(yr/BoxSize))+1)*BoxSize;
	if(zr<0) zr += (int(abs(zr/BoxSize))+1)*BoxSize;
	  
	int Nx = xr / BoxSize; int xc = (xr - BoxSize*Nx)/BoxSize * N;
	int Ny = yr / BoxSize; int yc = (yr - BoxSize*Ny)/BoxSize * N;
	int Nz = zr / BoxSize; int zc = (zr - BoxSize*Nz)/BoxSize * N;
	
	int xclocal = xc - Nxmin;

	if(xclocal>=0 && xclocal<Nlocal){

	  int index_zr = xclocal*N*N     + yc*N     + zc;
	  int index_dv = xclocal*N*(N+2) + yc*(N+2) + zc;

	  float Delta = delta[index_dv] * growthfactors[k];

	  if(zreion[index_zr]>redshift)
	    {
	      ionizedfraction[k] = 1;
	    }
	  else
	    {
	      ionizedfraction[k] = 0;
	    }
	  
	  if(clParameters.uniformu==1){
	    Delta = 0;
	    ionizedfraction[k] = meanhistory[k];
	  }
	  else if(clParameters.uniformx==1){
	    ionizedfraction[k] = meanhistory[k];
	  }
	  else if(clParameters.uniformd==1){
	    Delta = 0;
	  }
	  	  
	  dkszdtau[k] = kszfac * ionizedfraction[k] * (1 + Delta) *
	    ( x*xdv_x[index_dv] + y*xdv_y[index_dv] + z*xdv_z[index_dv] ) ;

	  dtbcur[k] = (1-ionizedfraction[k]) * (1 + Delta) ;

	  ionizedfraction[k] *= (1+Delta);
	  
	}
      }
      
      int mapindex = i*NPixels + j;

      taumapl[mapindex] = tau(redshifts,ionizedfraction,Nza,Omegam,Omegab,h);
      kszmapl[mapindex] = ksz(redshifts,       dkszdtau,Nza,Omegam,Omegab,h);
      dtbmapl[mapindex] = dtb(redshifts,         dtbcur,Nza,Omegam,Omegab,h);

    }
  }

  // sum up all maps

  MPI_Allreduce(taumapl, taumap, mapsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(kszmapl, kszmap, mapsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(dtbmapl, dtbmap, mapsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

  float taumean=0;
  float kszmean=0;
  float dtbmean=0;

  for(i=0;i<mapsize;i++){
    taumap[i]+=taul;
    taumean+=taumap[i];
    kszmean+=kszmap[i];
    dtbmean+=dtbmap[i];
  }

  taumean/=mapsize;
  kszmean/=mapsize;
  dtbmean/=mapsize;

  float taurms=0;
  float kszrms=0;
  float dtbrms=0;

  for(i=0;i<mapsize;i++){
    taurms+=pow((taumean-taumap[i]),2);
    kszrms+=pow((kszmean-kszmap[i]),2);
    dtbrms+=pow((dtbmean-dtbmap[i]),2);
  }

  taurms/=mapsize;
  kszrms/=mapsize;
  dtbrms/=mapsize;

  taurms=sqrt(taurms);
  kszrms=sqrt(kszrms);
  dtbrms=sqrt(dtbrms);

  if(myid==0) printf("\n\n tau mean = %f",taumean);
  if(myid==0) printf("\n tau rms  = %f",taurms);

  if(myid==0) printf("\n ksz mean = %f",kszmean);
  if(myid==0) printf("\n ksz rms  = %f",kszrms);

  if(myid==0) printf("\n 21-cm mean = %f",dtbmean);
  if(myid==0) printf("\n 21-cm rms  = %f",dtbrms);

}
