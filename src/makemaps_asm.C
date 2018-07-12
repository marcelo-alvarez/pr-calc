#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto_asm.h"
#include "allvars_asm.h"

void MakeMapsRayTrace()
{
  
  int i,j,k,index;
  
  float xm,ym,zm;    // position in map coordinates
  float x,y,z;       // position in cube coordinates
  
  double theta,phi;   // angle in cube coordinates
  float deg2rad=2.*M_PI/360.;

  // Local copies of parameters

  float BoxSize=Parameters.BoxSize;

  int N=Parameters.N;
  int NPixels=Parameters.NPixels;  

  float z0=Parameters.InitialRedshift;
  float z1=Parameters.FinalRedshift;
  int Nz=Parameters.NRedshifts;

  float dzl = z0/Nz;

  int Nxmin = Nlocal * myid;

  float observer_offset = 1.5e4;

  float Omegam = Parameters.Omegam;
  float Omegab = Parameters.Omegab;
  float Omegal = 1-Omegam;
  float      h = Parameters.h;
  float  zInit = Parameters.zInit;

  float Tcmb = 2.726e6; // Tcmb in muK
  float hoc = h * 1e2 / 3e5; // H0/c in units of 1/Mpc

  float thompson = 6.65e-25;
  float YHe = 0.25;
  float hubble0 = 100./3.086e19*h;
  float ne0=Omegab*h*h*1.88e-29/1.67e-24;
  float c=3e10;
  float nu0=1.4e3;

  float deltanu = clParameters.bandwidth;
  float nubin = clParameters.frequency;

  float dtau0=thompson*ne0*c/hubble0;

  float dkszdtau_ray, dtau_ray, dkap_ray;
  float ksz_pixel, tau_pixel, kap_pixel;
  float NHe;

  float los_epsilon = 0.5; // sets number of samples per cell along los

  // local copies of maps

  taumapl = (float *)malloc(mapsize*sizeof(float));
  kszmapl = (float *)malloc(mapsize*sizeof(float));
  dtbmapl = (float *)malloc(mapsize*sizeof(float));
  kapmapl = (float *)malloc(mapsize*sizeof(float));

  ReportMemory("before map projection",total_local_size,ovrt,oram);

  // Set Redshift to Radius Table
  
  Redshift2RadiusTable = new double[NZTABLE];
  SetRedshift2RadiusTable(h, Omegam, Omegal, Redshift2RadiusTable);

  // Set Redshift to Lensing Kernel Table

  Redshift2WKappaTable = new double[NZTABLE];
  SetRedshift2WKappaTable(h, Omegam, Omegal, Redshift2WKappaTable, Redshift2RadiusTable);

  // make redshift spacings such that dz = dz/dr * cell size

  float CellSize = BoxSize / N ;
  
  // first determine number of spacings in redshift
  int Nza = 0;
  float zcur = z0;
  float ra = Redshift2Float(zcur, Redshift2RadiusTable);	
  while(zcur < z1){
    float dza = hoc * sqrt(Omegam * pow(1+zcur,3) + 1 - Omegam) * CellSize * los_epsilon;
    zcur+=dza;
    Nza++;
  }

  float *redshiftsl       = new float[Nz]();
  float *growthfactorsl   = new float[Nz]();
  float *ionizedfractionl = new float[Nz]();

  float *redshifts        = new float[Nza]();
  float *meanhistory      = new float[Nza]();
  float *radii            = new float[Nza]();
  float *growthfactors    = new float[Nza]();
  float *ionizedfraction  = new float[Nza]();
  float *dkszdtau         = new float[Nza]();
  float *kszfac           = new float[Nza]();
  float *dz               = new float[Nza]();
  float *taufac           = new float[Nza]();
  float *kapfac           = new float[Nza]();

  redshifts[0]=z0;
  ionizedfraction[0]=1;
  zcur=z0;
  radii[0]=Redshift2Float(z0,Redshift2RadiusTable);

  for(k=1;k<Nza-1;k++){

    float hofz = hoc * sqrt(Omegam * pow(1+zcur,3) + 1 - Omegam);
    float dza = hofz * CellSize * los_epsilon ;
    zcur+=dza;

    // Precomputed quantities along rays

    redshifts[k] = zcur;
    dz[k] = dza;

    growthfactors[k] = growth(redshifts[k],Parameters.Omegam,
			      Parameters.Omegal, Parameters.w)/DInit;
    meanhistory[k]   = Redshift2Float(zcur,Redshift2HistoryTable);
    radii[k]         = Redshift2Float(zcur,Redshift2RadiusTable);
    kapfac[k]        = Redshift2Float(zcur,Redshift2WKappaTable);
    kapfac[k]       *= CellSize * los_epsilon ;

    // Note the negative sign here
    // dksz/dtau = - Tcmb * v = kszfac * vi
    // v = H(z) * D(z) / D(zi) / (1+z) * vi
    // kszfac = -Tcmb * v / vi
    // kszfac = -Tcmb * H(z) * D(z) / D(zi) / (1+z)
    
    kszfac[k] = - Tcmb * hofz * growthfactors[k] / (1+zcur) ; 

    NHe=1;
    if(zcur<3) NHe=2;
    taufac[k] = dtau0*pow((1+zcur),2) / sqrt(Omegam*pow((1+zcur),3)+1-Omegam) * 
                (1-YHe+NHe*YHe/4)*dz[k];
  }

  // Identify half-ionized redshift
  
  k=1;
  while(meanhistory[k++]>0.5);
  float zhalf = redshifts[k];
  float rhalf = Redshift2Float(zhalf, Redshift2RadiusTable);	
  float Dhalf = growth(zhalf,Parameters.Omegam,Parameters.Omegal, Parameters.w);
  float Tb0   = 2.3e4*(Omegab*h*h/0.02)*sqrt((0.15/Omegam/h/h)*(1+zhalf)/10.); // Tb in muK
  float nuhalf = nu0 / (1+zhalf);

  for(k=0;k<Nz;k++){
    redshiftsl[k] = (k-0.5)*dzl;
    growthfactorsl[k]=growth(redshiftsl[k],Parameters.Omegam,Parameters.Omegal,
    	       Parameters.w);
    ionizedfractionl[k] = 1;
  }

  float taul = tau(redshiftsl,ionizedfractionl,Nz,Omegam,Omegab,h);

  MPI_Barrier(MPI_COMM_WORLD);
  double t1 = MPI_Wtime();

  float numn = nubin - deltanu / 2;
  float numx = nubin + deltanu / 2;
  float zmn  = nu0 / numn - 1;
  float zmx  = nu0 / numx - 1;
  float xmn  = Redshift2Float(zmn,Redshift2HistoryTable);
  float xmx  = Redshift2Float(zmx,Redshift2HistoryTable);

  if(myid==0) printf("\n universe half ionized at: z = %f, nu = %f MHz\n",zhalf,nuhalf);
  if(myid==0) printf(" ionized fraction at numin = %f, zmin = %f: %f\n",numn,zmn,xmn);
  if(myid==0) printf(" ionized fraction at numax = %f, zmax = %f: %f\n",numx,zmx,xmx);

  // All processes loop over all lines of sight, but only the segments of lines of 
  // sight that overlap the domain belonging to the process are accumulated.
  // At the end, the maps from all the processes are added to obtain the final
  // map

  if(clParameters.verbose==1 && myid==0) printf("\n");
  i=0;
  for(int ii=0; ii<Parameters.NSide; ii++){
    if(myid==0 && (ii+1)%512==0) printf("\n side %d out of %d done",ii,Parameters.NSide);
  for(int jj=0; jj<Parameters.NSide; jj++){
  for(int kk=0; kk<12; kk++){
    
    pix2ang_nest(Parameters.NSide, (long)i, &theta, &phi);
    
    // Cartesian unit vector

    x = sin(theta)*cos(phi);
    y = sin(theta)*sin(phi);
    z = cos(theta);

    // Determine the number of points along line of sight to skip the ray leaves the slab

    int nskip = (N-Nlocal)/fabs(x)-3; // subtract 3 for roundoff errors to be on safe side
    if(x==0 || nskip < 0 ) nskip=0;    
    nskip = 0;

    // Loop over all points along the line of sight for kSZ

    float Tb = 0; // going to do l.o.s. average for 21-cm
    int InSlab = 0;
    tau_pixel=0;
    ksz_pixel=0;
    kap_pixel=0;
    for(k=1; k<Nza-1; k++){
      
      float redshift = redshifts[k];
      float r = radii[k];
      
      // The point (xr,yr,zr) represents the comoving position along the line
      // of sight with respect to the observer, obtained by multiplying the
      // unit vector, (x,y,z), with the radius, r.
      
      float xr = x*r + observer_offset;
      int Nx = xr / BoxSize; int xc = (xr - BoxSize*Nx)/BoxSize * N;
      int xclocal = xc - Nxmin;
      
      if(xclocal>=0 && xclocal<Nlocal){

	float yr = y*r + observer_offset;
	int Ny = yr / BoxSize; int yc = (yr - BoxSize*Ny)/BoxSize * N;
	
	float zr = z*r + observer_offset;
	int Nz = zr / BoxSize; int zc = (zr - BoxSize*Nz)/BoxSize * N;
	
	InSlab = 1;
	
	int index_zr = xclocal*N*N     + yc*N     + zc;
	int index_dv = xclocal*N*(N+2) + yc*(N+2) + zc;
	
	float DeltaInit = delta[index_dv];
	float xcur = meanhistory[k];
	float dcur = growthfactors[k];

	if(DeltaInit>900) {
	  DeltaInit -= 1000;	
	  if(DeltaInit*dcur > tau_pixel &&
	     abs(xcur-0.5)<0.05) tau_pixel = DeltaInit*dcur;
	}

	float Delta = DeltaInit * growthfactors[k];

	// 21-cm average over dx21
	float nu  = nu0 / ( 1 + redshift ) ;
	float dnu = abs(nubin - nu);

	if(zreion[index_zr]>redshift)
	  {
	    ionizedfraction[k] = 1;
	  }
	else
	  {
	    ionizedfraction[k] = 0;
            if(dnu < deltanu) Tb += 1+Delta;
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
	
	dkszdtau_ray = kszfac[k] * 
	  ( x*xdv_x[index_dv] + y*xdv_y[index_dv] + z*xdv_z[index_dv] ) ;
	dtau_ray = (1+Delta) * taufac[k] * ionizedfraction[k];
	dkap_ray =    Delta  * kapfac[k];

	tau_pixel += dtau_ray ; 
	ksz_pixel += dkszdtau_ray*dtau_ray ; 	
	kap_pixel += dkap_ray ;
	
      } else {
	
	// Here we will skip to the next slab,
	// unless we were not previously in the slab
	
	//	if(InSlab == 1) k+=nskip;
	InSlab = 0;
	
      }
      
    }

    int mapindex = i*NPixels + j;
    
    taumapl[i] = tau_pixel;
    kszmapl[i] = ksz_pixel;
    kapmapl[i] = kap_pixel;
    dtbmapl[i] = Tb;

    i++;
    
  }
  }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  double dt = MPI_Wtime() - t1;

  if(myid==0) printf(" Projection took %le seconds\n",dt);

  // sum up all maps
  
  MPI_Allreduce(taumapl, taumap, mapsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(kszmapl, kszmap, mapsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(dtbmapl, dtbmap, mapsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(kapmapl, kapmap, mapsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

  double taumean=0;
  double kszmean=0;
  double dtbmean=0;
  double kapmean=0;

  for(i=0;i<mapsize;i++){
    taumap[i]+=taul;
    taumean+=taumap[i];
    kszmean+=kszmap[i];
    dtbmean+=dtbmap[i];
    kapmean+=kapmap[i];
  }

  taumean/=mapsize;
  kszmean/=mapsize;
  dtbmean/=mapsize;
  kapmean/=mapsize;

  double taurms=0;
  double kszrms=0;
  double dtbrms=0;
  double kaprms=0;

  for(i=0;i<mapsize;i++){
    taurms+=pow((taumean-taumap[i]),2);
    kszrms+=pow((kszmean-kszmap[i]),2);
    dtbrms+=pow((dtbmean-dtbmap[i]),2);
    kaprms+=pow((kapmean-kapmap[i]),2);
  }

  taurms/=mapsize;
  kszrms/=mapsize;
  dtbrms/=mapsize;
  kaprms/=mapsize;

  taurms=sqrt(taurms);
  kszrms=sqrt(kszrms);
  dtbrms=sqrt(dtbrms);
  kaprms=sqrt(kaprms);

  if(myid==0) printf("\n tau mean = %e",taumean);
  if(myid==0) printf("\n tau rms  = %e",taurms);

  if(myid==0) printf("\n ksz mean = %e",kszmean);
  if(myid==0) printf("\n ksz rms  = %e",kszrms);

  if(myid==0) printf("\n 21-cm mean = %e",dtbmean);
  if(myid==0) printf("\n 21-cm rms  = %e",dtbrms);

  if(myid==0) printf("\n kappa mean = %e",kapmean);
  if(myid==0) printf("\n kappa rms  = %e",kaprms);

}


void MakeMapsBinCells()
{
  
  int i,j,k,index;
  
  float xm,ym,zm;    // position in map coordinates
  float x,y,z;       // position in cube coordinates
  
  double theta,phi;   // angle in cube coordinates
  float deg2rad=2.*M_PI/360.;

  // Local copies of parameters

  float BoxSize=Parameters.BoxSize;

  int N=Parameters.N;

  float zmin=Parameters.InitialRedshift;
  float zmax=Parameters.FinalRedshift;
  int Nz=Parameters.NRedshifts;

  float dzl = zmin/Nz;

  int Nxmin = Nlocal * myid;
  float slabsize = BoxSize / nproc;
  float offset   = slabsize * myid;

  float Omegam = Parameters.Omegam;
  float Omegab = Parameters.Omegab;
  float Omegal = 1-Omegam;
  float      h = Parameters.h;
  float  zInit = Parameters.zInit;

  float Tcmb = 2.726e6; // Tcmb in muK
  float hoc = h * 1e2 / 3e5; // H0/c in units of 1/Mpc

  float thompson = 6.65e-25;
  float YHe = 0.25;
  float hubble0 = 100./3.086e19*h;
  float ne0=Omegab*h*h*1.88e-29/1.67e-24;
  float c=3e10;

  float dtau0=thompson*ne0*c/hubble0;

  float NHe;

  // local copies of maps

  float *taumapl = new float[mapsize];
  float *kszmapl = new float[mapsize];
  float *dtbmapl = new float[mapsize];

  ReportMemory("before map projection",total_local_size,ovrt,oram);

  // Set Redshift to Radius Table
  
  Redshift2RadiusTable = new double[NZTABLE];
  SetRedshift2RadiusTable(h, Omegam, Omegal, Redshift2RadiusTable);

  float rmin=Redshift2Float(zmin,Redshift2RadiusTable);
  float rmax=Redshift2Float(zmax,Redshift2RadiusTable);

  float CellSize = BoxSize / N ;

  // Set Radius to Redshift Table

  /*  
  Radius2RedshiftTable = new double[NZTABLE];
  SetRadius2RedshiftTable(h, Omegam, Omegal, Radius2RedshiftTable);
  
  k=1;
  while(meanhistory[k++]>0.5);
  float zhalf = redshifts[k];
  float rhalf = Redshift2Float(zhalf, Redshift2RadiusTable);	
  float Dhalf = growth(zhalf,Parameters.Omegam,Parameters.Omegal, Parameters.w);
  float Tb0   = 2.3e4*(Omegab*h*h/0.02)*sqrt((0.15/Omegam/h/h)*(1+zhalf)/10.); // Tb in muK
  */

  float *redshifts       = new float[Nz]();
  float *ionizedfraction = new float[Nz]();

  for(k=0;k<Nz;k++){
    redshifts[k] = (k-0.5)*dzl;
    ionizedfraction[k] = 1;
  }

  float taul = tau(redshifts,ionizedfraction,Nz,Omegam,Omegab,h);

  // before looping over periodic slabs, find out maximum number of images in each dimension
  // use 15 Gpc as largest possible radius

  int nperiodic = (int)(1.5e4/BoxSize)+1;

  // slab corners

  float *xs = new float[8];
  float *ys = new float[8];
  float *zs = new float[8];

  MPI_Barrier(MPI_COMM_WORLD);
  double t1 = MPI_Wtime();

  // Each process loops over all periodic slab images overlapping with integration region
  // At the end, the maps from all the processes are added to obtain the final map

  for(int ip=-nperiodic;ip<=nperiodic;ip++){
  for(int jp=-nperiodic;jp<=nperiodic;jp++){
  for(int kp=-nperiodic;kp<=nperiodic;kp++){

    // set corners of slab

    xs[0]=xs[1]=xs[2]=xs[3]=offset+ip*BoxSize;
    xs[4]=xs[5]=xs[6]=xs[7]=offset+ip*BoxSize+offset;

    ys[0]=ys[2]=ys[4]=ys[6]=jp*BoxSize;
    ys[1]=ys[3]=ys[5]=ys[7]=(jp+1)*BoxSize;

    zs[0]=zs[1]=zs[4]=zs[5]=kp*BoxSize;
    zs[2]=zs[3]=zs[6]=zs[7]=(kp+1)*BoxSize;

    if(BoxOutsideOfShell(xs,ys,zs,rmin,rmax)) break;

    // now that we know the box intersects the shell loop over all cells

    for(int ic=0;ic<Nlocal;ic++){
      float x = xs[0] + (i-0.5)*CellSize;
    for(int jc=0;jc<N;jc++){
      float y = ys[0] + (j-0.5)*CellSize;
    for(int kc=0;kc<N;kc++){
      float z = zs[0] + (k-0.5)*CellSize;
      
      float r = sqrt(x*x + y*y + z*z);      
      float redshift = Radius2Float(r,Radius2RedshiftTable);      


    }
    }
    }    

  }
  }
  }    

  MPI_Barrier(MPI_COMM_WORLD);
  double dt = MPI_Wtime() - t1;

  if(myid==0) printf("\n Projection took %le seconds\n",dt);

  // sum up all maps
  
  MPI_Allreduce(taumapl, taumap, mapsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(kszmapl, kszmap, mapsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(dtbmapl, dtbmap, mapsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

  // find mean

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

  // find rms

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

  // report results

  if(myid==0) printf("\n\n tau mean = %e",taumean);
  if(myid==0) printf("\n\n tau low  = %e",taul);
  if(myid==0) printf("\n tau rms  = %e",taurms);

  if(myid==0) printf("\n ksz mean = %e",kszmean);
  if(myid==0) printf("\n ksz rms  = %e",kszrms);

  if(myid==0) printf("\n 21-cm mean = %e",dtbmean);
  if(myid==0) printf("\n 21-cm rms  = %e",dtbrms);

}


