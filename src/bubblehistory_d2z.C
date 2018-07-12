#include "allvars_d2z.h"
#include "proto_d2z.h"
#include <math.h>

float GetHistory(){

  history   = new float[NHISTORY];
  historyl  = new float[NHISTORY];
  history_z = new float[NHISTORY];
  
  dzhistory=(HISTORY_FINAL-HISTORY_INITIAL)/(NHISTORY-1);

  int N=Parameters.N;

  for(int i=0;i<NHISTORY;i++) historyl[i]=0;

  for(int i=0;i<Nlocal;i++){
    for(int j=0;j<N;j++){
      for(int k=0;k<N;k++){
	
	int index=i*(N+2)*N+j*(N+2)+k;
	float zcur=zreion[index];
	int bin=(int)((zcur-HISTORY_INITIAL)/dzhistory);
	if(bin>=0 && bin<NHISTORY) historyl[bin]++;
	
      }
    }
  }

  MPI_Allreduce(historyl, history, NHISTORY, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);  

  for(int i=0;i<NHISTORY;i++) history[i]/=pow(N,3);
  for(int i=1;i<NHISTORY;i++) history[i]+=history[i-1];
  for(int i=1;i<NHISTORY;i++) history_z[i]=HISTORY_INITIAL+i*dzhistory;

  float tau_local;
  if(myid==0) tau_local = tau(history_z, history, NHISTORY,
			Parameters.Omegam, Parameters.Omegab,
			Parameters.h);

  return tau_local;
  
}

void GetCHistory(){

  chistory   = new float[NCHISTORY];
  chistory_z = new float[NCHISTORY];
  
  dzchistory=(CHISTORY_FINAL-CHISTORY_INITIAL)/(NCHISTORY-1);

  float rmfp   = clParameters.Rmax;
  float zeta   = clParameters.zeta;

  float zInit  = Parameters.zInit;
  float omegam = Parameters.Omegam;
  float omegal = Parameters.Omegal;
  float w      = Parameters.w;
  float smin   = sigma_min;

  chistory[0]=0.;
  for(int i=1;i<NCHISTORY;i++){

    float x1 = chistory[i-1];
    
    float z1 = CHISTORY_INITIAL+    i*dzchistory;
    float z2 = CHISTORY_INITIAL+(i+1)*dzchistory;

    float f1 = fcoll(z1, zInit, omegam, omegal, w, smin);
    float f2 = fcoll(z2, zInit, omegam, omegal, w, smin);
    float df = f2-f1;

    float rb = BubbleMFP(x1);

    chistory[i] = chistory[i-1] + df/(1+rb/rmfp)*zeta;

    if(chistory[i]>=1) chistory[i]=1;

  }
  
  return;
  
}

void GetBubbleMFP(){
  
  int i,j,k,m,imfp,istart,jstart,kstart,index,ionized,neutral, bin;
  
  float zcur, fint, mean_distance, current_distance,dmfp;
  float skewer_startx, skewer_starty, skewer_startz;
  float zreioncur,zlast,frac,fpow,clow,chigh;

  mfp     = new float[NX_MFP];
  mfp_x   = new float[NX_MFP];
  dfdmfp  = new float[NMFP*NX_MFP]();
  dfdmfpl = new float[NMFP*NX_MFP](); 

  int         N  = Parameters.N;
  float BoxSize  = Parameters.BoxSize;
  float CellSize = BoxSize / N;

  float MFP_FINAL   = log(MFP_FINAL_BOX    *  BoxSize);
  float MFP_INITIAL = log(MFP_INITIAL_CELL * CellSize);

  dmfp  = (  MFP_FINAL -   MFP_INITIAL) / (    NMFP);
  dxmfp = (X_MFP_FINAL - X_MFP_INITIAL) / (NX_MFP-1);
  
  // Loop over all ionized fractions
  for(int ix=0;ix<NX_MFP;ix++){

    float xcur=X_MFP_INITIAL+ix*dxmfp;

    fpow = 1.5 + 1.5*(xcur-0.5);

    // First find the redshift corresponding to the current ionized fraction
    if(xcur<history[0]){
      zcur=HISTORY_INITIAL;
    }else if(xcur>history[NHISTORY-1]){
      zcur=HISTORY_FINAL;
    } else{
      for(int i=0;i<NHISTORY-1;i++)
	if(xcur>history[i] && xcur<history[i+1]){
	  fint=(xcur-history[i])/(history[i+1]-history[i]);
	  zcur=(1-fint)*(HISTORY_INITIAL+i*dzhistory)+fint*(HISTORY_INITIAL+(i+1)*dzhistory);
	  break;
	}
    }
      
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Shoot skewers through the redshift field
    float mean_distancel=0;
    int nskewers=NSKEWERS/nproc;
    for(int iskewer=0;iskewer<nskewers;iskewer++){
      
      // Pick random cells until an ionized one is found
      neutral=1;
      while(neutral){
	skewer_startx = ((float)rand()/RAND_MAX);
	skewer_starty = ((float)rand()/RAND_MAX);
	skewer_startz = ((float)rand()/RAND_MAX);

	istart = (int)(skewer_startx*Nlocal);
	jstart = (int)(skewer_starty*N);
	kstart = (int)(skewer_startz*N);

	if(istart>=Nlocal) istart -= Nlocal;
	if(jstart>=N)      jstart -= N;
	if(kstart>=N)      kstart -= N;

	if(istart<0) istart += Nlocal;
	if(jstart<0) jstart += N;
	if(kstart<0) kstart += N;

	index = istart*N*(N+2)+jstart*(N+2)+kstart;

	if(index<0 || index >= size_fftw) printf("Error: index = %d\n",index);

	if(zreion[index]>zcur && jstart < N-2 && kstart < N-2) neutral=0;
      }
 
      // Move along y-direction until neutral and record distance travelled
      i = istart;
      j = jstart;
      k = kstart;
      
      ionized=1;
      m = 0;
      index = i*N*(N+2)+j*(N+2)+k;
      zlast = zreion[index];
      while(ionized){
	j++;
	if(j>=N) j-=N;
	index = i*N*(N+2)+j*(N+2)+k;
	zreioncur = zreion[index];
	if(zreioncur<zcur) {
	  ionized=0;
	  frac = (zcur - zlast) / (zreioncur - zlast);	  
	}
	if(j==jstart) {
	  ionized = 0;
	  frac = 0;
	}
	zlast = zreioncur;
	m++;
      }
      
      frac = (float)rand()/RAND_MAX ;
      clow = log(m-0.5);
      if(m>1) clow = log(m-1.);
      chigh = log(m*1.0);
      current_distance = exp(clow + (chigh-clow)*frac) * CellSize;
      mean_distancel += current_distance ;

      imfp = (int)((log(current_distance)-MFP_INITIAL)/dmfp);
      bin = ix * NMFP + imfp;
      if(imfp>=0 && imfp<NMFP && log(current_distance) < MFP_FINAL &&
	 log(current_distance) > MFP_INITIAL) dfdmfpl[bin]++;
	
      // Move along z-direction until neutral and record distance travelled
      i = istart;
      j = jstart;
      k = kstart;

      index = i*N*(N+2)+j*(N+2)+k;
      ionized=1;
      m = 0;
      zlast = zreion[index];
      while(ionized){
	k++;
	if(k>=N) k-=N;
	index = i*N*(N+2)+j*(N+2)+k;
	zreioncur = zreion[index];
	if(zreioncur<zcur) {
	  ionized=0;
	  frac = (zcur - zlast) / (zreioncur - zlast);
	}
	if(k==kstart) {
	  ionized = 0;
	  frac = 0;
	}
	zlast = zreioncur;
	m++;
      }
      
      frac = (float)rand()/RAND_MAX ;
      clow = log(m-0.5);
      if(m>1) clow = log(m-1.);
      chigh = log(m*1.0);
      current_distance = exp(clow + (chigh-clow)*frac) * CellSize;
      mean_distancel += current_distance;      
      
      imfp = (int)((log(current_distance)-MFP_INITIAL)/dmfp);
      bin = ix * NMFP + imfp;
      if(imfp>=0 && imfp<NMFP && log(current_distance) < MFP_FINAL &&
	 log(current_distance) > MFP_INITIAL) dfdmfpl[bin]++;
      
      }
      
    MPI_Allreduce(&mean_distancel, &mean_distance, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    mean_distance/=(2.*nskewers*nproc);

    mfp[ix] = mean_distance;

  }

  for(int i=1;i<NX_MFP;i++) mfp_x[i]=X_MFP_INITIAL+i*dxmfp;

  MPI_Allreduce(dfdmfpl, dfdmfp, NX_MFP*NMFP, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

  WriteMFP();
  WriteMFPPDF();

}

float BubbleMFP(float x){
  
  float f;
  int itab;

  float rb;

  itab = (int)((x-X_MFP_INITIAL)/dxmfp);
  f    = (x-itab*dxmfp-X_MFP_INITIAL)/dxmfp;

  if(itab<0){
    rb = mfp[0];
    return rb;
  } else if(itab>=NX_MFP-1){
    rb = mfp[NX_MFP-1];
    return rb;
  }

  rb = (1-f)*mfp[itab]+f*mfp[itab+1];

  return rb;

}
 
 float zreioncorrect(float x){

   float f, z;

  if(x<chistory[0]){
    z=CHISTORY_INITIAL;
    return z;
  } else if(x>chistory[NCHISTORY-1]){
    z=CHISTORY_FINAL;
    return z;
  }

  for(int i=0;i<NCHISTORY-1;i++){

    if(x>=chistory[i] && x<=chistory[i+1]){
      f = (x-chistory[i])/(chistory[i+1]-chistory[i]);
      z = (1-f)*(CHISTORY_INITIAL+    i*dzchistory)+
	      f*(CHISTORY_INITIAL+(i+1)*dzchistory);
      return z;
    }

  }

  if(myid==0) printf("did not find correct zreion value, exiting ...\n");
  MPI_Finalize();
  exit(0);

}

void CorrectHistory(){

  float x,z;
  int N;

  N=Parameters.N;

  float *zvz = new float[NHISTORY];

  for(int i=0;i<NHISTORY;i++) zvz[i]=-1;

  for(int i=1;i<NHISTORY;i++){
    x = history[i];
    if(x!=0){
      z = zreioncorrect(x);
      zvz[i] = z;
    }
  }

  for(int i=NHISTORY-2;i>=0;i--) if(zvz[i]<0) zvz[i]=zvz[i+1];

  for(int i=0;i<Nlocal;i++){
    for(int j=0;j<N;j++){
      for(int k=0;k<N;k++){

	int index = i*N*(N+2)+j*(N+2)+k;
	z = zreion[index];

	int itab = (int)((z-HISTORY_INITIAL)/dzhistory);

	zreion[index] = zvz[itab];

      }
    }
  }

}





