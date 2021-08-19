#include "math.h"
#include "geometry.h"

bool BoxOutsideOfShell(float *x, float *y, float *z, float r1, float r2) {

  float d0=pow(x[0],2)+pow(y[0],2)+pow(z[0],2);
  float d1=pow(x[1],2)+pow(y[1],2)+pow(z[1],2);
  float d2=pow(x[2],2)+pow(y[2],2)+pow(z[2],2);
  float d3=pow(x[3],2)+pow(y[3],2)+pow(z[3],2);
  float d4=pow(x[4],2)+pow(y[4],2)+pow(z[4],2);
  float d5=pow(x[5],2)+pow(y[5],2)+pow(z[5],2);
  float d6=pow(x[6],2)+pow(y[6],2)+pow(z[6],2);
  float d7=pow(x[7],2)+pow(y[7],2)+pow(z[7],2);
  
  r1*=r1;
  r2*=r2;
  
  return (
	  ((d0<r1)&&(d1<r1)&&(d2<r1)&&(d3<r1)&&(d4<r1)&&(d5<r1)&&(d6<r1)&&(d7<r1)) ||
	  ((d0>r2)&&(d1>r2)&&(d2>r2)&&(d3>r2)&&(d4>r2)&&(d5>r2)&&(d6>r2)&&(d7>r2))
	  );
  
}
