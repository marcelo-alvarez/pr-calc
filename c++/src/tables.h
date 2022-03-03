#define NZTABLE 150000
#define ZTABLE_INITIAL 0
#define ZTABLE_FINAL 25

#define NRTABLE 150000
#define RTABLE_INITIAL 0
#define RTABLE_FINAL 1.3e4

#ifndef TABLES_H
#define TABLES_H

float Redshift2Float(float, double *);
float Radius2Float(float, double *);
float Wavenumber2Float(int, float, double *, double *);
float Wavenumber2FloatLogSpace(int, float, float, float, double *);
void SetRedshift2WKappaTable(float, float, float, double *, double *);
void SetRedshift2RadiusTable(float, float, float, double *);
void SetRadius2RedshiftTable(float, float, float, double *);
void SetRedshift2HistoryTable(int, int, float *, double *);
void SetRedshift2DgTable(double *, double *);
void SetRedshift2TauTable(float, float, float, double *, double *);
void SetWavenumber2P3DTable(int, double *, double *, double *, double *, double *, double *);
void SetWavenumber2P1DTable(int, double *, double *, double *, double *);

class Float2FloatTable
{
public: 
  float* table; 
  float  minval, maxval, delta;
  int    N;
  
  Float2FloatTable(char* fname);
  float Float2Float(float inval);

};

#endif