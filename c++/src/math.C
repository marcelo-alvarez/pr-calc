#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

/////////////////////////////////////////////////////////
// My very simple math stuff
/////////////////////////////////////////////////////////

struct myparams
{
  double a;
};

double myerfc(double x, void *params)
{
  struct myparams *p = (struct myparams *) params;
  double a = p->a;
  return erfc(x)-a;
}
		
float inverse_erfc(float x)
{

  int status;
  int iter = 0, max_iter = 100;

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;

  double x_lo = 0.0, x_hi = 100.0;

  gsl_function MYERFC;
  struct myparams params = {x};
     
  MYERFC.function = &myerfc;
  MYERFC.params   = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &MYERFC, x_lo, x_hi);
     
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
				       0, 0.0001);      
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  
  gsl_root_fsolver_free (s);
  
  if(status != GSL_SUCCESS) printf("WARNING, inverse error function didn't work!!!\n");
    
  return x_lo;
}


