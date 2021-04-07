
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

// Uses gsl-integ for computing the integral of ln(x)/sqrt(x) from x = 0 to x = 1.

double f (double x, void * params) {
  double alpha = *(double *) params;
  double f = log(alpha*x) / sqrt(x);
  return f;
}

double err(double x, void* params){
	double alpha = *(double *) params;
	double erf = 2/sqrt(alpha*M_PI)* exp(-pow(x,2));
	return erf;
}


int
main (void)
{
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

  double result, error;
  double expected = -4.0;
  double alpha = 1.0;

  gsl_function F;
// f
  F.function = &f;
  F.params = &alpha;

  gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
                        w, &result, &error);

  printf ("result          = % .18f\n", result);
  printf ("exact result    = % .18f\n", expected);
  printf ("estimated error = % .18f\n", error);
  printf ("actual error    = % .18f\n", result - expected);
  printf ("intervals       = %zu\n", w->size);

  gsl_integration_workspace_free (w);


 return 0;
}
