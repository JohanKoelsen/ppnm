#include <stdio.h>
#include <math.h>
#include <complex.h>

complex plainmc(int dim, double f(double* x), double* a, double* b, int N);

double f(double * x){
	return 1./(M_PI*M_PI*M_PI) * 1./(1. - cos(x[0]) * cos(x[1]) * cos(x[2]));
}

double g(double* x){
	return 1./sqrt(x[0]);
}

int main(){
	//Defining integration boarders and number of iterations, N.
	double a[] = {0,0,0};
	double b[] = {M_PI, M_PI, M_PI};
	int N = 1e7;
	//Integrating f over dx, dy, dz from 0 to pi

	complex result = plainmc(3,f,a,b,N);

	//Printing results and compare them with analytical expression found in the homework page.
	printf ("result = %g +/- %g*i \n", creal(result), cimag(result));
	printf ("analytical result = %g \n", pow(tgamma(1./4), 4)/(4*M_PI*M_PI*M_PI));

	//Trying f = 1/sqrt(x) for 0 to pi
	double a_test[] = {0};
	double b_test[] = {M_PI};
	double result_test = plainmc(1,g,a_test,b_test,N);
	printf("result from 1/sqrt(x) from 0 to pi = %g +/- %g*i\n", creal(result_test),cimag(result_test));

	return 0;
}
