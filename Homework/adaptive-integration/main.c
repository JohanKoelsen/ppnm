#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>

//Defining functions from ai.c
double recursive_integrate(double f(double), double a, double b, double acc, double eps, double f2, double f3, int limit);
double integrate(double f(double), double a, double b, double acc, double eps);
double clenshaw_curtis(double f(double), double a, double b, double acc, double eps);
double CC24(double f(double), double a, double b, double acc, double eps, double f2, double f3, int nrec);



int calls;

double g1(double x){
	calls++;
	return sqrt(x);
};
double g2(double x){
	calls++;
	return 4*sqrt(1-x*x);
}

double g2_GSL(double x, void* params){
	calls++;
	double alpha = *(double*) params;
	return alpha*4*sqrt(1-x*x);}

int main(){
	//Interval from [0,1] with acc = 0.001 and eps = 0.001.
	double a = 0;
	double b = 1;
	double acc = 0.0001;
	double eps = 0.0001;
	calls = 0;
	double Q1 = integrate(g1,a,b,acc,eps);
	double exact = 2./3;
	//Printing results in a with f = 1/sqrt(x).
	printf("Part A: recursive adaptive integration. In this part the code integrate two expressions, f = sqrt(x) and f = 4*(1-x^2) over some intervals.\n\n");
	printf("f = sqrt(x) from %g to %g\n",a,b);
	printf("Q = %g\n",Q1);
	printf("exact = %g\n",exact);
	printf("calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q1)*eps);
	printf("actual error = %g\n",fabs(Q1-exact));

	printf("\n");
	//Printing results with f = 4*(1-x^2).
	calls = 0;
	exact = M_PI;
	double Q2 = integrate(g2,a,b,acc,eps);
	printf("f = 4*sqrt(1-x^2) from %g to %g\n",a,b);
	printf("Q = %g\n",Q2);
	printf("exact = %g\n",exact);
	printf("calls = %d\n",calls);
	printf("estimated error = %g\n",acc + fabs(Q2)*eps);
	printf("actual error = %g\n",fabs(Q2-exact));



	//---B---
	calls = 0;

	double Q3 = clenshaw_curtis(g2,a,b,acc,eps);
	exact = M_PI;
	printf("\nPart B: an  adaptive integrator with the Clenshaw–Curtis variable transformation is implemented, tested and specifically calculating the\n");
	printf("integral of f(x) = 4*srqt(1-x^2) from 0 to 1\n\n");
	printf("f from %g to %g : CLENSHAW_CURTIS\n",a,b);
	printf("Q = %g\n",Q3);
	printf("exact = %g\n",exact);
	printf("calls = %d\n",calls);
	printf("estimated error = %g\n",acc+fabs(Q3)*eps);
	printf("actual error = %g\n",fabs(Q3-exact));
	printf("By comparing with the previous integration from part A, we can see that the result is more correct, however the number of calls and the error is larger.\n");

	//Comparing with GSL
	double alpha = 1.0;
	int limit = 10000;
	double result;
	double error;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
	calls = 0;
	gsl_function G;
	G.function = &g2_GSL;
	G.params = &alpha;
	gsl_integration_qags(&G, a, b, acc, eps, limit, w, &result, &error);
	printf("\nThe GSL's integration routines provide = %f with %d number of calls. The error is %g.\n", result, calls, error);
	printf("While providing the correct result, the number of calls is larger.\n");

	gsl_integration_workspace_free(w);


	return 0 ;
}

