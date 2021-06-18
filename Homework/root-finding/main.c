#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double eps);
void driver(void f(double t, gsl_vector* y, gsl_vector* dydt), double a, gsl_vector* ya, double b, gsl_vector* yb, double h, double acc, double eps);

//Rosenbrock function
void Rosenbrock_grad(gsl_vector* r, gsl_vector* Rr){
	double x = gsl_vector_get(r,0);
	double y = gsl_vector_get(r,1);
	double Rx = 2*(x-1) + 400 * (x*x - y)*x;
	double Ry = 200*(y-x*x);
	gsl_vector_set(Rr, 0, Rx);
	gsl_vector_set(Rr,1,Ry);
}


double e;
void schroedinger(double x, gsl_vector* y, gsl_vector* dydx){
	gsl_vector_set(dydx,0,gsl_vector_get(y,1));

	double dydx_val = -2*e*gsl_vector_get(y,0) - 2.0/x * gsl_vector_get(y,0);
	gsl_vector_set(dydx,1,dydx_val);

}


void SE_solved(gsl_vector* x, gsl_vector* M){
	e = gsl_vector_get(x,0);
	int n = 2; //Since it is a second order diff. equation
	//Defining variables for the ode45 solver.
	double a = 1e-3;
	double b = 8; //rmax

	gsl_vector* ya = gsl_vector_alloc(n);
	gsl_vector* yb = gsl_vector_alloc(n);
	double h = 0.01;
	double acc = 1e-3;
	double eps = 1e-3;


	//initial values
	gsl_vector_set(ya,0,a - a*a);
	gsl_vector_set(ya,1,1 - 2*a);


	driver(schroedinger, a, ya, b, yb, h, acc, eps);

	gsl_vector_set(M,0,gsl_vector_get(yb,0));

	gsl_vector_free(ya);
	gsl_vector_free(yb);
}

int main(){
	//------Part A---------
	printf("Part A: finding extremum of Rosenbrock valley\n");
	double eps_root = 0.01;
	gsl_vector* r = gsl_vector_alloc(2); //2d
	double x0 = -2;
	double y0 = 8; //Guesses for roots.
	gsl_vector_set(r,0,x0);
	gsl_vector_set(r,1,y0);
	newton(Rosenbrock_grad, r, eps_root);
	printf("The extremums of the valley is at: %g and %g.\n",gsl_vector_get(r,0),gsl_vector_get(r,1));

	/*

	//------Part B-------
	gsl_vector* x = gsl_vector_alloc(1); //1d
	gsl_vector_set(x,0,x0);
	newton(SE_solved, x, eps_root);
	printf("Part B: solving and finding extremum of SE\n");
	printf("The extremum is %g.\n", gsl_vector_get(r,0));
	*/
	//Cleaning
	gsl_vector_free(r);
	return 0;
}
