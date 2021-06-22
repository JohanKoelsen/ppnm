#include <math.h>
#include <gsl/gsl_vector.h>

//Defining functions
double norm(gsl_vector* v);

void driver(void f(double t, gsl_vector* b, gsl_vector* dydt), double a, gsl_vector* ya, double b, gsl_vector* yb, double h, double acc, double eps, char* path);

void rkstep12(void f(double t, gsl_vector* y, gsl_vector* dydt), double t, gsl_vector* yt, double h, gsl_vector* yh, gsl_vector* err);


//Defining the differential equation
void fun(double t, gsl_vector* y, gsl_vector* dydt){

	double y0 = gsl_vector_get(y,0);
	double y1 = gsl_vector_get(y,1);

	gsl_vector_set(dydt, 0, y1);
	gsl_vector_set(dydt, 1, -y0);
}


double N;
double T_c;
double T_r;

void fSIR(double x, gsl_vector* y, gsl_vector* dydx){
	gsl_vector_set(dydx,0,-gsl_vector_get(y,0)*gsl_vector_get(y,1)/(N*T_c));
	gsl_vector_set(dydx,1,gsl_vector_get(y,0)*gsl_vector_get(y,1)/(N*T_c)-gsl_vector_get(y,1)/T_r);
	gsl_vector_set(dydx,2,gsl_vector_get(y,1)/T_r);
}


int main() {
	//Running from a = 0 to b = 2*pi with step (h) = 0.1, accuracy goal (acc) = 1e-2 and relative accuracy goal (eps) = 1e-2
	double a = 0;
	double b = 2*M_PI;
	double h = 0.1;
	double acc = 1e-2;
	double eps = 1e-2;
	int n=2;

	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* res = gsl_vector_alloc(n);

	//Part A: solving the harmonics

	gsl_vector_set(y,0,0);
	gsl_vector_set(y,1,1);
	char* path = "harmonics.txt";
	driver(fun,a,y,b,res,h,acc,eps,path);

	printf("Part A: ode45 has solved the harmonic differential equation. The plot is shown in harmonics.png\n");

	printf("Part A 2): Investigating the effect of increasing T_c is done by solving the fSIR multiple times with increasing T_c.\n\n");
	//Part A: solving the fSIR -- defining parameters:
	n = 3;
	a = 0;
	b = 75;
	gsl_vector* ya = gsl_vector_alloc(n);
	gsl_vector* yb = gsl_vector_alloc(n);
	path = "SIR1.txt";
	N = 6e6; //An estimation of the population

	gsl_vector_set(ya,0,N);
	gsl_vector_set(ya,1,50);
	gsl_vector_set(ya,2,0);

	T_c = 0.5;
	T_r = 10;

	//RUNNING SIR 3 TIMES FOR DIFFERENT T_C.
	driver(fSIR, a, ya, b, yb, h, acc, eps, path);


	//Next T_c
	T_c = 1.5;
	T_r = 10;
	gsl_vector_set(ya,0,N);
	gsl_vector_set(ya,1,50);
	gsl_vector_set(ya,2,0);
	path = "SIR2.txt";
	driver(fSIR, a, ya,b ,yb ,h ,acc ,eps ,path);

	//Last T_c
	T_c = 3;
	T_r = 10;
	gsl_vector_set(ya,0,N);
	gsl_vector_set(ya,1,50);
	gsl_vector_set(ya,2,0);
	path = "SIR3.txt";
	driver(fSIR,a, ya, b, yb, h, acc, eps, path);

	printf("We can see that increasing T_C decreases the maximum amount of infected as well as shifting the solution to the left.\n");
	printf("The solutions are plotted in sir1.png, sir2.png and sir3.png.\n");
	printf("\nPart B: {ti, y(ti)} is stored in the given path in the driver function\n");
	//Cleaning
	gsl_vector_free(y);
	gsl_vector_free(res);
	gsl_vector_free(ya);
	gsl_vector_free(yb);
return 0;
}
