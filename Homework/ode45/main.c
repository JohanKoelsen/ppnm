#include <math.h>
#include <gsl/gsl_vector.h>

//Defining functions
double norm(gsl_vector* v);

void driver(void f(double t, gsl_vector* b, gsl_vector* dydt), double a, gsl_vector* ya, double b, gsl_vector* yb, double h, double acc, double eps);

void rkstep12(void f(double t, gsl_vector* y, gsl_vector* dydt), double t, gsl_vector* yt, double h, gsl_vector* yh, gsl_vector* err);


//Defining the differential equation
void fun(double t, gsl_vector* y, gsl_vector* dydt){

	double y0 = gsl_vector_get(y,0);
	double y1 = gsl_vector_get(y,1);

	gsl_vector_set(dydt, 0, y1);
	gsl_vector_set(dydt, 1, -y0);
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

	gsl_vector_set(y,0,0);
	gsl_vector_set(y,1,1);

	driver(fun,a,y,b,res,h,acc,eps);

return 0;
}
