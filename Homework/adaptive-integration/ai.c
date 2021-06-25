#include <math.h>
#include <assert.h>
#include <stdio.h>



#define SQR2 1.41421356237309504880

double recursive_integrate(double f(double),double a, double b,double acc, double eps, double f2, double f3, int limit){
	double f1 = f(a+(b-a)/6);
	double f4 = f(a+5*(b-a)/6);

	//From the trapedium rule wi = [2/6, 1/6, 1/6, 2/6] and vi = [1/4, 1/4, 1/4, 1/4] we can define Q and q.
	double Q = (2 * f1 + f2 + f3 + 2 * f4)/6*(b-a);
	double q = (f1 + f4 + f2 + f3)/4 * (b-a);

	//Defining tolerance and error from the chapter
	double tolerance = acc + eps * fabs(Q);
	double error = fabs(Q-q)/3;


	/*Recursion conditions. If limit == 0 (999 iterations for a limit = 999 in integrate func) or error < tolerance, the recursion is done.
	However, during each recursion the interval is split into two bits which are then recursively split.
	*/
	if(limit == 0){
		fprintf(stderr,"recursion limit reached\n");
		return Q;
		}
	if(error < tolerance) return Q;
	else {
		double Q1 = recursive_integrate(f,a,(a+b)/2, acc/sqrt(2), eps, f1, f2, limit-1);
		double Q2 = recursive_integrate(f,(a+b)/2,b, acc/sqrt(2), eps, f3, f4 ,limit-1);
		return Q1 + Q2;
	}
}

double integrate(double f(double),double a,double b, double acc, double eps){
	//From equation 51 we define xi -> a + xi(b-a) where xi = [1/6, 2/6, 4/6, 5/6].
	double f2 = f(a+2*(b-a)/6);
	double f3 = f(a+4*(b-a)/6);
	//Set limit = 999.
	int limit = 999;
	return recursive_integrate(f,a,b,acc,eps,f2,f3,limit);
}



static double A,B; // not accessible from other files
double F(double f(double),double t){ // auxilliary function for Clenshaw-Curtis
	return f( (A+B)/2+(A-B)/2*cos(t) )*sin(t)*(B-A)/2;
	}


double CC24(double f(double),double a, double b, double acc, double eps, double f2, double f3, int nrec){
	assert(nrec<99);
	//Using xi = [1/6, 2/6, 4/6, 5/6] again in order to define f1 and f4, eq 51 and 48 in the chapter.
	double f1 = F(f,a+(b-a)/6);
	double f4 = F(f,a+5*(b-a)/6);

	//defining Q and q from eq. 44 & 45 in the chapter.
	double Q = (2*f1+f2+f3+2*f4)/6*(b-a);
	double q = (f1+f4+f2+f3)/4*(b-a);

	double tolerance = acc + eps * fabs(Q);
	double error = fabs(Q-q)/3;
	//Using recursion again. The end of recursion happens when error < tolerance.
	if(error < tolerance) return Q;
	else {
		double Q1 = CC24(f, a, (a+b)/2, acc/sqrt(2), eps, f1, f2, nrec+1);
		double Q2 = CC24(f, (a+b)/2, b, acc/sqrt(2), eps, f3, f4, nrec+1);
		return Q1+Q2; }
}
double clenshaw_curtis(double f(double), double a, double b, double acc, double eps){
	A = a;
	B = b;
	a = 0;
	b = M_PI;
	double f2 = F(f, a+2*(b-a)/6);
	double f3 = F(f, a+4*(b-a)/6);
	int nrec = 0;
	return CC24(f,a,b,2*acc,2*eps,f2,f3,nrec);
}

// gcc nested function
double clenshaw_curtis_nested( double f(double),double a,double b,double acc,double eps){
	double g(double t){return f( (a+b)/2+(a-b)/2*cos(t) )*sin(t)*(b-a)/2;}

	return integrate(g,0,M_PI,2*acc,2*eps);
}



//Part C ---- infitnite boundaries
static double A; // the left integration limit

static double F_c(double f(double),double t){// variable transformation formula
	return f( A+(1-t)/t )/t/t;
	}


double wrap24( double f(double),double a, double b, double acc, double eps, double f2, double f3, int nrec){
	assert(nrec<99);
	double f1 = F_c(f,a+(b-a)/6), f4 = F_c(f,a+5*(b-a)/6);
	double Q = (2*f1 + f2 + f3 + 2*f4)/6*(b-a), q = (f1 + f4 + f2 + f3)/4*(b-a);
	double tolerance = acc + eps*fabs(Q), error = fabs(Q-q);
	if(error < tolerance) return Q;
	else {
		double Q1 = wrap24(f, a, (a+b)/2, acc/SQR2, eps, f1, f2, nrec + 1);
		double Q2 = wrap24(f, (a+b)/2, b, acc/SQR2, eps, f3, f4, nrec + 1);
		return Q1+Q2; }
}
//New func with infinite boundaries
double integrate_infinite(double f(double),double a,double acc,double eps ){
	A = a;
	a = 0;
	double b = 1;
	double f2 = F_c(f,a+2*(b-a)/6);
	double f3 = F_c(f,a+4*(b-a)/6);
	int nrec = 0;
	return wrap24(f,a,b,2*acc,2*eps,f2,f3,nrec);
}


