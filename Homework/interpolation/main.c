#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
typedef struct {int n; double *x, *y, *b, *c;} qspline;
typedef struct {gsl_vector* x, *y, *b, *c;} qinterp;

double linterp(gsl_vector* x, gsl_vector* y, double z);
double linterp_integ(gsl_vector* x, gsl_vector* y, double z);
int binsearch(gsl_vector* x, double z);
qspline* qspline_alloc(int n, double* x, double* y);
double qspline_eval(qspline *s, double z);
void qspline_free(qspline* s);
double qspline_der(qspline *s, double z);
double qspline_integ(qspline *s, double z);


void printVec(gsl_vector* x){
	for(int i = 0; i < x->size; i++){
		printf("%g ",gsl_vector_get(x,i));
	}
	printf("\n");

}

void func(gsl_vector* x){
	int d = x -> size;
	double val;
	for(int i = 0; i < d; i++){
		val = gsl_vector_get(x,i);
		gsl_vector_set(x,i,val * val + 2.3 *val + 1);
	}
}
void linfunc(gsl_vector* x){
	int d = x -> size;
	for(int i = 0; i < d; i++){
		gsl_vector_set(x,i,i);

	}

}

int main(){

	//Defining params
	int n = 7;
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	linfunc(x);
	gsl_vector_memcpy(y,x);
	func(y);
	double z = 3.5;
	double val = linterp(x, y, z);
	printf("Part A: First we test the linear interpolation function at value z.\\ With the x and y list:\n");
	printf("x: "); printVec(x);
	printf("y: "); printVec(y);
	printf("Now we test our linear interpolator and integrator vs GSL.\n");
	printf("The interpolated value at z = %g is %g.\n",z,val);

	//linterp_integ
	printf("We now calculate the integral of the linear spline of the tabulated data from the point x[0] = %g to the point z = %g.\n",gsl_vector_get(x,0), z);
	double integ_val = linterp_integ(x,y,z);
	printf("The value is %g.\n",integ_val);

	//GSL_interp
	double xa[x->size];
	double ya[x->size];
	for(int i = 0; i < x->size; i++){
		xa[i] = gsl_vector_get(x,i);
		ya[i] = gsl_vector_get(y,i);
	}
	gsl_interp* cspline = gsl_interp_alloc(gsl_interp_cspline,n);
	gsl_interp_init(cspline, xa, ya, n);
	double interp_eval_gsl = gsl_interp_eval(cspline, xa, ya, z, NULL);
	printf("\nThe GSL interpolation  at z = %g yield %g.\n",z,interp_eval_gsl);
	//plots linear interpolation
	FILE* functions = fopen("functions.txt", "w");
	for(int i = 0; i < x -> size; i++){
		fprintf(functions,"%20g %20g\n",gsl_vector_get(x,i), gsl_vector_get(y,i));
	}
	FILE* points = fopen("points.txt","w");
	for(int i = 0; i < n; i++){
		val = linterp(x,y,i);
		interp_eval_gsl = gsl_interp_eval(cspline, xa, ya, i, NULL);
		fprintf(points,"%20d %20g %20g\n", i, val, interp_eval_gsl);
	}


	//plots integral

	double interp_integ_gsl = gsl_interp_eval_integ(cspline, xa, ya, xa[0], z, NULL);
	printf("The integrated value from GSL from x[0] = %g to z = %g yield %g.\n",xa[0],z,interp_integ_gsl);
	FILE* integration_points = fopen("integration_points.txt","w");
	for(int i = 0; i < n; i++){
		integ_val = linterp_integ(x,y,i);
		interp_integ_gsl = gsl_interp_eval_integ(cspline,xa,ya,xa[0],i,NULL);
		fprintf(integration_points,"%20d %20g %20g\n",i,integ_val, interp_integ_gsl);
	}
	printf("The comparison between the interpolations and integrations are shown in the linplot.png and integplot.png.\n");
	int bin_val = binsearch(x,z);
	printf("\nBinsearch at z = %g yields %d.\n",5.21,bin_val);

	//-------- Part B---------

	qspline* Q = qspline_alloc(n,xa,ya);
	FILE* qspline_file = fopen("qspline.txt", "w");
	double di = 0.1;
	for(double i = xa[0]; i <= xa[n-1]; i+= di){
		double qz = qspline_eval(Q,i);
		double qd = qspline_der(Q,i);
		double qi = qspline_integ(Q,i);
		fprintf(qspline_file,"%20g %20g %20g %20g\n",i, qz,qd,qi );
	}


	//Cleaning
	qspline_free(Q);
	fclose(qspline_file);
	fclose(functions);
	fclose(points);
	fclose(integration_points);
	gsl_interp_free(cspline);

	return 0;
}
