#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

double linterp(gsl_vector* x, gsl_vector* y, double z);
double linterp_integ(gsl_vector* x, gsl_vector* y, double z);


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
		gsl_vector_set(x,i,val * 2.3 + 1);
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
	printVec(x);
	printVec(y);
	double val = linterp(x, y, z);
	printf("Part A: First we test the linear interpolation function at value z.\\ With the x and y list:\n");
	printf("x: "); printVec(x);
	printf("y: "); printVec(y);
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
	//plots
	FILE* functions = fopen("functions.txt", "w");
	for(int i = 0; i < x -> size; i++){
		fprintf(functions,"%20g %20g\n",gsl_vector_get(x,i), gsl_vector_get(y,i));
	}
	FILE* points = fopen("points.txt","w");
	fprintf(points,"%20g %20g\n",val, interp_eval_gsl);

	//Cleaning
	fclose(functions);
	fclose(points);
	gsl_interp_free(cspline);

	return 0;
}
