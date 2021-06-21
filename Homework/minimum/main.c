#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


//Defining functions in main
void gradient(double F(gsl_vector*), gsl_vector* x, gsl_vector* grad);
int qnewton(double F(gsl_vector*), gsl_vector* x, double acc);


//Rosenbrock function
double Rosenbrock(gsl_vector* r){
	double x = gsl_vector_get(r,0);
	double y = gsl_vector_get(r,1);
	return (1-x)*(1-x)+100*(y-x*x)*(y-x*x);
}

//Himmelblau's function
double Himmelblau(gsl_vector* r){
	double x = gsl_vector_get(r,0);
	double y = gsl_vector_get(r,1);
	return (x*x+y-11)*(x*x+y-11)+(x+y*y-7)*(x+y*y-7);
}


double BW(double Energy, gsl_vector* p){
	double m = gsl_vector_get(p, 0);
	double Gamma = gsl_vector_get(p, 1);
	double A = gsl_vector_get(p, 2);
	return A/((Energy-m)*(Energy-m)+Gamma*Gamma/4);
}

// Deviation function (B)

int ndata = 30;
gsl_vector* energy;
gsl_vector* sigma;
gsl_vector* dsigma;

double D(gsl_vector* p){
	double sum = 0;

	for(int i=0; i<ndata; i++){
		double Ei = gsl_vector_get(energy,i);
		double sigmai = gsl_vector_get(sigma, i);
		double dsigmai = gsl_vector_get(dsigma,i);

		double val = gsl_vector_get(p,2)/((Ei - gsl_vector_get(p,0))*(Ei-gsl_vector_get(p,0)) + dsigmai*dsigmai);
		double Di = (val-sigmai)*(val-sigmai)/dsigmai/dsigmai;
		sum += Di;
	}
	return sum;
}



int main(){
	//-----Part A-----
	int n = 2;
	gsl_vector* x0 = gsl_vector_alloc(n);
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y0 = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);

	gsl_vector_set(x0, 0, -2);
	gsl_vector_set(x0, 1, 8); // initial guess

	double acc = 1e-4;
	gsl_vector_memcpy(x,x0);
	int steps = qnewton(Rosenbrock, x, acc);
	printf("Minimum part A: finding the minimum of Rosenbrock and Himmelblau functions:\n\n");
	printf("Rosenbrock function:\n");
	printf("Rosenbrock minimum found in %d steps.\n",steps);
	printf("Rosenbrock init guess: (%g ,%g)\n", gsl_vector_get(x0,0),gsl_vector_get(x0,1) );
	printf("Rosenbrock found val: (%g,%g)\n", gsl_vector_get(x,0), gsl_vector_get(x,1) );

	acc = 1e-5;
	gsl_vector_set(y0,0,4);
	gsl_vector_set(y0,1,2);
	gsl_vector_memcpy(y,y0);
	steps = qnewton(Himmelblau,y,acc);
	printf("\nThe Himmelblau:\n");
	printf("Himmelblau min found in %d steps.\n",steps);
	printf("Himmelblau init guess: (%g,%g)\n",gsl_vector_get(y0,0),gsl_vector_get(y0,1));
	printf("Himmelblau found val: (%g,%g)\n",gsl_vector_get(y,0),gsl_vector_get(y,1));


	gsl_vector_free(x0);
	gsl_vector_free(x);
	gsl_vector_free(y0);
	gsl_vector_free(y);

	//------Part B-------
	int ndata = 30;
	int N = 3;
	gsl_matrix* data = gsl_matrix_alloc(ndata,N);

	FILE* Higgs_file = fopen("Higgs.txt", "r");
	gsl_matrix_fscanf(Higgs_file, data);
	fclose(Higgs_file);


	//defining vectors:
	gsl_vector* param0 = gsl_vector_alloc(N);
	gsl_vector* param = gsl_vector_alloc(N);

	energy = gsl_vector_alloc(ndata);
	sigma = gsl_vector_alloc(ndata);
	dsigma = gsl_vector_alloc(ndata);

	//Defining init values
	gsl_vector_set(param0,0,125);
	gsl_vector_set(param0,1,0.5);
	gsl_vector_set(param0,2,7);

	gsl_vector_memcpy(param,param0);

	//Setting values from Higgs_file into Energy, sigma and dsigma
	gsl_matrix_get_col(energy,data,0);
	gsl_matrix_get_col(sigma,data,1);
	gsl_matrix_get_col(dsigma,data,2);

	acc = 1e-5;
	printf("\n\nPart B: The tabulated data are scanned and evaluated.\n");
	steps = qnewton(D,param,acc);
	printf("Initial values: (%g,%g,%g)\n",gsl_vector_get(param0,0),gsl_vector_get(param0,1),gsl_vector_get(param0,2));
	printf("Minima found: (%g,%g,%g)\n",gsl_vector_get(param,0),gsl_vector_get(param,1),gsl_vector_get(param,2));
	printf("solution found in %d steps.\n",steps);

	//Cleaning
	gsl_vector_free(energy);
	gsl_vector_free(sigma);
	gsl_vector_free(dsigma);
	gsl_vector_free(param0);
	gsl_vector_free(param);
	gsl_matrix_free(data);

	return 0;
}
