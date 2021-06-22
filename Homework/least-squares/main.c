#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

double linfun(int i, double x){
	if(i==0) return 1;
	else if(i==1) return x;
	else return NAN;
}


//Defining functions that are used in this script. lsfit is from lsfit.c and printVector from gs.c
void lsfit(gsl_vector* x, gsl_vector* y, gsl_vector* dy, double funs(int, double), int m, gsl_vector* c, gsl_vector* dc, gsl_matrix* Cov);
void printVector(gsl_vector* v);

int main(){
	//Task A) We fit the data ith exponential function in the usual logarithmic way.
	//We define vector t and y from the data given in the exercise.
	printf("Part A: we fit the data with an exponential function in the usual logarithmic way.\n");
	int N = 9; //size of data set
	gsl_vector* t = gsl_vector_alloc(N);
	gsl_vector* y = gsl_vector_alloc(N);

	//Insering data points into vector t and y.
	double time[9] = {1,  2,  3, 4, 6, 9,   10,  13,  15};
	double activity[9] = {117,100,88,72,53,29.5,25.2,15.2,11.1};

	for(int k = 0; k < 9; k++){
		gsl_vector_set(t,k,time[k]);
		gsl_vector_set(y,k,activity[k]);
	}
	printf("The data points are inserted:\n");
	printVector(t);
	printVector(y);
	//Since we fit with exponential function in the logarithmic way of ln(y) = ln(a) - lambda, we make the vector logy and dlogy.
	gsl_vector* logy = gsl_vector_alloc(N);
	gsl_vector* dlogy = gsl_vector_alloc(N);
	for(int i = 0; i < N; i++){
		double yi = gsl_vector_get(y,i);

		double dyi = yi/20; //this error is  defined in the exercise
		double log_yi = log(yi);
		double dlog_yi = dyi/yi;
		gsl_vector_set(logy,i,log_yi);
		gsl_vector_set(dlogy,i,dlog_yi);
	}

	//Solving the system
	int M = 2; //since this is a linear fit
	gsl_vector *c = gsl_vector_alloc(M); //The coefficients
	gsl_vector* dc = gsl_vector_alloc(M); //The uncertainties of coefficients
	gsl_matrix* Cov = gsl_matrix_alloc(M,M);
	lsfit(t, logy, dlogy, linfun, 2, c, dc, Cov);


	printf("The coefficients of the linear fit is\n");
	printVector(c);
	printf("And the uncertainties of c\n");
	printVector(dc);

	//Saving c0, c1 and dc0, dc1 in vars
	double c0 = gsl_vector_get(c,0);
	double dc0 = gsl_vector_get(dc,0);
	double c1 = gsl_vector_get(c,1);
	double dc1 = gsl_vector_get(dc,1);

	//Finding half life
	double half_life = log(2)/(-c1);
	printf("The half life of 224Ra is found to be tau = %g\n",half_life);

	//Saving values for the fit and the linear experimental values in a txt file, so that I can make a plot.
	FILE* linfit = fopen("linfit.txt","w");
	FILE* exp_vals = fopen("exp_vals.txt","w");
	double max = time[8];
	double ts = 0;
	double step_size = 0.10;
	while(ts <= max){
		fprintf(linfit,"%20g %20g %20g %20g\n",ts, c0 + c1*ts, (c0 - dc0) + (c1 - dc1)*ts, (c0 + dc0) + (c1 + dc1)*ts);
		ts = ts + step_size;
	}
	for(int i = 0; i < 9; i++){
		fprintf(exp_vals,"%20g %20g %20g\n",gsl_vector_get(t,i),gsl_vector_get(logy,i), gsl_vector_get(dlogy,i));

	}


	//Finding uncertainty of the lifetime.
	double err_half_life = log(2)/(c1*c1)*dc1;
	printf("The uncertainty of the half life is found by error-propagation law to be %g\n",err_half_life);

	printf("The solution is shown in fit.png with uncertainties of the fitting coefficients (Part B) included aswell.\n");

	//Clearing vars and closing files
	fclose(linfit);
	fclose(exp_vals);
	gsl_vector_free(t);
	gsl_vector_free(y);
	gsl_vector_free(logy);
	gsl_vector_free(dlogy);
	gsl_vector_free(c);
	gsl_vector_free(dc);
	gsl_matrix_free(Cov);
	return 0;
}
