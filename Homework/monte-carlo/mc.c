#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#define RND (double ) rand ( ) /RAND_MAX


complex plainmc(int dim,double f(double* x),double* a,double* b,int N){
        //Function from the homework page

	double V = 1;
	for(int i = 0; i < dim; i++){
		V*= b[i] - a[i];
	}
        double sum = 0;
	double sum2 = 0;
	double x[dim];
        for(int i = 0; i < N; i++){
                for(int j = 0; j < dim; j++){
			x[j] = a[j] + RND * (b[j] - a[j]);
		}
                double fx = f(x);
		sum += fx;
		sum2 += fx*fx;
        }
        double mean = sum/N;
	double sigma = sqrt(sum2/N - mean*mean);
        complex result = mean * V + I * sigma * V/sqrt(N);


        return result;
}
/*

#define fracl(x) ((x) − floorl(x))
#define real long double

void lattice (int d, double *x ){
	static int dim = 0;
	static int n = 0;
	static real *alpha;
	int i;

	if (d < 0){
		dim = -d;
		n = 0;
		alpha = (real*)realloc(alpha,dim*sizeof(real));
		for(i = 0; i < dim; i++){
			alpha[i] = fracl(sqrtl(M_PI + i));
		}
	}
	else if (d > 0){
		n++;
		assert(d==dim && n > 0);
		for(i = 0; i < dim; i++){
			x[i] = fracl(n*alpha[i]);
		}
	}
	else if(alpha != NULL){
		free(alpha);
	}
	return ;
}
*/







