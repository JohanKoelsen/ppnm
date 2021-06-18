
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <stdlib.h>
#include <float.h>

double dot(gsl_vector* x, gsl_vector* y);
double norm(gsl_vector* x);
void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

//x = start conditions.
void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps){
	//Defining parameters:
	int n = x -> size;
	int limit = 10000;
	int steps = 0;


	//Allocating space to all matrices
	gsl_matrix* J = gsl_matrix_alloc(n,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_vector* Dx = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* fy = gsl_vector_alloc(n);
	gsl_vector* fx = gsl_vector_alloc(n);
	while (steps < limit){
		f(x,fx); //f calculates f(x) and stores the value in fx.

		//The Jacobian is defined in the following:
		for (int j = 0; j < n; j++){
			double xj = gsl_vector_get(x,j);

			gsl_vector_set(x,j,xj + sqrt(DBL_EPSILON));
			f(x,df); //Since f(x) --> f(x+sqrt(epsilon)), we store df = f(x+sqrt(eps)), and set df = f(x + sqrt(eps)) - f(x)
			gsl_vector_sub(df,fx);


			for (int i = 0; i < n; i++){
				double Jij = (gsl_vector_get(df,i)/sqrt(DBL_EPSILON));
				gsl_matrix_set(J,i,j,Jij);

			}
			gsl_vector_set(x,j,xj); // Setting x[j] -= dx

		}
		//We use GS decomposition from the linear-equation homework to solve J*Dx = -f(x).
		GS_decomp(J,R);
		GS_solve(J,R,fx,Dx); //The solution is stored in Dx.
		gsl_vector_scale(Dx,-1.0);

		double s = 1.0;
		while(s>1./64){
			gsl_vector_memcpy(y,x);
			gsl_vector_add(y,Dx); //y = x + Dx*s
			f(y,fy); //fy = f(y)

			if(norm(fy) < (1-s*0.5)*norm(fx)) break;
			s *= 0.5;
			gsl_vector_scale(Dx,0.5); //Such that Dx -> Dx*s
			}

		gsl_vector_memcpy(x,y);
		gsl_vector_memcpy(fx,fy);
		if (norm(Dx) <sqrt(DBL_EPSILON) || norm(fx) < eps) break;
		steps++;

	}

	//Cleaning up
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(Dx);
	gsl_vector_free(y);
	gsl_vector_free(fy);
	gsl_vector_free(df);
	gsl_vector_free(fx);
}
