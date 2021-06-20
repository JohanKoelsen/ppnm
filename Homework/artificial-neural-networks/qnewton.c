#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <float.h>
static const double DELTA = sqrt(DBL_EPSILON);

void numeric_gradient(double F(gsl_vector* ), gsl_vector* x, gsl_vector* gradient){
	//defining dim, d
	int d = x -> size;
	double fx;
	fx = F(x);

	for(int i = 0; i < d; i++){
		double dx = DELTA;
		double xi = gsl_vector_get(x,i);
		if (fabs(xi) > sqrt(DELTA)) dx = fabs(xi)*DELTA;

		//Setting x[i] to xi + dx and then back to xi again.
		gsl_vector_set(x, i, xi + dx);
		gsl_vector_set(gradient, i, (F(x) - fx)/(dx));
		gsl_vector_set(x,i,xi);
	}

}

int qnewton(double F(gsl_vector* x), gsl_vector* x, double acc){
	int nsteps = 0;
	int n = x -> size;
	//Defining matrices
	gsl_matrix* B = gsl_matrix_alloc(n,n);
	gsl_vector* gx = gsl_vector_alloc(n);
	gsl_vector* steps = gsl_vector_alloc(n);
	gsl_vector* z  = gsl_vector_alloc(n);
	gsl_vector* gz = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* u = gsl_vector_alloc(n);
	gsl_vector* a = gsl_vector_alloc(n);

	//Following table 1 in the chapter:
	gsl_matrix_set_identity(B);
	numeric_gradient(F, x, gx);
	double fx = F(x);
	double fz;
	double stepsi;
	double ui;
	double sTg;
	double sTy;
	double norm2gx;
	double lambda;
	double uTy;
	double gamma;
	while(nsteps < 1000){
		norm2gx = 0;

		for(int i = 0; i < n; i++){
			norm2gx += gsl_vector_get(gx,i)*gsl_vector_get(gx,i);
		}// i for-loop bracket


		if (norm2gx < acc*acc) break;
		nsteps++;
		for(int i = 0; i < n; i++){
			gsl_vector_set(steps,i,0);
			for(int j = 0; j < n; j++){
				stepsi = gsl_vector_get(steps,i);
				gsl_vector_set(steps, i, stepsi - gsl_matrix_get(B,i,j) * gsl_vector_get(gx,j));
			} //j for-loop bracket
		}//i for-loop bracket

		lambda = 1;
		while(1){
			for(int i = 0; i < n; i++){
				gsl_vector_set(z,i,gsl_vector_get(x,i) + gsl_vector_get(steps,i));
			}//i for-loop bracket

			fz = F(z);
			sTg = 0;

			for(int i = 0; i < n; i++){
			sTg += gsl_vector_get(steps,i) * gsl_vector_get(gx,i);
			} //i for-loop bracket


			if (fz < fx + 0.01 * sTg) break;
			if (lambda < DELTA){
				gsl_matrix_set_identity(B);
				break;
			}//if bracket
			lambda /= 2;
			for(int i = 0; i < n; i++){
				stepsi = gsl_vector_get(steps,i);
				gsl_vector_set(steps,i,stepsi * 0.5);
			} //i for-loop bracket

		} //while(1) bracket

		numeric_gradient(F,z,gz);
		for(int i = 0; i < n; i++){
			gsl_vector_set(y,i,gsl_vector_get(gz,i) - gsl_vector_get(gx,i));
		} //i for-loop bracket
		for(int i = 0; i < n; i++){
			gsl_vector_set(u,i,gsl_vector_get(steps,i));
			for(int j = 0; j < n; j++){
				ui = gsl_vector_get(u,i);
				gsl_vector_set(u,i,ui - gsl_matrix_get(B,i,j) * gsl_vector_get(y,j));
			}//j for-loop bracket
		}// i for-loop bracket
		sTy = 0;
		for(int i = 0; i < n; i++){
			sTy += gsl_vector_get(steps,i) * gsl_vector_get(y,i);
		}//i for-loop bracket

		if(fabs(sTy) > DELTA){
			uTy = 0;
			for(int i = 0; i < n; i++){
				uTy += gsl_vector_get(u,i)*gsl_vector_get(y,i);
			}//i for-loop bracket
			gamma = uTy/2/sTy;
			for(int i = 0; i < n; i++){
				gsl_vector_set(a,i, gsl_vector_get(u,i) - gamma * gsl_vector_get(steps,i));
			}//i for-loop bracket

			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					gsl_matrix_set(B,i,j, gsl_matrix_get(B,i,j) + gsl_vector_get(a,i)*gsl_vector_get(steps,j)/sTy);
				}//j for-loop
			}// i for-loop

			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					gsl_matrix_set(B,i,j, gsl_matrix_get(B,i,j) + gsl_vector_get(steps,i) * gsl_vector_get(a,j)/sTy);
				}//j for-loop
			}//i for-loop
		}//if bracket
		for(int i = 0; i < n; i++){
			gsl_vector_set(x,i,gsl_vector_get(z,i));
		}
		for(int i = 0; i < n; i++){
			gsl_vector_set(gx,i,gsl_vector_get(gz,i));
		}
		fx = fz;
	} //first while bracket
	return nsteps;
	//Cleaning
	gsl_matrix_free(B);
	gsl_vector_free(gx);
	gsl_vector_free(steps);
	gsl_vector_free(z);
	gsl_vector_free(gz);
	gsl_vector_free(y);
	gsl_vector_free(u);
	gsl_vector_free(a);
}//func closes bracket


