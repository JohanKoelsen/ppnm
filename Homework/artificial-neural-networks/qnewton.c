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

	double fx = F(x);

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

	//Cleaning
	gsl_matrix_free(B);
	gsl_vector_free(gx);
	gsl_vector_free(steps);
	gsl_vector_free(z);
	gsl_vector_free(gz);
	gsl_vector_free(y);
	gsl_vector_free(u);
	gsl_vector_free(a);

        return nsteps;
}//func closes bracket




#if 0
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#define TINY (1.0/524288)
//#define TINY (1.0/1048576)
//#define TINY (1.0/2097152)
//#define TINY 1.5e-8

void numeric_gradient
B(double beta(gsl_vector*), gsl_vector*x, gsl_vector*grad){
double fx=beta(x);
for(int i=0;i<x->size;i++){
	double xi=gsl_vector_get(x,i);
	double dx=fabs(xi)*TINY;
	if(fabs(xi)<sqrt(TINY)) dx=TINY;
	gsl_vector_set(x,i,xi+dx);
	gsl_vector_set(grad,i,(beta(x)-fx)/dx);
	gsl_vector_set(x,i,xi);
	}
}

int qnewton(double beta(gsl_vector*), gsl_vector*x, double acc) {
int n=x->size,nsteps=0,nbad=0,ngood=0;
gsl_matrix* B=gsl_matrix_alloc(n,n);
gsl_vector* gx=gsl_vector_alloc(n);
gsl_vector* Dx=gsl_vector_alloc(n);
gsl_vector* z=gsl_vector_alloc(n);
gsl_vector* gz=gsl_vector_alloc(n);
gsl_vector* y=gsl_vector_alloc(n);
gsl_vector* u=gsl_vector_alloc(n);
gsl_matrix_set_identity(B);
numeric_gradient(beta,x,gx);
double fx=beta(x),fz;

while(nsteps<2000){
	nsteps++;
if(fx<acc){fprintf(stderr,"qnewton:converged: fx<acc=%g\n",acc); break;}
	gsl_blas_dgemv(CblasNoTrans,-1,B,gx,0,Dx);
//	if(gsl_blas_dnrm2(Dx)<TINY*gsl_blas_dnrm2(x))
//		{fprintf(stderr,"qnewton: |Dx|<TINY*|x|\n"); break;}
	if(gsl_blas_dnrm2(gx)<acc)
		{fprintf(stderr,"qnewton: |grad|<acc\n"); break;}
	double lambda=1;
	while(1){
		gsl_vector_memcpy(z,x);
		gsl_vector_add(z,Dx);
		fz=beta(z);
		double sTg; gsl_blas_ddot(Dx,gx,&sTg);
		if(fz<fx+0.01*sTg){
			ngood++;
			break;
			}
		if(lambda<1.0/32){
			nbad++;
			break;
			}
		lambda*=0.5;
		gsl_vector_scale(Dx,0.5);
		}
	numeric_gradient(beta,z,gz);
	gsl_vector_memcpy(y,gz);
	gsl_blas_daxpy(-1,gx,y); /* y=grad(z)-grad(x) */
	gsl_vector_memcpy(u,Dx); /* u=s */
	gsl_blas_dgemv(CblasNoTrans,-1,B,y,1,u); /* u=s-By */
	double sTy,uTy;
	gsl_blas_ddot(Dx,y,&sTy);
	if(fabs(sTy)>1e-6){
		gsl_blas_ddot(u,y,&uTy);
		double gamma=uTy/2/sTy;
		gsl_blas_daxpy(-gamma,Dx,u); /* u=u-gamma*s */
		gsl_blas_dger(1.0/sTy,u,Dx,B);
		gsl_blas_dger(1.0/sTy,Dx,u,B);
		}
	gsl_vector_memcpy(x,z);
	gsl_vector_memcpy(gx,gz);
	fx=fz;
	}
gsl_matrix_free(B);
gsl_vector_free(gx);
gsl_vector_free(Dx);
gsl_vector_free(z);
gsl_vector_free(gz);
gsl_vector_free(y);
gsl_vector_free(u);
fprintf(stderr,"qnewton: nsteps=%i ngood=%i nbad=%i fx=%.1e\n"
		,nsteps,ngood,nbad,fx);
return nsteps;
}
#endif
