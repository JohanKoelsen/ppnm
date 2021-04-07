#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

//Getting GS-functions from gs.c
void GS_decomp(gsl_matrix* A, gsl_matrix* B);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* X);

//Create a least-square fit from QR decomposition.
void lsfit(gsl_vector* x, gsl_vector* y, gsl_vector* dy, double funs(int, double), int m, gsl_vector* c, gsl_vector* dc, gsl_matrix* Cov){
	//This function follows the steps from the notes, which finds the coefficients c from eq 8.
	//First we need to define the matrix A and vector b from eq 7 in the notes.

	int n = x -> size;
	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_vector* b = gsl_vector_alloc(n);
	for(int i = 0; i < n; i++){
		double xi = gsl_vector_get(x,i);
		double yi = gsl_vector_get(y,i);
		double dyi = gsl_vector_get(dy,i);
		double bi = yi/dyi;

		//Defining b
		gsl_vector_set(b,i,bi);

		//Defining A
		for(int k = 0; k < m; k++){
			gsl_matrix_set(A,i,k,funs(k,xi)/dyi);
		}

	}
	//Now we can can use QR-decomposition on A defined in gs.c
	gsl_matrix* R = gsl_matrix_alloc(m,m);
	GS_decomp(A,R);
	//And then solve by using GS_solve defined in gs.c
	GS_solve(A,R,b,c);


	//And the uncertainties can now be computed
	gsl_matrix* I = gsl_matrix_alloc(m,m);
	gsl_matrix_set_identity(I);
	gsl_matrix* Rinv = gsl_matrix_alloc(m,m);
	GS_inverse(I,R,Rinv); //This function is from gs.c. It calculates the inverse of R and saves it in Rinv.
	//Now we can compute the CoV matrix which is given as Cov = R^-1*(R^-1)^T (from eq. 15 in the notes)
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, Rinv, Rinv, 0, Cov);
	
	//Now the uncertainties of c can be found using eq. 16.
	for(int k=0; k<m; k++){
		double Covkk = gsl_matrix_get(Cov, k, k);
		double dck = sqrt(Covkk);
		gsl_vector_set(dc, k, dck);
	}

	//Clearing vars
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_vector_free(b);
	gsl_matrix_free(I);
	gsl_matrix_free(Rinv);
}
