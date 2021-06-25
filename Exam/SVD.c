#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <assert.h>
#include <gsl/gsl_cblas.h>


double dot(gsl_vector* x, gsl_vector* y){
	double x_dot_y;
	gsl_blas_ddot(x, y, &x_dot_y); //GSL function that calculates dot product
	return x_dot_y;
}

double norm(gsl_vector* x){
	return sqrt(dot(x,x));
}

void printM(gsl_matrix *A){
	for(int i = 0; i < A->size1; i++){
		for(int j = 0; j < A-> size2; j++){
			printf("%7.3f",gsl_matrix_get(A,i,j));
		}
		printf("\n");
	}
}

void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c = cos(theta);
	double s = sin(theta);
	for(int i = 0; i < A->size1; i++){
		double new_aip = c * gsl_matrix_get(A,i,p) - s * gsl_matrix_get(A,i,q);
		double new_aiq = s * gsl_matrix_get(A,i,p) + c * gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}


void SVD(gsl_matrix* A, gsl_matrix* U, gsl_matrix* D, gsl_matrix* V){
	//Checking that n = m
	assert(A->size1 == A->size2);

	//Allocating space for vectors and matrices
	int n = A -> size1;
	gsl_vector* ap = gsl_vector_alloc(n);
	gsl_vector* aq = gsl_vector_alloc(n);
	gsl_vector* A_col = gsl_vector_alloc(n);

	//The one-sided Singular Value Decomposition. Only small changes from Jacobi diagonalization with cyclic sweeps were made.
	//Theta is computed from the p, q columns of A, and we only the one-side A -> AJ.
	int changed;
	do{
		changed=0;
  		int m = A->size2;
		for(int p = 0; p < m - 1; p++){
			for(int q = p + 1; q < m; q++){
  				gsl_matrix_get_col(aq,A,q);
  				gsl_matrix_get_col(ap,A,p);


				double apq = dot(ap,aq);
				double aqq = dot(aq,aq);
				double app = dot(ap,ap);
				double theta = 0.5 * atan2(2*apq,aqq-app);
				double c = cos(theta);
				double s = sin(theta);

				double new_app = c * app - s * apq;
				double new_aqq = s * apq + c * aqq;

				if(new_app!=app || new_aqq!=aqq){
					changed=1;
					timesJ(A,p,q, theta);
					timesJ(V,p,q, theta); // V←V*J
				}
			}
		}
	}while(changed!=0);

	//Creating D = ||a'i||, ui=a'i/||a'i||
	for(int i = 0; i< n; i++){
		gsl_matrix_get_col(A_col,A,i);
		double norm_a = norm(A_col);
		gsl_matrix_set(D,i,i,norm_a);

		gsl_vector_scale(A_col,1/norm_a);
		gsl_matrix_set_col(U,i,A_col);
	}

	//Cleaning allocations
	gsl_vector_free(A_col);
	gsl_vector_free(ap);
	gsl_vector_free(aq);
}

