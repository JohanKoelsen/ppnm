#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void GS_decomp(gsl_matrix* A, gsl_matrix* R){
	// Initializing m, the size of matrix R is (m, m)
	int m = A-> size2;
	// Defining a double for later use for the dot product value
	double dotProduct;
	// Gram-Schmidt decomposition
	for (int i = 0; i < m; i++) {
		// Create a view of the i-th column of the matrix
		gsl_vector_view col = gsl_matrix_column(A, i);
		// Compute the norm of this column (sum of entries)
		gsl_matrix_set(R, i, i, gsl_blas_dnrm2(&col.vector));
		// Divide all entries of this column by its norm
		gsl_vector_scale(&col.vector, 1 / gsl_matrix_get(R, i, i));
		//
		for (int j = i + 1; j < m; j++) {
			// Create a view of the j-th colum of the matrix
			gsl_vector_view colJ = gsl_matrix_column(A, j);
			// Compute the dot product between A_i and A_j
			gsl_blas_ddot(&col.vector, &colJ.vector, &dotProduct);
			gsl_matrix_set(R, i, j, dotProduct);
			// Set A_j to be A_j - q_i * R_{ij}

			gsl_blas_daxpy(-dotProduct,&col.vector,&colJ.vector);
		}
	}
}

void printMatrix(gsl_matrix* M){
	for (int i = 0; i < M->size1; i++) {
		for (int j = 0; j < M->size2; j++) {
			printf("%f	", gsl_matrix_get(M, i, j));
		}
		printf("\n");
	}
}
void printVector(gsl_vector* v){
	for (int i = 0; i < v->size; i++){
		printf("%f       ",gsl_vector_get(v,i));
	}
	printf("\n");
}
double randomBetweenPlusMinus1(){
	// Initialize
	double randomDouble;
	// For a random seed (random due to dependence on computer's internal clock).
	//srand(time(NULL));
	// Generating random double value from random int
	randomDouble = (double)rand()/RAND_MAX*2.0-1.0; // Double in range -1 to 1
	// Returning the random double value between -1 and 1
	return randomDouble;
}

void backsub(gsl_matrix* A,gsl_vector* v){
	int n = v -> size;
	for(int i = n - 1; i >= 0; i--){
		double vi = gsl_vector_get(v,i);
		for(int j = i + 1; j < n; j++){
			vi -= gsl_matrix_get(A,i,j) * gsl_vector_get(v,j);
		}
		gsl_vector_set(v,i,vi/gsl_matrix_get(A,i,i));

	}

}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	gsl_blas_dgemv(CblasTrans, 1, Q, b, 0 ,x);
	backsub(R,x);

}

void QtimesQ(gsl_matrix* A){
	int n = A->size2;
	gsl_matrix* Q = gsl_matrix_alloc(n,n);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans,1, A, A, 0, Q);
	printMatrix(Q);
	gsl_matrix_free(Q);

}

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
	//Finding the inverse of A. This is done by solving QR = A, QR * xi = ei. The solution will be stored in B.
	gsl_vector* ei = gsl_vector_alloc(Q->size2);
	for(int i = 0; i < Q->size2; i++){
		gsl_vector_set_basis(ei,i);
		gsl_vector_view col = gsl_matrix_column(B,i);
		GS_solve(Q,R,ei,&col.vector);
	}
	gsl_vector_free(ei);

}


