#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void printMatrix(gsl_matrix* M);
double randomBetweenPlusMinus1();
void GS_solve(gsl_matrix* Q, gsl_matrix* R,gsl_vector* b,gsl_vector* x);
void printVector(gsl_vector* v);
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);


int main(void){
	// Part A,1

	int n = 5; // dimension
	int m = 3; // dimension
	gsl_matrix* A = gsl_matrix_alloc(n, m);
	// Assigning entries in A with floats of random values ranging from -1 to 1.
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			gsl_matrix_set(A, i, j, randomBetweenPlusMinus1());
		}
	}
	printf("Part A: We are checking that the QR-decomposition is working properly.\n\n");
	printf("Matrix A before QR-decomposition\n");
	printMatrix(A);
	gsl_matrix* R = gsl_matrix_alloc(m, m);
	GS_decomp(A, R);
	printf("Matrix A after QR-decomposition\n");
	printMatrix(A);
	printf("Matrix R after QR-decomposition\n");
	printMatrix(R);
	printf("Q^T*Q\n");
	//Checking for QTQ:
	gsl_matrix* Q = gsl_matrix_alloc(m,m);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,A,0,Q);
	printMatrix(Q);
	gsl_matrix_free(Q);


	//Checking that Q*R = A
	printf("Checking that Q*R = A again:\n");
	gsl_matrix* QR = gsl_matrix_alloc(n,m);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, A, R, 0, QR);
	printMatrix(QR);
	gsl_matrix_free(QR);


	/*
	Part A,2
	*/

	//Making GS_solve. Generating n x n matrix, R and A, with A being filled with random numbers. Also creating a vec b with random numbers.
	n = 5;
	gsl_matrix * A2 = gsl_matrix_alloc(n,n);
	gsl_matrix* R2 = gsl_matrix_alloc(n,n);
	gsl_vector* b = gsl_vector_alloc(n);
	for (int i = 0; i < n; i++){
		gsl_vector_set(b,i,randomBetweenPlusMinus1());
		for( int j = 0; j < n; j++){
			gsl_matrix_set(A2,i,j,randomBetweenPlusMinus1());
		}
	}
	//Computing GS_solve with empty vec, x
	GS_decomp(A2,R2);
	gsl_vector* x = gsl_vector_alloc(n);
	GS_solve(A2,R2,b,x);
	printf("\nPrinting vector x = \n");
	printVector(x);
	//Checking that Ax = b yields the same b.
	gsl_matrix * QR2 = gsl_matrix_alloc(n,n);
	gsl_vector * QRx = gsl_vector_alloc(n);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,A2,R2,0,QR2);
	gsl_blas_dgemv(CblasNoTrans,1,QR2,x,0,QRx);

	printf("\nChecking that Q*R*x. It should be equal to b.\n");
	printVector(QRx);
	printf("vector b is:\n");
	printVector(b);

	/*
	Part B
	*/
	printf("\nPart B: checking that the Gram-Schmidt QR factorization is working.\n");
	//Creating a square n x n matrix, A, with random inputs.
	gsl_matrix* A3 = gsl_matrix_alloc(n,n);
	gsl_matrix* R3 = gsl_matrix_alloc(n,n);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			gsl_matrix_set(A3,i,j,randomBetweenPlusMinus1());
		}
	}
	printf("Generated matrix, A:\n");
	printMatrix(A3);
	GS_decomp(A3,R3);
	gsl_matrix* B3 = gsl_matrix_alloc(n,n);
	GS_inverse(A3,R3,B3);

	//Checking that B3 is the inversed.
	gsl_matrix* QR3 = gsl_matrix_alloc(n,n);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A3, R3, 0, QR3);
	gsl_matrix* AB = gsl_matrix_alloc(n,n);
	gsl_matrix* BA = gsl_matrix_alloc(n,n);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, QR3, B3, 0, AB);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B3, QR3, 0, BA);
	printf("\nChecking that the matrices, BA = B*A = A*B = AB = I\n");
	printf("A*B =\n");
	printMatrix(AB);
	printf("B*A = \n");
	printMatrix(BA);


	//Clearing allocations.
	gsl_matrix_free(QR2);
	gsl_vector_free(QRx);
	gsl_matrix_free(A2);
	gsl_matrix_free(R2);
	gsl_vector_free(b);
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(A3);
	gsl_matrix_free(B3);
	gsl_matrix_free(R3);
	gsl_matrix_free(QR3);
	gsl_matrix_free(AB);
	gsl_matrix_free(BA);
	return 0;
}
