#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>

//Defining functions made in SVD.c
void printM(gsl_matrix* A);
void SVD(gsl_matrix* A,gsl_matrix* U, gsl_matrix* D, gsl_matrix* V);

int main(){
	//Generating a random symmetric matrix, A
	int N = 5; //the dimensions of the square matrix A.

	//Allocating place for matrix
	gsl_matrix* A = gsl_matrix_alloc(N,N);
	gsl_matrix* V = gsl_matrix_alloc(N,N);
	gsl_matrix* D = gsl_matrix_alloc(N,N);
	gsl_matrix* U = gsl_matrix_alloc(N,N);
	gsl_matrix* UD = gsl_matrix_alloc(N,N);
	gsl_matrix* UDVT = gsl_matrix_alloc(N,N);
	gsl_matrix* UU = gsl_matrix_alloc(N,N);
	gsl_matrix* VV = gsl_matrix_alloc(N,N);

	//Generating A
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){

			double Aij = (double) rand()/RAND_MAX*10;
			gsl_matrix_set(A,i,j,Aij);
		}
	}


	//Setting V, D, U to identity matrices.
	gsl_matrix_set_identity(V);
  	gsl_matrix_set_identity(U);
  	gsl_matrix_set_identity(D);
	printf("The One-sided Jacobi algorithm for Singular Value Decomposition, exam question 8\n\n");
	printf("Generating random n x m (with m = n) matrix A with n = %d. We make sure in the SVD function that m = n.\n\nMatrix A:\n",N);
	printM(A);

	printf("Now running the one-sided Jacobi SVD algorithm.\n");
	SVD(A,U,D,V);
	printf("The one-sided Jacobi SVD algorithm decomposes A -> SDV^T.\nWe not print the three produced matrices:\n\n");
	printf("The matrix U, ui = a'i/||a'i||, with a' being the i'th column of A' = A*R where R is the accumulation of the individual rotations.\n");
	printM(U);
	printf("The matrix V = R.\n");
	printM(V);
	printf("The matrix D, Dii=||a'i||.\n");
	printM(D);

	printf("\nWe now check that A = UDV^T.\n");

	//Creating U*D*V^T
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,U,D,0.0,UD);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,UD,V,0.0,UDVT);
	printM(UDVT);
	printf("We can see that A is equal to UDVT.\n");


	//-------- testing orthogonality of U and V
	printf("\n\nAt last we test the orthogonality of VV^T and UU^T:\n");
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,V,V,0,VV);
	printf("V*V^T:\n");
	printM(VV);

	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,U,U,0,UU);
	printf("U*U^T:\n");
	printM(UU);
	printf("We can see that they're orthogonal as anticipated.\n");

	//Clearing allocations
	gsl_matrix_free(V);
	gsl_matrix_free(A);
	gsl_matrix_free(U);
	gsl_matrix_free(D);
	gsl_matrix_free(UD);
	gsl_matrix_free(UDVT);
	gsl_matrix_free(UU);
	gsl_matrix_free(VV);
	return 0;
}

