#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>

//Defining functions made in Jd.c
void printMatrix(gsl_matrix* A);
void jacobi_diag(gsl_matrix* A, gsl_matrix* V);


int main(){
	//Generating a random symmetric matrix, A
	int N = 5; //the dimensions
	gsl_matrix* A = gsl_matrix_alloc(N,N);
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){

			double Aij = (double) rand()/RAND_MAX*10;
			gsl_matrix_set(A,i,j,Aij);
			gsl_matrix_set(A,j,i,Aij);
		}
	}
	//Printing A:
	printf("The matrix A:");
	printMatrix(A);

	//Creating new matrix V
	gsl_matrix* V = gsl_matrix_alloc(N,N);
	gsl_matrix_set_identity(V); //Setting V equal to an identity matrix. We are now able to use the jacobi_diag function.
	jacobi_diag(A,V);
	//After  applying the jacobi eigenvalue algorithm. A is now the matrix where the diagonal elements are the eigenvalues
	printf("The eigenvalues of matrix A:");
	printMatrix(A);

	printf("Calculing VDV^T to check that this is equal to the original matrix, A.\n");
	//Checking that A=VDV^T
	gsl_matrix* VD = gsl_matrix_alloc(N,N);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,A,0,VD);
	printf("VDV^T is equal to the original matrix A:");
	gsl_matrix* VDVT = gsl_matrix_alloc(N,N);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,VD,V,0,VDVT);
	printMatrix(VDVT);



	//Question B
	//Creating the hamiltonial matrix
	int n = 20;
	double s=1.0/(n+1);
	gsl_matrix* H = gsl_matrix_alloc(n,n);
	for(int i=0;i<n-1;i++){
  		gsl_matrix_set(H,i,i,-2);
  		gsl_matrix_set(H,i,i+1,1);
  		gsl_matrix_set(H,i+1,i,1);
	}
	gsl_matrix_set(H,n-1,n-1,-2);
	gsl_matrix_scale(H,-1/s/s);
	printf("Question B\n\n");
	printf("The Hamiltonian Matrix");
	printMatrix(H);
	//Now it is time to diagonalize it by our jacobi_diag function.
	gsl_matrix* VB = gsl_matrix_alloc(n,n);
	gsl_matrix_set_identity(VB);
	jacobi_diag(H,VB);
	printf("The energies are the following\n");
	//The energies are printed
	for (int k=0; k < n/3; k++){
 		double exact = M_PI*M_PI*(k+1)*(k+1);
		double calculated = gsl_matrix_get(H,k,k);
		printf("n = %i, calculated eigenvalues = %g, exact eigenvalues =  %g\n",k,calculated,exact);
	}

	//Comparing them to the analytical solution, which is sqrt(2/L) sin(m*pi*ksi)
	printf("At last we compare them to the analytical solution, which is sqrt(2/L) sin(m*pi*ksi)\n");
	double Coeff = sqrt(2.0/(n+1)); // L = n + 1
	for(int i = 0; i < 3; i++){
		for(int k = 0; k < n; k++){
			double ksi = (k + 1.0)/(n + 1);
			printf("Analytical solution: %g, our eigenfunction: %g\n",Coeff*sin((i+1)*M_PI*ksi),gsl_matrix_get(VB,k,i)*pow(-1,k));
		}
	}


	//Clearing vars
	gsl_matrix_free(A);
	gsl_matrix_free(VD);
	gsl_matrix_free(VDVT);
	gsl_matrix_free(VB);

	return 0;
}
