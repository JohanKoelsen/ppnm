//This code will solve the matrix equation A*x = b

#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>


void vec_print(char s[],gsl_vector* v){
	printf("%s\n",s);
	for(int i = 0;i <v->size;i++)printf("%f, ",gsl_vector_get(v,i));
	printf("\n");
}

void check_ans(char s[], gsl_matrix* M, gsl_vector* v){
	gsl_vector* y = gsl_vector_alloc(3);

	gsl_blas_dgemv(CblasNoTrans,1,M,v,0,y);

	vec_print(s,y);

	gsl_vector_free(y);
}

int main(){
	int n = 3;
	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* Acopy = gsl_matrix_alloc(n,n);
	gsl_vector* b = gsl_vector_alloc(n);
	gsl_vector* x = gsl_vector_alloc(n);

	//Inserting values in matrix A:
	double val[3][3] ={
	{6.13, -2.90, 5.86},
	{8.08, -6.31, -3.89},
	{ -4.36, 1.00, 0.19}
	};
	for(int i=0; i < A->size1; i++)
		for(int j=0; j<A->size2;j++){
			gsl_matrix_set(A,i,j,val[i][j]);
		}
	//Defining vector b:
	double vec_val[3] = {6.23, 5.37, 2.29};
	for(int i=0;i<b->size;i++){
		gsl_vector_set(b,i,vec_val[i]);
	}

	//solving for x
	gsl_matrix_memcpy(Acopy,A);
	gsl_linalg_HH_solve(Acopy,b,x); //løser Ax = b, hvor x er ukendt

	//Plotting answers
	vec_print("answer of equation - printing vector x in A*x = b:",x);

	//Checking A*x = b is correct
	check_ans("found b by A*x",A,x);
	vec_print("printing values of b",b);






gsl_matrix_free(A);
gsl_matrix_free(Acopy);
gsl_vector_free(b);
gsl_vector_free(x);
return 0;
}
