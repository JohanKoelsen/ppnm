//hello.c
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>

char hello[] = "file scope";

void bar(void){
	char hello[] = "function scope";
	printf("hello (function scope?):%s\n",hello);
	{
		char hello[] = "block scope";
		printf("hello (block scope?):%s\n",hello);
	}
}
	


int main() {
	printf("hello world\n");
	bar();



	return 0;
}
