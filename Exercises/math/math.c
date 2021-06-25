#include <stdio.h>
#include <math.h>
#include <tgmath.h>
#include <complex.h>

/* 
TØ uge 2. Beregne bessel og gamma, samt undersøge matematiske funktioner
*/


#define M_E 2.7182818284 /* e */
#define M_pi 3.141592 /* pi */
void printDash(void){
	printf("-------------\n");
}
void print(char funName[], double fun, double realVal){
	printDash();
	printf("%s found value: %g\n", funName, fun);
	printf("%s real value: %g\n",funName, realVal);
}

void printcomplex(char funName[], double complex fun, double complex realVal){
	printDash();
	printf("%s found value: %g+i%g\n",funName,crealf(fun),cimagf(fun));
	printf("%s real value: %g+i%g\n",funName,crealf(realVal),cimagf(realVal));
}



int main(void){
	print("Gamma(5)",gamma(5),3.17);
	print("J1(0.5)",j1(0.5),0.5);
	printcomplex("sqrt(-2)",csqrt(-2),1.44*I);
	printcomplex("exp(i*pi)",cexp(I*M_pi),-1);
	printcomplex("exp(i)",cexp(I),-1);
	printcomplex("i^e",pow(I,M_E),-0.428219);
	printcomplex("i^i",pow(I,I),0.2078);
	return 0;
}


