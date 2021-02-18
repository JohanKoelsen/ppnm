#include <stdio.h>
#include <math.h>


int equal(double a, double b, double tau, double epsilon){
	
	if (fabs(a-b) < tau || fabs(a-b)/(fabs(a)+fabs(b)) < epsilon/2) {return 1;}
	else {return 0;}

}

void printdash(){
	printf("--------------------------\n");
}


int main(void){
	printdash();
	printf("a = 1, b = 2, tau = 1, epsilon = 1/4 results in %i\n",equal(2,1,1,1/4));
	printf("a = 2,b = 1, tau = 2, epsilon = 3 results in %i\n",equal(1,2,2,3));

	return 0;
}
