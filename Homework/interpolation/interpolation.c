#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_interp.h>
#include<assert.h>
#include "my_h_file.h"

//The main file computes different interpolations from the my_h_file.h.


int main(void){
	printf("Part A: interpolates the list shown below, and calculates the the linear spline\n");
	int n = 7;
	double x_list[7] = {1,3,4,7,10,12,16};
	double y_list[7] = {2,3,5,7,8,11,12};
	double z = 6;
	double answer = linterp(n, x_list, y_list, z);
	printf("a.1 test = %f\n",answer);
	printf("a.2 test = %f\n",linterp_integ(x_list,y_list,5));


	return 0;
}
