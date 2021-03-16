

#include<assert.h>
#include<stdio.h>
#include<gsl/gsl_integration.h>
#include<math.h>


/*
linterp returns an interpolated value at z from two lists using a linear spline.
linterp_integ returns the value of the integral of two lists, x and y from 0 to z.
*/


double linterp(int n,double x[],double y[], double z){
	printf("The values are n = %d, x[3] = %f, y[4] = %f, z = %f\n",n,x[3],y[4],z);
	assert(n>1 && z>=x[0] && z<=x[n-1]);

	int i = 0, j = n - 1;

	while(j-i > 1){int m = (i+j)/2; if(z>x[m]) i=m; else j=m;}
	assert(x[i+j]>x[i]);

	return y[i] + (y[i+1] - y[i]) / (x[i+1]-x[i]) * (z-x[i]);

}


double linterp_integ(double x[], double y[], double z){
	assert(z >= x[0]);
	// defining integral value which the function returns.
	double integral = 0;
	int i = 0;
	while(z>x[i+1]){
		// Since it is p_i * (x - x_i)^2   = (y2 - y1)/(x2 - x1)  = (y[i+1] - y[i])*(x[i+1] - x[i])
		double sum = y[i] * (x[i + 1] - x[i]) + 0.5 * (y[i + 1] - y[i]) * (x[i + 1] - x[i]);
		integral += sum;
		i++;
	}
	integral += y[i] * (z - x[i]) + 0.5 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i])) * (z - x[i]) * (z - x[i]);
	return integral;
}
