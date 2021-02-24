
#include <math.h>
#include <stdio.h>
// from xs list this script calculates the Erf(x) values and prints them.


double Erf(double x){
	if (x<0) return -Erf(-x);
	double a[] = {0.254829592, -0.284496736,1.421413741, -1.453152027, 1.061405429};
	double t = 1/(1+0.3275911*x);
	double sum = t* (a[0]+t*(a[1]+t*(a[2] + t*(a[3] + t*a[4]))));
	return 1 - sum*exp(-x*x);
}




int main(){
	double xs[] = {0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
	int i;
	double ys;
	for (i=0;i<15;i++){
		ys = Erf(xs[i]);
		printf("%g, %g\n",ys,xs[i]);

}
	return 0;
}
