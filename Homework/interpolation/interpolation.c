

#include<assert.h>
#include<stdio.h>
#include<gsl/gsl_integration.h>
#include<math.h>
#include <gsl/gsl_vector.h>

/*
linterp returns an interpolated value at z from two lists using a linear spline.
linterp_integ returns the value of the integral of two lists, x and y from 0 to z.
*/

double linterp(gsl_vector* x, gsl_vector* y, double z){
	int n = x -> size;
	assert(n > 1 && z >= gsl_vector_get(x,0) && z <= gsl_vector_get(x,n-1));
	int i = 0;
	int j = n - 1;
	while(j - i > 1){
		int m = (i + j)/2;
		if(z > gsl_vector_get(x,m)) i = m;
		else j = m;
	}
	return gsl_vector_get(y,i) + (gsl_vector_get(y,i+1) - gsl_vector_get(y,i)) / (gsl_vector_get(x,i+1) - gsl_vector_get(x,i)) * (z - gsl_vector_get(x,i));

}


double linterp_integ(gsl_vector* x, gsl_vector* y, double z){
	assert(z >= gsl_vector_get(x,0));
	double integral = 0;
	int i = 0;
	while(z > gsl_vector_get(x,i + 1)){
		double sum = gsl_vector_get(y,i) * (gsl_vector_get(x,i + 1) - gsl_vector_get(x,i)) + 0.5 * (gsl_vector_get(y, i + 1) - gsl_vector_get(y,i)) * (gsl_vector_get(x,i+1) - gsl_vector_get(x,i));
		integral += sum;
		i++;
	}
	integral += gsl_vector_get(y,i) * (z - gsl_vector_get(x,i)) + 0.5* ((gsl_vector_get(y,i+1) - gsl_vector_get(y,i))) / (gsl_vector_get(x,i+1) - gsl_vector_get(x,i)) * (z - gsl_vector_get(x,i))*(z - gsl_vector_get(x,i));
	return integral;
}

int binsearch(gsl_vector* x, double z){/* locates the interval for z by bisection */ 
	int n = x -> size;
	assert(n>1 && gsl_vector_get(x,0) <= z && z <= gsl_vector_get(x,n-1));
	int i = 0;
	int j = n-1;
	while(j - i > 1){
		int mid = (i + j)/2;
		if(z > gsl_vector_get(x,mid)) i = mid; else j = mid;
		}
	return i;
	}

int binsearch_pointer(int n, double* x, double z){/* locates the interval for z by bisection */ 
	assert(n>1 && x[0]<=z && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
}

typedef struct {int n; double *x, *y, *b, *c;} qspline;
//Create qspline
qspline* qspline_alloc(int n,double* x,double* y){
	qspline* s = malloc(sizeof(qspline));//spline
	s -> b = malloc((n-1)*sizeof(double));  // b_i
	s -> c = malloc((n-1)*sizeof(double));  // c_i
	s -> x = malloc(n*sizeof(double));      // x_i
	s -> y = malloc(n*sizeof(double));      // y_i
	s -> n = n;
	for(int i=0;i<n;i++){
		s->x[i]=x[i];
		s->y[i]=y[i];
	}
	int i;
	double p[n-1];
	double h[n-1];
	for(i = 0; i < n - 1; i++){
		h[i] = x[i+1] - x[i];
		p[i] = (y[i+1] - y[i])/h[i];
	}
	s -> c[0] = 0;
	for(i = 0; i < n - 2; i++)
		s -> c[i+1] = (p[i+1] - p[i] - s->c[i] * h[i])/h[i+1];
	s->c[n-2]/=2;
	for(i = n-3; i >= 0; i--)
		s -> c[i] = (p[i+1] - p[i] - s->c[i+1] * h[i+1])/h[i];
	for(i = 0; i<n - 1; i++)
		s -> b[i] = p[i] - s->c[i]*h[i];
	return s;
}

double qspline_eval(qspline *s, double z){     //evaluates s(z)
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i = 0;
	int j = s->n - 1;                     //binary search: implemented explicitely here.
	while(j-i>1){
		int m = (i+j)/2;
		if(z > s-> x[m]) i = m;
		else j = m;
	}
	double h= z - s->x[i];
	return s->y[i]+h*(s->b[i]+h*s->c[i]);
}//inerpolating polynomial

void qspline_free(qspline *s){ //free the allocated memory
	free(s->x); free(s->y); free(s->b); free(s->c); free(s);
}


// Function that evaluates the derivative in point z
double qspline_der(qspline *s, double z){
	//Binary search
	int i=0, j=s->n-1;
		while(j-i>1){
		int m=(i+j)/2;
		if(z>s->x[m]) i=m;
		else j=m;
	}
	double h = z- s -> x[i];
	double bi = s -> b[i];
	double ci = s -> c[i];
	double qder = bi+2*ci*h;
	return qder;
	}

double qspline_integ(qspline *s, double z){
	//Binary search
	int i=0, j=s->n-1;
	while(j-i>1){
		int m=(i+j)/2;
		if(z>s->x[m]) i=m;
		else j=m;
	}
	double area = 0;
	for(int k=0; k<j; k++){ // summing over all intervals except the last one
		double h = s ->x[k+1] - s -> x[k];
		double yi = s -> y[k];
		double bi = s -> b[k];
		double ci = s -> c[k];
		area += yi*h+1./2*bi*pow(h,2)+1./3*ci*pow(h,3);
	}
	double hj = z - s -> x[j]; // and then the last one
	double yj = s -> y[j];
	double bj = s -> b[j];
	double cj = s -> c[j];
	area += yj*hj+1./2*bj*pow(hj,2)+1./3*cj*pow(hj,3);
	return area;
	}

