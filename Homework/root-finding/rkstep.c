#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<assert.h>


// Calculates the norm of a vector v
double norm(gsl_vector* v) {
	double res=0;
	for(int i=0; i<v->size; i++) {
		res += gsl_vector_get(v,i)*gsl_vector_get(v,i);
	}
	return sqrt(res);
}

void rkstep12(void f(double t, gsl_vector* y, gsl_vector* dydt),double t, gsl_vector* yt, double h, gsl_vector* yh, gsl_vector* err){
	int n = yt->size;

	//Defining vectors
	gsl_vector* k0 = gsl_vector_alloc(n);
	gsl_vector* k12 = gsl_vector_alloc(n);
	gsl_vector* y_hold = gsl_vector_alloc(n);


	//Using the second order Runge-Kutta method in the rkstep

	f(t, yt, k0);

	//Copying elements of k0 onto y_hold
	gsl_vector_memcpy(y_hold, yt);

	//Adding k0 * h/2 to elements of y_hold. Such that y_hold = yt + h/2 * k0
	gsl_blas_daxpy(h/2,k0,y_hold);

	f(t + h/2, y_hold, k12);

	//Calculating yh = yt + h * k12. First by setting yh = yt, and then adding h*k12 to the elements
	gsl_vector_memcpy(yh,yt);
	gsl_blas_daxpy(h,k12,yh);

	// Computing uncertainty
	for(int i=0; i<n; i++) {
		double k0i = gsl_vector_get(k0,i);
		double k12i = gsl_vector_get(k12,i);
		gsl_vector_set(err,i, (k0i-k12i)*h/2);
	}


	//Freeing vectors
	gsl_vector_free(k0);
	gsl_vector_free(k12);
	gsl_vector_free(y_hold);
}


void driver(void f(double t, gsl_vector* b, gsl_vector* dydt), double a, gsl_vector* ya, double b, gsl_vector* yb, double h, double acc, double eps, char* path){
	//Defining vectors and values
	FILE* list = fopen(path, "w");
	assert(ya->size == yb->size && b > a && h < (b-a));
	int n = ya->size;
	double tau;
	double e;
	double x = a;

	//yh is the value each step and dy is the error for each step
	gsl_vector* yh 	= gsl_vector_alloc(n);
	gsl_vector* dy = gsl_vector_alloc(n);

	fprintf(list,"%20g ", x);
	for(int i = 0; i < n; i++){
		fprintf(list,"%20g ", gsl_vector_get(ya, i));
	}

	fprintf(list,"\n");

	while(x < b) {
		if(x + h > b) h = b-x;
		rkstep12(f, x, ya, h, yh, dy);
		tau = (eps*norm(ya) + acc) * sqrt(h/(b-a)); //calculate absolute tolerance tau

		e = norm(dy); //calculate local tolerance e

		if(e < tau) {
			x += h;
			gsl_vector_memcpy(ya,yh);


			fprintf(list,"%20g ", x); // Print to the path log
			for(int i=0; i<n; i++) {
				fprintf(list,"%20g" , gsl_vector_get(ya, i));
			}
			fprintf(list,"\n");

		}
		else if(e > 0) h *= pow(tau/e,0.25) * 0.95;
		else h *= 2;
	}
	gsl_vector_free(yh);
	gsl_vector_free(dy);
	fclose(list);
}
