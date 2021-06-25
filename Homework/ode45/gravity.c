#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_vector.h>

#define FMT "%9.3f %9.3f "
#define TIME 6
#define M1 1.0
#define M2 1.0
#define M3 1.0
#define GRAVITY 1


void f(double t,gsl_vector *u,gsl_vector *dudt){
	int n = u -> size;
	assert(n==12);
	double x1 = gsl_vector_get(u,0), y1 = gsl_vector_get(u,1);
	double x2 = gsl_vector_get(u,2), y2 = gsl_vector_get(u,3);
	double x3 = gsl_vector_get(u,4), y3 = gsl_vector_get(u,5);
	double vx1 = gsl_vector_get(u,6), vy1 = gsl_vector_get(u,7);
	double vx2 = gsl_vector_get(u,8), vy2 = gsl_vector_get(u,9);
	double vx3 = gsl_vector_get(u,10), vy3 = gsl_vector_get(u,11);
	double dx1dt = vx1, dy1dt = vy1;
	double dx2dt = vx2, dy2dt = vy2;
	double dx3dt = vx3, dy3dt = vy3;
	double dx12 = x2 - x1, dy12 = y2 - y1;
	double dx13 = x3 - x1, dy13 = y3 - y1;
	double dx23 = x3 - x2, dy23 = y3 - y2;
	double r12 = sqrt(dx12*dx12 + dy12*dy12);
	double r13 = sqrt(dx13*dx13 + dy13*dy13);
	double r23 = sqrt(dx23*dx23 + dy23*dy23);
	double f12 = M1*M2*GRAVITY/r12/r12;
	double f13 = M1*M3*GRAVITY/r13/r13;
	double f23 = M2*M3*GRAVITY/r23/r23;
	double dvx1dt = 1./M1*( f12*dx12/r12 + f13*dx13/r13);
	double dvy1dt = 1./M1*( f12*dy12/r12 + f13*dy13/r13);
	double dvx2dt = 1./M2*(-f12*dx12/r12 + f23*dx23/r23);
	double dvy2dt = 1./M2*(-f12*dy12/r12 + f23*dy23/r23);
	double dvx3dt = 1./M3*(-f13*dx13/r13 - f23*dx23/r23);
	double dvy3dt = 1./M3*(-f13*dy13/r13 - f23*dy23/r23);
	gsl_vector_set(dudt,0,dx1dt); gsl_vector_set(dudt,1,dy1dt);
	gsl_vector_set(dudt,2,dx2dt); gsl_vector_set(dudt,3,dy2dt);
	gsl_vector_set(dudt,4,dx3dt); gsl_vector_set(dudt,5,dy3dt);
	gsl_vector_set(dudt,6,dvx1dt); gsl_vector_set(dudt,7,dvy1dt);
	gsl_vector_set(dudt,8,dvx2dt); gsl_vector_set(dudt,9,dvy2dt);
	gsl_vector_set(dudt,10,dvx3dt);gsl_vector_set(dudt,11,dvy3dt);
	}
