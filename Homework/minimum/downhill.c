#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void reflection(gsl_vector* highest, gsl_vector* centroid, int d, gsl_vector* reflected){
	for(int i = 0; i < d; i++){
		gsl_vector_set(reflected,i, 2*gsl_vector_get(centroid,2) - gsl_vector_get(highest,i));
	}
}

void expansion(gsl_vector* highest, gsl_vector* centroid, int d, gsl_vector* expanded){
	for(int i = 0; i < d; i++){
		gsl_vector_set(expanded, i, 3*gsl_vector_get(centroid,i) - 2*gsl_vector_get(highest,i));
	}
}
void contraction(gsl_vector* highest, gsl_vector* centroid, int d, gsl_vector* contracted){
	for(int i = 0; i < d; i++){
		gsl_vector_set(contracted,i, 0.5*gsl_vector_get(centroid, i) + 0.5 * gsl_vector_get(highest,i));
	}
}

void reduction(gsl_matrix* simplex, int d, int lo){
	for(int k = 0; k < d + 1; k++){
		if(k != lo){
			for(int i = 0; i < d; i++){
				gsl_matrix_set(simplex, k, i, 0.5 * (gsl_matrix_get(simplex, k, i) + gsl_matrix_get(simplex, lo, i)));
			}
		}
	}

}

double distance(gsl_vector* a, gsl_vector* b, int d){
	double s = 0;
	for(int i = 0; i < d; i++){
		s += pow(gsl_vector_get(b,i) - gsl_vector_get(a,i), 2);
	}
	return sqrt(s);
}

double size(gsl_matrix* simplex, int d){
	double s = 0;
	for(int k = 0; k < d + 1;k++){
		gsl_vector* simplex0 = gsl_vector_alloc(simplex -> size1);
		gsl_vector* simplexk = gsl_vector_alloc(simplex -> size1);
		gsl_matrix_get_col(simplex0, simplex, 0);
		gsl_matrix_get_col(simplexk, simplex, k);
		double dist = distance(simplex0, simplexk,d);
		if(dist > s) s = dist;
		gsl_vector_free(simplex0);
		gsl_vector_free(simplexk);
	}
	return s;
}

void simplex_update(gsl_matrix* simplex, gsl_vector* f_values, int d, int* hi, int* lo, gsl_vector* centroid){
	*hi = 0;
	*lo = 0;
	double highest = gsl_vector_get(f_values,0);
	double lowest = gsl_vector_get(f_values,0);
	for(int k = 1; k < d + 1; k++){
		double next = gsl_vector_get(f_values, k);
		if(next > highest){
			highest = next;
			*hi = k;
		}
		if(next < lowest){
			lowest = next;
			*lo = k;
		}
	}
	for(int i = 0; i < d; i++){
		double s = 0;
		for(int k = 0; k < d + 1; k++)if(k != *hi) s+= gsl_matrix_get(simplex,k,i);
	gsl_vector_set(centroid, i, s/d);
	}
}


void simplex_initiate(double fun(gsl_vector*), gsl_matrix* simplex, gsl_vector* f_values, int d, int* hi, int* lo, gsl_vector* centroid){
	for(int k = 0; k < d + 1; k++){
		gsl_vector* simplexk = gsl_vector_alloc(simplex -> size1);
		gsl_matrix_get_col(simplexk, simplex, k);
		gsl_vector_set(f_values, k, fun(simplexk));
		gsl_vector_free(simplexk);
	}
	simplex_update(simplex, f_values, d, hi, lo, centroid);
}


int downhill_simplex(double F(gsl_vector*), gsl_matrix* simplex, int d, double simplex_size_goal){
	int hi, lo, k = 0;
	gsl_vector* simplexhi = gsl_vector_alloc(simplex -> size1);
	gsl_vector* centroid = gsl_vector_alloc(d);
	gsl_vector* F_value = gsl_vector_alloc(d + 1);
	gsl_vector* p1 = gsl_vector_alloc(d);
	gsl_vector* p2 = gsl_vector_alloc(d);

	simplex_initiate(F,simplex, F_value, d,&hi,&lo, centroid);

	while(size(simplex, d) > simplex_size_goal){
        simplex_update(simplex, F_value, d, &hi, &lo, centroid);

	gsl_matrix_get_col(simplexhi,simplex,hi);
        reflection(simplexhi, centroid, d, p1);
        double f_re = F(p1);

        if(f_re < gsl_vector_get(F_value,lo)){
            //  Reflection looks good, try expansion
		gsl_matrix_get_col(simplexhi, simplex, hi);
            expansion(simplexhi, centroid, d, p2);
            double f_ex = F(p2);

            if(f_ex < f_re){
                // Accept expansion
                for(int i = 0; i < d; ++i) {
                    gsl_matrix_set(simplex,hi,i, gsl_vector_get(p2,i));
                }
                gsl_vector_set(F_value, hi, f_ex);
            }
            else{
                // Reject expansion and accept reflection
                for(int i = 0; i < d; ++i){
                    gsl_matrix_set(simplex, hi, i, gsl_vector_get(p1,i));
                }
                gsl_vector_set(F_value, hi, f_re);
            }
        }
        else{
            // Reflection wasn’t good
            if(f_re < gsl_vector_get(F_value,hi)){
                // Ok, accept reflection
                for(int i = 0; i < d; ++i) {
			gsl_matrix_set(simplex,hi,i,gsl_vector_get(p1,i));
                }
		gsl_vector_set(F_value, hi, f_re);
            }
            else{
                // Try  contraction
		gsl_matrix_get_col(simplexhi, simplex,hi);
		contraction(simplexhi, centroid, hi, p1);
                double f_co = F(p1);

                if(f_co < gsl_vector_get(F_value,hi)){
                    // Accept contraction
                    for(int i = 0; i < d; ++i){
			gsl_matrix_set(simplex, hi, i, gsl_vector_get(p1,i));
                    }
		gsl_vector_set(F_value,hi, f_co);
                }
                else{
                    // Do reduction
                    reduction(simplex, d, lo);
                    simplex_initiate(F, simplex, F_value, d, &hi, &lo, centroid);
                }
            }
        }
        k++;
    	}
	return k ;

	//Cleaning
	gsl_vector_free(centroid);
	gsl_vector_free(F_value);
	gsl_vector_free(p1);
	gsl_vector_free(p2);
	gsl_vector_free(simplexhi);

}
