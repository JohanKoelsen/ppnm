#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

void printdash(){
	printf("----------\n");
}


void printMaxMinInt(char maxOrMin[], char loop[], int integer){
	printf("My %s integer for %s is %i.\n", maxOrMin, loop, integer);
}

int (*intOperation)(int x, int y);
int lessThan(int x,int y){
	return x < y;
}
int greaterThan(int x,int y){
	return x > y;
}

void calculateMaxMinInt(int posNeg){
	char maxMin[10];
	char intMaxMinName[10];
	int intMaxMin;
	if(posNeg == 1) {
		strcpy(maxMin,"max");
		strcpy(intMaxMinName, "INT_MAX");
		intMaxMin = INT_MAX;
		intOperation = greaterThan;
	} else if (posNeg ==-1){
		strcpy(maxMin,"min");
		strcpy(intMaxMinName,"INT_MIN");
		intOperation = lessThan;
	} else {
		exit(EXIT_FAILURE);
	}
	printMaxMinInt(maxMin,intMaxMinName,intMaxMin);
	int i;
	for(i = posNeg; intOperation(i+(posNeg*1),i); i = i + posNeg) {}
	printMaxMinInt(maxMin,"for loop",i);

	while(intOperation(i+(posNeg*1),i)) {i = i + posNeg;}
	printMaxMinInt(maxMin,"while loop",i);

	i=posNeg;
	do {i = i+ posNeg;
	} while (intOperation(i+(posNeg*1),i));
	printMaxMinInt(maxMin,"do while loop",i);
}


void calcsum(){
	int max = INT_MAX/3;
	int i;
	float sum_up_float = 0;
	float sum_down_float = 0;
	for(i = 1;i<=max;i = i + 1){
		sum_up_float = sum_up_float + 1.0f/i;
		sum_down_float = sum_down_float +1.0f/(max - i+1);
	}
	printdash();
	printf("sum_up_float is %f\n",sum_up_float);
	printf("sum_down_float is %f\n",sum_down_float);
	//There seems to be a difference between the two sums.

	double sum_up_double = 0;
	double sum_down_double = 0;
	for(i = 1;i<=max;i = i + 1){
		sum_up_double = sum_up_double + 1.0f/i;
		sum_down_double = sum_down_double + 1.0f/(max-i+1);
	}
	printdash();
	printf("sum_up_double is %g\n",sum_up_double);
	printf("sum_down_double if %g\n",sum_down_double);


}
void TheMachineEpsilon(){
	printdash();
	//while
	float x_float=1;
	while(1+x_float!=1){x_float/=2;}x_float*=2;
	printf("float x for while loop is %f\n",x_float);
	double x_double=1;
	while(1+x_double!=1){x_double/=2;} x_double*=2;
	printf("double x for while loop is %g\n",x_double);
	long double x_longdouble=1;
	while(1+x_longdouble!=1){x_longdouble/=2;}x_longdouble*=2;
	printf("long double x for while loop is %Lg\n",x_longdouble);
	
	//for loop
	float e_float; double e_double; long double e_longdouble;
	for(e_float=1;1+e_float!=1;e_float/=2){} e_float*=2;
	for(e_double=1;1+e_double!=1;e_double/=2){} e_double*=2;
	for(e_longdouble=1;1+e_longdouble!=1;e_longdouble/=2){ }e_longdouble*=2;

	printf("float e for for loop is %f\n",e_float);
	printf("double e for for loop is %g\n",e_double);
	printf("long double e for for loop is %Lg\n",e_longdouble);
	
	printdash();
	printf("float check %f\n",FLT_EPSILON);
	printf("double check %g\n",DBL_EPSILON);
	printf("long double check %Lg\n",LDBL_EPSILON);
	printdash();

	//while do
	float  y_f=1;double y_d=1; long double y_ld=1;
	do {y_f/=2;} while(1+y_f!=1); y_f*=2;
	do {y_d/=2;} while(1+y_d!=1); y_d*=2;
	do {y_ld/=2;} while(1+y_ld!=1); y_ld*=2;
	printf("float y for do while is %f\n",y_f);
	printf("double y for do while is %g\n",y_d);
	printf("long double y for do while is %Lg\n",y_ld);
}



int main(void){

	calculateMaxMinInt(1);
	calculateMaxMinInt(-1);
	calcsum();
	TheMachineEpsilon();
	return 0;
}
