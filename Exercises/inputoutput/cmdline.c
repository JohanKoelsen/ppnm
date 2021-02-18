#include <stdio.h>
#include <stdlib.h>
#include <math.h>




int main(int argc, char** argv){
	printf("program name is %s\n",argv[0]);
	if (argc<2)printf("there were no arguments\n");
	else{
		for (int i=1;i<argc;i = i + 1){
			double x = atof(argv[i]);
			printf("x=%g, sin(x) = %g and cos(x) = %g\n",x,sin(x),cos(x));
		}
	}
	return 0;
}
