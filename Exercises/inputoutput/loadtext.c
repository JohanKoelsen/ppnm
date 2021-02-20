#include <stdio.h>
#include <math.h>

int main(int argc,char** argv){

	double x;
	int items;
	if(argc < 3){
	printf("no argument, klaphat\n");
	return 0;}

	FILE* my_input_stream = fopen(argv[1],"r");
	FILE* my_out_stream = fopen(argv[2],"w");

	do  {
		items = fscanf(my_input_stream, "%lg",&x);
		if(items != EOF){
			printf("x = %lg, cos(x) = %lg\n",x, cos(x));
			fprintf(my_out_stream,"x=%lg, cos(x) = %lg\n",x,cos(x));
			}
	} while(items != EOF);
	fclose(my_input_stream);
	fclose(my_out_stream);

	return 0;

}
