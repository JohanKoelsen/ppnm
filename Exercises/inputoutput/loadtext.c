#include <stdio.h>
#include <math.h>

int main(void){

	double x;
	int items;
	FILE* my_input_stream = fopen("inputtext.txt","r");
	FILE* my_out_stream = fopen("my_out_stream.txt","w");
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
