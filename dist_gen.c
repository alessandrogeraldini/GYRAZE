// Author: Robbie Ewart 
// MODIFIED: 15/08/2019
/* this code generates the distribution function (currently just a maxwellian) that will be read by a different file*/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
void main()
{
	FILE *fptr;
	//the pointer to the file that will contain all the distribution mallarce
	double *F, v_max, v_cut;
	v_max = 7.0;
	printf("enter the cut off velocity in electron thermal velocities: ");
	scanf("%lf", &v_cut);

	//an array that will contain the values of the distribution before they are written to the file
	int v_pts,i;
	printf("Enter the number of points in velocity space: ");
	scanf("%d", &v_pts);
	//specifies the maxium velocity (in thermal velocities) and the number of grid points between zero and the max velocity
	F = malloc(v_pts*sizeof(double)); //sets the number of points we want in F
	
	fptr = fopen("dist_file.txt","w");// open the distribution file in write mode
	
	if(fptr == NULL)
	{
		printf("Crap");
	}
	for(i=0; i<v_pts; ++i)
	{
		F[i] = exp(-0.5*pow(((i*v_max)/(v_pts-1)),2.0))/sqrt(2.0*M_PI);
	}
	for(i=0; i<v_pts; ++i)
	{
		fprintf(fptr,"%.15e\n",F[i]);
	}
	fclose(fptr);

	//stores relevant parameters for later use
	fptr = fopen("dist_parameters.txt", "w");
	if (fptr == NULL)
	{
		printf("Crap");
	}
	fprintf(fptr, "%f\n", v_max);
	fprintf(fptr, "%f\n", v_cut);
	fclose(fptr);
}
