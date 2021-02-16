#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "Headerfiles/iondens.h"
#include "Headerfiles/newguess.h"
#include "Headerfiles/cspline.h"
#include "Headerfiles/makelookup.h"
#include "Headerfiles/linetodata.h"
#include "Headerfiles/cutlookup.h"

//struct spline
//{
//	double* a;
//	double* b;
//	double* c;
//	double* d;
//
//};

int main() {
clock_t begin_it = clock(); // Finds the start time of the computation
int convergence = 0, problem = 0, p_size = 1000,i,cols = 0, cut_points = 10000; 
/* convergence = 0 tells newguess.c to set smoothing for the next iteration, convergence = 1 tells it to not smooth, and convergence = 2 means that the convergence criterion has been satisfied and the iteration is finished */
int N_iterations = 1000, N , use_new = 0;
/* Number of maximum iterations is set to some large number e.g. 100 but iterations should converge in about 20. If they don't, then there is a problem. */
double tot_time, v_max, v_cut;

double* ne_grid, * phi_grid, *ne_cut_grid, *phi_cut_grid, *ne_new_grid, *phi_new_grid;

//struct spline spl_phitone, spl_netophi, spl_new_phitone, spl_new_netophi;

char line_pot1[20];
char* string_pot1;
double* storevalspot1; //temporarily holds values to be moved around

FILE* fptr1; 

//getting some of the parameters
fptr1 = fopen("dist_parameters.txt", "r");
if (fptr1 == NULL)
{
	printf("Crap");
}
i = 0;
while (fgets(line_pot1, 20, fptr1) != NULL)
{
	printf("YOLO\n");
	string_pot1 = malloc(strlen(line_pot1) * sizeof(char));
	string_pot1 = line_pot1;
	storevalspot1 = linetodata(string_pot1, strlen(string_pot1), &cols);
	if (i == 0)
	{
		v_max = *storevalspot1;
	}
	else
	{
		if (i == 1)
		{
			v_cut = *storevalspot1;
		}
	}
	i++;
}
phi_grid = malloc(p_size * sizeof(double));
ne_grid = malloc(p_size * sizeof(double));
ne_new_grid = malloc(p_size * sizeof(double));
phi_new_grid = malloc(p_size * sizeof(double));
phi_cut_grid = malloc(cut_points * sizeof(double));
ne_cut_grid = malloc(cut_points * sizeof(double));

makelookup(phi_grid, ne_grid, p_size, v_cut);
//generate_cspline(&spl_netophi, ne_grid, phi_grid, p_size);
//generate_cspline(&spl_phitone, phi_grid, ne_grid, p_size);

cutlookup(phi_cut_grid, ne_cut_grid, cut_points);

//ne_new_grid = ne_grid;
//phi_new_grid = phi_grid;
//spl_new_netophi = spl_netophi;
//spl_new_phitone = spl_phitone;


srand ( time ( NULL));
for (N=0;N<N_iterations;N++)
{	
	printf("Iteration number is %d\n", N);
	if (use_new == 0)
	{
		problem = iondens(0, ne_grid, phi_grid, p_size); //The argument of iondens set to zero makes the module not compute the ion distribution function
	}
	else if (use_new == 1)
	{
		problem = iondens(0, ne_new_grid, phi_new_grid, p_size); //The argument of iondens set to zero makes the module not compute the ion distribution function
	}
	convergence = newguess(convergence, problem, ne_grid, phi_grid, p_size, ne_new_grid, phi_new_grid, ne_cut_grid, phi_cut_grid, cut_points,&use_new);
	if (convergence == 2)
	{
		if (use_new == 0)
		{
			iondens(1,ne_grid,phi_grid,p_size); //The argument of iondens set to one makes the module compute the ion distribution function at x=0
		}
		else if (use_new == 1)
		{
			iondens(1,ne_new_grid,phi_new_grid,p_size); //The argument of iondens set to one makes the module compute the ion distribution function at x=0
		}
		clock_t end_it = clock(); // finds end time of last iteration
		tot_time = (double) (end_it - begin_it) / CLOCKS_PER_SEC;
		printf("At %dth iteration converged successfully in %f seconds\n", N, tot_time);
		exit(0); 
	} 
}
if (N_iterations == 0)
{	
	problem = iondens(1, ne_grid, phi_grid, p_size);
}
printf("No convergence after %d iterations\n", N_iterations);
exit(0);
}
