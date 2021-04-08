#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mps.h"

int main() {
double rho_over_grid=1.0, temp_over_grid=1.0, ioncharge= 1.0;
double **dist_i_GK, *mu_i, *U_i;
double v_max, **dist_e_DK, *vpar_e, *mu_e;
int size_mu_i, size_U_i, size_mu_e, size_vpar_e;
int nrows_distfile = 0, ncols=0;
char line_million[1000000];
//nrows_distfile counts the number of rows in the input file with the distribution function; ncols_distfile does the same for columns; ncols_Umufile  finds the number of columns for both rows or Umufile.txt. The first rows contains values of U, the second of mu. The columns are the possible values; sizeUU and sizemumu stores the final number of columns */
clock_t begin_it = clock(); // Finds the start time of the computation
int convergence = 0, problem = 0, i;
/* convergence = 0 tells newguess.c to set smoothing for the next iteration, convergence = 1 tells it to not smooth, and convergence = 2 means that the convergence criterion has been satisfied and the iteration is finished */
int N_iterations = 1000, N;
/* Number of maximum iterations is set to some large number e.g. 100 but iterations should converge in about 20. If they don't, then there is a problem. */
double tot_time, v_cut, u_e;

double *ne_grid, *ni_grid, v_cutDS = 0.5;

//struct spline spl_phitone, spl_netophi, spl_new_phitone, spl_new_netophi;

double *storevals;
// Pointers used to read the files
double *x_grid, *phi_grid;
double Te, alpha;
int size_phigrid, size_ngrid = 0;

char line_hundred[100];

//getting some of the parameters
//fptr1 = fopen("dist_parameters.txt", "r");
//if (fptr1 == NULL)
//{
//	printf("Crap");
//}
//i = 0;
//while (fgets(line_hundred, 100, fptr1) != NULL)
//{
//	string_pot1 = malloc(strlen(line_hundred) * sizeof(char));
//	string_pot1 = line_hundred;
//	storevalspot1 = linetodata(string_pot1, strlen(string_pot1), &ncols);
//	if (i == 0)
//	{
//		v_max = *storevalspot1;
//	}
//	else
//	{
//		if (i == 1)
//		{
//			v_cut = *storevalspot1;
//		}
//	}
//	i++;
//}

/* READ INPUT FILE 
This contains the values of alpha and tau */

FILE *input;
if ((input = fopen("inputfile.txt", "r")) == NULL) { 
	printf("Cannot open %s\n", "inputfile.txt");
	exit(EXIT_FAILURE); 
}
i=0;
ncols = 0;
while (fgets(line_hundred, 100, input) != NULL) {	
	storevals = linetodata(line_hundred, strlen(line_hundred), &ncols);
	if (i==0)
		alpha = *storevals;
	else if (i==1)
		Te = *storevals; 
	else if (i==2)
		v_cut = *storevals; 
	i += 1;	
}
v_cutDS = v_cut;
fclose(input);
i=0;
printf("alpha = %f\n", alpha);
/* READ POTENTIAL FILE */
// Now we read the file containing the potential profile and the grid of values of g = sqrt(x) on which is is defined
FILE *fp;
if ((fp = fopen("phidata.txt", "r")) == NULL) {
	printf("Cannot open %s\n", "fp.txt"); 
	exit(EXIT_FAILURE);
}
/* The while loop counts the lines in the file to determine the size of the arrays to be created */
i=0;
while (fgets(line_hundred, 100, fp) != NULL)
	i += 1;
size_phigrid = i; // n is the array size
rewind(fp); // go back to beginning of file
phi_grid = malloc(size_phigrid*sizeof(double));
x_grid = malloc(size_phigrid*sizeof(double));
ne_grid = malloc(size_phigrid*sizeof(double));
ni_grid = malloc(size_phigrid*sizeof(double));
free(storevals);
i=0;
ncols = 0;
while (fgets(line_hundred, 100, fp) != NULL) {
	storevals = linetodata(line_hundred, strlen(line_hundred), &ncols);
	x_grid[i] = storevals[0]*storevals[0];
	phi_grid[i] = storevals[1];
	//printf("gg = %f\tphi = %f\n", storevals[0], phi[i]);
	i += 1;
}
fclose(fp); // Close file
//free(storevals);


/* This part of the code extracts the ion distribution function from a file. We import the distribution function into a 2 dimensional array, F(mu,U). */
FILE *distfile, *Umufile;
if ((distfile = fopen("distfuncin.txt", "r")) == NULL)
{	
	printf("cannot open file %s\n", "distfuncin.txt");
	exit(-1); 
}
/* Count the number of rows in the distfuncin.txt file */
while(fgets(line_million, 1000000, distfile) != NULL) {	
	nrows_distfile += 1; 
}
/* Allocate the right amount of memory to the distribution function pointer */
dist_i_GK = (double**) calloc(nrows_distfile,sizeof(double*));
/* The number of rows is also the the size of the array in U_i */
nrows_distfile = 0; // Set number of rows counter to zero again
fclose(distfile); // Close file
// Check file exists and can be opened
if ((distfile = fopen("distfuncin.txt", "r")) == NULL) {
	printf("cannot open file %s\n", "distfuncin.txt");
	exit(-1);
}
/* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and counting the columns. Once the columns are counted we allocate memory to dist_i_GK and assign each number to a different column of dist_i_GK (for a fixed row). The while loop then jumps to a new line and repeats the process */
while (fgets(line_million, 1000000, distfile) != NULL) {
	storevals = linetodata(line_million, strlen(line_million), &ncols);
	//dist_i_GK[nrows_distfile] = (double*)calloc(ncols,sizeof(double));
	//mu_i = (double*)calloc(ncols,sizeof(double));
	dist_i_GK[nrows_distfile] = storevals;
	nrows_distfile +=1; 
}
fclose(distfile);
printf("~~~~~The second element of dist_i_GK is %f~~~~~\n",dist_i_GK[0][1]);
/* Now we extract the two 1D arrays representing the grid points on which dist_i_GK (the distribution function array) is defined. */
// first open file and check for error
if ((Umufile = fopen("Umufile.txt", "r")) == NULL) {	
	printf("cannot open file %s\n", "Umufile.txt");
	exit(-1); 
}
/* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and assigning the values of the first line to mu_i, second to U_i */
i=0;
while (fgets(line_million, 1000000, Umufile) != NULL) {	
	storevals = linetodata(line_million, strlen(line_million), &ncols);
	if (i == 0) {
		size_mu_i = ncols;
		mu_i = storevals;
	}
	else {
		size_U_i = ncols;
		U_i = storevals;
	}
	i += 1; 
}
fclose(Umufile);


/* This part of the code extracts the electron distribution function from a file. We import the distribution function into a 2 dimensional array, F(mu,U). */
distfile = fopen("dist_file.txt", "r");
if (distfile == NULL)
{	
	printf("cannot open file %s\n", "dist_file.txt");
	exit(-1); 
}
/* Count the number of rows in the distfuncin.txt file */
while(fgets(line_million, 1000000, distfile) != NULL) {	
	nrows_distfile += 1; 
}
/* Allocate the right amount of memory to the distribution function pointer */
dist_e_DK = (double**) calloc(nrows_distfile,sizeof(double*));
/* The number of rows is also the the size of the array in mu_i */
nrows_distfile = 0; // Set number of rows counter to zero again
rewind(distfile); // rewind to first line of file
/* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and counting the columns. Once the columns are counted we allocate memory to dist_i_GK and assign each number to a different column of dist_i_GK (for a fixed row). The while loop then jumps to a new line and repeats the process */
while (fgets(line_million, 1000000, distfile) != NULL) {
	storevals = linetodata(line_million, strlen(line_million), &ncols);
	//dist_e_DK[nrows_distfile] = (double*)calloc(ncols,sizeof(double));
	dist_e_DK[nrows_distfile] = linetodata(line_million, strlen(line_million), &ncols);
	mu_e = (double*)calloc(ncols,sizeof(double));
	//dist_e_DK[nrows_distfile] = storevals;
	nrows_distfile +=1; 
}
fclose(distfile);
printf("~~~~~The second element of dist_e_DK is %f~~~~~\n",dist_e_DK[0][1]);
/* Now we extract the two 1D arrays representing the grid points on which dist_i_GK (the distribution function array) is defined. */
// first open file and check for error
if ((Umufile = fopen("dist_file_arguments.txt", "r")) == NULL) {	
	printf("cannot open file %s\n", "dist_file_arguments.txt");
	exit(-1); 
}
/* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and assigning the values of the first line to mu_e, second to vpar_e */
i=0;
while (fgets(line_million, 1000000, Umufile) != NULL) {	
	storevals = linetodata(line_million, strlen(line_million), &ncols);
	if (i == 0) {
		size_mu_e = ncols;
		mu_e = storevals;
	}
	else {
		size_vpar_e = ncols;
		vpar_e = storevals;
	}
	i += 1; 
}
fclose(Umufile);
v_max = vpar_e[size_vpar_e-1];



//// Now get electron distribution function
////getting some of the parameters
//fptr1 = fopen("dist_parameters.txt", "r");
//if (fptr1 == NULL) {
//	printf("Crap");
//}
//i = 0;
//while (fgets(line_hundred, 100, fptr1) != NULL) {
//	storevalspot1 = linetodata(line_hundred, strlen(line_hundred), &ncols);
//	if (i == 0)
//	{
//		v_max = *storevalspot1;
//	}
//	else
//	{
//		if (i == 1)
//		{
//			//v_cut = *storevalspot1;
//		}
//	}
//	i++;
//}
//
//
////Opening the distribution file to determine the the number of values of the distribution function
//fptr1 = fopen("dist_file.txt", "r");
//if (fptr1 == NULL) {
//	printf("Crap");
//}
////first reading through the file to find the values of the distribution function
//i = 0;
//while (fgets(line_hundred, 100, fptr1) != NULL)
//{
//	i++;
//}
//len_F = i;
//size_vpar_e = len_F;
//size_mu_e = 1;
//dist_e_DK = malloc(size_mu_e * sizeof(double));
//dist_e_DK = malloc(len_F * sizeof(double));
//rewind(fptr1); // rewind through the files
////reread throught the files now storing the values of the Distribution function
//i = 0;
//while (fgets(line_hundred, 100, fptr1) != NULL)
//{
//	storevalspot1 = linetodata(line_hundred, strlen(line_hundred), &ncols);
//	dist_e_DK[0][i] = *storevalspot1;
//	i++;
//	free(storevalspot1);
//}
//fclose(fptr1);

srand ( time ( NULL));
for (N=0;N<N_iterations;N++)
{	
	printf("Iteration number is %d\n", N);
	problem = densfinorb(0, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, rho_over_grid, temp_over_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i); 
	printf("size_phigrid = %d\n", size_phigrid);
	printf("size_ngrid = %d\n", size_ngrid);
	printf("vcut = %f\nphi[0] = %f\n", v_cut, phi_grid[0]);
	if (v_cut*v_cut > - 2.0*phi_grid[0])
		v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
	else v_cutDS = 0.001;
	makelookup(phi_grid, ne_grid, size_phigrid, v_cut, v_cutDS, &u_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e);
	convergence = newguess(convergence, problem, Te, alpha, x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, size_ngrid, u_e, &v_cut, 1, dist_i_GK, U_i, mu_i, size_U_i, size_mu_i);
	if (convergence == 2)
	{
		densfinorb(1, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, rho_over_grid, temp_over_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i); 
		//The argument of densfinorb set to one makes the module compute the ion distribution function at x=0
		clock_t end_it = clock(); // finds end time of last iteration
		tot_time = (double) (end_it - begin_it) / CLOCKS_PER_SEC;
		printf("At %dth iteration converged successfully in %f seconds\n", N, tot_time);
		exit(0); 
	} 
}
if (N_iterations == 0)
{	
	problem = densfinorb(1, Te, alpha, size_phigrid, &size_ngrid, ne_grid, x_grid, phi_grid, rho_over_grid, temp_over_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i);
}
printf("No convergence after %d iterations\n", N_iterations);
exit(0);
}
