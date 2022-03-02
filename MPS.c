//Authors: Alessandro Geraldini and Robbie Ewart 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <sys/stat.h>
#include "mps.h"
//#include "simparams.h"
#define STOP_MP 0.99
#define STOP_DS 0.99
#define INITIAL_GRID_PARAMETER 1.0
#define TEST_EL 0
#define MAX_IT 200
// number of maximum iterations set to some large number but iterations should converge in about 20-100. If they don't, then there is a problem
#define SYS_SIZ 20.0
#define ZOOM_DS 1
#define ZOOM_MP 3
#define GRIDSIZE_DS 0.1
#define GRIDSIZE_MP 0.3
#define WEIGHT 0.4
#define CAUTIOUSWEIGHT 0.1

const char *strqty[7] = {"alpha=","gamma=","nspec=","Ti:Te=","mi:me=","jwall=","pwall="};
const int lenstrqty = 7;

void remake_MPgrid(double *x_grid, double *phi_grid, int *psize_phigrid, double delta) {
	int i, phithenx=1, size_phigrid = *psize_phigrid;
	double xi, *new_phi, *new_x;
	//FILE *fp;
	new_phi = malloc(size_phigrid*sizeof(double));
	new_x = malloc(size_phigrid*sizeof(double));
	//fp = fopen("OUTPUT/phidata.txt", "w");
	//if (fp == NULL) {
	//	printf("Error opening new phidata.txt file\n");
	//	exit(EXIT_FAILURE);
	//}
	printf("in remake_MPgrid: size_phigrid = %d\n", size_phigrid);
	gsl_interp_accel *acc
	= gsl_interp_accel_alloc ();
	gsl_spline *spline
	= gsl_spline_alloc (gsl_interp_cspline, size_phigrid);

	gsl_spline_init (spline, phi_grid, x_grid, size_phigrid);

	new_phi[0] = phi_grid[0];
	new_x[0] = x_grid[0];
	//for (i=1; i < size_phigrid; i++) {
	i=1;
	while (new_x[i-1] < x_grid[size_phigrid-1] - delta) {
		if (phithenx == 1) {
			xi = i*delta;
			new_phi[i] = new_phi[i-1] + delta;
			new_x[i] = gsl_spline_eval (spline, new_phi[i], acc);
		}
		if ( (phithenx == 0) || (new_x[i] - new_x[i-1] > delta) ) {
			new_x[i] = new_x[i-1] + delta;
			if (phithenx == 1)  {
				gsl_spline_free (spline);
				gsl_interp_accel_free (acc);

				gsl_interp_accel *acc
				= gsl_interp_accel_alloc ();
				gsl_spline *spline
				= gsl_spline_alloc (gsl_interp_cspline, size_phigrid);
				gsl_spline_init (spline, x_grid, phi_grid, size_phigrid);
				phithenx = 0;
			}
		printf(" i = %d\n", i);
			new_phi[i] = gsl_spline_eval (spline, new_x[i], acc);
		printf(" i = %d\n", i);
			printf("x phi = (%f %f)\n", new_x[i], new_phi[i]);
		printf(" i = %d\n", i);
		}
		//g = sqrt(new_x[i]);
		printf("x phi = (%f %f)\n", new_x[i], new_phi[i]);
		//fprintf(fp, "%f %f\n", g, new_phi[i]);
		i++;
	}
	*psize_phigrid = i;
	for (i=0; i < size_phigrid; i++) {
		phi_grid[i] = new_phi[i];
		x_grid[i] = new_x[i];
	}
	free(new_phi);
	free(new_x);
	return;
}

void make_phigrid(double *x_grid, double *phi_grid, int size_phigrid, double grid_parameter, double deltax, int initial, double phi_jump, double len_scale) {
	int i;
	double g, xi, *ff, *new_phi, *new_x, power = 0.5;
	//FILE *fp;
	new_phi = malloc(size_phigrid*sizeof(double));
	new_x = malloc(size_phigrid*sizeof(double));
  	ff = malloc(size_phigrid*sizeof(double));
	//fp = fopen("OUTPUT/phidata.txt", "w");
	//if (fp == NULL) {
	//	printf("Error opening new phidata.txt file\n");
	//	exit(EXIT_FAILURE);
	//}
	printf("in make_phigrid: size_phigrid = %d\n", size_phigrid);
	if (initial != 0) {
			printf("grid_parameter = %f\n", grid_parameter);
		for (i=0; i<size_phigrid; i++) {
		    //ff[i] = pow(sqrt(grid_parameter)+sqrt(x_grid[i]), 2.0) - grid_parameter;
		    ff[i] = pow(pow(grid_parameter, power)+pow(x_grid[i], power), 1.0/power) - grid_parameter;
			//printf("ff[%d] = %f\n", i, ff[i]);
		}
		gsl_interp_accel *acc
		= gsl_interp_accel_alloc ();
		gsl_spline *spline
		= gsl_spline_alloc (gsl_interp_cspline, size_phigrid);

		gsl_spline_init (spline, ff, phi_grid, size_phigrid);

		for (i=0; i < size_phigrid; i++) {
			xi = i*deltax;
			new_phi[i] = phi_grid[i];
			new_x[i] = pow( pow(grid_parameter+xi, power) - pow(grid_parameter, power), 1.0/power);
			g = sqrt(new_x[i]);
			//fprintf(fp, "%f %f\n", g, new_phi[i]);
		}
		gsl_spline_free (spline);
		gsl_interp_accel_free (acc);
		for (i=0; i < size_phigrid; i++) {
			phi_grid[i] = new_phi[i];
			x_grid[i] = new_x[i];
		}
	}
	else {
		for (i=0; i < size_phigrid; i++) {
			xi = i*deltax;
			ff[i] = xi;
			//new_x[i] = pow( pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), 2.0);
			new_x[i] = pow( pow(grid_parameter+xi, power) - pow(grid_parameter, power), 2.0);
			g = sqrt(new_x[i]);
			new_phi[i] = phi_jump*pow(len_scale, 2.0)/pow(new_x[i] + len_scale, 2.0);
			//if (new_x[i] < 10.0/sqrt(2)) 
			//	new_phi[i] = 3.0*pow((new_x[i]*sqrt(2.0)/10.0 - 1.0), 5.0);
			//else 
			//	new_phi[i] = 0.0;
			//fprintf(fp, "%f %f\n", g, new_phi[i]);
		}
		for (i=0; i < size_phigrid; i++) {
			phi_grid[i] = new_phi[i];
			x_grid[i] = new_x[i];
		}
	}
	free(new_phi);
	free(new_x);
	free(ff);
	//fclose(fp);
	return;
}

// this function is no longer used, but keep it just in case
void rescale_array(double *array, int size_array, double jump) {
	int i;
	double old_jump = array[0];
	for (i=0; i < size_array; i++) {
		array[i] *= (jump/old_jump);
	}
	return;
}

/* the set of functions below are necessary to model electron reflection from the infinitely thin Debye sheath obtain parallel velocity cutoff as a function of magnetic moment with assumption rho_e >> lambda_D
*/
double FF(double beta) {
	double integrand, integral=0.0, beta_s = 0.005, betaval;
	int ind, num_beta = (int) (beta/beta_s);
	for (ind = 0; ind < num_beta+1; ind++) {
		betaval = beta*ind/num_beta;
		integrand = 0.5*(1.0-cos(2.0*betaval)) / ( M_PI - betaval + 0.5*sin(2.0*betaval) );
		integral += integrand*beta/num_beta;
	}
        return integral;
}

double FFp(double beta) {
	double integrand = 0.5*(1.0-cos(2.0*beta)) / ( M_PI - beta + 0.5*sin(2.0*beta) );
        return integrand;
}

double vparcut(double beta, double vcut) {
	double vparcutval = sqrt((1.0-exp(-2.0*FF(beta))))*vcut/sin(beta);
	return vparcutval;
}

double mucut(double beta, double vcut) {
	double mucutval = 0.5*exp(-2.0*FF(beta))*pow(vcut/sin(beta), 2.0);
	return mucutval;
}

// The main function of MAGSHEATH
int main() {

double power = 1.0;
clock_t begin_it = clock(); // Finds the start time of the computation
int zoomfactor_DS;
int size_mu_i, size_U_i, size_mu_e, size_vpar_e, number_species;
int i, j, nrows_distfile = 0, ncols=0, ndirname, fix_current=0, ind; //s
double error_current = 0.0005;
double deltaphi=0.0, weight_j=1.0, weight=WEIGHT;
double ioncharge= 1.0, *TioverTe;
double **dist_i_GK, *mu_i, *U_i, zero = 0.0;
double **dist_e_DK, **dist_e_GK, *vpar_e, *mu_e, *U_e_DS, *mu_e_DS;
double *fx_i_DS, *vx_i_DS, **F_i_DS;
double dvpar = 0.2/sqrt(2.0), Uminmu_MPE;
char line_million[1000000], dirname[200];
//nrows_distfile counts the number of rows in the input file with the distribution function; ncols_distfile does the same for columns; sizeUU and sizemumu stores the final number of columns */
int convergence = 0, convergence_DS = 0, convergence_j = 0;
/* convergence = 0 tells newguess.c to set smoothing for the next iteration, convergence = 1 tells it to not smooth, and convergence = 2 means that the convergence criterion has been satisfied and the iteration is finished */
int N=0, N_DS=0;
double tot_time, v_cut , current, flux_e=0.0, flux_i=0.0, flux_i_DS=0.0, flux_eDS = 0.0, mioverme;
double gamma_ref;
int size_cut=100, sizevxopen;
double *mue_cut_lookup, *vpar_cut_lookup;//, *vpar_cut_lookup_old;
double *ne_grid, *ni_grid, v_cutDS = 0.5, *ne_DSgrid;
double target_current;
double *storevals;
// Pointers used to read the files
double *x_grid, *phi_grid, *phi_DSgrid, *x_DSgrid, *ni_DSgrid;
double Te, alpha, alpha_deg, deltaxDS, factor_small_grid_parameter=2.0, system_size, DS_size;
int size_phigrid, size_ngrid = 0, size_neDSgrid = 0, size_phiDSgrid;
char line_hundred[100];
double grid_parameter, deltax;

//gsl_permutation *p;
//gsl_matrix *m;
//	const char *strqty[6];
//	strqty[0] = "alpha";
//	//const char strgamma = "gamma";
//	//const char strnspec = "nspec";
//	//const char strtempi = "tempi";
//	//const char strmassi = "tempi";
//	//const char strfixj = "fixj";
//	//const char strjorphi = "value";

// set the parameter g that controls the grid spacing in x 
// transition from evenly spaced in sqrt(x) near x=0 to evenly spaced in x at x=infty
grid_parameter = INITIAL_GRID_PARAMETER;
// set the grid spacing
deltax = GRIDSIZE_MP;

mkdir("OUTPUT", S_IRWXU);
FILE *fout = fopen("OUTPUT/output_MAGSHEATH.txt", "w");
if (fout == NULL) {
	printf("Problem opening output file for writing\n");
	exit(1);
}
/* READ INPUT FILE 
Contains the values of: 
     alpha
     number of ion species
     Ti/Te for all species (line to array)
     mi/me for all ion species (line to array)
     fix_current = 1 or 0 (fixes potential),
     fixed current/potential 
*/
// Note: multiple species not yet implemented
FILE *input;
if ((input = fopen("inputfile.txt", "r")) == NULL) { 
	printf("Cannot open %s\n", "inputfile.txt");
	fprintf(fout, "Cannot open %s\n", "inputfile.txt");
	exit(EXIT_FAILURE); 
}

i=0; // for counting rows (lines) of file
j=0;
ncols = 0; // for counting numbers in each row (line) of file
ndirname=0;
char *outputstr = "OUTPUT/";
for (ndirname=0; ndirname<strlen(outputstr);ndirname++) 
	dirname[ndirname] = outputstr[ndirname];
while (fgets(line_hundred, 100, input) != NULL) {	
	storevals = linetodata(line_hundred, strlen(line_hundred), &ncols);
	//printf("ndirname = %d\n", ndirname);
	if (i!= 5) {
		for (j=0; j < lenstrqty-1; j++) {
			dirname[ndirname+j] = strqty[i-fix_current][j];
		}
		ndirname += (lenstrqty -1);
		//dirname[0] = 'a'; dirname[1] = 'l'; dirname[2] 
		for (j=0;j<strlen(line_hundred)-1;j++) {
			//printf("%c\n", line_hundred[j]);
			dirname[j+ndirname] = line_hundred[j];
			//printf("%s\n", dirname);
		}
	}
	if (i==0)
		// input in degrees only for convenience
		// should be below 5 degrees (0.1 rad) for asymptotic theory to be valid
		alpha_deg = *storevals; 
	else if (i==1)
		gamma_ref = *storevals; 
	else if (i==2)
		number_species = (int) (*storevals); 
	else if (i==3)
		TioverTe = linetodata(line_hundred, strlen(line_hundred), &ncols); 
	else if (i==4)
		mioverme = *storevals; 
	else if (i==5)
		fix_current = (int) (*storevals); 
	else if (i==6) {
		if (fix_current != 0) {
			target_current = current = *storevals; 
			v_cut = 3.01; // set to a reasonable value
			// we will set v_cut to some value later
		}
		else {
			v_cut = *storevals; 
			target_current = current = 0.0; // gets calculated afterwords
		}
	}
	if (i<5) {
		ndirname += strlen(line_hundred)-1;
		dirname[ndirname] = '_';
		ndirname += 1;
	}
	i += 1; // count the rows in the file
	//printf("dirname = %s\n", dirname);
}
fclose(input);
printf("directory where output will be stored is %s\n", dirname);
mkdir(dirname, S_IRWXU);
Te = 1/TioverTe[0]; // Te/Ti is used and labelled Te in code
alpha = alpha_deg*M_PI/180; // alpha used in radians from now on
// initial iteration assumes flat potential profile in magnetic presheath
// therefore, the parameter v_cutDS is equal to v_cut
v_cutDS = v_cut;
if (gamma_ref > 1.0)
	deltaxDS = GRIDSIZE_DS/gamma_ref;
else
	deltaxDS = GRIDSIZE_DS;
printf("deltaxDS = %f\n\n\n", deltaxDS);
if (Te < 1.0) system_size = SYS_SIZ;
else system_size = SYS_SIZ*sqrt(Te);

if (gamma_ref < TINY) DS_size = SYS_SIZ;
else if (gamma_ref < 1.0) DS_size = SYS_SIZ/gamma_ref;
else DS_size = SYS_SIZ;
i=0;
printf("INPUT PARAMETERS:\n\tmagnetic field angle = α = %f deg (%f rad)\n\telectron gyroradius over Debye length (reference value at MP entrance) = ρ_e/λ_D = %f\n\tnumber of species = %d\n", alpha_deg, alpha, gamma_ref, number_species);
fprintf(fout, "INPUT PARAMETERS:\n\tmagnetic field angle = α = %f deg (%f rad)\n\telectron gyroradius over Debye length  (reference value at MP entrance) = ρ_e/λ_D = %f\n\tnumber of species = %d\n", alpha_deg, alpha, gamma_ref, number_species);
for (i=0; i<number_species; i++) {
	printf("\tion temperature = τ = ZT_i/T_e for species %d = %f\n\tmass ratio = m_i/m_e for species ? = %f\n", i+1, TioverTe[i], mioverme);
	fprintf(fout, "\tion temperature = τ = ZT_i/T_e for species %d = %f\n\tmass ratio = m_i/m_e for species ? = %f\n", i+1, TioverTe[i], mioverme);
}
printf("\tfix_current = %d\n", fix_current);
fprintf(fout, "\tfix_current = %d\n", fix_current);
if (fix_current != 0) {
	printf("\tcurrent = j/(e*n_MPE*v_t,i) = %f\n", target_current);
	fprintf(fout, "\tcurrent = j/(e*n_MPE*v_t,i) = %f\n", target_current);
}
else {
	printf("\twall potential = eφ_W/T_e = %f\n", -0.5*v_cut*v_cut);
	fprintf(fout, "\twall potential = eφ_W/T_e = %f\n", -0.5*v_cut*v_cut);
}

printf("magnetic presheath size in simulation (in ion gyroradii) = %f\n", system_size);
fprintf(fout, "magnetic presheath size in simulation (in ion gyroradii) = %f\n", system_size);
size_phigrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size), 2.0) - grid_parameter ) / deltax );
size_phiDSgrid = (int) ( DS_size/deltaxDS );
printf("size phigrid = %d\n", size_phiDSgrid);

zoomfactor_DS = ZOOM_DS;
printf("zoomfactor_DS = %d\n", zoomfactor_DS);
 
printf("size of coarse potential grid in Debye sheath = %d\n", size_phiDSgrid);
fprintf(fout, "size of coarse potential grid in Debye sheath = %d\n", size_phiDSgrid);
printf("size of coarse potential grid in magnetic presheath = %d\n", size_phigrid);
fprintf(fout, "size of coarse potential grid in magnetic presheath = %d\n", size_phigrid);

// form x, phi, ne and ni grids for magnetic presheath
x_grid = malloc(size_phigrid*sizeof(double));
phi_grid = malloc(size_phigrid*sizeof(double));
ne_grid = malloc(size_phigrid*sizeof(double));
ni_grid = malloc(size_phigrid*sizeof(double));

// form x, phi, ne and ni grids for Debye sheath
x_DSgrid = malloc(size_phiDSgrid*sizeof(double));
phi_DSgrid = malloc(size_phiDSgrid*sizeof(double));
ne_DSgrid = malloc(size_phiDSgrid*sizeof(double));
ni_DSgrid = malloc(size_phiDSgrid*sizeof(double));

/*
EXTRACT ION DISTRIBUTION FUNCTION
import the distribution function into a 2 dimensional array, F(mu,U)
*/
FILE *file;
if ((file = fopen("Fi_mpe.txt", "r")) == NULL)
{	
	printf("cannot open file %s\n", "Fi_mpe.txt");
	fprintf(fout, "cannot open file %s\n", "Fi_mpe.txt");
	exit(-1); 
}
/* Count the number of rows in the distfuncin.txt file */
while(fgets(line_million, 1000000, file) != NULL) {	
	nrows_distfile += 1; 
}
// Allocate the right amount of memory to the distribution function pointer
dist_i_GK = (double**) calloc(nrows_distfile,sizeof(double*));
// The number of rows is also the the size of the array in U_i 
nrows_distfile = 0; // Set number of rows counter to zero again
rewind(file); // rewind file
// Read each line of file, extract the data (ie the numbers in the line) and count the columns
// Once columns are counted, allocate memory to dist_i_GK and assign each number to a different column of dist_i_GK (for a fixed row)
while (fgets(line_million, 1000000, file) != NULL) {
	storevals = linetodata(line_million, strlen(line_million), &ncols);
	dist_i_GK[nrows_distfile] = storevals;
	nrows_distfile +=1; 
}
fclose(file);
//printf("~~~~~The second element of dist_i_GK is %f~~~~~\n",dist_i_GK[0][1]);
// Extract two 1D arrays representing the grid points on which dist_i_GK is defined
// first open file and check for error
if ((file = fopen("Fi_mpe_args.txt", "r")) == NULL) {	
	printf("cannot open file %s\n", "Fi_mpe_args.txt");
	fprintf(fout, "cannot open file %s\n", "Fi_mpe_args.txt");
	exit(-1); 
}
// Read each line of file, extract the data and assign the values of the first line to mu_i, second to U_i
i=0;
while (fgets(line_million, 1000000, file) != NULL) {	
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
fclose(file);

/*
EXTRACT ELECTRON DISTRIBUTION FUNCTION
Import the distribution function into a 2 dimensional array, F(mu,vpar). 
*/
file = fopen("Fe_mpe.txt", "r");
if (file == NULL)
{	
	printf("Cannot open file Fe_mpe.txt");
	fprintf(fout, "Cannot open file Fe_mpe.txt");
	exit(-1); 
}
/* Count the number of rows in the distfuncin.txt file */
while(fgets(line_million, 1000000, file) != NULL) {	
	nrows_distfile += 1; 
}
/* Allocate the right amount of memory to the distribution function pointer */
dist_e_DK = (double**) calloc(nrows_distfile,sizeof(double*));
dist_e_GK = (double**) calloc(nrows_distfile,sizeof(double*));
nrows_distfile = 0; // Set number of rows counter to zero again
rewind(file); // rewind to first line of file
while (fgets(line_million, 1000000, file) != NULL) {
	// assign each number to a different column of dist_e_DK (for a fixed row)
	dist_e_DK[nrows_distfile] = linetodata(line_million, strlen(line_million), &ncols);
	// allocate memory to dist_e_GK which might be used to solve Debye sheath
	// dist_e_GK depends on magnetic presheath and Debye sheath solutions
	// will assign different values to dist_e_GK at each iteration
	dist_e_GK[nrows_distfile] = (double*) calloc(ncols,sizeof(double*));
	nrows_distfile +=1; 
	// jump to a new line and repeat the process
}
fclose(file);
//printf("~~~~~The second element of dist_e_DK is %f~~~~~\n",dist_e_DK[0][1]);
// Extract the two 1D arrays representing the grid points on which dist_i_GK is defined
// first open file and check for error
if ((file = fopen("Fe_mpe_args.txt", "r")) == NULL) {	
	printf("cannot open file Fe_mpe_args.txt");
	fprintf(fout, "Cannot open file Fe_mpe_args.txt");
	exit(-1); 
}
// Read each line of the file, extract the data (ie numbers in a line) 
// Assign values of first line to mu_e, second to vpar_e
i=0;
while (fgets(line_million, 1000000, file) != NULL) {	
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
printf("size_mu_e = %d\nsize_vpar_e=%d\n", size_mu_e, size_vpar_e);
fclose(file);
F_i_DS = calloc(1,sizeof(double*));
sizevxopen = (int) 100*sqrt((1.0+Te)*(1.0+1.0/Te));
vx_i_DS = malloc(sizevxopen*sizeof(double));
fx_i_DS = malloc(sizevxopen*sizeof(double));
U_e_DS  = malloc(size_vpar_e*sizeof(double));
mu_e_DS  = malloc(size_mu_e*sizeof(double));

mue_cut_lookup = malloc((1+size_cut)*sizeof(double));
vpar_cut_lookup = malloc((1+size_cut)*sizeof(double)); 
//vpar_cut_lookup_old = malloc((1+size_cut)*sizeof(double)); 
for (i=0; i< size_cut+1; i++) 
	mue_cut_lookup[i] = i*(mu_e[size_mu_e-1]+TINY)/size_cut;

/*
CARRY OUT A SINGLE ION DENSITY CALCULATION
*/
if (MAX_IT == 0)
{	
	make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, 0, 0.0, 1.0);
	//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 5.0), 2.0) - grid_parameter ) / deltax );
	densfinorb(1, sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ne_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i, ZOOM_MP, STOP_MP);
}

if (TEST_EL==1) {
	//size_phiDSgrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(2.0*system_size), 2.0) - 0.3 ) / deltax );
	//printf("size_phiDSgrid = %d\n", size_phiDSgrid);
	F_i_DS[0] = malloc(sizevxopen*sizeof(double));
	make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, 0.0, deltaxDS, 0, log(0.001), 0.01);
	make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, 0, 0.0, 1.0);
	//phi_DSgrid[0] = -5.0;
	printf("v_cutDS = %f\n", v_cutDS);
	printf("now evaluate electron density in MPS\n");
	v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]) ;
	for (i=0; i< size_cut+1; i++) {
		vpar_cut_lookup[i] = v_cutDS;
		mue_cut_lookup[i] = i*(mu_e[size_mu_e-1]+TINY)/size_cut;
		//printf("%f %f\n", vpar_cut_lookup[i], mue_cut_lookup[i]);
	}
	printf("v_cutDS = %f\n", v_cutDS);
	denszeroorb(-1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, mue_cut_lookup, vpar_cut_lookup, size_cut);
	//v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);

	printf("now evaluate ion density in MPS\n");
	//size_ngrid = 0;
	//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 5.0), 2.0) - grid_parameter ) / deltax );
	densfinorb(1, sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i, ZOOM_MP, STOP_MP); 
	for (ncols=0; ncols<size_mu_e; ncols+=1) {
		for (ind=0; ind<size_vpar_e; ind+=1) {
			U_e_DS[ind] = ind*dvpar*ind*dvpar;
			Uminmu_MPE = sqrt(2.0*U_e_DS[ind] - 2.0*phi_grid[0]);
			dist_e_GK[ncols][ind] = 0.5*bilin_interp(mu_e[ncols], Uminmu_MPE, dist_e_DK, mu_e, vpar_e, size_mu_e, size_vpar_e, -1, -1)/ne_grid[0];
			//printf("%f ", dist_e_GK[ncols][ind]);
		}
		//printf("\n");
	}
	printf("size grid = %d\n", size_phiDSgrid);
	printf("size_neDSgrid = %d\n", size_neDSgrid);
	densfinorb(2, size_cut+1, 1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e, U_e_DS, size_mu_e, size_vpar_e, 0.0, mue_cut_lookup, vpar_cut_lookup, &flux_eDS, ZOOM_DS, 100.05); 

	exit(0);
}



N=0;
//ITERATE MAGNETIC PRESHEATH AND DEBYE SHEATH POTENTIAL TO FIND SOLUTION
//if ( (gamma_ref > 2.0) || (gamma_ref < 0.5) ) {
//if ( (gamma_ref > 10.0) && (gamma_ref < 0.1) ) 
// SIMPLIFIED ELECTRON MODEL USED TO CALCULATE DEBYE SHEATH POTENTIAL DROP
	while ( ( (convergence <= 1) || ( convergence_j <= 1 && fix_current == 1 ) ) && (N<MAX_IT) ) {
		printf("ITERATION # = %d\n", N);
		fprintf(fout, "ITERATION # = %d\n", N);
		current = target_current;

		if ( (factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut) < grid_parameter) )  {// && (phi_grid[0] + 0.5*v_cut*v_cut >1.0e-13 ) )   
		//	grid_parameter = 0.5;
			grid_parameter = factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut);
			weight = CAUTIOUSWEIGHT;
			//weight_j = CAUTIOUSWEIGHT;
		}
		//if (N != 0)
		//power = log((phi_grid[2] - phi_grid[0])/(phi_grid[1] - phi_grid[0]))/log(x_grid[2]/x_grid[1]);
		//if (power < 0.5) grid_parameter = 1.0;
		//else if (power > 1.0) grid_parameter = 0.0;
		//else grid_parameter = 2.0*pow((1.0-power), 2.0); 
		//if ( (phi_grid[2] - phi_grid[1])/(phi_grid[1] - phi_grid[0]) > 1.0 )
		//	grid_parameter -= 0.1;
		//else if ( (phi_grid[2] - phi_grid[1])/(phi_grid[1] - phi_grid[0]) < 0.95 )
		//	grid_parameter += 0.1;

		//if (grid_parameter <= 0.0001) grid_parameter = 0.0001;
		
		// MAKE ELECTROSTATIC POTENTIAL GRID
		//if (N==0) 
		make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, 0.0, 1.0);
		//else remake_MPgrid(x_grid, phi_grid, &size_phigrid, deltax);//
		printf("grid parameter = %f\n", grid_parameter);
		fprintf(fout, "grid parameter = %f\n", grid_parameter);
		printf("\t(phi_DSE, phi_wall) = (%f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut);
		fprintf(fout, "\t(phi_DSE, phi_wall) = (%f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut);
		printf("evaluate ion density in MPS\n");
		//size_ngrid = 0;
		//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 7.0), 2.0) - grid_parameter ) / deltax );
		densfinorb(0, 0, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i, ZOOM_MP, STOP_MP); 
		//printf("densfinorb module ran in %f seconds\n", jobtime);
		printf("size of ion density grid = %d\n", size_ngrid);
		fprintf(fout, "size of ion density grid = %d\n", size_ngrid);
		if (v_cut*v_cut < - 2.0*phi_grid[0]) v_cut = sqrt(-2.0*phi_grid[0]) + TINY;
		v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
		printf("v_cutDS = %f\n", v_cutDS);
		if (gamma_ref >= 0.1) {
			vpar_cut_lookup[0] = 1e10;
			mue_cut_lookup[0] = 0.0;
			for (i=1; i< size_cut; i++) {
				vpar_cut_lookup[i] = vparcut(M_PI - i*M_PI/size_cut, v_cutDS);
				mue_cut_lookup[i] = mucut(M_PI - i*M_PI/size_cut, v_cutDS);
				
				//printf("vpar = %f at mu = %f\n", vpar_cut_lookup[i], mue_cut_lookup[i]);
			}
			vpar_cut_lookup[size_cut] = 0.0;
			mue_cut_lookup[size_cut] = 1e10;
			//printf("%f %f\n", vpar_cut_lookup[size_cut], mue_cut_lookup[size_cut]);
		}
		else if (gamma_ref < 0.1) { // first solve MPS w/ simplified e- reflection
			for (i=0; i<= size_cut; i++) {
				vpar_cut_lookup[i] = v_cutDS;// + i;
				mue_cut_lookup[i] = i*(mu_e[size_mu_e-1]+TINY)/size_cut;;
				//printf("%f %f\n", vpar_cut_lookup[i], mue_cut_lookup[i]);
			}
		}
		denszeroorb(-1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, mue_cut_lookup, vpar_cut_lookup, size_cut);
		newguess(&convergence, Te, x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, &size_ngrid, 0.0, 0.0, 1.5, weight);//, p, m);
		current = flux_i - flux_e*sqrt(mioverme/2.0);
		if (fix_current == 1) {
			newvcut(&v_cut, v_cutDS, mioverme, flux_i, flux_e*sqrt(mioverme/2.0), target_current, error_current, &convergence_j, weight_j);
		}
		if (v_cut*v_cut < - 2.0*phi_grid[0]) {
			printf("WARNING: The total sheath and presheath potential drop is smaller than the presheath potential drop\n");
			printf("To avoid a non-monotonic potential, I will set the total potential drop to be just above the presheath potential drop\n");
			//v_cut = sqrt(-2.0*phi_grid[0] + TINY);
			//v_cutDS = sqrt(TINY);
			rescale_array(phi_grid, size_phigrid, -0.5*v_cut*v_cut + TINY );
			v_cutDS = sqrt(2.0*TINY);
		}
		printf("\tcurrent = %f (target = %f)\n", current, target_current);
		fprintf(fout, "\tcurrent = %f (target = %f)\n", current, target_current);
		printf("convergence (MP, j) = (%d, %d)\n", convergence, convergence_j);
		if ( (convergence > 1) && ( convergence_j > 1 || fix_current == 0) ) {	
			printf("ENTER ION DENSITY EVALUATION IN MPS\n");
			//size_ngrid = 0;
			//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 7.0), 2.0) - grid_parameter ) / deltax );
			densfinorb(1, sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i, ZOOM_MP, 0.999); 
			//The argument of densfinorb set to one makes the module compute the ion distribution function at x=0
			clock_t end_it = clock(); // finds end time of last iteration
			tot_time = (double) (end_it - begin_it) / CLOCKS_PER_SEC;
			printf("At %dth iteration MP iteration with simplified DS model converged successfully in %f seconds\n", N, tot_time);
			fprintf(fout, "At %dth iteration MP iteration with simplified DS model converged successfully in %f seconds\n", N, tot_time);
		} 
		N++;
	}
//}
if ( (gamma_ref <= 10.0) && (gamma_ref >= 0.1) ) {
// FULL DEBYE SHEATH SOLUTION CALCULATED WITH FINITE (DISTORTED) ELECTRON GYROORBITS
	F_i_DS[0] = malloc(sizevxopen*sizeof(double));
	if (gamma_ref < 1.0) 
		make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, 0.0, deltaxDS, 0, -phi_grid[0] - 0.5*v_cut*v_cut, 0.5/gamma_ref);
	else 
		make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, 0.0, deltaxDS, 0, -phi_grid[0] - 0.5*v_cut*v_cut, 0.5);
	N_DS = 0;
	while ( ( ( convergence <= 1) || (convergence_DS <= 1) || (convergence_j <= 1) ) && (N_DS < MAX_IT) ) {
		make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, 0.0, 1.0);
		//remake_MPgrid(x_grid, phi_grid, &size_phigrid, deltax);
		//for (i=0; i < size_cut+1; i++) 
		//	vpar_cut_lookup[i] = vpar_cut_lookup_old[i];	

		printf("in while loop for combined DS+MP iteration\n");
		v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
		printf("v_cutDS = %f\n", v_cutDS);
		//if (gamma_ref < 1.0) zoomfactor_DS = ZOOM_DS * (int) (1.0/gamma_ref);
		//else zoomfactor_DS = ZOOM_DS;
		//if (v_cutDS*v_cutDS < 1.0) zoomfactor_DS = (int) (2.0*ZOOM_DS);///(v_cutDS*v_cutDS));


		printf("evaluate electron density in MPS\n");
		denszeroorb(-1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, mue_cut_lookup, vpar_cut_lookup, size_cut);

		printf("evaluate ion density in MPS\n");
		densfinorb(1, sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i, ZOOM_MP, STOP_MP); 
		//clock_t gsl_t = clock(); 
		//m = gsl_matrix_alloc (size_neDSgrid-2, size_neDSgrid-2);
		//p = gsl_permutation_alloc (size_neDSgrid-2);
		//for (i = 0; i < size_neDSgrid-2; i++) {
		//	for (j = 0; j < size_neDSgrid-2; j++) {
		//		if (i == j) 
		//			gsl_matrix_set (m, i, j, -2.0);
		//      //gsl_matrix_set (m, i, j, -2.0 - deltaxsq*exp(phi_grid[i])/gammasq);
		//		else if ( (i==j+1) || (i==j-1) )
		//		      gsl_matrix_set (m, i, j, 1.0);
		//		else 
		//		      gsl_matrix_set (m, i, j, 0.0);
		//	}
		//}
		//gsl_linalg_LU_decomp (m, p, &s);
		//clock_t gsl_t_end = clock(); 
		//printf("time for LU decomp is %f\n",  (double) (gsl_t_end - gsl_t) / CLOCKS_PER_SEC);
		printf("flux_eDS = (%f, %f)\tflux_i = %f\n", flux_eDS*ne_grid[0]*sqrt(2.0*mioverme), flux_e*sqrt(mioverme), flux_i);
		if (factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut) < grid_parameter) { // && (phi_grid[0] + 0.5*v_cut*v_cut > 1.0e-13 ) )   {
			grid_parameter = factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut);
			weight = CAUTIOUSWEIGHT;
			weight_j = CAUTIOUSWEIGHT;
		}
		//power = log((phi_grid[2] - phi_grid[0])/(phi_grid[1] - phi_grid[0]))/log(x_grid[2]/x_grid[1]);
		//if ( (phi_grid[2] - phi_grid[1])/(phi_grid[1] - phi_grid[0]) > 0.95 )
		//	grid_parameter *= 0.9;
		//else if ( (phi_grid[2] - phi_grid[1])/(phi_grid[1] - phi_grid[0]) < 0.9 )
		//	grid_parameter *= 1.1;

		///printf("grid parameter = %f\n", grid_parameter);
		///fprintf(fout, "grid parameter = %f\n", grid_parameter);

		zero = 0.0;
		printf("evaluate ion density in DS\n");
		for (ncols=0; ncols<sizevxopen; ncols+=1) {
			F_i_DS[0][ncols] = fx_i_DS[ncols]/sqrt(2.0);
			vx_i_DS[ncols] *= sqrt(2.0);
		}
		denszeroorb(ioncharge, phi_DSgrid, ni_DSgrid, size_phiDSgrid, &flux_i_DS, F_i_DS, vx_i_DS, NULL, sizevxopen, 1, NULL, &zero, 0);



		printf("evaluate electron density in DS\n");
		for (ncols=0; ncols<size_mu_e; ncols+=1) {
			mu_e_DS[ncols] = mu_e[ncols];
			for (ind=0; ind<size_vpar_e; ind+=1) {
				U_e_DS[ind] = ind*dvpar*ind*dvpar;
				Uminmu_MPE = sqrt(2.0*U_e_DS[ind] - 2.0*phi_grid[0]);
				dist_e_GK[ncols][ind] = 0.5*bilin_interp(mu_e[ncols], Uminmu_MPE, dist_e_DK, mu_e, vpar_e, size_mu_e, size_vpar_e, -1, -1)/ne_grid[0];
				//printf("%f ", dist_e_GK[ncols][ind]);
			}
			//printf("\n");
		}
		flux_eDS = -100000.0;
		printf("size_phiDSgrid = %d\n", size_phiDSgrid);
		densfinorb(2, size_cut+1, 1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e_DS, U_e_DS, size_mu_e, size_vpar_e, 0.0, mue_cut_lookup, vpar_cut_lookup, &flux_eDS, zoomfactor_DS, STOP_DS); 

		printf("new MP potential guess\n");
		newguess(&convergence, Te, x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, &size_ngrid, 0.0, deltaphi, 1.5, weight);// p, m);

		printf("new WALL potential guess\n");
		flux_e = flux_eDS*ne_grid[0]*sqrt(2.0);
		current = flux_i - flux_e*sqrt(mioverme/2.0);
		if (fix_current == 1)  // && (convergence_DS*convergence >= 1) )
			newvcut(&v_cut, v_cutDS, mioverme, flux_i, flux_e*sqrt(mioverme/2.0), target_current, error_current, &convergence_j, weight_j/1.0);

		printf("new DS potential guess\n");
		if (v_cut*v_cut < - 2.0*phi_grid[0]) {
			printf("WARNING: The total sheath and presheath potential drop is smaller than the presheath potential drop\n");
			printf("To avoid a non-monotonic potential, I will set the total potential drop to be just above the presheath potential drop\n");
			v_cut = sqrt(-2.0*phi_grid[0] + 0.01);
			v_cutDS = sqrt(0.01);
		}
		newguess(&convergence_DS, Te, x_DSgrid, ne_DSgrid, ni_DSgrid, phi_DSgrid, size_phiDSgrid, &size_neDSgrid, 1/(gamma_ref*gamma_ref), sqrt(2.0*phi_grid[0] + v_cut*v_cut), 2.0, 1.0);//, p, m);
		//printf("before rescaling\n");
		//for (i=0; i < size_phiDSgrid; i++) {
		//	printf("x = %f\t phi = %f\n", x_DSgrid[i], phi_DSgrid[i]);
		//}
		//printf("-phi_grid[0] - 0.5*v_cut^2 = %f\n", - phi_grid[0] - 0.5*v_cut*v_cut);
		//rescale_array(phi_DSgrid, size_phiDSgrid, - phi_grid[0] - 0.5*v_cut*v_cut);
		//printf("after rescaling\n");
		//for (i=0; i < size_phiDSgrid; i++) {
		//	printf("x = %f\t phi = %f\n", x_DSgrid[i], phi_DSgrid[i]);
		//}
		printf("\tcurrent = %f (target = %f)\n", current, target_current);
		fprintf(fout, "\tcurrent = %f (target = %f)\n", current, target_current);
		printf("\t(phi_DSE, phi_wall) = (%f, %f)\n\tcurrent = %f\n", phi_grid[0], -0.5*v_cut*v_cut, current);
		fprintf(fout, "\t(phi_DSE, phi_wall) = (%f, %f)\n\t current = %f\n", phi_grid[0], -0.5*v_cut*v_cut, current);
		printf("ITERATION N = %d of which N_DS = %d\n\n", N, N_DS);
		fprintf(fout, "ITERATION N = %d of which N_DS = %d\n\n", N, N_DS);
		printf("convergence (MP,DS, j) = (%d, %d, %d)\n", convergence, convergence_DS, convergence_j);
		fprintf(fout, "convergence (MP, DS, j) = (%d, %d, %d)\n\n", convergence, convergence_DS, convergence_j);
		N_DS++; N++;
		//gsl_permutation_free (p);
		//gsl_matrix_free (m);
	}
	printf("FINISHED?\n");
	densfinorb(1, sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i, ZOOM_MP, 0.999); 
	densfinorb(2, size_cut+1, 1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e_DS, U_e_DS, size_mu_e, size_vpar_e, 0.0, mue_cut_lookup, vpar_cut_lookup, &flux_eDS, zoomfactor_DS, 0.999); 
}

if (N== MAX_IT) {
	printf("No convergence after %d iterations :(\n", MAX_IT);
	fprintf(fout, "No convergence after %d iterations :(\n", MAX_IT);
}
else {
	//densfinorb(sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i); 
	clock_t end_itDS = clock(); // finds end time of last iteration
	tot_time = (double) (end_itDS - begin_it) / CLOCKS_PER_SEC;
	printf("At %dth iteration MP+DS combined iteration converged successfully in %f seconds\n", N, tot_time);
	fprintf(fout, "At %dth iteration MP+DS combined iteration converged successfully in %f seconds\n", N, tot_time);
	printf("\t(phi_DSE, phi_wall) = (%f, %f)\n\tcurrent = %f\n", phi_grid[0], -0.5*v_cut*v_cut, current);
	fprintf(fout, "\t(phi_DSE, phi_wall) = (%f, %f)\n\t current = %f\n", phi_grid[0], -0.5*v_cut*v_cut, current);
}
fclose(fout);

// Print potential and density profiles to files
FILE *fp; 
char fpstr[150];
snprintf(fpstr, 150, "%s/phi_n_DS.txt", dirname);
//fp = fopen("OUTPUT/phi_n_DS.txt", "w");
fp = fopen(fpstr, "w");
if (fp == NULL)  
	printf("error when opening file\n");

for (i=0; i<size_phiDSgrid; i++) {
	fprintf(fp, "%f %f %f %f\n", x_DSgrid[i], phi_DSgrid[i], ni_DSgrid[i], ne_DSgrid[i]);
}
fclose(fp);
snprintf(fpstr, 150, "%s/phi_n_MP.txt", dirname);
fp = fopen(fpstr, "w");
if (fp == NULL) {
	printf("error when opening file\n");
}
for (i=0; i<size_phigrid; i++) {
	fprintf(fp, "%f %f %f %f\n", x_grid[i], phi_grid[i], ni_grid[i], ne_grid[i]);
}
fclose(fp);

free(mue_cut_lookup);
free(vpar_cut_lookup);
free(x_grid);
free(phi_grid);
free(ne_grid);
free(ni_grid);
free(x_DSgrid);
free(phi_DSgrid);
free(ne_DSgrid);
free(ni_DSgrid);
free(U_e_DS);
free(mu_e_DS);
free(fx_i_DS);
free(vx_i_DS);
free(F_i_DS[0]);
free(F_i_DS);
exit(0);

}
