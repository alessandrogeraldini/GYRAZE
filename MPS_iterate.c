#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include "mps.h"
#define INITIATE_GRID 1
#define INITIAL_GRID_PARAMETER 1.0
#define GRID_SPACING 0.2
#define TINY 1e-8
#define TEST_EL 0
#define MAX_IT 1000

void make_phigrid(double *x_grid, double *phi_grid, int size_phigrid, double grid_parameter, double deltax, int initial, double phi_jump, double len_scale) {
	int i;
	double g, xi, *ff, *new_phi, *new_x;
	double dtwophidsqrtxtwo, dphidfzero, dtwophidftwozero;
	FILE *fp;
	new_phi = malloc(size_phigrid*sizeof(double));
	new_x = malloc(size_phigrid*sizeof(double));
  	ff = malloc(size_phigrid*sizeof(double));
	fp = fopen("phidata.txt", "w");
	if (fp == NULL) {
		printf("Error opening new phidata.txt file\n");
		exit(EXIT_FAILURE);
	}

	dphidfzero = (phi_grid[1] - phi_grid[0])/deltax;
	dtwophidftwozero = (phi_grid[2] + phi_grid[0] - 2.0*phi_grid[1])/(deltax*deltax);
	dtwophidsqrtxtwo = pow(2.0*(sqrt(x_grid[1]) + sqrt(grid_parameter)), 2.0)*dtwophidftwozero + 2.0*dphidfzero; //*(phi_grid[2] + phi_grid[0] - 2.0*phi_grid[1])/(deltax*deltax);
	//printf("dphidfzero = %f, dtwophidftwozero=%f\n", dphidfzero*(-phi_grid[0]), dtwophidftwozero);
	//printf("ddtwophidsqrtxtwo = %f\n", dtwophidsqrtxtwo);
	printf("in make_phigrid: size_phigrid = %d\n", size_phigrid);
	if (initial != 0) {
			printf("grid_parameter = %f\n", grid_parameter);
		
		for (i=0; i<size_phigrid; i++) {
		    ff[i] = pow(sqrt(grid_parameter)+sqrt(x_grid[i]), 2.0) - grid_parameter;
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
			new_x[i] = pow( pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), 2.0);
			g = sqrt(new_x[i]);
			fprintf(fp, "%f %f\n", g, new_phi[i]);
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
			new_x[i] = pow( pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), 2.0);
			g = sqrt(new_x[i]);
			//if (new_x[i] < len_scale) {
			//	new_phi[i] = phi_jump*(1.0 - new_x[i]/len_scale); //pow(len_scale, 2.0)/pow(new_x[i] + len_scale, 2.0);
			//}	
			//else new_phi[i] = 0.0;
			new_phi[i] = phi_jump*pow(len_scale, 2.0)/pow(new_x[i] + len_scale, 2.0);
			fprintf(fp, "%f %f\n", g, new_phi[i]);
		}
		for (i=0; i < size_phigrid; i++) {
			phi_grid[i] = new_phi[i];
			x_grid[i] = new_x[i];
		}
	}
	free(new_phi);
	free(new_x);
	free(ff);
	fclose(fp);
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

// the set of functions below are necessary to model electron reflection from the infinitely thin Debye sheath 
// obtain parallel velocity cutoff as a function of magnetic moment with assumption rho_e >> lambda_D
double FF(double beta) {
	double integrand, integral=0.0, beta_s = 0.005, betaval;
	int ind, num_beta = (int) (beta/beta_s);
	for (ind = 0; ind < num_beta+1; ind++) {
		betaval = beta*ind/num_beta;
		integrand = 0.5*(1.0-cos(2.0*betaval)) / ( M_PI - betaval + 0.5*sin(2.0*betaval) );
		integral += integrand*beta/num_beta;
	}
        //#returnvalue = beta**3/(3*np.pi) 
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
	int zoomfactor = 2, zoomfactor_DS;
	int i; //, j, s;
	double error_current = 0.005;
	double ioncharge= 1.0, *TioverTe;
	double **dist_i_GK, *mu_i, *U_i, zero = 0.0;
	double v_max, **dist_e_DK, **dist_e_GK, *vpar_e, *mu_e, *U_e_DS;
	double *fx_i_DS, *vx_i_DS, **F_i_DS;
	double dvpar = 0.2/sqrt(2.0), Uminmu_MPE;
	int size_mu_i, size_U_i, size_mu_e, size_vpar_e, number_species;
	int nrows_distfile = 0, ncols=0, fix_current, ind;
	char line_million[1000000];
	//nrows_distfile counts the number of rows in the input file with the distribution function; ncols_distfile does the same for columns; ncols_Umufile  finds the number of columns for both rows or Umufile.txt. The first rows contains values of U, the second of mu. The columns are the possible values; sizeUU and sizemumu stores the final number of columns */
	clock_t begin_it = clock(); // Finds the start time of the computation
	int convergence = 0, problem = 0, problem_el = 0, convergence_DS = 0;
	/* convergence = 0 tells newguess.c to set smoothing for the next iteration, convergence = 1 tells it to not smooth, and convergence = 2 means that the convergence criterion has been satisfied and the iteration is finished */
	int N=0, N_DS=0;
	/* Number of maximum iterations is set to some large number e.g. 100 but iterations should converge in about 20. If they don't, then there is a problem. */
	double tot_time, v_cut , current, flux_e=0.0, flux_i=0.0, flux_i_DS=0.0, flux_eDS = 0.0, mioverme;
	double rhoeoverlambdaD_ref;
	int size_cut=200, sizevxopen;
	double *mue_cut_lookup=malloc((1+size_cut)*sizeof(double)), *vpar_cut_lookup=malloc((1+size_cut)*sizeof(double));

	double *ne_grid, *ni_grid, v_cutDS = 0.5, *ne_DSgrid;
	double target_current;

	//struct spline spl_phitone, spl_netophi, spl_new_phitone, spl_new_netophi;

	double *storevals;
	// Pointers used to read the files
	double *x_grid, *phi_grid, *phi_DSgrid, *x_DSgrid, *ni_DSgrid;
	double Te, alpha, alpha_deg, grid_parameter, deltax, deltaxDS, factor_small_grid_parameter=2.0, system_size;
	int size_phigrid, size_ngrid = 0, size_neDSgrid = 0, size_phiDSgrid;
	char line_hundred[100];
	//gsl_permutation *p;
	//gsl_matrix *m;

	// set the parameter g that controls the grid spacing in x 
	// transition from evenly spaced in sqrt(x) near x=0 to evenly spaced in x at x=infty
	grid_parameter = INITIAL_GRID_PARAMETER;
	// set the grid spacing
	deltax = GRID_SPACING;

	FILE *fout = fopen("output_MAGSHEATH.txt", "w");
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
	ncols = 0; // for counting numbers in each row (line) of file
	while (fgets(line_hundred, 100, input) != NULL) {	
		storevals = linetodata(line_hundred, strlen(line_hundred), &ncols);
		if (i==0)
			// input in degrees only for convenience
			// should be below 5 degrees (0.1 rad) for asymptotic theory to be valid
			alpha_deg = *storevals; 
		else if (i==1)
			rhoeoverlambdaD_ref = *storevals; 
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
				//v_cut = sqrt(log(mioverme/(2.0*M_PI))); //v_max*0.5; // some substantial fraction of v_max
				v_cut = 4.1; // set to a reasonable value
				// we will set v_cut to some value later
			}
			else {
				v_cut = *storevals; 
				target_current = current = 0.0; // gets calculated afterwords
			}
		}
		i += 1; // count the rows in the file
	}
	fclose(input);
	Te = 1/TioverTe[0]; // Te/Ti is used and labelled Te in code
	alpha = alpha_deg*M_PI/180; // alpha used in radians from now on
	// initial iteration assumes flat potential profile in magnetic presheath
	// therefore, the parameter v_cutDS is equal to v_cut
	v_cutDS = v_cut;
	//if (rhoeoverlambdaD_ref > 1.0) deltaxDS = GRID_SPACING/(rhoeoverlambdaD_ref);
	//else deltaxDS = GRID_SPACING;
	if (rhoeoverlambdaD_ref < 1.0) deltaxDS = GRID_SPACING/(rhoeoverlambdaD_ref);
	else deltaxDS = GRID_SPACING;
	printf("deltaxDS = %f\n\n\n", deltaxDS);
	if (Te < 1.0) system_size = 10.0;
	else system_size = 8.0*sqrt(Te);
	i=0;
	printf("INPUT PARAMETERS:\n\tmagnetic field angle = α = %f deg (%f rad)\n\telectron gyroradius over Debye length (reference value at MP entrance) = ρ_e/λ_D = %f\n\tnumber of species = %d\n", alpha_deg, alpha, rhoeoverlambdaD_ref, number_species);
	fprintf(fout, "INPUT PARAMETERS:\n\tmagnetic field angle = α = %f deg (%f rad)\n\telectron gyroradius over Debye length (reference value at MP entrance) = ρ_e/λ_D = %f\n\tnumber of species = %d\n", alpha_deg, alpha, rhoeoverlambdaD_ref, number_species);
	for (i=0; i<number_species; i++) {
	printf("\tion temperature = τ = ZT_i/T_e for species %d = %f\n\tmass ratio = m_i/m_e for species ? = %f\n", i+1, TioverTe[i], mioverme);
	fprintf(fout, "\tion temperature = τ = ZT_i/T_e for species %d = %f\n\tmass ratio = m_i/m_e for species ? = %f\n", i+1, TioverTe[i], mioverme);
	}
	printf("\tfix_current = %d\n", fix_current);
	fprintf(fout, "\tfix_current = %d\n", fix_current);
	if (fix_current != 0) {
		printf("\tcurrent = j/(e*n_MPE*v_t,i) = %f\n", target_current);
		fprintf(fout, "\tcurrent = j/(e*n_MPE*v_t,i) = %f\n", target_current);
		//printf("FIX CURRENT TO %f (in ion thermal velocity)\n", current);
		//fprintf(fout, "FIX CURRENT TO %f (in ion thermal velocity)\n", current);
	}
	else {
		printf("\twall potential = eφ_W/T_e = %f\n", -0.5*v_cut*v_cut);
		fprintf(fout, "\twall potential = eφ_W/T_e = %f\n", -0.5*v_cut*v_cut);
		//printf("FIX WALL POTENTIAL TO %f (in electron temperature units)\n", 0.5*v_cut*v_cut);
		//fprintf(fout, "FIX WALL POTENTIAL TO %f (in electron temperature units)\n", 0.5*v_cut*v_cut);
	}

	printf("magnetic presheath size in simulation (in ion gyroradii) = %f\n", system_size);
	fprintf(fout, "magnetic presheath size in simulation (in ion gyroradii) = %f\n", system_size);
	size_phigrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size), 2.0) - grid_parameter ) / deltax );
	size_phiDSgrid = (int) ( system_size*1.0*(1.0+1.0/rhoeoverlambdaD_ref)/deltaxDS );
	//size_phiDSgrid = (int) ( system_size/deltaxDS );
	printf("size phigrid = %d\n", size_phiDSgrid);
	//size_neDSgrid = (int) ( ( system_size*(1.0+1.0/rhoeoverlambdaD_ref) - 5.0)/deltaxDS );
	//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 5.0), 2.0) - grid_parameter ) / deltax );
	//size_ngrid = (int) ( ( system_size - 7.0)/deltaxDS );
	if (rhoeoverlambdaD_ref < 1.0) zoomfactor_DS = zoomfactor* (int) (1.0+1.0/rhoeoverlambdaD_ref);
	else zoomfactor_DS = zoomfactor;
	 
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
	FILE *distfile, *Umufile;
	if ((distfile = fopen("distfuncin.txt", "r")) == NULL)
	{	
		printf("cannot open file %s\n", "distfuncin.txt");
		fprintf(fout, "cannot open file %s\n", "distfuncin.txt");
		exit(-1); 
	}
	/* Count the number of rows in the distfuncin.txt file */
	while(fgets(line_million, 1000000, distfile) != NULL) {	
		nrows_distfile += 1; 
	}
	// Allocate the right amount of memory to the distribution function pointer
	dist_i_GK = (double**) calloc(nrows_distfile,sizeof(double*));
	// The number of rows is also the the size of the array in U_i 
	nrows_distfile = 0; // Set number of rows counter to zero again
	rewind(distfile); // rewind file
	// Read each line of file, extract the data (ie the numbers in the line) and count the columns
	// Once columns are counted, allocate memory to dist_i_GK and assign each number to a different column of dist_i_GK (for a fixed row)
	while (fgets(line_million, 1000000, distfile) != NULL) {
		storevals = linetodata(line_million, strlen(line_million), &ncols);
		dist_i_GK[nrows_distfile] = storevals;
		nrows_distfile +=1; 
	}
	fclose(distfile);
	//printf("~~~~~The second element of dist_i_GK is %f~~~~~\n",dist_i_GK[0][1]);
	// Extract two 1D arrays representing the grid points on which dist_i_GK is defined
	// first open file and check for error
	if ((Umufile = fopen("Umufile.txt", "r")) == NULL) {	
		printf("cannot open file %s\n", "distfuncin.txt");
		fprintf(fout, "cannot open file %s\n", "Umufile.txt");
		exit(-1); 
	}
	// Read each line of file, extract the data and assign the values of the first line to mu_i, second to U_i
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

	/*
	EXTRACT ELECTRON DISTRIBUTION FUNCTION
	Import the distribution function into a 2 dimensional array, F(mu,vpar). 
	*/
	distfile = fopen("dist_file.txt", "r");
	if (distfile == NULL)
	{	
		printf("Cannot open file %s\n", "dist_file.txt");
		fprintf(fout, "Cannot open file %s\n", "dist_file.txt");
		exit(-1); 
	}
	/* Count the number of rows in the distfuncin.txt file */
	while(fgets(line_million, 1000000, distfile) != NULL) {	
		nrows_distfile += 1; 
	}
	/* Allocate the right amount of memory to the distribution function pointer */
	dist_e_DK = (double**) calloc(nrows_distfile,sizeof(double*));
	dist_e_GK = (double**) calloc(nrows_distfile,sizeof(double*));
	nrows_distfile = 0; // Set number of rows counter to zero again
	rewind(distfile); // rewind to first line of file
	while (fgets(line_million, 1000000, distfile) != NULL) {
		// assign each number to a different column of dist_e_DK (for a fixed row)
		dist_e_DK[nrows_distfile] = linetodata(line_million, strlen(line_million), &ncols);
		// allocate memory to dist_e_GK which might be used to solve Debye sheath
		// dist_e_GK depends on magnetic presheath and Debye sheath solutions
		// will assign different values to dist_e_GK at each iteration
		dist_e_GK[nrows_distfile] = (double*) calloc(ncols,sizeof(double*));
		nrows_distfile +=1; 
		// jump to a new line and repeat the process
	}
	fclose(distfile);
	//printf("~~~~~The second element of dist_e_DK is %f~~~~~\n",dist_e_DK[0][1]);
	// Extract the two 1D arrays representing the grid points on which dist_i_GK is defined
	// first open file and check for error
	if ((Umufile = fopen("dist_file_arguments.txt", "r")) == NULL) {	
		printf("cannot open file %s\n", "dist_file_arguments.txt");
		fprintf(fout, "Cannot open file %s\n", "dist_file_arguments.txt");
		exit(-1); 
	}
	// Read each line of the file, extract the data (ie numbers in a line) 
	// Assign values of first line to mu_e, second to vpar_e
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
	F_i_DS = calloc(1,sizeof(double*));
	sizevxopen = (int) 100*sqrt((1.0+Te)*(1.0+1.0/Te));
	vx_i_DS = malloc(sizevxopen*sizeof(double));
	fx_i_DS = malloc(sizevxopen*sizeof(double));
	U_e_DS  = malloc(size_vpar_e*sizeof(double));

	/*
	CARRY OUT A SINGLE ION DENSITY CALCULATION
	*/
	if (MAX_IT == 0)
	{	
		make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, 0, 0.0, 1.0);
		//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 5.0), 2.0) - grid_parameter ) / deltax );
		problem = densfinorb(sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ne_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i, zoomfactor);
	}

	if (TEST_EL==1) {
		printf("size_phiDSgrid = %d\n", size_phiDSgrid);
		F_i_DS[0] = malloc(sizevxopen*sizeof(double));
		make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, zero, deltaxDS, 0, -3.02, 0.001);
		make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, 0, 0.0, 1.0);
		//phi_DSgrid[0] = -5.0;
		printf("in while loop for combined DS+MP iteration\n");
		printf("v_cutDS = %f\n", v_cutDS);
		printf("now evaluate electron density in MPS\n");
		v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
		for (i=0; i< size_cut+1; i++) {
			vpar_cut_lookup[i] = v_cutDS;
			mue_cut_lookup[i] = i*(mu_e[size_mu_e-1]+TINY)/size_cut;
			//printf("%f %f\n", vpar_cut_lookup[i], mue_cut_lookup[i]);
		}
		printf("v_cutDS = %f\n", v_cutDS);
		makelookup(-1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, mue_cut_lookup, vpar_cut_lookup, size_cut);
		//v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);

		printf("now evaluate ion density in MPS\n");
		//size_ngrid = 0;
		//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 5.0), 2.0) - grid_parameter ) / deltax );
		problem = densfinorb(sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i, zoomfactor); 
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
		//size_neDSgrid = (int) ( ( system_size*(1.0+1.0/rhoeoverlambdaD_ref) -5.0)/deltaxDS );
		printf("size_neDSgrid = %d\n", size_neDSgrid);
		problem_el = densfinorb(-size_cut-1, 1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e, U_e_DS, size_mu_e, size_vpar_e, 0.0, mue_cut_lookup, vpar_cut_lookup, &flux_eDS, 1); 

		exit(0);
	}


	/*
	ITERATE MAGNETIC PRESHEATH AND DEBYE SHEATH POTENTIAL TO FIND SOLUTION
	*/
	if ( (rhoeoverlambdaD_ref < 0.1) || (rhoeoverlambdaD_ref > 2.0) ) {
		/*
		A SIMPLIFIED ELECTRON MODEL IS USED TO CALCULATE THE DEBYE SHEATH POTENTIAL DROP
		*/
		N=0;
		while ( (convergence != 2) && (N<MAX_IT) ) {
			current = target_current;
			printf("ITERATION # = %d\n", N);
			fprintf(fout, "ITERATION # = %d\n", N);
			if ( (factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut) < grid_parameter) && (phi_grid[0] + 0.5*v_cut*v_cut >1.0e-13 ) )   {
				grid_parameter = factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut);
			}
			// MAKE ELECTROSTATIC POTENTIAL GRID
			if (problem == 0) {
				make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, 0.0, 1.0);
			}
			else  {
				make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, 0, 0.0, 1.0);
			}
			printf("grid parameter = %f\n", grid_parameter);
			fprintf(fout, "grid parameter = %f\n", grid_parameter);
			printf("\t(phi_DSE, phi_wall) = (%f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut);
			fprintf(fout, "\t(phi_DSE, phi_wall) = (%f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut);
			printf("ENTER ION DENSITY EVALUATION IN MPS\n");
			//size_ngrid = 0;
			//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 7.0), 2.0) - grid_parameter ) / deltax );
			problem = densfinorb(0, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i, zoomfactor); 
			//printf("densfinorb module ran in %f seconds\n", jobtime);
			printf("size of ion density grid = %d\n", size_ngrid);
			fprintf(fout, "size of ion density grid = %d\n", size_ngrid);
			if (v_cut*v_cut < - 2.0*phi_grid[0]) v_cut = sqrt(-2.0*phi_grid[0]) + TINY;
			v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
			printf("v_cutDS = %f\n", v_cutDS);
			if (rhoeoverlambdaD_ref >= 5.0) {
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
			else if (rhoeoverlambdaD_ref < 0.2) { // first solve MPS w/ simplified e- reflection
				for (i=0; i<= size_cut; i++) {
					vpar_cut_lookup[i] = v_cutDS;// + i;
					mue_cut_lookup[i] = i*(mu_e[size_mu_e-1]+TINY)/size_cut;;
					//printf("%f %f\n", vpar_cut_lookup[i], mue_cut_lookup[i]);
				}
			}
			makelookup(-1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, mue_cut_lookup, vpar_cut_lookup, size_cut);
			newguess(&convergence, problem, Te, alpha, x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, &size_ngrid, 0.0, dist_i_GK, U_i, mu_i, size_U_i, size_mu_i, 0.0);//, p, m);
			current = flux_i - flux_e*sqrt(mioverme/2.0);
			if ( (fix_current == 1) && (fabs(current - target_current) > error_current) ) {
				newvcut(&v_cut, v_cutDS, mioverme, flux_i, flux_e*sqrt(mioverme/2.0), target_current, &convergence);
			}
			printf("\tcurrent = %f (target = %f)\n", current, target_current);
			fprintf(fout, "\tcurrent = %f (target = %f)\n", current, target_current);
			//size_ngrid
			if (convergence == 2)
			{
				printf("ENTER ION DENSITY EVALUATION IN MPS\n");
				//size_ngrid = 0;
				//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 7.0), 2.0) - grid_parameter ) / deltax );
				densfinorb(sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i, zoomfactor); 
				//The argument of densfinorb set to one makes the module compute the ion distribution function at x=0
				clock_t end_it = clock(); // finds end time of last iteration
				tot_time = (double) (end_it - begin_it) / CLOCKS_PER_SEC;
				printf("At %dth iteration MP iteration with simplified DS model converged successfully in %f seconds\n", N, tot_time);
				fprintf(fout, "At %dth iteration MP iteration with simplified DS model converged successfully in %f seconds\n", N, tot_time);
			} 
			N++;
		}
	}
	else { 

		/*
		THE FULL DEBYE SHEATH SOLUTION IS CALCULATED WITH FINITE (DISTORTED) ELECTRON GYROORBITS
		*/
		F_i_DS[0] = malloc(sizevxopen*sizeof(double));
		//make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, zero, deltaxDS, 0, -phi_grid[0] - 0.5*v_cut*v_cut, 0.2+0.2/(rhoeoverlambdaD_ref));
		make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, zero, deltaxDS, 0, -phi_grid[0] - 0.5*v_cut*v_cut, 0.2+0.2/(rhoeoverlambdaD_ref));
		while ( (convergence != 2) || (convergence_DS != 2) ) {
			if (problem == 0) {
				make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, 0.0, 1.0);
			}
			else  {
				make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, 0, 0.0, 1.0);
			}
			printf("in while loop for combined DS+MP iteration\n");
			printf("v_cutDS = %f\n", v_cutDS);
			printf("now evaluate electron density in MPS\n");
			v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
			if (N_DS == 0) {
				for (i=0; i< size_cut+1; i++) {
					vpar_cut_lookup[i] = v_cutDS;
					mue_cut_lookup[i] = i*(mu_e[size_mu_e-1]+TINY)/size_cut;
					//printf("%f %f\n", vpar_cut_lookup[i], mue_cut_lookup[i]);
				}
			}
			printf("v_cutDS = %f\n", v_cutDS);
			makelookup(-1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, mue_cut_lookup, vpar_cut_lookup, size_cut);
			//v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);

			printf("now evaluate ion density in MPS\n");
			//size_ngrid = 0;
			//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 5.0), 2.0) - grid_parameter ) / deltax );
			problem = densfinorb(sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i, zoomfactor); 
			for (ncols=0; ncols<size_mu_e; ncols+=1) {
				for (ind=0; ind<size_vpar_e; ind+=1) {
					U_e_DS[ind] = ind*dvpar*ind*dvpar;
					Uminmu_MPE = sqrt(2.0*U_e_DS[ind] - 2.0*phi_grid[0]);
					dist_e_GK[ncols][ind] = 0.5*bilin_interp(mu_e[ncols], Uminmu_MPE, dist_e_DK, mu_e, vpar_e, size_mu_e, size_vpar_e, -1, -1)/ne_grid[0];
					//printf("%f ", dist_e_GK[ncols][ind]);
				}
				//printf("\n");
			}
			flux_eDS = -100000.0;
			printf("size_phiDSgrid = %d ????\n", size_phiDSgrid);
			printf("1st electron density evaluation in sheath\n");
			//i=0;
			//while (phi_DSgrid[i]/phi_DSgrid[0] < 0.05) {
			//	i+=1;
			//}
			//size_neDSgrid = i;
			//size_neDSgrid = (int) ( ( system_size*(1.0+1.0/rhoeoverlambdaD_ref) -5.0)/deltaxDS );
			problem_el = densfinorb(-size_cut-1, 1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e, U_e_DS, size_mu_e, size_vpar_e, 0.0, mue_cut_lookup, vpar_cut_lookup, &flux_eDS, zoomfactor_DS); 
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
			printf("size_neDSgrid = %d\n", size_neDSgrid);
			printf("flux_eDS = (%f, %f)\tflux_i = %f\n", flux_eDS*ne_grid[0]*sqrt(2.0*mioverme), flux_e*sqrt(mioverme), flux_i);
			printf("phi_DSE (%f) - phi_wall (%f) = %f\n", phi_grid[0], 0.5*v_cut*v_cut, phi_grid[0] +0.5*v_cut*v_cut);
			if ( (factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut) < grid_parameter) && (phi_grid[0] + 0.5*v_cut*v_cut > 1.0e-13 ) )   {
				grid_parameter = factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut);
			}
			printf("grid parameter = %f\n", grid_parameter);
			fprintf(fout, "grid parameter = %f\n", grid_parameter);
			zero = 0.0;
			printf("Evaluate ion density in sheath\n");
			for (ncols=0; ncols<sizevxopen; ncols+=1) {
				F_i_DS[0][ncols] = fx_i_DS[ncols]/sqrt(2.0);
				vx_i_DS[ncols] *= sqrt(2.0);
			}
			makelookup(ioncharge, phi_DSgrid, ni_DSgrid, size_phiDSgrid, &flux_i_DS, F_i_DS, vx_i_DS, NULL, sizevxopen, 1, NULL, &zero, 0);

			printf("size of ion density grid = %d\n", size_ngrid);
			fprintf(fout, "size of ion density grid = %d\n", size_ngrid);
			printf("now evaluate the new potential guess in the MPS\n");
			newguess(&convergence, problem, Te, alpha, x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, &size_ngrid, 0.0, dist_i_GK, U_i, mu_i, size_U_i, size_mu_i, 0.0);// p, m);
			flux_e = flux_eDS*ne_grid[0]*sqrt(2.0);
			current = flux_i - flux_e*sqrt(mioverme/2.0);
			if ( (fix_current == 1) && (fabs(current - target_current) > 0.01) ) {
				newvcut(&v_cut, v_cutDS, mioverme, flux_i, flux_e*sqrt(mioverme/2.0), target_current, &convergence);
				printf("changed v_cut\n\n\n");
			}
			printf("\tcurrent = %f (target = %f)\n", current, target_current);
			fprintf(fout, "\tcurrent = %f (target = %f)\n", current, target_current);
			if (convergence == 2) {
				clock_t end_it = clock(); // finds end time of last iteration
				tot_time = (double) (end_it - begin_it) / CLOCKS_PER_SEC;
				printf("MP converged, DS not yet\n");
				fprintf(fout, "MP converged, DS not yet\n");
			}
			printf("New sheath potential guess\n");
			//printf("%f ????\n", sqrt(2.0*phi_grid[0] + v_cut*v_cut));
			if (v_cut*v_cut < - 2.0*phi_grid[0]) {
				printf("WARNING: The total sheath and presheath potential drop is smaller than the presheath potential drop\n");
				printf("To avoid a non-monotonic potential, I will set the total potential drop to be just above the presheath potential drop\n");
				v_cut = sqrt(-2.0*phi_grid[0] + 0.01);
				v_cutDS = sqrt(0.01);
			}
			newguess(&convergence_DS, problem, Te, alpha, x_DSgrid, ne_DSgrid, ni_DSgrid, phi_DSgrid, size_phiDSgrid, &size_neDSgrid, 1.0/(rhoeoverlambdaD_ref*rhoeoverlambdaD_ref), dist_i_GK, U_i, mu_i, size_U_i, size_mu_i, sqrt(2.0*phi_grid[0] + v_cut*v_cut));//, p, m);
			printf("N iterations including DS = %d\n", N_DS);
			rescale_array(phi_DSgrid, size_phiDSgrid, - phi_grid[0] - 0.5*v_cut*v_cut);
			printf("\t(phi_DSE, phi_wall) = (%f, %f)\n\tcurrent = %f\n", phi_grid[0], -0.5*v_cut*v_cut, current);
			fprintf(fout, "\t(phi_DSE, phi_wall) = (%f, %f)\n\t current = %f\n", phi_grid[0], -0.5*v_cut*v_cut, current);
			printf("END OF ITERATION (N, N_DS) = (%d, %d)\n\n\n\n", N, N_DS);
			fprintf(fout, "END OF ITERATION (N, N_DS) = (%d, %d)\n\n\n\n", N, N_DS);
			//gsl_permutation_free (p);
			//gsl_matrix_free (m);
			N++; N_DS++;
		}
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
	
	// Print potential solutions to files
	FILE *potMP, *potDS;
	potMP = fopen("PostProcessing/phi_MP.txt", "w");
	potDS = fopen("PostProcessing/phi_DS.txt", "w");
	if ( (potMP == NULL) || (potDS == NULL) ) {
		printf("error when opening file\nmake sure PostProcessing directory is present\n");
	}
	for (i=0; i<size_phigrid; i++) {
		fprintf(potMP, "%f %f\n", x_grid[i], phi_grid[i]);
	}
	for (i=0; i<size_phiDSgrid; i++) {
		fprintf(potDS, "%f %f\n", x_DSgrid[i], phi_DSgrid[i]);
	}
	fclose(potMP);
	fclose(potDS);

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
	exit(0);
}
