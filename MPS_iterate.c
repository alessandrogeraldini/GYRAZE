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
#define INITIAL_GRID_PARAMETER 0.5
#define GRID_SPACING 0.2
#define TINY 1e-8
#define test 1
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
			//fprintf(fp, "%f %f\n", xi, phi[i]);
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
			//new_phi[i] = phi_jump*(-1.0  + new_x[i]/pow(pow(grid_parameter+(size_phigrid-1)*deltax + 1.0e-13, 0.5) - sqrt(grid_parameter), 2.0));
			new_phi[i] = phi_jump*pow(len_scale, 2.0)/pow(new_x[i] + len_scale, 2.0);
			//fprintf(fp, "%f %f\n", xi, phi[i]);
			fprintf(fp, "%f %f\n", g, new_phi[i]);
		}
		for (i=0; i < size_phigrid; i++) {
			phi_grid[i] = new_phi[i];
			x_grid[i] = new_x[i];
		}
	}
	fclose(fp);
	return;
}

void rescale_array(double *array, int size_array, double jump) {
	int i;
	double old_jump = array[0];
	for (i=0; i < size_array; i++) {
		array[i] *= (jump/old_jump);
	}
	return;
}

// the set of functions below are necessary to model electron reflection from the infinitely thin Debye sheath as a function of energy and magnetic moment under the assumption rho_e >> lambda_D
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
// end of model for rho_e >> lambda_D

int main() {
	double error_current = 0.01;
	double ioncharge= 1.0, *TioverTe;
	double **dist_i_GK, *mu_i, *U_i, zero = 0.0;
	double v_max, **dist_e_DK, **dist_e_GK, *vpar_e, *mu_e, *U_e_DS;
	double *fx_i_DS, *vx_i_DS, **F_i_DS;
	double dvpar = 0.2, Uminmu_MPE;
	int size_mu_i, size_U_i, size_mu_e, size_vpar_e, number_species;
	int nrows_distfile = 0, ncols=0, fix_current, ind;
	char line_million[1000000];
	//nrows_distfile counts the number of rows in the input file with the distribution function; ncols_distfile does the same for columns; ncols_Umufile  finds the number of columns for both rows or Umufile.txt. The first rows contains values of U, the second of mu. The columns are the possible values; sizeUU and sizemumu stores the final number of columns */
	clock_t begin_it = clock(); // Finds the start time of the computation
	int convergence = 0, problem = 0, problem_el = 0, i, convergence_DS = 0;
	/* convergence = 0 tells newguess.c to set smoothing for the next iteration, convergence = 1 tells it to not smooth, and convergence = 2 means that the convergence criterion has been satisfied and the iteration is finished */
	int N, N_DS;
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
	double Telarge, Tesmall;
	int size_phigrid, size_ngrid = 0, size_neDSgrid = 0, size_phiDSgrid;

	char line_hundred[100];
	FILE *fout = fopen("output_MAGSHEATH.txt", "w");
	if (fout == NULL) {
		printf("Problem opening output file for writing\n");
		exit(1);
	}

	grid_parameter = INITIAL_GRID_PARAMETER;
	deltax = GRID_SPACING;
	deltaxDS = GRID_SPACING;

	/* READ INPUT FILE 
	This contains the values of: 
	  alpha
	  number of ion species
	  Ti/Te for all species (array)
	  mi/me for all ion species (array)
	  fix_current = 1 or 0 (fixes potential),
	  fixed current/potential */
	FILE *input;
	if ((input = fopen("inputfile.txt", "r")) == NULL) { 
		printf("Cannot open %s\n", "inputfile.txt");
		fprintf(fout, "Cannot open %s\n", "inputfile.txt");
		exit(EXIT_FAILURE); 
	}
	i=0;
	ncols = 0;
	while (fgets(line_hundred, 100, input) != NULL) {	
		storevals = linetodata(line_hundred, strlen(line_hundred), &ncols);
		if (i==0)
			alpha_deg = *storevals; 
		// should be a number below 5 for asymptotic theory to be valid
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
				// we will set v_cut = v_max later
				// when v_max is assigned
			}
			else {
				v_cut = *storevals; 
				target_current = current = 0.0; // gets calculated afterwords
			}
		}
		i += 1;	
	}
	Te = 1/TioverTe[0];
	if (Te>1.0) {	
		Telarge = Te;
		Tesmall = 1.0;
	}
	else {
		Telarge = 1.0;
		Tesmall = Te;
	}
	alpha = alpha_deg*M_PI/180;
	v_cutDS = v_cut;
	fclose(input);
	system_size = 8.0*sqrt(0.5*(Te+1.0));
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
	size_phiDSgrid = (int) 1.0*system_size*(1.0+1.0/rhoeoverlambdaD_ref)/deltaxDS;
	//size_phiDSgrid = (int) ( 1.5*size_phigrid*(1.0+1.0/rhoeoverlambdaD_ref) ) ;
	printf("size of coarse potential grid in Debye sheath = %d\n", size_phiDSgrid);
	fprintf(fout, "size of coarse potential grid in Debye sheath = %d\n", size_phiDSgrid);
	printf("size of coarse potential grid in magnetic presheath = %d\n", size_phigrid);
	fprintf(fout, "size of coarse potential grid in magnetic presheath = %d\n", size_phigrid);

	x_grid = malloc(size_phigrid*sizeof(double));
	phi_grid = malloc(size_phigrid*sizeof(double));
	ne_grid = malloc(size_phigrid*sizeof(double));
	ni_grid = malloc(size_phigrid*sizeof(double));

	x_DSgrid = malloc(size_phiDSgrid*sizeof(double));
	phi_DSgrid = malloc(size_phiDSgrid*sizeof(double));
	ne_DSgrid = malloc(size_phiDSgrid*sizeof(double));
	ni_DSgrid = malloc(size_phiDSgrid*sizeof(double));

	// EXTRACT ION DISTRIBUTION FUNCTION
	/* This part of the code extracts the ion distribution function from a file. We import the distribution function into a 2 dimensional array, F(mu,U). */
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
	/* Allocate the right amount of memory to the distribution function pointer */
	dist_i_GK = (double**) calloc(nrows_distfile,sizeof(double*));
	/* The number of rows is also the the size of the array in U_i */
	nrows_distfile = 0; // Set number of rows counter to zero again
	fclose(distfile); // Close file
	// Check file exists and can be opened
	if ((distfile = fopen("distfuncin.txt", "r")) == NULL) {
		printf("cannot open file %s\n", "distfuncin.txt");
		fprintf(fout, "cannot open file %s\n", "distfuncin.txt");
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
	//printf("~~~~~The second element of dist_i_GK is %f~~~~~\n",dist_i_GK[0][1]);
	/* Now we extract the two 1D arrays representing the grid points on which dist_i_GK (the distribution function array) is defined. */
	// first open file and check for error
	if ((Umufile = fopen("Umufile.txt", "r")) == NULL) {	
		printf("cannot open file %s\n", "distfuncin.txt");
		fprintf(fout, "cannot open file %s\n", "Umufile.txt");
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


	// EXTRACT ELECTRON DISTRIBUTION FUNCTION
	/* This part of the code extracts the ion distribution function from a file. We import the distribution function into a 2 dimensional array, F(mu,U). */
	/* This part of the code extracts the electron distribution function from a file. We import the distribution function into a 2 dimensional array, F(mu,U). */
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
	/* The number of rows is also the the size of the array in mu_i */
	nrows_distfile = 0; // Set number of rows counter to zero again
	rewind(distfile); // rewind to first line of file
	/* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and counting the columns. Once the columns are counted we allocate memory to dist_i_GK and assign each number to a different column of dist_i_GK (for a fixed row). The while loop then jumps to a new line and repeats the process */
	while (fgets(line_million, 1000000, distfile) != NULL) {
		storevals = linetodata(line_million, strlen(line_million), &ncols);
		//dist_e_DK[nrows_distfile] = (double*)calloc(ncols,sizeof(double));
		dist_e_DK[nrows_distfile] = linetodata(line_million, strlen(line_million), &ncols);
		dist_e_GK[nrows_distfile] = (double*) calloc(ncols,sizeof(double*));
		//dist_e_DK[nrows_distfile] = storevals;
		nrows_distfile +=1; 
	}
	fclose(distfile);
	//printf("~~~~~The second element of dist_e_DK is %f~~~~~\n",dist_e_DK[0][1]);
	/* Now we extract the two 1D arrays representing the grid points on which dist_i_GK (the distribution function array) is defined. */
	// first open file and check for error
	if ((Umufile = fopen("dist_file_arguments.txt", "r")) == NULL) {	
		printf("cannot open file %s\n", "dist_file_arguments.txt");
		fprintf(fout, "Cannot open file %s\n", "dist_file_arguments.txt");
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
	if (fix_current != 0) v_cut = v_max*0.5; // some substantial fraction of v_max
	F_i_DS = calloc(1,sizeof(double*));
	//sizevxopen = (int) 100*sqrt(Telarge)/sqrt(Tesmall);
	sizevxopen = (int) 100*sqrt(Telarge)/sqrt(Tesmall);
	vx_i_DS = malloc(sizevxopen*sizeof(double));
	fx_i_DS = malloc(sizevxopen*sizeof(double));
	U_e_DS  = malloc(size_vpar_e*sizeof(double));

	srand ( time ( NULL));

	if (MAX_IT == 0)
	{	
		make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, 0, 0.0, 1.0);
		problem = densfinorb(sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ne_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i);
	}

	if ( (rhoeoverlambdaD_ref < 0.2) || (rhoeoverlambdaD_ref > 5.0) ) {
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
			problem = densfinorb(0, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i); 
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
			newguess(&convergence, problem, Te, alpha, x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, size_ngrid, 0.0, dist_i_GK, U_i, mu_i, size_U_i, size_mu_i, 0.0);
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
				densfinorb(sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i); 
				//The argument of densfinorb set to one makes the module compute the ion distribution function at x=0
				clock_t end_it = clock(); // finds end time of last iteration
				tot_time = (double) (end_it - begin_it) / CLOCKS_PER_SEC;
				printf("At %dth iteration MP iteration with simplified DS model converged successfully in %f seconds\n", N, tot_time);
				fprintf(fout, "At %dth iteration MP iteration with simplified DS model converged successfully in %f seconds\n", N, tot_time);
			} 
			N++;
		}
	}
	


	F_i_DS[0] = malloc(sizevxopen*sizeof(double));
	if ( (rhoeoverlambdaD_ref > 0.2) && (rhoeoverlambdaD_ref < 5.0) ) {

		current = target_current;
		N=0;
		printf("ITERATION # = %d\n", N);
		fprintf(fout, "ITERATION # = %d\n", N);
		if ( (factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut) < grid_parameter) && (phi_grid[0] + 0.5*v_cut*v_cut >1.0e-13 ) )   {
			grid_parameter = factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut);
		}
		printf("grid parameter = %f\n", grid_parameter);
		make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, 0.0, 1.0);
		make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, zero, deltaxDS, N, -phi_grid[0] - 0.5*v_cut*v_cut, 1.0+1.0/(rhoeoverlambdaD_ref));
		convergence = 0;
		N=1;
		N_DS = 0;
		if (test==0) {
		while ( (convergence != 2) || (convergence_DS != 2) ) {
			printf("in while loop for combined DS+MP iteration\n");
			v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
			printf("v_cutDS = %f\n", v_cutDS);
			printf("now evaluate electron density in MPS\n");
			if (N_DS == 0) {
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
			makelookup(-1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, mue_cut_lookup, vpar_cut_lookup, size_cut);

			printf("now evaluate ion density in MPS\n");
			problem = densfinorb(sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i); 
			//printf("densfinorb module ran in %f seconds\n", jobtime);

			// this is the part where DS scale is solved
			// fill in array corresponding to U_e_DS
			for (ncols=0; ncols<size_mu_e; ncols+=1) {
				for (ind=0; ind<size_vpar_e; ind+=1) {
					U_e_DS[ind] = ind*dvpar*ind*dvpar;
					Uminmu_MPE = sqrt(2.0*U_e_DS[ind] - 2.0*phi_grid[0]);
					dist_e_GK[ncols][ind] = 0.5*bilin_interp(mu_e[ncols], Uminmu_MPE, dist_e_DK, mu_e, vpar_e, size_mu_e, size_vpar_e, -1, -1)/ne_grid[0];
				}
			}
			flux_eDS = -100000.0;
			printf("1st electron density evaluation in sheath\n");
			problem_el = densfinorb(-size_cut-1, 1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e, U_e_DS, size_mu_e, size_vpar_e, 0.0, mue_cut_lookup, vpar_cut_lookup, &flux_eDS); 
			printf("flux_eDS = (%f, %f)\n", flux_eDS*ne_grid[0], flux_e);
			printf("electron density at infinity in sheath is ne_DSgrid[%d] = %f\n", size_neDSgrid-1, ne_DSgrid[size_neDSgrid-1]);
			printf("phi_DSE (%f) - phi_wall (%f) = %f\n", phi_grid[0], 0.5*v_cut*v_cut, phi_grid[0] + 0.5*v_cut*v_cut);
			/* This part rescales the electron density in the Debye sheath such that it 
			   is unity at infinity
			*/
			//while ( (ne_DSgrid[size_neDSgrid-1] < 0.98) || (ne_DSgrid[size_neDSgrid-1] > 0.995) ) {
				//printf("Rescaling electron density to normalize correctly in sheath\n");
				//verify this
				//rescale_array(phi_DSgrid, size_phiDSgrid, phi_grid[0] + 0.5*v_cut*v_cut);
				//problem_el = densfinorb(-size_cut-1, 1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, 1.0, dist_e_GK, mu_e, U_e_DS, size_mu_e, size_vpar_e, 0.0, mue_cut_lookup, vpar_cut_lookup, &flux_eDS); 
				//printf("electron density at infinity in sheath is ne_DSgrid[%d] = %f\n", size_neDSgrid-1, ne_DSgrid[size_neDSgrid-1]);
				//printf("phi_DSE (%f) - phi_wall (%f) = %f\n", phi_grid[0], 0.5*v_cut*v_cut, phi_grid[0] +0.5*v_cut*v_cut);
				//if (ne_DSgrid[size_neDSgrid-1] > 1.0)
				//	v_cut -= 0.05;
				//else if (ne_DSgrid[size_neDSgrid-1] < 0.99)
				//	v_cut += 0.05;
				
				//printf("phi_cut changed to %f\n", 0.5*v_cut*v_cut);
				//v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
				//printf("v_cutDS = %f\n", v_cutDS);
			//} 	
			//printf("After electron density rescaling:\n\t(phi_DSE, phi_wall) = (%f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut);
			//fprintf(fout, "After electron density rescaling:\n\t(phi_DSE, phi_wall) = (%f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut);
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
			//problem = 0;
			//if (problem == 0) {
			//	make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, 0.0, 1.0);
			//}
			//else  {
			//	make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, 0, 0.0, 1.0);
			//}

//ion density was here

			printf("size of ion density grid = %d\n", size_ngrid);
			fprintf(fout, "size of ion density grid = %d\n", size_ngrid);
			if (v_cut*v_cut < - 2.0*phi_grid[0]) v_cut = sqrt(-2.0*phi_grid[0]) + TINY;
			printf("now evaluate the new potential guess in the MPS\n");
			newguess(&convergence, problem, Te, alpha, x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, size_ngrid, 0.0, dist_i_GK, U_i, mu_i, size_U_i, size_mu_i, 0.0);
			current = flux_i - flux_e*sqrt(mioverme/2.0);
			if ( (fix_current == 1) && (fabs(current - target_current) > error_current) ) {
				newvcut(&v_cut, v_cutDS, mioverme, flux_i, flux_e*sqrt(mioverme/2.0), target_current, &convergence);
				printf("changed v_cut\n\n\n");
			}
			printf("\tcurrent = %f (target = %f)\n", current, target_current);
			fprintf(fout, "\tcurrent = %f (target = %f)\n", current, target_current);
			if (convergence == 2)
			{
				//densfinorb(sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i); 
				//The argument of densfinorb set to one makes the module compute the ion distribution function at x=0
				clock_t end_it = clock(); // finds end time of last iteration
				tot_time = (double) (end_it - begin_it) / CLOCKS_PER_SEC;
				printf("MP converged, DS not yet\n");
				fprintf(fout, "MP converged, DS not yet\n");
			} 
			printf("New sheath potential guess\n");
			newguess(&convergence_DS, problem, Te, alpha, x_DSgrid, ne_DSgrid, ni_DSgrid, phi_DSgrid, size_phiDSgrid, size_neDSgrid, 1.0/(rhoeoverlambdaD_ref*rhoeoverlambdaD_ref), dist_i_GK, U_i, mu_i, size_U_i, size_mu_i, 0.0);
			printf("N iterations including DS = %d\n", N_DS);
			rescale_array(phi_DSgrid, size_phiDSgrid, - phi_grid[0] - 0.5*v_cut*v_cut);
			printf("\t(phi_DSE, phi_wall) = (%f, %f)\n\tcurrent = %f\n", phi_grid[0], -0.5*v_cut*v_cut, current);
			fprintf(fout, "\t(phi_DSE, phi_wall) = (%f, %f)\n\t current = %f\n", phi_grid[0], -0.5*v_cut*v_cut, current);
			printf("\n\n\n\n\n\n\n");
			N++;
			N_DS++;
		}
	}
	}




///////////////////////////////////////////////////////////
// test alternative method here
///////////////////////////////////////////////////////////

	if (test == 1) {
		if ( (rhoeoverlambdaD_ref > 0.2) && (rhoeoverlambdaD_ref < 5.0) ) {
		while ( (convergence != 2) || (convergence_DS != 2) ) {
			printf("in while loop for combined DS+MP iteration\n");
			printf("v_cutDS = %f\n", v_cutDS);
			printf("now evaluate electron density in MPS\n");
			v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
			if (N_DS == 0) {
			//else {
			//	//if (0.05*( -2.0*log(2.0*flux_eDS*ne_grid[0]*sqrt(M_PI)) ) + 0.95*v_cut*v_cut + 2.0*phi_grid[0] < 0.0)  {
			//	//	printf("WARNING: avoid blow up\n");
			//	//	v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
			//	//}
			//	//else {
			//	v_cutDS = sqrt( 0.1*( -2.0*log(2.0*flux_eDS*ne_grid[0]*sqrt(M_PI)) ) + 0.9*v_cut*v_cut + 2.0*phi_grid[0] );
			//	//}
			//}
			//printf("v_cutDS fictitious= %f\n", v_cutDS);
			for (i=0; i< size_cut+1; i++) {
				vpar_cut_lookup[i] = v_cutDS;//
				//vpar_cut_lookup[i] = 0.1;//;
				mue_cut_lookup[i] = i*(mu_e[size_mu_e-1]+TINY)/size_cut;
				//printf("%f %f\n", vpar_cut_lookup[i], mue_cut_lookup[i]);
			}
			}
			printf("v_cutDS = %f\n", v_cutDS);
			makelookup(-1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, mue_cut_lookup, vpar_cut_lookup, size_cut);
			//v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);

			printf("now evaluate ion density in MPS\n");
			problem = densfinorb(sizevxopen, Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, vx_i_DS, fx_i_DS, &flux_i); 
			for (ncols=0; ncols<size_mu_e; ncols+=1) {
				for (ind=0; ind<size_vpar_e; ind+=1) {
					U_e_DS[ind] = ind*dvpar*ind*dvpar;
					Uminmu_MPE = sqrt(2.0*U_e_DS[ind] - 2.0*phi_grid[0]);
					dist_e_GK[ncols][ind] = 0.5*bilin_interp(mu_e[ncols], Uminmu_MPE, dist_e_DK, mu_e, vpar_e, size_mu_e, size_vpar_e, -1, -1)/ne_grid[0];
				}
			}
			flux_eDS = -100000.0;
			printf("1st electron density evaluation in sheath\n");
			problem_el = densfinorb(-size_cut-1, 1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e, U_e_DS, size_mu_e, size_vpar_e, 0.0, mue_cut_lookup, vpar_cut_lookup, &flux_eDS); 
			printf("flux_eDS = (%f, %f)\tflux_i = %f\n", flux_eDS*ne_grid[0]*sqrt(2.0*mioverme), flux_e*sqrt(mioverme), flux_i);
			printf("electron density at infinity in sheath is ne_DSgrid[%d] = %f\n", size_neDSgrid-1, ne_DSgrid[size_neDSgrid-1]);
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
			newguess(&convergence, problem, Te, alpha, x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, size_ngrid, 0.0, dist_i_GK, U_i, mu_i, size_U_i, size_mu_i, 0.0);
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
			//(1.0/sqrt(ne_grid[0]))
			newguess(&convergence_DS, problem, Te, alpha, x_DSgrid, ne_DSgrid, ni_DSgrid, phi_DSgrid, size_phiDSgrid, size_neDSgrid, 1.0/(rhoeoverlambdaD_ref*rhoeoverlambdaD_ref), dist_i_GK, U_i, mu_i, size_U_i, size_mu_i, sqrt(2.0*phi_grid[0] + v_cut*v_cut));
			printf("N iterations including DS = %d\n", N_DS);
			rescale_array(phi_DSgrid, size_phiDSgrid, - phi_grid[0] - 0.5*v_cut*v_cut);
			printf("\t(phi_DSE, phi_wall) = (%f, %f)\n\tcurrent = %f\n", phi_grid[0], -0.5*v_cut*v_cut, current);
			fprintf(fout, "\t(phi_DSE, phi_wall) = (%f, %f)\n\t current = %f\n", phi_grid[0], -0.5*v_cut*v_cut, current);
			printf("\n\n\n\n\n\n\n");
			N++;
			N_DS++;
		}
		}
	}

//end of test










	if (N== MAX_IT) {
		printf("No convergence after %d iterations\n", MAX_IT);
		fprintf(fout, "No convergence after %d iterations\n", MAX_IT);
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
	FILE *potMP, *potDS;
	potMP = fopen("phi_MP.txt", "w");
	potDS = fopen("phi_DS.txt", "w");
	if ( (potMP == NULL) || (potDS == NULL) ) {
		printf("error when opening file\n");
	}
	for (i=0; i<size_phigrid; i++) {
		fprintf(potMP, "%f %f\n", x_grid[i], phi_grid[i]);
	}
	for (i=0; i<size_phiDSgrid; i++) {
		fprintf(potDS, "%f %f\n", x_DSgrid[i], phi_DSgrid[i]);
	}
	exit(0);
}
