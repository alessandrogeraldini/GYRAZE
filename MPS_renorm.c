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
#include <sys/types.h>
#include "mps_renorm.h"
#define MAX_IT 500
// number of maximum iterations set to some large number but iterations should converge in about 20-100. If they don't, then there is a problem
#define INITIAL_GRID_PARAMETER 1.4
#define TEST_EL 0
#define SYS_SIZ 20.0
#define ADHOCFS 1
#define MAXV 6.0
#define DV 0.14

const double tol_MP[2] = {0.0005, 0.004}, tol_DS[2] = {0.001, 0.004}, tol_current = 0.01;
const double WEIGHT_j = 0.1, WEIGHT_MP = 0.3, CAUTIOUSWEIGHT = 0.15;
const double STOP_MP = 0.98, STOP_DS = 0.98;
const double GRIDSIZE_MP = 0.4, GRIDSIZE_DS = 0.2;
const int ZOOM_DS = 1, ZOOM_MP = 2;

const char *strqty[7] = {"alpha=","gamma=","nspec=","Ti:Te=","mi:me=","jwall=","pwall="};
const int lenstrqty = 7;

//void remake_MPgrid(double *x_grid, double *phi_grid, int *psize_phigrid, double delta) {
//	int i, phithenx=1, size_phigrid = *psize_phigrid;
//	double xi, *new_phi, *new_x;
//	//FILE *fp;
//	new_phi = malloc(size_phigrid*sizeof(double));
//	new_x = malloc(size_phigrid*sizeof(double));
//	//fp = fopen("OUTPUT/phidata.txt", "w");
//	//if (fp == NULL) {
//	//	printf("Error opening new phidata.txt file\n");
//	//	exit(EXIT_FAILURE);
//	//}
//	printf("in remake_MPgrid: size_phigrid = %d\n", size_phigrid);
//	gsl_interp_accel *acc
//	= gsl_interp_accel_alloc ();
//	gsl_spline *spline
//	= gsl_spline_alloc (gsl_interp_cspline, size_phigrid);
//
//	gsl_spline_init (spline, phi_grid, x_grid, size_phigrid);
//
//	new_phi[0] = phi_grid[0];
//	new_x[0] = x_grid[0];
//	//for (i=1; i < size_phigrid; i++) {
//	i=1;
//	while (new_x[i-1] < x_grid[size_phigrid-1] - delta) {
//		if (phithenx == 1) {
//			xi = i*delta;
//			new_phi[i] = new_phi[i-1] + delta;
//			new_x[i] = gsl_spline_eval (spline, new_phi[i], acc);
//		}
//		if ( (phithenx == 0) || (new_x[i] - new_x[i-1] > delta) ) {
//			new_x[i] = new_x[i-1] + delta;
//			if (phithenx == 1)  {
//				gsl_spline_free (spline);
//				gsl_interp_accel_free (acc);
//
//				gsl_interp_accel *acc
//				= gsl_interp_accel_alloc ();
//				gsl_spline *spline
//				= gsl_spline_alloc (gsl_interp_cspline, size_phigrid);
//				gsl_spline_init (spline, x_grid, phi_grid, size_phigrid);
//				phithenx = 0;
//			}
//			printf(" i = %d\n", i);
//				new_phi[i] = gsl_spline_eval (spline, new_x[i], acc);
//			printf(" i = %d\n", i);
//				printf("x phi = (%f %f)\n", new_x[i], new_phi[i]);
//			printf(" i = %d\n", i);
//		}
//		//g = sqrt(new_x[i]);
//		printf("x phi = (%f %f)\n", new_x[i], new_phi[i]);
//		//fprintf(fp, "%f %f\n", g, new_phi[i]);
//		i++;
//	}
//	*psize_phigrid = i;
//	for (i=0; i < size_phigrid; i++) {
//		phi_grid[i] = new_phi[i];
//		x_grid[i] = new_x[i];
//	}
//	free(new_phi);
//	free(new_x);
//	return;
//}

void densionDS(double alpha, double *ni_DS, double *phi_DS, double phi0, double **FF, double *mu, double *Uminmu, double *vy, double *chiM, double *twopidmudvy, int size_phi, int size_mu, int size_U) {
	int i, j, k, count;
	double deltaUperp, halfVx0sq, intgrd, intgrdold, vzk, vzkm, n_inf;
	n_inf = 0.0;
	intgrd = 0.0;
	for (j=0; j< size_mu; j++) {
		intgrdold = intgrd;
		intgrd=0.0;
		for (k=1; k < size_U; k++) {
			halfVx0sq = chiM[j] - 0.5*vy[j]*vy[j] - phi0;
			deltaUperp = mu[j] - chiM[j] ;
			//if ( (deltaUperp < 0.0) && (deltaUperp > -0.01) ) deltaUperp = 0.0; // printf("YOLO\n");
			if (phi0 == 0.0) {
				halfVx0sq = 0.0; 
				deltaUperp = 0.0;
			}
			vzk = sqrt(2.0*(deltaUperp + Uminmu[k]));
			vzkm = sqrt(2.0*(deltaUperp + Uminmu[k-1]));
			intgrd += ( (sqrt(2.0*(halfVx0sq + alpha*vzk*twopidmudvy[j])) - sqrt(2.0*halfVx0sq)) * FF[j][k] + (sqrt(2.0*(halfVx0sq + alpha*vzkm*twopidmudvy[j])) - sqrt(2.0*halfVx0sq)) * FF[j][k-1] ) * 0.5 * ( vzk - vzkm );
		}
		if (j != 0) 
			n_inf += (intgrd + intgrdold)*0.5*(mu[j] - mu[j-1]);
		//printf("n_inf = %f\n", n_inf);
	}
	for (i=0; i < size_phi; i++) { 
		ni_DS[i] = 0.0;
		for (j=0; j< size_mu; j++) {
			intgrdold = intgrd;
			intgrd = 0.0;
			count = 0;
			for (k=1; k < size_U; k++) {
				halfVx0sq = chiM[j] - 0.5*vy[j]*vy[j] - phi_DS[i] - phi0;
				deltaUperp = mu[j] - chiM[j];
				if (phi0 == 0.0) {
					halfVx0sq = -phi_DS[i]; 
					deltaUperp = 0.0;
				}
				vzk = sqrt(2.0*(deltaUperp + Uminmu[k]));
				vzkm = sqrt(2.0*(deltaUperp + Uminmu[k-1]));
				intgrd += ( (sqrt(2.0*(halfVx0sq + alpha*vzk*twopidmudvy[j])) - sqrt(2.0*halfVx0sq)) * FF[j][k] + (sqrt(2.0*(halfVx0sq + alpha*vzkm*twopidmudvy[j])) - sqrt(2.0*halfVx0sq)) * FF[j][k-1] ) * 0.5 * ( vzk - vzkm );
				if ((count == 0) && (intgrd != intgrd) ) {
					count = 1;
					//printf("???%f %f %f %f\n", halfVx0sq, deltaUperp, vzk, vzkm);
				}
			}
			if (j != 0) 
				ni_DS[i] += (intgrd + intgrdold)*0.5*(mu[j] - mu[j-1]);
		}
		ni_DS[i] /= n_inf;
		//printf("ni_DS = %f\tphi_DS = %f\n", ni_DS[i], phi_DS[i]);
	}
	return;
}

void Figen(double **ffarr, double *Uminmuarr, double *muarr, int number_species, double mioverme, double *TioverTe, int sizevpar, int sizevperp, double dvpar, double dvperp) {
	int n, i, j, coldelectrons;
	double mu, Uminmu, ff, Te;;
	double u, condition, chodura, normalization, flow=0.0;

	for (n=0; n<number_species; n++) {
		Te = 1.0/TioverTe[n];
		if (Te < 1e-9) coldelectrons = 1;
		else coldelectrons = 0;
		printf("are electrons cold (1=yes, 0=no): %d\n", coldelectrons);
		u = 0.0;
		condition = 1.0/Te; 
		if (coldelectrons == 0) {
			if (Te>=1.0) {
				chodura = 9999999.0;
				u = 0.0;
				while ( (chodura > condition) || (chodura < condition - 0.05) ) {
					if (chodura > condition)
						u += 0.0001;
					else 
						u -= 0.0001;
					normalization =  (1 + erf(u/sqrt(2.0)))*(1.0+u*u) + sqrt(2.0/M_PI)*u*exp(-0.5*u*u);
					printf("normalization = %f\n", normalization);
					chodura = (1+erf(u/sqrt(2.0)))/(normalization);
					flow = ( 0.5*(1+0.5*u*u)*exp(-0.5*u*u) + (1.5 + 0.5*u*u)*(sqrt(M_PI/2.0)*u/2.0)*(1 + erf(u/sqrt(2.0))))*4.0/sqrt(2.0*M_PI)/normalization;
					printf("%f, %f, %f, %f %f\n", normalization, u, flow, chodura, condition);
				}
			}
			else {
				chodura = 0.0;
				while  ( (chodura > condition) || (chodura < condition - 0.05) ) {
					if (chodura > condition)
						u -= 0.01;
					else
						u += 0.01;
					normalization =  2.0*(sqrt(M_PI*u) - M_PI*exp(1/u)*(1-erf(1.0/sqrt(u))))/(2.0*u*sqrt(u));
					printf("normalization = %f\n", normalization);
					chodura = 2.0*M_PI*exp(1.0/u)*(1.0-erf(1.0/sqrt(u)))/(2.0*sqrt(u)*normalization);
					flow = -9999.9;
					printf("%f, %f, %f, %f %f\n", normalization, u, flow, chodura, condition);
				}
			}
		}

		if (coldelectrons == 0) {
			if ( Te >= 1.0 ) 
				normalization =  (1 + erf(u/sqrt(2.0)))*(1.0+u*u) + sqrt(2.0/M_PI)*u*exp(-0.5*u*u);
				//normalization =  (1 + erf(u))*(1.0+2.0*u*u) + (2.0/sqrt(M_PI))*u*exp(-u*u);
			else
				normalization =  2.0*(sqrt(M_PI*u) - M_PI*exp(1/u)*(1-erf(1/sqrt(u))))/(2.0*u*sqrt(u));
		}

		if (coldelectrons == 0) {
			if (Te>=1.0) {
				printf("Te = %f >= 1\n", Te);
				for (i=0; i<sizevperp; i++) {
					mu = 0.5*i*dvperp*i*dvperp;
					muarr[i] = mu;
					for (j=0; j<sizevpar; j++) { 
						Uminmu = 0.5*j*dvpar*j*dvpar;
						if (i==0)
							Uminmuarr[j] = Uminmu;
						ff = (1.0/(M_PI*sqrt(M_PI)))*(1.0/normalization)*(Uminmu)*exp(- Uminmu - mu + u*sqrt(2.0*Uminmu) - 0.5*u*u);
						ffarr[i][j] = ff;
					}
				}
			}
			else {
				printf("Te = %f < 1\n", Te);
				for (i=0; i<sizevperp; i++) {
					mu = 0.5*i*dvperp*i*dvperp;
					muarr[i] = mu;
					for (j=0; j<sizevpar; j++) { 
						Uminmu = 0.5*j*dvpar*j*dvpar;
						if (i==0)
							Uminmuarr[j] = Uminmu;

						ff = (1.0/(4.0*M_PI))*(1.0/normalization)*((Uminmu)/(1+u*(Uminmu)))*exp(- Uminmu - mu );
						ffarr[i][j] = ff;
					}
				}
			}
		}
		else { // coldelectrons == 1:
			printf("Te = infinity\n");
			normalization = 1.0;
			for (i=0; i<sizevperp; i++) {
				mu = 0.5*i*dvperp*i*dvperp;
				muarr[i] = mu;
				for (j=0; j<sizevpar; j++) { 
					Uminmu = 0.5*j*dvpar*j*dvpar;
					if (i==0)
						Uminmuarr[j] = Uminmu;
					ff = (1.0/(2.0*(M_PI)*sqrt(M_PI)))*exp(- Uminmu - mu );
					ffarr[i][j] = ff;
				}
			}
		}
	}
	return;
}

void Fegen(double **ffarr, double *mue, double *vpar_e, int sizemu_e, int sizevpar_e, double dvperp, double dvpar) {
	double mu, vpar, ff, Uminmu;
	int i,j;
	for (i=0; i<sizemu_e; i++) { 
		mu = 0.5*i*dvperp*i*dvperp;
		mue[i] = mu;
		for (j=0; j<sizevpar_e; j++) { 
			vpar = j*dvpar;
			Uminmu = 0.5*vpar*vpar;
			if (i==0) vpar_e[j] = vpar;
			ff = (1.0/(2.0*(M_PI*sqrt(M_PI))))*exp(- Uminmu - mu );
			ffarr[i][j] = ff;
		}
	}
	return;
}

//double xifp(double phi) {
//	double y;
//	y = 1.0/sqrt( -3.0 -2.0*phi + 4.0*exp(phi) - exp(2.0*phi) );
//	return y;
//}

void phi_fluidMPS(double *phigrid, double *xgrid, int size_x, double alpha) {
	double *x_for_interp, *phi_for_interp, *psi_for_interp, phi0=log(alpha), dxdphi=0.0, dxdphiold;
	int i, size_for_interp = 500;
	double dphi = -phi0/size_for_interp;
	x_for_interp = malloc(size_for_interp*sizeof(double));
	phi_for_interp = malloc(size_for_interp*sizeof(double));
	psi_for_interp = malloc(size_for_interp*sizeof(double));
	for (i=0; i<size_for_interp; i++) {
		dxdphiold = dxdphi;
		phi_for_interp[i] = log(alpha) + i*dphi;
		psi_for_interp[i] = phi_for_interp[i] + 0.5*alpha*alpha*(exp(-2.0*phi_for_interp[i])-1.0);
		dxdphi = 1.0/sqrt( -3.0 -2.0*psi_for_interp[i] + 4.0*exp(psi_for_interp[i]) - exp(2.0*psi_for_interp[i]) ) ; //xifp(psi_for_interp[i]);
		if (i!=0) x_for_interp[i] = x_for_interp[i-1] + (dxdphi + dxdphiold)*0.5*dphi;
		else x_for_interp[0] = 0.0;
		printf(" x,phi = %f, %f\n", x_for_interp[i], phi_for_interp[i]);
	}
	for (i=0; i<size_x; i++) {
		phigrid[i] = lin_interp(x_for_interp, phi_for_interp, xgrid[i], size_for_interp, 1228);
	}
	return;
}


void make_phigrid(double *x_grid, double *phi_grid, int size_phigrid, double grid_parameter, double deltax, int initial, double phi_jump, double len_scale, double alpha) {
	// improved = 0 is only working case
	int i, improved = 0;//, size_phigrid_improved;
	double xi, *ff, *new_phi, *new_x, power = 0.5, deltaxmin = 0.2, deltaphi=0.1;
	new_phi = malloc(size_phigrid*sizeof(double));
	new_x = malloc(size_phigrid*sizeof(double));
  	ff = malloc(size_phigrid*sizeof(double));
	printf("in make_phigrid: size_phigrid = %d\n", size_phigrid);
	if (initial != 0) {
		if (improved == 0) {
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
				//g = sqrt(new_x[i]);
			}
			gsl_spline_free (spline);
			gsl_interp_accel_free (acc);
			for (i=0; i < size_phigrid; i++) {
				phi_grid[i] = new_phi[i];
				x_grid[i] = new_x[i];
			}
		}
		else { // improved =1 STILL NOT WORKING
			gsl_interp_accel *acc
			= gsl_interp_accel_alloc ();
			gsl_spline *spline
			= gsl_spline_alloc (gsl_interp_cspline, size_phigrid);

			gsl_spline_init (spline, phi_grid, x_grid, size_phigrid);

			gsl_interp_accel *acc2
			= gsl_interp_accel_alloc ();
			gsl_spline *spline2
			= gsl_spline_alloc (gsl_interp_cspline, size_phigrid);

			gsl_spline_init (spline2, x_grid, phi_grid, size_phigrid);

			new_x[0] = 0.0;
			new_phi[0] = phi_grid[0];
			i=0;
			while (i < size_phigrid) {
			//for (i=1; i < size_phigrid; i++) 
				if (phi_grid[0] + i*deltaphi < phi_grid[size_phigrid-1])
					xi = gsl_spline_eval (spline, phi_grid[0] + deltaphi*i, acc);
				else xi = x_grid[i-1] + deltaxmin + TINY;
				printf("xi = %f\n", xi);
				if (xi < x_grid[i-1] + deltaxmin) { 
					new_x[i] = xi;
					new_phi[i] = phi_grid[0] + deltaphi*i;
					printf("new_phi = %f\n", new_phi[i]);
				}
				else if (new_x[i-1] + deltaxmin > x_grid[size_phigrid-1]) {
					//i = size_phigrid-1;
					new_x[i] = new_x[i-1] + deltaxmin; 
					new_phi[i] = new_phi[i-1]  + TINY;
				}
				else {
					new_x[i] = new_x[i-1] + deltaxmin;
					new_phi[i] = gsl_spline_eval (spline2, new_x[i], acc);
				}
				printf("i=%d\tnew_x = %f\tnew_phi = %f\n", i, new_x[i], new_phi[i]);
				i++;
			}
			gsl_spline_free (spline);
			gsl_interp_accel_free (acc);
			gsl_spline_free (spline2);
			gsl_interp_accel_free (acc2);
			for (i=0; i < size_phigrid; i++) {
				phi_grid[i] = new_phi[i];
				x_grid[i] = new_x[i];
			}
		}
	}
	else {
		for (i=0; i < size_phigrid; i++) {
			xi = i*deltax;
			ff[i] = xi;
			//new_x[i] = pow( pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), 2.0);
			new_x[i] = pow( pow(grid_parameter+xi, power) - pow(grid_parameter, power), 2.0);
			//g = sqrt(new_x[i]);
			new_phi[i] = phi_jump*pow(len_scale, 2.0)/pow(new_x[i] + len_scale, 2.0);
			//if (new_x[i] < 10.0/sqrt(2)) 
			//	new_phi[i] = 3.0*pow((new_x[i]*sqrt(2.0)/10.0 - 1.0), 5.0);
			//else 
			//	new_phi[i] = 0.0;
		}
		//phi_fluidMPS(new_phi, new_x, size_phigrid, alpha);
		for (i=0; i < size_phigrid; i++) {
			phi_grid[i] = new_phi[i];
			x_grid[i] = new_x[i];
		}
	}
	free(new_phi);
	free(new_x);
	free(ff);
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

double vparcut_mu(double mu, double vcut) {
	double beta, musearch=-10000.0, vparcutn;
	double betalow = 0.0, betaup = M_PI;
	do {
		beta = (betaup + betalow)/2.0;
		musearch = 0.5*exp(-2.0*FF(beta))*pow(vcut/sin(beta), 2.0);
		if (musearch < mu) {
			betaup = beta;
		}
		else {
			betalow = beta;
		}
		printf("%f %f\n", musearch, mu);
	} while (fabs(musearch-mu) > TINY) ;
	vparcutn = sqrt((1.0-exp(-2.0*FF(beta))))*vcut/sin(beta);
	return vparcutn;
}

// The main function of MAGSHEATH
int main() {
	//double power = 1.0;
	clock_t begin_it = clock(); // Finds the start time of the computation
	double *vy_e_wall, *vy_i_wall, *chiM_e, *chiM_i, *twopidmudvy_e, *twopidmudvy_i;
	int zoomfactor_DS;
	int size_mu_i, size_U_i, size_mu_e, size_vpar_e, number_species;
	int i, j, nrows_distfile = 0, ncols=0, ndirname, fix_current=0, ind; //s
	double error_MP[2], error_DS[2], olderror_MP, olderror_DS, Bohmshouldbe;
	double deltaphi=0.0, weight_j=WEIGHT_j, weight_MP=WEIGHT_MP, phiDSbump=-999.0;
	double ioncharge= 1.0, *TioverTe, phi0_init_MP = -0.0;
	double **dist_i_GK, *mu_i, *U_i, zero = 0.0;
	double **dist_e_DK, **dist_e_GK, *vpar_e, *mu_e, *U_e_DS, *vpar_e_DS;
	double dvpar = DV, Uminmu_MPE;
	char line_million[1000000], dirname[200];
	//nrows_distfile counts the number of rows in the input file with the distribution function; ncols_distfile does the same for columns; sizeUU and sizemumu stores the final number of columns */
	int convergence_MP = 0, convergence_DS = 0, convergence_j = 0;
	int N=0, N_DS=0;
	double tot_time, v_cut , current, flux_e=0.0, flux_i=0.0, flux_i_DS=0.0, flux_eDS = 0.0, mioverme;
	double gamma_ref;
	int sizevxopen;
	double *vpar_e_cut;
	double *ne_grid, *ni_grid, v_cutDS = 0.5, *ne_DSgrid;
	double target_current;
	double *storevals;
	// Pointers used to read the files
	double *x_grid, *phi_grid, *phi_DSgrid, *x_DSgrid, *ni_DSgrid;
	double Te, alpha, alpha_deg, deltaxDS, factor_small_grid_parameter=1.0, system_size, DS_size;
	int size_phigrid, size_ngrid = 0, size_neDSgrid = 0, size_phiDSgrid;
	char line_hundred[100];
	double grid_parameter, deltax;
	int n_dims;

	//printf("%f %f \n", vparcut(0.01, 3.0), mucut(0.01, 3.0)); 
	//printf("%f %f \n", vparcut_mu(mucut(1.0, 3.0), 3.0), vparcut(1.0, 3.0));
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
			TioverTe = storevals; //linetodata(line_hundred, strlen(line_hundred), &ncols); 
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
		else if (i==7) 
			n_dims = (int) (*storevals); 
		if (i<5) {
			ndirname += strlen(line_hundred)-1;
			dirname[ndirname] = '_';
			ndirname += 1;
		}
		i += 1; // count the rows in the file
		//printf("dirname = %s\n", dirname);
	}
	dirname[ndirname+1] = '\0';
	ndirname+=1;
	fclose(input);
	if (alpha_deg < 2.0) weight_MP = WEIGHT_MP/2.0;
	if (alpha_deg < 1.0) weight_MP = WEIGHT_MP/3.0;
	if (TioverTe[0] < 0.6) weight_MP = WEIGHT_MP/3.0;
	//if (TioverTe[0] < 0.4) weight_MP = WEIGHT_MP/6.0;
	printf("directory where output will be stored is %s\n", dirname);
	mkdir(dirname, S_IRWXU);
	Te = 1/TioverTe[0]; // Te/Ti is used and labelled Te in code
	printf("Te = %f\n", Te);
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


	if (ADHOCFS == 1) {
		size_mu_i = (int) (MAXV/DV);
		size_U_i = (int)  (MAXV/DV);
		size_mu_e = (int) (MAXV/DV);
		size_vpar_e = (int) (MAXV/DV);
		mu_i = malloc(size_mu_i*sizeof(double));
		mu_e = malloc(size_mu_e*sizeof(double));
		U_i = malloc(size_U_i*sizeof(double));
		vpar_e = malloc(size_mu_e*sizeof(double));
		dist_i_GK = malloc(size_mu_i*sizeof(double));
		dist_e_DK = malloc(size_mu_e*sizeof(double));
		dist_e_GK = malloc(size_mu_e*sizeof(double));
		for (i=0; i < size_mu_e; i++) {
			dist_e_DK[i] = malloc(size_vpar_e*sizeof(double));
			dist_e_GK[i] = malloc(size_vpar_e*sizeof(double));
		}
		for (i=0; i < size_mu_i; i++) 
			dist_i_GK[i] = malloc(size_U_i*sizeof(double));
		Figen(dist_i_GK, mu_i, U_i, number_species, mioverme, TioverTe, size_U_i, size_mu_i, DV, DV);
		Fegen(dist_e_DK, mu_e, vpar_e, size_vpar_e, size_mu_e, DV, DV);
	}
	else { 
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
			else if (i==1) {
				size_U_i = ncols;
				//printf("size_U_i = %d\n", size_U_i);
				U_i = storevals;
				//for (i=0; i< size_U_i; i++) printf("U_i[%d] = %f\n", i, U_i[i]);
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
	}
	sizevxopen = (int) 50*sqrt((1.0+Te)*(1.0+1.0/Te));
	printf("sizevxopen = %d\n", sizevxopen);
	U_e_DS  = malloc(size_vpar_e*sizeof(double));
	vpar_e_DS  = malloc(size_mu_e*sizeof(double));
	vpar_e_cut  = malloc(size_mu_e*sizeof(double));
	vy_e_wall  = malloc(size_mu_e*sizeof(double));
	chiM_e  = malloc(size_mu_e*sizeof(double));
	twopidmudvy_e  = malloc(size_mu_e*sizeof(double));
	vy_i_wall  = malloc(size_mu_i*sizeof(double));
	chiM_i  = malloc(size_mu_i*sizeof(double));
	twopidmudvy_i  = malloc(size_mu_i*sizeof(double));


	/*
	CARRY OUT A SINGLE ION DENSITY CALCULATION
	*/
	if (MAX_IT == 0) {	
		make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, 0, phi0_init_MP, 1.0, alpha);
		//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 5.0), 2.0) - grid_parameter ) / deltax );
		densfinorb(Te, alpha, size_phigrid, &size_ngrid, ne_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, &flux_i, ZOOM_MP, STOP_MP, -999.9, vy_i_wall, chiM_i, twopidmudvy_i);
	}

	if (TEST_EL==1) {
		//size_phiDSgrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(2.0*system_size), 2.0) - 0.3 ) / deltax );
		//printf("size_phiDSgrid = %d\n", size_phiDSgrid);
		v_cut = 2.0;
		make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, 0, -0.0, 1.0, alpha);
		//phi_DSgrid[0] = -5.0;
		v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]) ;
		make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, 0.0, deltaxDS, 0, -0.5*v_cutDS*v_cutDS, 1.0/gamma_ref, alpha);
		for (i=0; i< size_mu_e; i++) 
			vpar_e_cut[i] = v_cutDS;
		printf("v_cutDS = %f\nphi_DSgrid[0] = %f\n", v_cutDS, phi_DSgrid[0]);
		printf("now evaluate electron density in MPS\n");
		denszeroorb(-1.0, 1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, vpar_e_cut, 0.0, x_grid);
		//v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);

		printf("now evaluate ion density in MPS\n");
		//size_ngrid = 0;
		//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 5.0), 2.0) - grid_parameter ) / deltax );
		densfinorb(Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, &flux_i, ZOOM_MP, STOP_MP, -999.9, vy_i_wall, chiM_i, twopidmudvy_i); 
		for (ncols=0; ncols<size_mu_e; ncols+=1) {
			for (ind=0; ind<size_vpar_e; ind+=1) {
				U_e_DS[ind] = 0.5*ind*dvpar*ind*dvpar;
				Uminmu_MPE = sqrt(2.0*U_e_DS[ind] - 2.0*phi_grid[0]);
				if (ncols == 0) vpar_e_DS[ind] = sqrt(2.0*U_e_DS[ind] - 2.0*phi_grid[0]);
				dist_e_GK[ncols][ind] = bilin_interp(mu_e[ncols], Uminmu_MPE, dist_e_DK, mu_e, vpar_e, size_mu_e, size_vpar_e, -1, -1)/ne_grid[0];
				//printf("%f ", dist_e_GK[ncols][ind]);
			}
			//printf("\n");
		}
		printf("size grid = %d\n", size_phiDSgrid);
		printf("size_neDSgrid = %d\n", size_neDSgrid);
		densfinorb(1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e, U_e_DS, size_mu_e, size_vpar_e, 0.0, &flux_eDS, ZOOM_DS, 1.9, -999.9, vy_e_wall, chiM_e, twopidmudvy_e); 
		for (i=0; i< size_mu_e; i++) 
			vpar_e_cut[i] = 0.0; //sqrt(-2.0*lin_interp(x_DSgrid, phi_DSgrid, sqrt(mue_cut_lookup[i]), size_phiDSgrid, 1234)); 
		for (i=0; i<80; i++) printf("x=%f\tanalytical density for negligible gyroradius = %f and phi = %f\n", x_DSgrid[i], (1+erf(sqrt(-2.0*(phi_DSgrid[0]-phi_DSgrid[i]))))*exp(phi_DSgrid[i])/((1+erf(sqrt(-2.0*(phi_DSgrid[0]))))), phi_DSgrid[i]);
		denszeroorb(-1.0, 1.0, phi_DSgrid, ne_DSgrid, size_phiDSgrid, &flux_eDS, dist_e_DK, vpar_e_DS, mu_e, size_vpar_e, size_mu_e, vpar_e_cut, gamma_ref, x_DSgrid);
		exit(0);
	}



	N=0;
	//ITERATE MAGNETIC PRESHEATH AND DEBYE SHEATH POTENTIAL TO FIND SOLUTION
	//if ( (gamma_ref > 5.0) || (gamma_ref < 0.05) ) { 
		weight_j=WEIGHT_j;
	// SIMPLIFIED ELECTRON MODEL USED TO CALCULATE DEBYE SHEATH POTENTIAL DROP
		error_MP[0] = error_DS[0] = 10000.0;
		while ( ( (convergence_MP <= 1) || ( convergence_j <= 1 && fix_current == 1 ) ) && (N<MAX_IT) ) {
			olderror_MP = error_MP[0];
			olderror_DS = error_DS[0];
			printf("weight_MP = %f\n", weight_MP);
			printf("ITERATION # = %d\n", N);
			fprintf(fout, "ITERATION # = %d\n", N);
			current = target_current;

			make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, phi0_init_MP, 1.0, alpha);
			if ( (factor_small_grid_parameter*sqrt(phi_grid[0] + 0.5*v_cut*v_cut) < INITIAL_GRID_PARAMETER) )  {// && (phi_grid[0] + 0.5*v_cut*v_cut >1.0e-13 ) )   
				//grid_parameter = 0.0001 + factor_small_grid_parameter*sqrt(phi_grid[0] + 0.5*v_cut*v_cut);
				//weight = CAUTIOUSWEIGHT;
				//weight_j = CAUTIOUSWEIGHT/2.0;
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
			//else remake_MPgrid(x_grid, phi_grid, &size_phigrid, deltax);//
			printf("grid parameter = %f\n", grid_parameter);
			fprintf(fout, "grid parameter = %f\n", grid_parameter);
			printf("\t(phi_mp0, phi_ds0, phi_wall) = (%f, %f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut - phi_grid[0], -0.5*v_cut*v_cut);
			fprintf(fout, "\t(phi_mp0, phi_ds0, phi_wall) = (%f, %f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut - phi_grid[0], -0.5*v_cut*v_cut);

			printf("evaluate ion density in MPS\n");
			//size_ngrid = 0;
			//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 7.0), 2.0) - grid_parameter ) / deltax );
			phiDSbump = - 0.5*v_cut*v_cut - phi_grid[0];  
			//phiDSbump = -999.9;
			densfinorb(Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, &flux_i, ZOOM_MP, STOP_MP, phiDSbump, vy_i_wall, chiM_i, twopidmudvy_i); 
			printf("densfinorb module ran\n");
			printf("size of ion density grid = %d\n", size_ngrid);
			fprintf(fout, "size of ion density grid = %d\n", size_ngrid);
			//if (v_cut*v_cut < - 2.0*phi_grid[0]) v_cut = sqrt(-2.0*phi_grid[0]) + TINY;
			v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
			if (v_cut*v_cut < - 2.0*phi_grid[0]) v_cutDS = 0.0;
			printf("v_cutDS = %f\n", v_cutDS);
			if (gamma_ref >= TINY) {
				//vpar_cut_lookup[0] = 1e10;
				//mue_cut_lookup[0] = 0.0;
				//for (i=1; i< size_cut; i++) {
				//	vpar_cut_lookup[i] = vparcut(M_PI - i*M_PI/size_cut, v_cutDS);
				//	mue_cut_lookup[i] = mucut(M_PI - i*M_PI/size_cut, v_cutDS); }
				//vpar_cut_lookup[size_cut] = 0.0;
				//mue_cut_lookup[size_cut] = 1e10;
				for (i=0; i< size_mu_e; i++) 
					vpar_e_cut[i] = vparcut_mu(mu_e[i], v_cutDS);
			}
			else if (gamma_ref < TINY) { // first solve MPS w/ simplified e- reflection
			      for (i=0; i< size_mu_e; i++) 
			      	vpar_e_cut[i] = v_cutDS;
			}
			denszeroorb(-1.0, 1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, vpar_e_cut, 0.0, x_grid);
			Bohmshouldbe = (1.0/Te)*(ne_grid[1] - ne_grid[0])/((phi_grid[1] - phi_grid[0])*ni_grid[0]);
			printf("Bohm integral should be %f\n", Bohmshouldbe);
			fprintf(fout, "Bohm integral of converged MP solution should be %f\n", Bohmshouldbe);
			//newguess(x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, &size_ngrid, 0.0, 0.0, 1.5, weight);//, p, m);
			error_Poisson(error_MP, x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, size_ngrid, 0.0);
			printf("error_av = %f\terror_max = %f\n", error_MP[0], error_MP[1]);
			if ( (error_MP[0] < tol_MP[0]) && (error_MP[1] < tol_MP[1]) ) convergence_MP += 1 ;
			else convergence_MP = 0;
			if (convergence_MP == 0) {
				printf("MP not converged --> calculate new MP potential guess\n");
				fprintf(fout, "MP not converged --> calculate new MP potential guess\n");
				newguess(x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, &size_ngrid, 0.0, deltaphi, 1.5, weight_MP);// p, m);
			}
			else {
				printf("MP converged --> no iteration needed\n");
				fprintf(fout, "MP converged --> no iteration needed\n");
			}
			current = flux_i - flux_e*sqrt(Te*mioverme);
			if (fix_current == 1) {
				if (fabs((target_current - current)/flux_i) > tol_current) convergence_j = 0;
				else convergence_j += 1;
				printf("target current = %f +/- %f\n", target_current, tol_current*flux_i);
				fprintf(fout, "target current = %f +/- %f\n", target_current, tol_current*flux_i);
				printf("current = %f = ion_current (%f) - electron_current (%f) = %f x ion_current\n", current, flux_i, flux_e*sqrt(Te*mioverme), current/flux_i);
				fprintf(fout, "current = %f = ion_current (%f) - electron_current (%f) = %f x ion_current\n", current, flux_i, flux_e*sqrt(Te*mioverme), current/flux_i);
				if (convergence_j == 0) {
					fprintf(fout, "new WALL potential guess\n");
					printf("new WALL potential guess\n");
					newvcut(&v_cut, v_cutDS, mioverme, TioverTe[0], flux_i, flux_e*sqrt(Te*mioverme), target_current, tol_current, weight_j);
				}
				else {
					printf("current converged, NO wall potential iteration\n");
					fprintf(fout, "current converged, NO wall potential iteration\n");
				}
			}
			else {
				printf("current = %f = ion_current (%f) - electron_current (%f) = %f x ion_current\n", current, flux_i, flux_e*sqrt(Te*mioverme), current/flux_i);
				fprintf(fout, "current = %f = ion_current (%f) - electron_current (%f) = %f x ion_current\n", current, flux_i, flux_e*sqrt(Te*mioverme), current/flux_i);
			}
			//if (v_cut*v_cut < - 2.0*phi_grid[0]) {
			//	printf("WARNING: The total sheath and presheath potential drop is smaller than the presheath potential drop\n");
			//	printf("To avoid a non-monotonic potential, I will set the total potential drop to be just above the presheath potential drop\n");
			//	//v_cut = sqrt(-2.0*phi_grid[0] + TINY);
			//	//v_cutDS = sqrt(TINY);
			//	rescale_array(phi_grid, size_phigrid, -0.5*v_cut*v_cut + TINY );
			//	v_cutDS = sqrt(2.0*TINY);
			//}
			printf("\tcurrent = %f (target = %f)\n", current, target_current);
			fprintf(fout, "\tcurrent = %f (target = %f)\n", current, target_current);
			printf("convergence (MP, j) = (%d, %d)\n", convergence_MP, convergence_j);
			if ( (convergence_MP > 1) && ( convergence_j > 1 || fix_current == 0) ) {	
				//printf("ENTER ION DENSITY EVALUATION IN MPS\n");
				//The argument of densfinorb set to one makes the module compute the ion distribution function at x=0
				densfinorb(Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, &flux_i, ZOOM_MP, 1.000, -999.9, vy_i_wall, chiM_i, twopidmudvy_i); 
				error_Poisson(error_MP, x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, size_ngrid, 0.0);
				printf("error_av = %f\terror_max = %f\n", error_MP[0], error_MP[1]);
				//if ( error_MP[1] > 2.0*tol_MP[1] ) {
				//	convergence_MP = 0;
				//	N++;
				//}
			} 
			else N++;
			//if (error_MP[0] > olderror_MP) weight_MP /= 1.05;
		}
		clock_t end_it = clock(); // finds end time of last iteration
		tot_time = (double) (end_it - begin_it) / CLOCKS_PER_SEC;
		printf("At %dth iteration MP iteration with simplified DS model converged successfully in %f seconds\n", N, tot_time);
		fprintf(fout, "At %dth iteration MP iteration with simplified DS model converged successfully in %f seconds\n", N, tot_time);
	//}
	if ( (gamma_ref <= 5.0) && (gamma_ref >= 0.05) ) {
		if (alpha_deg < 2.0) weight_j = WEIGHT_j/2.5;
		if (alpha_deg < 1.0) weight_j = WEIGHT_j/4.0;
	// FULL DEBYE SHEATH SOLUTION CALCULATED WITH FINITE (DISTORTED) ELECTRON GYROORBITS
		make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, phi0_init_MP, 1.0, alpha);
		if (gamma_ref < 1.0) 
			make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, 0.0, deltaxDS, 0, -phi_grid[0] - 0.5*v_cut*v_cut, 1.0/gamma_ref, alpha);
		else 
			make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, 0.0, deltaxDS, 0, -phi_grid[0] - 0.5*v_cut*v_cut, 1.0, alpha);
		//vpar_cut_lookup[0] = 1e10;
		//mue_cut_lookup[0] = 0.0;
		//for (i=1; i< size_cut; i++) {
		//	vpar_cut_lookup[i] = vparcut(M_PI - i*M_PI/size_cut, v_cutDS);
		//	mue_cut_lookup[i] = mucut(M_PI - i*M_PI/size_cut, v_cutDS);
		//	//printf("vpar = %f at mu = %f\n", vpar_cut_lookup[i], mue_cut_lookup[i]);
		//}
		//vpar_cut_lookup[size_cut] = 0.0;
		//mue_cut_lookup[size_cut] = 1e10;
		v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
		if (gamma_ref >= TINY) {
			//vpar_cut_lookup[0] = 1e10;
			//mue_cut_lookup[0] = 0.0;
			//for (i=1; i< size_cut; i++) {
			//	vpar_cut_lookup[i] = vparcut(M_PI - i*M_PI/size_cut, v_cutDS);
			//	mue_cut_lookup[i] = mucut(M_PI - i*M_PI/size_cut, v_cutDS); }
			//vpar_cut_lookup[size_cut] = 0.0;
			//mue_cut_lookup[size_cut] = 1e10;
			for (i=0; i< size_mu_e; i++) 
				vpar_e_cut[i] = vparcut_mu(mu_e[i], v_cutDS);
		}
		else if (gamma_ref < TINY) { // first solve MPS w/ simplified e- reflection
		      for (i=0; i< size_mu_e; i++) 
			vpar_e_cut[i] = v_cutDS;
		}
		N_DS = 0;
		convergence_MP = convergence_j = 0;
		error_MP[0] = error_DS[0] = 100000.0;
		if (fix_current == 0) convergence_j=2;
		while ( ( ( convergence_MP <= 1) || (convergence_DS <= 1) || (convergence_j <= 1) ) && (N_DS < MAX_IT) ) {
			printf("weight_MP = %f\n", weight_MP);
			olderror_MP = error_MP[0];	
			olderror_DS = error_DS[0];	
			make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, phi0_init_MP, 1.0, alpha);
			printf("in while loop for combined DS+MP iteration\n");
			printf("v_cutDS = %f\n", v_cutDS);
			//if (gamma_ref < 1.0) zoomfactor_DS = ZOOM_DS * (int) (1.0/gamma_ref);
			//else zoomfactor_DS = ZOOM_DS;
			//if (v_cutDS*v_cutDS < 1.0) zoomfactor_DS = (int) (2.0*ZOOM_DS);///(v_cutDS*v_cutDS));
			printf("evaluate electron density in MPS\n");
			denszeroorb(-1.0, 1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, vpar_e_cut, 0.0, x_grid);
			Bohmshouldbe = (1.0/Te)*(ne_grid[1] - ne_grid[0])/((phi_grid[1] - phi_grid[0])*ne_grid[0]);
			printf("Bohm integral should be %f\n", Bohmshouldbe);
			fprintf(fout, "Bohm integral of converged MP solution should be %f\n", Bohmshouldbe);
			//for (i=1; i< size_cut; i++) 
			//	vpar_cut_lookup[i] = sqrt(pow(vpar_cut_lookup[i], 2.0) - 2.0*phi_grid[i]);
			
			printf("evaluate ion density in MP\n");
			densfinorb(Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, &flux_i, ZOOM_MP, STOP_MP, -999.9, vy_i_wall, chiM_i, twopidmudvy_i); 
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
			printf("flux_eDS = (%f, %f)\tflux_i = %f\n", flux_eDS*ne_grid[0]*sqrt(Te*mioverme), flux_e*sqrt(Te*mioverme), flux_i);
			if (factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut) < INITIAL_GRID_PARAMETER) { // && (phi_grid[0] + 0.5*v_cut*v_cut > 1.0e-13 ) )   {
				//grid_parameter = factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut);
				//weight = CAUTIOUSWEIGHT;
				//weight_j = CAUTIOUSWEIGHT;
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
			//denszeroorb(ioncharge, Te, phi_DSgrid, ni_DSgrid, size_phiDSgrid, &flux_i_DS, F_i_DS, vx_i_DS, NULL, sizevxopen, 1, NULL, &zero, 0, 0.0, x_DSgrid);
			densionDS(alpha, ni_DSgrid, phi_DSgrid, phi_grid[0], dist_i_GK, mu_i, U_i, vy_i_wall, chiM_i, twopidmudvy_i, size_phiDSgrid, size_mu_i, size_U_i);

			printf("evaluate electron density in DS\n");
			for (ncols=0; ncols<size_mu_e; ncols+=1) {
				for (ind=0; ind<size_vpar_e; ind+=1) {
					U_e_DS[ind] = 0.5*ind*dvpar*ind*dvpar;
					Uminmu_MPE = sqrt(2.0*U_e_DS[ind] - 2.0*phi_grid[0]);
					if (ncols == 0) vpar_e_DS[ind] = Uminmu_MPE;
					dist_e_GK[ncols][ind] = bilin_interp(mu_e[ncols], Uminmu_MPE, dist_e_DK, mu_e, vpar_e, size_mu_e, size_vpar_e, -1, -1)/ne_grid[0];
					//printf("%f ", dist_e_GK[ncols][ind]);
				}
				//printf("\n");
			}
			flux_eDS = -100000.0;
			printf("size_phiDSgrid = %d\n", size_phiDSgrid);

			if (gamma_ref > SMALLGAMMA) {
				densfinorb(1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e, U_e_DS, size_mu_e, size_vpar_e, 0.0, &flux_eDS, zoomfactor_DS, STOP_DS, -999.9, vy_e_wall, chiM_e, twopidmudvy_e); 
				//printf("chiM mu\n");
				for (i=0; i<size_mu_e; i++) {
					vpar_e_cut[i] = sqrt(2.0*(chiM_e[i] - mu_e[i]));
					//printf("%f %f \n", chiM_e[i], mu_e[i]);
				} 
				printf("flux_eDS = %f\n", flux_eDS);
				//exit(1);
			}
			else {
				for (i=0; i< size_mu_e; i++) 
					vpar_e_cut[i] = 0.0; //sqrt(-2.0*lin_interp(x_DSgrid, phi_DSgrid, sqrt(mue_cut_lookup[i]), size_phiDSgrid, 1234)); 
					//vpar_e_cut[i] = sqrt(-2.0*lin_interp(x_DSgrid, phi_DSgrid, sqrt(2.0*mu_e[i]), size_phiDSgrid, 1234)); 
					//printf("vpar = %f at mu = %f\n", vpar_cut_lookup[i], mue_cut_lookup[i]);
				denszeroorb(-1.0, 1.0, phi_DSgrid, ne_DSgrid, size_phiDSgrid, &flux_eDS, dist_e_DK, vpar_e_DS, mu_e, size_vpar_e, size_mu_e, vpar_e_cut, gamma_ref, x_DSgrid);
				for (i=0; i< size_mu_e; i++) 
					vpar_e_cut[i] = sqrt(-2.0*lin_interp(x_DSgrid, phi_DSgrid, sqrt(2.0*mu_e[i]), size_phiDSgrid, 1234)); 
					//vpar_e_cut[i] = sqrt(-2.0*phi_DSgrid[0]); 
					//vpar_e_cut[i] = v_cutDS;
				i=0;
				while (ne_DSgrid[i] < 0.98) i++;
				size_neDSgrid = i;
				//printf("flux_eDS = %f\n", flux_eDS);
				//printf("%f\n", flux_eDS*ne_grid[0]*sqrt(2.0));
				//printf("flux_e = %f\n", flux_e);
				//exit(1);
			}


			error_Poisson(error_MP, x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, size_ngrid, 0.0);
			printf("error_av = %f\terror_max = %f\n", error_MP[0], error_MP[1]);
			if ( (error_MP[0] < tol_MP[0]) && (error_MP[1] < tol_MP[1]) ) convergence_MP += 1 ;
			else convergence_MP = 0;
			if (convergence_MP != -1) {
				printf("MP not converged --> calculate new MP potential guess\n");
				fprintf(fout, "MP not converged --> calculate new MP potential guess\n");
				newguess(x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, &size_ngrid, 0.0, deltaphi, 1.5, weight_MP);// p, m);
			}
			else {
				printf("MP converged --> no iteration needed\n");
				fprintf(fout, "MP converged --> no iteration needed\n");
			}

			//flux_e = flux_eDS*ne_grid[0];
			current = flux_i - flux_e*sqrt(Te*mioverme);
			if (fix_current == 1) {
				if (fabs((target_current - current)/flux_i) > tol_current) convergence_j = 0;
				else convergence_j += 1;
				printf("target current = %f +/- %f\n", target_current, tol_current*flux_i);
				fprintf(fout, "target current = %f +/- %f\n", target_current, tol_current*flux_i);
				printf("current = %f = ion_current (%f) - electron_current (%f) = %f x ion_current\n", current, flux_i, flux_e*sqrt(Te*mioverme), current/flux_i);
				fprintf(fout, "current = %f = ion_current (%f) - electron_current (%f) = %f x ion_current\n", current, flux_i, flux_e*sqrt(Te*mioverme), current/flux_i);
				if (convergence_j != -1) {
					fprintf(fout, "new WALL potential guess\n");
					newvcut(&v_cut, v_cutDS, mioverme, TioverTe[0], flux_i, flux_e*sqrt(Te*mioverme), target_current, tol_current, weight_j);
				}
				else {
					printf("current converged, NO wall potential iteration\n");
					fprintf(fout, "current converged, NO wall potential iteration\n");
				}
			}
			else {
				printf("current = %f = %f x ion_current", current, current/flux_i);
				fprintf(fout, "current = %f = %f x ion_current", current, current/flux_i);
			}

			printf("new DS potential guess\n");
			if (v_cut*v_cut < - 2.0*phi_grid[0]) {
				printf("WARNING: The total sheath and presheath potential drop is smaller than the presheath potential drop\n");
				printf("\tTo avoid a non-monotonic potential, setting the total potential drop to be just above the presheath potential drop\n");
				v_cut = sqrt(-2.0*phi_grid[0] + 0.01);
				v_cutDS = sqrt(0.01);
			}
			else v_cutDS = sqrt(2.0*phi_grid[0] + v_cut*v_cut);
			error_Poisson(error_DS, x_DSgrid, ne_DSgrid, ni_DSgrid, phi_DSgrid, size_phiDSgrid, size_neDSgrid, 1.0/(gamma_ref*gamma_ref));
			printf("error_av = %f\terror_max = %f\n", error_DS[0], error_DS[1]);
			if ( (error_DS[0] < tol_DS[0]) && (error_DS[1] < tol_DS[1]) ) convergence_DS += 1 ;
			else convergence_DS = 0;
			if (convergence_DS != -1) {
				printf("DS not converged --> calculate new DS potential guess\n");
				fprintf(fout, "DS not converged --> calculate new DS potential guess\n");
				newguess(x_DSgrid, ne_DSgrid, ni_DSgrid, phi_DSgrid, size_phiDSgrid, &size_neDSgrid, 1.0/(gamma_ref*gamma_ref), v_cutDS, 2.0, weight_MP);// p, m);
			}
			else {
				printf("DS converged --> no iteration needed\n");
				fprintf(fout, "DS converged --> no iteration needed\n");
			}

			//newguess(x_DSgrid, ne_DSgrid, ni_DSgrid, phi_DSgrid, size_phiDSgrid, &size_neDSgrid, 1/(gamma_ref*gamma_ref), sqrt(2.0*phi_grid[0] + v_cut*v_cut), 2.0, 1.0);//, p, m);
			//printf("\tcurrent = %f (target = %f)\n", current, target_current);
			//fprintf(fout, "\tcurrent = %f (target = %f)\n", current, target_current);
			printf("\t(phi_mp0, phi_ds0, phi_wall) = (%f, %f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut - phi_grid[0], -0.5*v_cut*v_cut);
			fprintf(fout, "\t(phi_mp0, phi_ds0, phi_wall) = (%f, %f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut - phi_grid[0], -0.5*v_cut*v_cut);
			printf("ITERATION N = %d of which N_DS = %d\n\n", N, N_DS);
			fprintf(fout, "ITERATION N = %d of which N_DS = %d\n\n", N, N_DS);
			printf("convergence (MP,DS, j) = (%d, %d, %d)\n\tfull convergence achieved when all of these numbers are 2 or above\n", convergence_MP, convergence_DS, convergence_j);
			fprintf(fout, "convergence (MP,DS, j) = (%d, %d, %d)\n\tfull convergence achieved when all of these numbers are 2 or above\n", convergence_MP, convergence_DS, convergence_j);
			//if (error_MP[0] > olderror_MP) weight_MP /= 1.05;
			N_DS++; N++;
			//gsl_permutation_free (p);
			//gsl_matrix_free (m);
		}
		printf("FINAL CHECK on accuracy of electrostatic potential solution\n");
		fprintf(fout, "FINAL CHECK on accuracy of electrostatic potential solution\n");
		densfinorb(Te, alpha, size_phigrid, &size_ngrid, ni_grid, x_grid, phi_grid, ioncharge, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, &flux_i, ZOOM_MP, 0.999, -999.9, vy_i_wall, chiM_i, twopidmudvy_i); 
		densfinorb(1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e, U_e_DS, size_mu_e, size_vpar_e, 0.0, &flux_eDS, zoomfactor_DS, 0.999, -999.9, vy_e_wall, chiM_e, twopidmudvy_e); 
		error_Poisson(error_MP, x_grid, ne_grid, ni_grid, phi_grid, size_phigrid, size_ngrid, 0.0);
		printf("error_av = %f\terror_max = %f\n", error_MP[0], error_MP[1]);
		if ( error_MP[1] > tol_MP[1] ) {
			printf("ERROR: MP solution rejected because it does not satisfy Poisson's equation accurately enough on the extended domain\n");
			fprintf(fout, "ERROR: MP solution rejected because it does not satisfy Poisson's equation accurately enough on the extended domain\n");
			//exit(-1);
		}
		error_Poisson(error_DS, x_DSgrid, ne_DSgrid, ni_DSgrid, phi_DSgrid, size_phiDSgrid, size_neDSgrid, 1.0/(gamma_ref*gamma_ref));
		printf("error_av = %f\terror_max = %f\n", error_DS[0], error_DS[1]);
		i=1;
		//while (i < size_neDSgrid) {
		//	if (ni_DSgrid[i] - ne_DSgrid[i] - ni_DSgrid[i-1] + ne_DSgrid[i-1] > tol_DS[1]) {
		//		printf("non-monotonicity = %f\n", ni_DSgrid[i] - ne_DSgrid[i] - ni_DSgrid[i-1] + ne_DSgrid[i-1]);
		//		printf("ERROR: DS solution rejected because it does not satisfy Poisson's equation accurately enough on the extended domain #2\n");
		//		fprintf(fout, "ERROR: DS solution rejected because it does not satisfy Poisson's equation accurately enough on the extended domain #2\n");
		//		i = size_neDSgrid;
		//	}
		//	i++;
		//}
		if ( error_DS[1] > tol_DS[1]) {
			// IDEA: check instead for non-monotonicity in n_i-n_e to accept or reject a Debye sheath solution
			printf("ERROR: DS solution rejected because it does not satisfy Poisson's equation accurately enough on the extended domain\n");
			fprintf(fout, "ERROR: DS solution rejected because it does not satisfy Poisson's equation accurately enough on the extended domain\n");
			//exit(-1);
		}
		printf("FINAL CHECK passed. HURRAY!\n");
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
	snprintf(fpstr, 150, "%s/vparcut.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file\n");
	}
	for (i=0; i<=size_mu_e; i++) {
		fprintf(fp, "%f %f\n", mu_e[i], vpar_e_cut[i]);
	}
	fclose(fp);
	snprintf(fpstr, 150, "%s/inputfile.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file\n");
	}
	fprintf(fp, "%f\n", alpha);
	fprintf(fp, "%f\n", gamma_ref);
	fprintf(fp, "%d\n", number_species);
	fprintf(fp, "%f\n", TioverTe[0]);
	fprintf(fp, "%f\n", mioverme);
	fprintf(fp, "%d\n", fix_current);
	if (fix_current == 1) fprintf(fp, "%f\n", target_current);
	else fprintf(fp, "%f\n", v_cut);
	fclose(fp);

	free(vpar_e_cut);
	free(x_grid);
	free(phi_grid);
	free(ne_grid);
	free(ni_grid);
	free(x_DSgrid);
	free(phi_DSgrid);
	free(ne_DSgrid);
	free(ni_DSgrid);
	free(U_e_DS);
	free(vpar_e_DS);
	exit(0);
}



	//if (output_type == 1)
	//{
	//	dvxopen = (5.0/dim_output)*sqrt(Telarge);
	//	Ucap = 18.0 + Te;
	//	/* The part below performs the integration over the open orbit distribution function at x=0 in a way that allows to extract the distribution function in v_x (because v_x is the only velocity component that matters in the Debye sheath).*/
	//	FILE *output, *outputyz, *outputyzvy, *outputyzvz;
	//	if ((output = fopen("OUTPUT/f0x.txt", "w")) == NULL)
	//	{        // Check for presence of file
	//		printf("Cannot open %s\n", "outputfile.txt");
	//		exit(EXIT_FAILURE);
	//	}
	//	if ((outputyz = fopen("OUTPUT/f0yz.txt", "w")) == NULL)
	//	{
	//		// Check for presence of file
	//		printf("Cannot open %s\n", "outputyz.txt");
	//		exit(EXIT_FAILURE);
	//	}
	//	if ((outputyzvy = fopen("OUTPUT/vy0.txt", "w")) == NULL)
	//	{
	//		// Check for presence of file
	//		printf("Cannot open %s\n", "outputyzvy.txt");
	//		exit(EXIT_FAILURE);
	//	}
	//	if ((outputyzvz = fopen("OUTPUT/vz0.txt", "w")) == NULL)
	//	{
	//		// Check for presence of file
	//		printf("Cannot open %s\n", "outputyzvz.txt");
	//		exit(EXIT_FAILURE);
	//	}
	//	intdxbaropen = 0.0;
	//	intdvxopen = 0.0;
	//	intdUopen = 0.0;
	//	F = 0.0;
	//	dvzopen = 0.1;
	//	sizeU = sqrt(2.0*Ucap)/dvzopen;
	//	if (DEBUG == 1) printf("maxj = %d, sizeU = %d\n", maxj, sizeU);
	//	for (i=0; i<dim_output; i++) {
	//		vxopen = i*dvxopen;
	//		intdxbaropenold = intdxbaropen;
	//		intdxbaropen = 0.0;
	//		for (j=0; j < maxj; j++)
	//		{	
	//			if (i==0)
	//			{
	//				fprintf(outputyzvy, "%f\n", xbar[j]);
	//			}
	//			intdUopenold = intdUopen;	
	//			intdUopen = 0.0;
	//			F = 0.0;
	//			//sizeU = (int) ( sqrt(Ucap - pow(vxopen, 2.0) - chi[j][0])/dvzopen );
	//			for (l=0; l<sizeU; l++)
	//			{	
	//				vz = dvzopen*l;
	//				if (j == 0 && i == 0)
	//					fprintf(outputyzvz, "%f\n", vz);
	//				Fold = F;
	//				vx0open = sqrt(2.0*(chiMax[j] - chi[j][0]));
	//				U = 0.5*vz*vz + 0.5*vx0open*vx0open + chi[j][0];
	//				Deltavx = sqrt(vx0open*vx0open + alpha*vz*openorbit[j]) -  vx0open;
	//				aa = tophat(vx0open, vx0open + Deltavx, vxopen);
	//				//Fopen = bilin_interp(mu[j][0], U, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
	//				Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1);
	//				//if (vx0open < dvxopen) 
	//				//printf("vx0open = %f\t Fopen = %f\n", vx0open, Fopen);
	//				F = Fopen*aa;
	//				if (i==0) // just because these have to be printed only once, and are independent of i
	//				{
	//					if (l==sizeU - 1 && i == 0)
	//					{	
	//						fprintf(outputyz, "%f\n", Fopen*Deltavx); 
	//					}
	//					else 
	//					{	
	//						fprintf(outputyz, "%f ", Fopen*Deltavx); 
	//					}
	//				}
	//				if (DEBUG == 1)
	//					printf("tophat = %f when vx0open = %f, Deltavx = %f and vxopen = %f\n, F = %f\n", aa, vx0open, Deltavx, vxopen, F);
	//				if (l!=0) {	
	//					intdUopen += 2.0*(F+Fold)*dvzopen; 
	//				}
	//			}
	//			if (j != 0) {	
	//				intdxbaropen += 0.5*(intdUopen + intdUopenold)*(xbar[j] - xbar[j-1]); 
	//			} 
	//		}
	//		if (DEBUG == 1)
	//			printf("FBohm = %f\n", intdxbaropen);
	//		if (i!=0)
	//		{
	//			intdvxopen += 0.5*(intdxbaropen + intdxbaropenold)*dvxopen;
	//			intdvxopenfluidBohm += 0.5*(intdxbaropen*vxopen + intdxbaropenold*(vxopen-dvxopen))*dvxopen;
	//			if (i>1)
	//			{
	//				intdvxopenBohm += 0.5*(intdxbaropen/pow(vxopen, 2.0) + intdxbaropenold/pow(vxopen-dvxopen, 2.0))*dvxopen;
	//			}
	//			else if (i==1)
	//			{
	//				intdvxopenBohm += 0.5*intdxbaropen/pow(vxopen, 2.0)*dvxopen;
	//			}
	//		}
	//		x_ax[i] = vxopen;
	//		y_ax[i] = intdxbaropen;
	//	}
	//	for (i=0; i<dim_output; i++) {	
	//		y_ax[i] /= intdvxopen;
	//		//printf("%f %f\n", x_ax[i], y_ax[i]);//n_grid[0]);
	//		fprintf(output, "%f %f\n", x_ax[i], y_ax[i]);
	//	}
	//
	//	fclose(output);
	//	fclose(outputyz);
	//	fclose(outputyzvy);
	//	fclose(outputyzvz);
	//
	//	// Calculating correction due to non boltzmann electrons
	//	// AG: to the Bohm condition
	//	/*the gradient of the function */
	//	//double ne_p0=0.0;
	//	//ne_p0 = (ne_grid[1] - ne_grid[0]) / (phi[1] - phi[0]);
	//	//printf("ne_p0 = %f (%f)\n", ne_p0, (ne_grid[1] - ne_grid[0]) / (phi_grid[1] - phi_grid[0]));
	//
	//	printf("flux calculated = %f\n", *flux);
	//	printf("in densfinorb: The density at x=0 obtained from the extracted distribution function is %f\nThe density at x=0 obtained from the ion density integral (more direct) is %f\n", intdvxopen, n_grid[0]);
	//	printf("Bohm integral = %f (obtained from the extracted distribution function at x=0)\nBohm integral = %f (obtained directly from the distribution function at infinity)\nThe Bohm condition is: Bohm integral = 2*dn_e/dphi/(Te*ni0)\n", intdvxopenBohm/intdvxopen, Bohm);//, (2.0*ne_p0)/(Te*n_grid[0])); 
	//	printf("in densfinorb: The flow at x=0 obtained from the extracted distribution function is %f\nThe flow at x=0 obtained from the distribution function at infinity is %f\nThe Bohm speed is %f\n", intdvxopenfluidBohm/intdvxopen, *flux, 1.0/sqrt(2.0)); 
	//	//printf("in densfinorb: kBohmsquared= %f (should be >0)\n", kBohmsquared);
	//
	//}
	//else if (output_type ==2) {
	//	FILE *output;
	//	if ((output = fopen("OUTPUT/Umucutfile.txt", "w")) == NULL) {       
	//		printf("Cannot open %s\n", "outputfile.txt");
	//		exit(EXIT_FAILURE);
	//	}
	//	if (DEBUG == 0) printf("muopen Ucritf\n");
	//	for (j=0;j<dim_output;j++) {
	//		if (muopen[maxj] < x_ax[j]) y_ax[j] = sqrt(2.0*Ucritf[maxj]);
	//		else y_ax[j] = sqrt(2.0*lin_interp(muopen, Ucritf, x_ax[j], maxj, 986));
	//		//else y_ax[j] = sqrt(2.0)*lin_interp(muopen, Ucritf, x_ax[j], maxj, 986);
	//		fprintf(output, "%f %f\n", x_ax[j], y_ax[j]);
	//		if (DEBUG == 1) printf("mu_crit = %f\tvpar_crit = %f\n", x_ax[j], y_ax[j]);
	//		if (DEBUG == 1) printf("mu_crit = %f\tU_crit = %f\n", x_ax[j], 0.5*y_ax[j]*y_ax[j]);
	//	}
	//	fclose(output);
	//}
