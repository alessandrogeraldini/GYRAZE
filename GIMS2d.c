// 2d version of Magnetic presheath / Chodura sheath code
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
#include "mps2d.h"
#define MAX_IT 500
// number of maximum iterations set to some large number but iterations should converge in about 20-100. If they don't, then there is a problem
#define INITIAL_GRID_PARAMETER 1.4
#define TEST 0
#define SYS_SIZ 18.0
#define ADHOCFS 1
#define MAXV 6.0
#define DV 0.2
#define EGYRO = 0
#define SIZEY 31

const double tol_MP[2] = {0.002, 0.01}, tol_current = 0.01;
const double WEIGHT_j = 0.1, WEIGHT_MP = 0.3, CAUTIOUSWEIGHT = 0.15;
const double STOP_MP = 0.95, GRIDSIZE_MP = 0.45;
const int ZOOM_MP = 3;

const char *strqty[6] = {"alpha=", "delta=", "Ti:Te=","mi:me=","jwall=","pwall="};
const int lenstrqty = 6;

void Figen(double mioverme, double TioverTe, double ***ffarr, double *phiinf, double *y_grid, double *muarr, double *vzinfarr, int sizey, int sizeU, int sizemu, double dvpar, double dvperp, double dy) {
	int k, i, j, coldelectrons;
	double mu, vzinf, ff, Te;
	double u, condition, chodura, normalization, flow=0.0;
	double *dphiinfdy = malloc(sizey*sizeof(double));
	dphiinfdy[0] = 0.0;
	dphiinfdy[sizey-1] = 0.0;
	for (k=1; k < sizey-1; k++) {
		dphiinfdy[k] = ( phiinf[k+1] - phiinf[k-1] ) / (2.0*y_grid[1]);
	}
	for (k=0; k<sizey; k++) {
		Te = 1.0/TioverTe;
		if (Te < 1e-9) coldelectrons = 1;
		else coldelectrons = 0;
		//printf("are electrons cold (1=yes, 0=no): %d\n", coldelectrons);
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
					//printf("normalization = %f\n", normalization);
					chodura = (1+erf(u/sqrt(2.0)))/(normalization);
					flow = ( 0.5*(1+0.5*u*u)*exp(-0.5*u*u) + (1.5 + 0.5*u*u)*(sqrt(M_PI/2.0)*u/2.0)*(1 + erf(u/sqrt(2.0))))*4.0/sqrt(2.0*M_PI)/normalization;
					//printf("%f, %f, %f, %f %f\n", normalization, u, flow, chodura, condition);
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
					//printf("normalization = %f\n", normalization);
					chodura = 2.0*M_PI*exp(1.0/u)*(1.0-erf(1.0/sqrt(u)))/(2.0*sqrt(u)*normalization);
					flow = -9999.9;
					//printf("%f, %f, %f, %f %f\n", normalization, u, flow, chodura, condition);
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
		printf("k=%d\tnormalization = %f\n", k, normalization);

		if (coldelectrons == 0) {
			if (Te>=1.0) {
				//printf("Te = %f >= 1\n", Te);
				for (i=0; i<sizemu; i++) {
					mu = i*dvperp*i*dvperp;
					muarr[i] = mu;
					for (j=0; j<sizeU; j++) { 
						vzinf = j*dvpar;
						if (i==0)
							vzinfarr[j] = vzinf;
// 						solve for vz given U, mu and ystar
						//while (fabs(vz - vzold) > TINY) {
						//	vzold = vz;
						//	vz = sqrt(2.0*(U-mu-phiinf[k]));
						//}
						ff = exp(phiinf[k])*(1.0/(M_PI*sqrt(M_PI)))*(1.0/normalization)*vzinf*vzinf*exp(- 0.5*vzinf*vzinf - mu + u*vzinf - 0.5*u*u);
						if ((i==0) && (j==1)) printf("k=%d\tvzinf=%f\tmu=%f\tff=%f\n", k, vzinf, mu, ff);
						ffarr[k][i][j] = ff;
					}
				}
			}
			else {
				//printf("Te = %f < 1\n", Te);
				for (i=0; i<sizemu; i++) {
					mu = i*dvperp*i*dvperp;
					muarr[i] = mu;
					for (j=0; j<sizeU; j++) { 
						vzinf = j*dvpar;
						if (i==0)
							vzinfarr[j] = vzinf;
						ff = exp(phiinf[k])*(1.0/(4.0*M_PI))*(1.0/normalization)*(0.5*vzinf*vzinf/(1+u*0.5*vzinf*vzinf))*exp(-0.5*vzinf*vzinf - mu);
						ffarr[k][i][j] = ff;
					}
				}
			}
		}
		else { // coldelectrons == 1:   // for the moment not active in 2d code
			//printf("Te = infinity\n");
			normalization = 1.0;
			for (i=0; i<sizemu; i++) {
				mu = i*dvperp*i*dvperp;
				muarr[i] = mu;
				for (j=0; j<sizeU; j++) { 
					vzinf = j*dvpar;
					if (i==0)
						vzinfarr[j] = vzinf;
					ff = (1.0/(2.0*(M_PI)*sqrt(M_PI)))*exp(- 0.5*vzinf*vzinf - mu );
					ffarr[k][i][j] = ff;
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
	
void make_phiinfgrid(int shape, double phiamp, double delta, double *y_grid, double *phiinf_grid, int size_ygrid) {
	int i;
	double ymid=(y_grid[0] + y_grid[size_ygrid-1])*0.5;
	if (shape == 0) {
		// this would be a single positive or negative perturbation (eg a Gaussian shape)
		for (i=0; i<size_ygrid; i++) {
			phiinf_grid[i] = phiamp*exp(-4.0*pow(delta*(y_grid[i] - ymid), 2.0));
			printf("phiinf[%d] = %f\n", i, phiinf_grid[i]);
		}
	}
	else if (shape == 1) {
		// an oscillatory perturbation
		for (i=0; i<size_ygrid; i++) {
			phiinf_grid[i] = phiamp*sin(delta*y_grid[i]);
		}
	}
	return;
}

void make_phigrid(double *x_grid, double *y_grid, double **phi_grid, double *phiinf, int size_xgrid, double deltax, int size_ygrid, double grid_parameter) {
	int i, j;
	double xi, *ff, *new_phi, *new_x, power = 0.5;
	//new_phi = malloc(size_ygrid*sizeof(double));
	//for (i=0; i<size_ygrid; i++) 
	new_phi = malloc(size_xgrid*sizeof(double)); 
	new_x = malloc(size_xgrid*sizeof(double)); 
  	ff = malloc(size_xgrid*sizeof(double)); 
	printf("in make_phigrid: size_xgrid = %d\n", size_xgrid);
	for (i=0; i < size_xgrid; i++) {
		xi = i*deltax;
		ff[i] = xi;
		new_x[i] = pow( pow(grid_parameter+xi, power) - pow(grid_parameter, power), 2.0);
		new_phi[i] = 0.0;
	}
	//phi_fluidMPS(new_phi, new_x, size_xgrid, alpha);
	for (i=0; i < size_xgrid; i++) {
		x_grid[i] = new_x[i];
		//printf("size_phigrid = %d\n", size_ygrid);
		for (j=0; j < size_ygrid; j++) {
			//printf("phi_grid = %f\n", phi_grid[j][i]);
			phi_grid[j][i] = new_phi[i] + phiinf[j]; //*exp(-1.0/(x_grid[i] + TINY))*;
		}
	}
	free(new_phi);
	free(new_x);
	free(ff);
	printf("exit make_phigrid");
	return;
}

// this function is no longer used, but keep it just in case
//void rescale_array(double *array, int size_array, double jump) {
//	int i;
//	double old_jump = array[0];
//	for (i=0; i < size_array; i++) 
//		array[i] *= (jump/old_jump);
//	return;
//}
//

void iondens2d(double Te, double alpha, int size_ygrid, int size_phigrid, int *size_ngrid, double **n_grid, double *y_grid, double *x_grid, double **phi_grid, double *phiinf,  double ***FF, double *mumu, double *UU, int sizemumu, int sizeUU, double grid_parameter, double *flux, int zoomfactor, double stopdens, double **vy_op, double **chiMax_op, double **dmudvy_op) {
	// declare variables
	clock_t begin = clock(); // Finds start time of computation
	double limit_rho = 8.0, *n_inf, n_norm, ystar = 0.0, vzmin ;
	double deltax, deltay, deltax_inf, xi, fluxinf1, densinf1, densinf1old, twopidmudxbarantycal, densinfnorm;
	int p=0, i=0, ic=0, j=0, k=0, l=0;
	int s, w, ind, stop = 0, sizexbar, j_inf, sizeU, size_finegrid;
	double oorbintgrd, oorbintgrdantycal, oorbintgrdBohm, square;
	double vz, U, dvz = 0.2, dvzopen = 0.2, dvx, dxbar, intdU=0.0, intdUopen=0.0, intdUopenBohm = 0.0, intdUflow=0.0, intdUflowold=0.0;
	double intdUold=0.0, intdvx=0.0, intdvxold = 0.0, intdvxflow=0.0, intdvxflowold=0.0, intdxbar=0.0, intdxbarflow = 0.0, intdxbaropen=0.0, intdxbaropenBohm = 0.0, F, Fold=0.0, Ucap, Bohm;
	double intdUopenflow = 0.0, intdUopenflowold = 0.0, intdxbaropenflow = 0.0, oorbintgrdflow = 0.0, oorbintgrdflowold = 0.0;
	double oorbintgrdold=0.0, oorbintgrdBohmold=0.0, Fopen=0.0, intdUopenold=0.0, intdUopenBohmold=0.0;
	double intdUantycal=0.0, intdvxantycal=0.0, vxnew=0.0, vxold = 0.0, Uperpnew = 0.0, **xtop, intdUopenantycal=0.0;
	double openorbitnew, yopenorbitnew, chinew, munew = 0.0, vx0open;
	double oorbintgrdxbar, oorbintgrdsquare, intdUopenxbar = 0.0, intdUopensquare = 0.0, intdxbarxbar=0.0, intdxbarsquare=0.0;
	double oorbintgrdxbarold, oorbintgrdsquareold, intdUopenxbarold, intdUopensquareold;
	double du, fluxinf1old, fluxinfintgrdold, fluxinfintgrd, u, Chodura2, Chodura2old, Chodura1old, Chodura1, *Chodura;
	// arrays to be allocated
	double **chiinf, **muinf, ***vxinf;
	double *gg, *ff, *xx, **phi, **phi_x, **phi_y, ***chi;
	double **twopidmudxbar, **muopen, **chiMopen, **xbaropen, **twopidmudxbaropen, **twopidmudy, **twopidmudyopen, **xbar, **xifunction;
	double  *xbarcrit, *chiMcrit;
	int *maxj, *icrit;
	int **jmclosed, **jmopen;
	int **crossed_max, **crossed_min, **lowerlimit, **upperlimit, **imax, **imin;
	int ***upper;
	double **chiMax, **chimpp, **chimin;
	double ***Uperp, ***mu;
	double ****vx;
	double *fluxinf, *fluxinf2, *fluxzero, *densinf, intfluxinf=0.0, intfluxzero=0.0; 
	double phiinfval, *dphiinfdy, dphiinfdyval, vzinf, yinf;

	dphiinfdy = malloc(size_ygrid*sizeof(double));
	dphiinfdy[0] = 0.0;
	dphiinfdy[size_ygrid-1] = 0.0;
	for (p=1; p<size_ygrid-1; p++)  {
		dphiinfdy[p] = (phiinf[p+1] - phiinf[p-1])/(y_grid[1]*2.0);
		printf("dphiinfdy[%d] = %f\n", p, dphiinfdy[p]);
	}

/* 
 * Note: first index of FF and FFprime is y, second one is mu, third one is U-mu
 * xx is the array containing the position (distance from the wall) on the fine grid; phi contains the values of phi(x), extracted from a file. phi_x is phi prime, first derivative of phi. n; newphi is the new electrostatic potential guess; chi is the effective potential
 * n is the domain of the position x (the largest value of the index i, so e.g. x_grid[n] = L_2 in the paper); sizexbar is the domain of the position xbar (the largest value of the index j);
 * Uperp stores the possible values of Uperp associated with closed orbits, and so does vx; chiMax and chimin store the local maxima and minima of the effective potential maximum
 * twopidmudxbar is the Delta_M = 2*pi*dmu/dxbar; twopidmudxbarantycal is the analytical value of openorbit for a flat potential; mu is the array containing values of mu(xbar, Uperp), index j for values of xbar, k for values of Uperp; xbar is the grid of values used in the closed orbit integral; FF contains the distribution function, read from the file distfile.txt. UU and mumu contain the values of U and mu corresponding to the function FF (which is F(mu, U)); FFprime is the numerical first derivative of F with respect to U; xifunction is the function of x which defines the grid of values of xbar by finding a chi whose minimum lies exactly at each grid point x. 
 * crossed_max and crossed_min is non-zero when a minimum (or maximum) of chi is found for some xbar[j]
 * jmclosed and jmopen represent minimum values of xbar above which we integrate open and closed orbit density integrals respectively (xbar_m,o and xbar_m in the paper); i is an index usually representing the positin x  j is an index usually representing the orbit position xbar; k is an index usually representing the energy Uperp (or velocity vx). It's always used in conjunction with j (and sometimes i); l is an index (used in for loops) usually representing the total energy (or velocity vz). It's used only in the DENSITY INTEGRALS part of the code; sizeU is the size of the integration range over U (or velocity vz). It is set later on in the code
 * lowerlimit represents the lower limit of k in the integrals over Uperp (or vx). It's needed because some of the earlies energies; (which are the largest because thy are values of chi stored after the maximum is found); may be so large that they are associated with very small values of the distribution function. This avoids integrating in an empty portion of phase space; upperlimit[j] represents the largest value of k (the smallest stored energy Uperp = chi_minimum) associated with some value of j; upper[j][i] represents the value of k associated with the smallest value of vx when integrating over Uperp. Going above upperlimit[j][i] makes Uperp < chi so velocities imaginary; imax/imin[j] stores the position of the maximum/minimum of the effective potential chi (It's x_M/x_m in the paper, which depends on xbar).
 * vz used in the density integral; U is the total energy, used in the density integral; dvz is the thickness of the vz grid used to take the integral over U (which is taken over vz in practice), dvzopen is the same for the open orbit piece; dvx is the thickness of the vx grid used to take the integral over Uperp ( which is taken over vx in practice). It must be evaluated because it depends on stored values of vx[j][i][k]; dxbar is the thickness of the xbar grid; intdU is the value of the integral over U in the closed orbit density integration process; intdUopen is the same as above, for the open orbit integral; intdUopenBohm same, for Bohm integral
 * oorbintgrd is the value of the integrand in the first open orbit integral (oorbintgrdantycal is the analytical result for flat potential); oorbintgrdBohm is the value of the integrand in the `Bohm' integral 
 * intdUold is a variable which stores the old intdU, so that the trapezium rule of integration can be applied (intdUold + intdU)*dvz; intdvx stores the integral over Uperp (hence over vx) in the closed orbit integral; intdxbar stores the value of the closed orbitintegral over xbar (which is the final result!), intdxbaropen does the same in the open orbit density integral; intdxbaropenBohm does the same for the Bohm integral; idealBohm is what the Bohm integral shoult be if Bohm condition is marginally satisfied; F is the value of the distribution function evaluated in the density integrals by interpolating FF, and Fold is the `old' needed to apply the trapezium rule; Fprime is the trilinearly interpolated value of FFprime, and Fprimeold is the same at the previous grid point (needed for trapezium rule); used in INTEGRALS OF DISTRIBUTION FUNCTION AT INFINITY; Ucap is the topmost total energy integrated to 
 * intdUantycal is the integral over U (or v_z) for a flat potential profile (phi =0) for some value of xbar and Uperp; intdvxantycal  is the integral over Uperp (or vx) for a flat potential profile for some value of xbar; vxnew is the value of vx at the 'new' grid point, used in the vx integral (taken using the trapezium rule); vxold is the value of vx at the 'old' grid point, used in the vx integral; Uperpnew is the value of Uperp (used in the closed orbit density integral); munew is the valye of mu (used in the closed orbit density integral); xtop is the top bounce point x_t of the last closed orbit; intdUopenantycal is the analytical value of the integral over U  in the open orbit density integral 
 * oorbintgrdxbar is an integral over the open orbit distribution function at x=0 which is needed to evaluate a coefficient that appears when; expanding quasineutrality near x=0. It was just for playing around and at the moment plays no role in the code; similarly with all other integrals here 
*/

	//make 2d phi grid which is same resolution or finer in x than input phi_grid

	// Make grids in g = sqrt(x) and the equidistant variable f
	gg = malloc(size_phigrid*sizeof(double)); 
	ff = malloc(size_phigrid*sizeof(double));
	for (i=0; i<size_phigrid; i++) {
		gg[i] = sqrt(x_grid[i]);
		ff[i] = pow(sqrt(grid_parameter)+gg[i], 2.0) - grid_parameter;
	}
	deltax = ff[1];
	deltay = y_grid[1] - y_grid[0];
	size_finegrid = zoomfactor*size_phigrid-zoomfactor;
	phi = malloc(size_ygrid*sizeof(double));
	for (p=0; p<size_ygrid; p++) 
		phi[p] = malloc(size_finegrid*sizeof(double)); // phi now has correct size
	printf("size_finegrid = %d\t it is %d x %d - %d\n", size_finegrid, zoomfactor, size_phigrid, zoomfactor);
	xx = malloc(size_finegrid*sizeof(double)); // Make finer x grid than input x_grid
	if (zoomfactor != 1) {
		for (p=0; p<size_ygrid;p++) {
			gsl_interp_accel *acc = gsl_interp_accel_alloc ();
			gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, size_phigrid);
			gsl_spline_init (spline, ff, phi_grid[p], size_phigrid);
			//gsl_spline_init (spline, gg, phi_grid, size_phigrid);
			for (i = 0; i < size_finegrid; i += 1) {
				xi = i*deltax/zoomfactor;
				if (i == 0) xi += 0.00001; //TINY;
				if (i == size_finegrid-1) xi -= 0.00001; //TINY;
				xx[i] = pow( pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), 2.0);
				//printf("i=%d/%d\txx = %f\n", i, size_finegrid, xx[i]);
				phi[p][i] = gsl_spline_eval (spline, xi, acc);
				phi[p][i] *= Te;
			}
			gsl_spline_free (spline);
			gsl_interp_accel_free (acc);
		}

	}
	else {
		for (p=0; p<size_ygrid;p++) {
			for (i = 0; i < size_finegrid; i += 1) {
				xi = i*deltax/zoomfactor;
				xx[i] = x_grid[i];
				phi[p][i] = phi_grid[p][i];
				phi[p][i] *= (Te);
			}
		}
	}
	printf("Te= %f\t\n", Te);
	//for (p=0; p<size_ygrid;p++)
	//	printf("phi[%d][0] = %f\tphi_grid[%d][0] = %f\n", p, phi[p][0], p, phi_grid[p][0]);

	// Introduce a cap in energy (U, Uperp) high enough that we can safely assume F = 0
	Ucap = 12.0 + 6.0*Te;
	for (p=0; p<size_ygrid; p++) {
		if ( (trilin_interp(y_grid[p], 0.0, sqrt(2.0*Ucap), FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1) > 1e-10) || ( trilin_interp(y_grid[p], Ucap, 0.0, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1) > 1e-10 ) ) {
			printf("ERROR in iondens2d_renorm.c: increase Ucap please because F values are %f and %f\n", trilin_interp(y_grid[p], 0.0, sqrt(2.0*Ucap), FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1), trilin_interp(y_grid[p], Ucap, 0.0, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1)) ;
			exit(1);
		}
	}
	printf("Ucap has been checked to be large enough\n");

	/* We initialize all arrays that contain functions of position x with the correct size n */
	/* FORM XBAR GRIDS 
	Take derivatives of phi and use them to obtain two grids for xbar, one to be used for closed orbits and one to be used for open orbits. */

	jmclosed    = (int**)      malloc(size_ygrid*sizeof(int*));
	jmopen      = (int**)      malloc(size_ygrid*sizeof(int*));
	upper       = (int***)     malloc(size_ygrid*sizeof(int**)); 
	crossed_min = (int**)      malloc(size_ygrid*sizeof(int*)); 
	crossed_max = (int**)      malloc(size_ygrid*sizeof(int*));
	upperlimit  = (int**)      malloc(size_ygrid*sizeof(int*));
	lowerlimit  = (int**)      malloc(size_ygrid*sizeof(int*));
	imax        = (int**)      malloc(size_ygrid*sizeof(int*));
	imin        = (int**)      malloc(size_ygrid*sizeof(int*));
	icrit       = (int*)       malloc(size_ygrid*sizeof(int));
	maxj        = (int*)       malloc(size_ygrid*sizeof(int));

	phi_x       = (double**)   malloc(size_ygrid*sizeof(double*)); 
	phi_y       = (double**)   malloc(size_ygrid*sizeof(double*)); 
	xifunction  = (double**)   malloc(size_ygrid*sizeof(double*));
	xbar        = (double**)   malloc(size_ygrid*sizeof(double));
	chi         = (double***)  malloc(size_ygrid*sizeof(double**)); // chi(xbar, x) indices j and i 
	Uperp       = (double***)  malloc(size_ygrid*sizeof(double**)); 
	mu          = (double***)  calloc(size_ygrid,sizeof(double**)); 
	vx          = (double****) malloc(size_ygrid*sizeof(double***)); 
	chimin      = (double**)   malloc(size_ygrid*sizeof(double*)); 
	chiMax      = (double**)   malloc(size_ygrid*sizeof(double*));
	chimpp      = (double**)   malloc(size_ygrid*sizeof(double*));
	twopidmudxbar   = (double**)   malloc(size_ygrid*sizeof(double*));
	twopidmudy   = (double**)   malloc(size_ygrid*sizeof(double*));
	xtop        = (double**)   malloc(size_ygrid*sizeof(double*));
	n_inf       = (double*)    malloc(size_ygrid*sizeof(double));
	chiMcrit    = (double*)    malloc(size_ygrid*sizeof(double));
	xbarcrit    = (double*)    malloc(size_ygrid*sizeof(double));
	muopen      = (double**)   malloc(size_ygrid*sizeof(double*));
	chiMopen    = (double**)   malloc(size_ygrid*sizeof(double*));
	xbaropen    = (double**)   malloc(size_ygrid*sizeof(double*));
	twopidmudxbaropen = (double**)   malloc(size_ygrid*sizeof(double*));
	twopidmudyopen = (double**)   malloc(size_ygrid*sizeof(double*));
	vxinf       = (double***)   malloc(size_ygrid*sizeof(double**));
	chiinf      = (double**)   malloc(size_ygrid*sizeof(double*));
	muinf       = (double**)   malloc(size_ygrid*sizeof(double*));

	densinf = malloc(size_ygrid*sizeof(double));
	fluxinf = malloc(size_ygrid*sizeof(double));
	fluxinf2 = malloc(size_ygrid*sizeof(double));
	fluxzero = malloc(size_ygrid*sizeof(double));
	Chodura = malloc(size_ygrid*sizeof(double));
	// Evaluate integrals of distribution function at infinity
	fluxinfintgrd = fluxinf1 = densinf1 = Chodura2 = Chodura1 = 0.0;
	du = 0.01;
	printf("sizemumu = %d\tsizeUU = %d\n", sizemumu, sizeUU);
	for (p=0; p<size_ygrid; p++) {
		densinf[p] = Chodura[p] = fluxinf[p] = 0.0;
		for (i=0; i<sizemumu; i++) {
			fluxinf1old = fluxinf1;
			Chodura1old = Chodura1;
			Chodura1 = 0.0;
			fluxinf1 = 0.0;
			densinf1old = densinf1;
			densinf1 = 0.0;
			//Fprimeold = Fprime = 0.0;
			Fold = F = 0.0;
			fluxinfintgrdold = fluxinfintgrd = 0.0;
			for (j=0; j< sizeUU; j++) {
				u = sqrt(2.0*UU[j]);
				fluxinfintgrdold = fluxinfintgrd;
				Chodura2old = Chodura2;
				Fold = F;
				F = FF[p][i][j];
				fluxinfintgrd = alpha*F*u;
				if (j>0) Chodura2 = F/(u*u);
				else Chodura2 = 0.0;
				if (j != 0) {
					du = u - sqrt(2.0*UU[j-1]); 
					fluxinf1 += 0.5*(fluxinfintgrd + fluxinfintgrdold)*du;
					Chodura1 += 0.5*(Chodura2 + Chodura2old)*du;
					densinf1 += 0.5*(F + Fold)*du;
				}
			}
			if (i!=0) {
				fluxinf[p] += 0.5*(mumu[i]-mumu[i-1])*(fluxinf1 + fluxinf1old);
				Chodura[p] += 0.5*(mumu[i]-mumu[i-1])*(Chodura1 + Chodura1old);
				densinf[p] += 0.5*(mumu[i]-mumu[i-1])*(densinf1 + densinf1old);
			}
		}
		densinf[p] *= (4.0*M_PI);
		fluxinf[p] *= (4.0*M_PI);
		Chodura[p] *= (4.0*M_PI);
		if (p==0) densinfnorm = densinf[0];
		fluxinf[p] /= densinfnorm;
		Chodura[p] /= densinfnorm;
		printf("at y index = %d, densinf = %f\tfluxinf = %f\tChodura = %f\n", p, densinf[p], fluxinf[p], Chodura[p]);

		xbar[p] = malloc(size_finegrid*sizeof(double));	
		phi_x[p] = malloc(size_finegrid*sizeof(double));
		phi_y[p] = malloc(size_finegrid*sizeof(double));
		xifunction[p]   = malloc(size_finegrid*sizeof(double));
		jmclosed[p] = (int*) malloc(size_finegrid*sizeof(int));
		jmopen[p]   = (int*) malloc(size_finegrid*sizeof(int));
		for (i=0; i<size_finegrid; i++) {
			// Evaluate derivative of phi
			jmopen[p][i] = jmclosed[p][i] = 0; //*
			if (i == 0) 	
				phi_x[p][0] = (phi[p][1] - phi[p][0])/(xx[1]-xx[0]);
			else if (i == size_finegrid-1) 
				phi_x[p][size_finegrid-1] = (phi[p][size_finegrid-1] - phi[p][size_finegrid-2])/(xx[size_finegrid-1]-xx[size_finegrid-2]);
			else 
				phi_x[p][i] = ((xx[i] - xx[i-1])/(xx[i+1]- xx[i-1]))*(phi[p][i+1] - phi[p][i])/(xx[i+1] - xx[i]) + ((xx[i+1] - xx[i])/(xx[i+1]- xx[i-1]))*(phi[p][i] - phi[p][i-1])/(xx[i] - xx[i-1]); 

			if (p>1) 
				phi_y[p-1][i] = (phi[p][i] - phi[p-2][i])/(2*deltay);
			// xifunction is xbar corresponding to given position x being a stationary point
			xifunction[p][i] = xx[i] + phi_x[p][i];
			//if (p==0)
			//	printf("p=%d/%d\ti=%d/%d\txifunction = %f\n", p, size_ygrid, i, size_finegrid, xifunction[p][i]);
			if (i == 1) {	
				if  (xifunction[p][i] > xifunction[p][i-1]) { // immediately found that xi is increasing at x=0 telling us x_c = 0
					icrit[p] = i-1;
				}
			}
			else if (i > 1) {	
				if ( (xifunction[p][i] > xifunction[p][i-1]) && (xifunction[p][i-1] < xifunction[p][i-2] ) ) 
					icrit[p] = i-1;
				else if ( (xifunction[p][i] < xifunction[p][i-1]) && (xifunction[p][i-1]  > xifunction[p][i-2]) ) {
				// found maximum of xi: this only happens when phi has noise in second derivative
					printf("ERROR in iondens2d: too much noise in second derivative\n");
					printf("xifunction = %f\n", xifunction[p][i-1]);
					printf("xx[%d] = %f\tphi_x = %f\n", i-1, xx[i-1], phi_x[p][i-1]);
				}
			} 
		}
		//printf("icrit = %d\n", icrit[p]);
	}
	//impose periodic boundary conditions in y
	for (i=0; i<size_finegrid; i++) {	
		phi_y[0][i] = (phi[1][i] - phi[size_ygrid-1][i])/(2*deltay);
		phi_y[size_ygrid-1][i] = (phi[0][i] - phi[size_ygrid-2][i])/(2*deltay);
	}

	//for (i=0; i<size_finegrid; i++) 
	//	printf("p=0/%d\ti=%d/%d\txifunction = %f\n", size_ygrid, i, size_finegrid, xifunction[0][i]);
	//FILE *fout; 
	//if ((fout = fopen("OUTPUT/iondens2d_out.txt", "w")) == NULL) {	
	//	printf("Cannot open iondens2d_out.txt");
	//	exit(EXIT_FAILURE);
	//}

	for (p=0; p<size_ygrid; p++) {
		//printf("p=%d/%d\n", p, size_ygrid);
		//printf("icrit=%d\n", icrit[p]);
		j=0; // set counting index to zero
		for (k=icrit[p]+1; k<size_finegrid; k++) {
			xbar[p][j] = xifunction[p][k]; 
			//printf("j=%d, k=%d, xbar = %f\n", j, k, xbar[p][j]);
			//if (p==0)
			//printf("p=%d/%d\tk=%d/%d\txifunction = %f\n", p, size_ygrid, k, size_finegrid, xifunction[p][k]);
			j++;
		}
		xbarcrit[p] = xifunction[p][icrit[p]];
		chiMcrit[p] = 0.5*phi_x[p][icrit[p]]*phi_x[p][icrit[p]] + phi[p][icrit[p]];
		if (p==0)
			sizexbar = j;
		maxj[p] = sizexbar + 1;
		if (DEBUG == 1) for (j=0;j<sizexbar; j++) printf("xbar[0][%d] = %f\n", j, xbar[0][j]); 
		//printf("sizexbar = %d, size_finegrid (x) =%d\n", sizexbar, size_finegrid);
		//
		// Lots of array allocations now that size of xbar (vy + x)
		chi[p] = (double **) calloc(sizexbar,sizeof(double*)); // chi(xbar, x) indices j and i 
		Uperp[p] = (double**)  calloc(sizexbar,sizeof(double*)); 
		mu[p] = (double**)  calloc(sizexbar,sizeof(double*)); 
		vx[p] = (double***) calloc(sizexbar,sizeof(double**)); 
		upper[p] = (int**)     calloc(sizexbar,sizeof(int*)); 
		chimin[p] = (double*)   calloc(sizexbar,sizeof(double)); 
		chiMax[p] = (double*)   calloc(sizexbar,sizeof(double));
		chimpp[p] = (double*)   calloc(sizexbar,sizeof(double));
		crossed_min[p] = (int*)      malloc(sizexbar*sizeof(int)); 
		crossed_max[p] = (int*)      malloc(sizexbar*sizeof(int));
		twopidmudxbar[p] = (double*)   calloc(sizexbar,sizeof(double));
		twopidmudy[p] = (double*)   calloc(sizexbar,sizeof(double));
		upperlimit[p] = (int*)      calloc(sizexbar,sizeof(int));
		lowerlimit[p] = (int*)      calloc(sizexbar,sizeof(int));
		xtop[p] = (double*)   calloc(sizexbar,sizeof(double));
		imax[p] = (int*)      calloc(sizexbar,sizeof(int));
		imin[p] = (int*)      calloc(sizexbar,sizeof(int));

		/////////////////////////////////////////////////////
		/* CLOSED ORBIT ARRAY FILLING
		Set up the grid in xbar and also initialize all arrays that contain a different number at different values of xbar, indexed j */
		/* Loop below initializes all 2d arrays which are functions of xbar and x. It allocates the right amount of memory to arrays of pointers of size xbar. The result is a 2D array indexed j (size sizexbar) and i or k (size n, see above) */
		for (j=0;j<sizexbar;j++)  {
			imax[p][j] = imin[p][j] = -1;
			twopidmudxbar[p][j] = 0.0;
			twopidmudy[p][j] = 0.0;
			crossed_min[p][j] = 0;
			crossed_max[p][j] = 0;
			xtop[p][j] = 0.0;
			chimin[p][j] = 0.0;
			chiMax[p][j] = 0.0;
			upperlimit[p][j] = -1;
			lowerlimit[p][j] = 0; 
			chi[p][j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing chi(xbar, x) indexes j and i 
			chi[p][j][0] =  0.5*pow((xx[0] - xbar[p][j]), 2.0) + phi[p][0];
			Uperp[p][j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing Uperp for closed orbits indexes j and k
			mu[p][j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing adiab invariant mu(xbar, Uperp), indexes j and k
			upper[p][j] = (int*)calloc(size_finegrid,sizeof(int)); // array containing the index of energy Uperp corresponding to chi(x) 
			vx[p][j] = (double**)calloc(size_finegrid,sizeof(double*)); // see below
			/* The loop below initializes the 3d array containing the velocity of a particle at a given orbit position xbar, with particle position x and with energy Uperp, indexed i, j and k. */
			for (i = 0; i < size_finegrid; i++) 
				vx[p][j][i] = (double*)calloc(size_finegrid,sizeof(double)); // 4D array containing value of vx at different y, xbar, x and Uperp
		}
		/* This for loop fills in the arrays, calculating the integrals where necessary */
		for (i=1; i<size_finegrid; i++) {
			for (j=0; j<sizexbar; j++) {
				chi[p][j][i] =  0.5*pow((xx[i] - xbar[p][j]), 2.0) + phi[p][i];
		//FINDING MAXIMA/MINIMA
		/* Below, we use the array elements to define arrays for the effective potential maxima and minima that exist for every xbar (index j). We search through the chi curve from x=0 to the largest value of x but we search through every chi curve first (so first scan in xbar at fixed x, then move to the next x). We create arrays of mu, Uperp and vx.
		Because I am comparing neighbouring values of x to find a maximum, I need to consider two distinct cases.The first one is treated in the if loop below, which the program should enter only if chi is decreasing at x=0. Implying that chi(x=0) is an effective potential maximum. The second one is the else if loop after that, which finds maxima of chi that are stationary points by comparing the value of the function before and after each point. Note that because we compare the function at a point with the function at points before and after, the point we consider at every iteration step in the index i is indexed (i-1), compared with (i-2) and i.  */
				if ( (i==1) && (chi[p][j][i] < chi[p][j][i-1]) ) {
					crossed_max[p][j] += 1;		
					imax[p][j] = 0;
					chiMax[p][j] = chi[p][j][0];
				} 
				else if ( (i > 1) && ((chi[p][j][i] < chi[p][j][i-1]) && (chi[p][j][i-1] > chi[p][j][i-2])) ) {
					crossed_max[p][j] += 1;
		// crossed_min and crossed_max let the code know whether we go up or down the effective potential well. We are starting to go down!
					if (crossed_max[p][j] > 1) {	
						printf("***WARNING*** There is more than one maximum!\n");
						printf("j is %d and maxima at %d and %d for first and second\n",  j, imax[p][j], i-1);
						for (k=imax[p][j]-1;k<i+1;k++) 
						printf("chi[%d][%d] = %f\n", j, k, chi[p][j][k]);
						crossed_max[p][j] = 1; 
						printf("exit code now\n");
						exit(-1);
					}
					imax[p][j] = i-1; // Store the index of the maximum 
					chiMax[p][j] = chi[p][j][i-1]; //Store chi Maximum itself, a function of xbar (index j)
				}
		/* We store the index corresponding to the position of a minimum for a given value of xbar.*/
				//watch out HERE
				else if ( (i > 1) && ((chi[p][j][i] > chi[p][j][i-1]) && (chi[p][j][i-1] < chi[p][j][i-2])) ) {
					crossed_min[p][j] += 1;
					//printf("i-1 = %d, imax[p][%d] = %d, crossed_min[p][%d] = %d, crossed_max[%d][%d] = %d\n", i-1, j,imax[p][j], j, crossed_min[p][j], p, j, crossed_max[p][j]);
					if (crossed_min[p][j] > 1) {	
						printf("***WARNING*** There is more than one minimum!\n");
						printf("j is %d and minima at %d and %d for first and second\n",  j, imin[p][j], i-1);
						for (k=imin[p][j]-1;k<i+1;k++) 
						printf("chi[%d][%d] = %f\n", j, k, chi[p][j][k]);
						crossed_min[p][j] = 1; 
					}
					imin[p][j] = i-1; // Store the index of the maximum 
					upperlimit[p][j] = imin[p][j] - imax[p][j];
					chimpp[p][j] = ( ( chi[p][j][i] - chi[p][j][i-1] ) / (xx[i] - xx[i-1]) - (chi[p][j][i-1] - chi[p][j][i-2]) / (xx[i-1] - xx[i-2]) ) *2.0/ (xx[i] - xx[i-2] ) ;  
					//printf("chimpp[%d] = %f\n", chimpp[j]);
					chimin[p][j] = chi[p][j][i-1];  //Store chi minimum itself, a function of xbar (index j)
					//crossed_max[p][j] = 1;
				}
				//printf("crossed_min[%d][%d] = %d\n", p, j, crossed_min[p][j]);
				//printf("i=%d/%d\tcrossed_max[%d][%d] = %d\n",i, size_finegrid, p, j, crossed_max[p][j]);

/* We will now start going up the effective potential well! The temporary flag below is used because we still have to store the effective potential minimum as a possible value of Uperp. If we don't have this flag we miss the region near the minimum of chi. */
/* FILLING IN ARRAYS */
/* Once we cross a maximum, we start filling in arrays for mu, Uperp and vx related to the orbit we are in. As we go down the maximum we store the values of Uperp = chi(x) we encounter which will form a grid of (unevenly spaced) allowed values of Uperp. We also store the value of the small deltamu associated with every value of Uperp above the current value of chi(x), and add it to the previous values */                                                                               	
				if ( ( (crossed_max[p][j] == 1  && crossed_min[p][j] == 0) || (i-1 == imin[p][j]) )  ) {
					Uperp[p][j][i-1-imax[p][j]] = chi[p][j][i-1];
					if (Uperp[p][j][i-1-imax[p][j]] > Ucap && lowerlimit[p][j] != 0)
						lowerlimit[p][j] = i-1-imax[p][j]; 
					mu[p][j][i-1-imax[p][j]] = 0.0;
					upper[p][j][i-1] = i-1-imax[p][j];
		/* Note that the size of the dimension of the array with values of Uperp is set to n, which is larger than the size it will turn out to be. This is because in C there is no way to append elements to arrays as I go along, enough memory has to be given to the array from the start. n is the largest possible size the array could have. */
					for (k=0;k<=upper[p][j][i-1]; k++) 
					{	
						//printf("k=%d/%d\n", k, upper[j][i-1]);
						// replaced k with upper below
						if ( (upper[p][j][i-1] == 0) || (i-1 == 0) ) {
							vx[p][j][i-1][k] = sqrt((2.0*Uperp[p][j][k] - chi[p][j][i-1])); // equiv 0
							mu[p][j][k] += 0.0; 

						}
						else if ( (k == imin[p][j] - imax[p][j]) && (crossed_min[p][j] == 1) ) {
							vx[p][j][i-1][k] = 0.0;
							mu[p][j][k] = 0.0; 
						}
						else if ( (k == imin[p][j] - imax[p][j] - 1) && (crossed_min[p][j] == 1) ) {
							vx[p][j][i-1][k] = sqrt(2.0*(Uperp[p][j][k] - chi[p][j][i-1])); // equiv 0
							mu[p][j][k] = 0.5*pow(xx[imin[p][j]-1] - xx[imin[p][j]], 2.0)*pow(chimpp[p][j], 0.5);
							//printf("1..%f\n", sqrt(0.5*chimpp[j])*pow(xx[imin[j]-1] - xx[imin[j]], 2.0));
						}
						else if ( k == upper[p][j][i-1] - 1 )  {	
							vx[p][j][i-1][k] = sqrt(2.0*(Uperp[p][j][k] - chi[p][j][i-1])); // equiv 0
							mu[p][j][k] += (sqrt(2.0)/M_PI)*sqrt(chi[p][j][i-2]-chi[p][j][i-1])*(2.0/3.0)*(xx[i-1] - xx[i-2]);
						}
						else if (k == upper[p][j][i-1]) {
							vx[p][j][i-1][k] = 0.0;
						}
						else {
							vx[p][j][i-1][k] = sqrt(2.0*(Uperp[p][j][k] - chi[p][j][i-1]));
							mu[p][j][k] += (1.0/M_PI)*0.5*(vx[p][j][i-1][k] + vx[p][j][i-2][k])*(xx[i-1] - xx[i-2]); 
						}
						if (mu[p][j][k] != mu[p][j][k]) {
							printf("BEFORE: mu[%d][%d] is NAN\n", j, k); 
							exit(-1);
						}  
					}
				}
		/* Once we cross the minimum, we stop creating array elements with values of Uperp. However, we keep storing the value of vx associated with any given point x on an effective potential curve with xbar, with energy Uperp and using this value to finish performing the mu integral. This should happen as long the effective potential at the point under consideration is smaller than the effective potential maximum. */
				else if ( ( crossed_min[p][j] == 1 && crossed_max[p][j] == 1 && chi[p][j][i-1] < chiMax[p][j] && ( i-1 != imin[p][j] ) ) ) {
					for (k=0;k <= upperlimit[p][j] ;k++) {	
						if ( (chi[p][j][i-1] < Uperp[p][j][k]) && (chi[p][j][i-2] < Uperp[p][j][k]) ) {	
							vx[p][j][i-1][k] = sqrt(2.0*(Uperp[p][j][k] - chi[p][j][i-1]));
							mu[p][j][k] += (1.0/M_PI)*0.5*(vx[p][j][i-1][k] + vx[p][j][i-2][k])*(xx[i-1] - xx[i-2]); 
						}
						else if (Uperp[p][j][k] <= chi[p][j][i-1] && Uperp[p][j][k-1] > chi[p][j][i-1]) {
							upper[p][j][i-1] = k;
					//mu[j][k] += (2.0/M_PI)*(vx[j][i-2][k])*(xx[i-1] - xx[i-2])*(Uperp[j][k] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]);
							ind = 0;
							while (Uperp[p][j][k] < chi[p][j][i-2-ind]) 
								ind++;
							
							mu[p][j][k] += (sqrt(2.0)/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[p][j][k] - chi[p][j][i-2-ind], 1.5)/(chi[p][j][i-1-ind] - chi[p][j][i-2-ind]); // double check normalization
				
						}
						else if (Uperp[p][j][k] <= chi[p][j][i-1] && Uperp[p][j][k] > chi[p][j][i-2]) 
							mu[p][j][k] += (sqrt(2.0)/M_PI)*(2.0/3.0)*(xx[i-1] - xx[i-2])*pow(Uperp[p][j][k] - chi[p][j][i-2], 1.5)/(chi[p][j][i-1] - chi[p][j][i-2]); 
						if (mu[p][j][k] != mu[p][j][k]) {
							printf("mu[%d][%d] is NAN\n", j, k); 
							exit(-1);
						}  
					}
				}
		/* When the effective potential at the iteration point (i-1) under consideration becomes larger than the effective potential maximum, we finish performing the mu integral. We also store the position of the top of the orbit which has chi = chiMax, in order to perform the open orbit integral. If the loop below is accessed, a switch it turned off to signify the no more closed orbits can be present */
				//else if ( ( crossed_min[j] == 1 ) && ( crossed_max[j] == 1) && ( chi[j][i-1] > chiMax[j] - TINY) ) 	
				else if ( ( crossed_min[p][j] == 1 ) && ( crossed_max[p][j] == 1) && ( chi[p][j][i-1] > chiMax[p][j] ) ) {
					//itop[j] = i-2;
					xtop[p][j] = xx[i-2] + ((chiMax[p][j] - chi[p][j][i-2])/(chi[p][j][i-1] - chi[p][j][i-2]))*(xx[i-1] - xx[i-2]);
					for (k=0; k<upper[p][j][i-2]; k++) {
						ind = 0;
						//while (Uperp[j][0] < chi[j][i-2-ind]) 
						//	ind++; // find top bounce point index i-2-ind for a given xbar[j] and for Uperp[j][0] = chiMax
						//mu[j][k] += (sqrt(2.0)/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[j][0] - chi[j][i-2-ind], 1.5)/(chi[j][i-1-ind] - chi[j][i-2-ind]);  // CHECK NORMALIZATION
						while (Uperp[p][j][k] < chi[p][j][i-2-ind]) 
							ind++; // find top bounce point index i-2-ind, with i the smallest index such that chi[j][i-1] > chiMax[j] for a given xbar[j], for orbits with Uperp[j][k] > chi[j][i-2] 
						mu[p][j][k] += (sqrt(2.0)/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[p][j][k] - chi[p][j][i-2-ind], 1.5)/(chi[p][j][i-1-ind] - chi[p][j][i-2-ind]);  // CHECK NORMALIZATION
						//printf("mu[%d][%d] = %f\n", j, k, mu[j][k]);
						//mu[j][k] += (1.0/M_PI)*(vx[j][i-2][k])*(xx[i-1] - xx[i-2])*(Uperp[j][k] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]);
					}
					crossed_max[p][j] = 0; 
				}
				if (j!=0) {
					if ( ((chiMax[p][j-1] < chi[p][j-1][i-1] ) && (chiMax[p][j] > chi[p][j][i-1])) ) { // || (imax[p][j-1] == -1 && imax[p][j] != -1) )
						jmclosed[p][i-1] = j-1;
						if (i-1>icrit[p]) { 	
							jmopen[p][i-1] = j-1; 
							if (DEBUG == 1) {
								printf("i = %d, j = %d, icrit = %d\n", i-1, j-1, icrit[p]); 
								printf("jmopen[%d] = %d\n", i-1, jmopen[p][i-1]); 
							}
						} 
					}
				} 
				//printf("upper = %d, upperlimit = %d\n", upper[j][i-1], upperlimit[j]);
			}
		} 

		upperlimit[p][sizexbar-1] = upperlimit[p][sizexbar-2] +1;
	}

	for (p=0; p<size_ygrid; p++) {
		// OPEN ORBIT INTEGRAL
		/* Now we perform the open orbit integral. We use a change of variables which makes the integrand smooth at the top bounce point. The change of variables is to some var = sqrt(x_t - x) */
		//maxj[p] = sizexbar +1 ; //*
		printf("p=%d/%d\n", p, size_ygrid);
		muopen[p] = malloc(maxj[p]*sizeof(double)); //*
		chiMopen[p] = malloc(maxj[p]*sizeof(double));//*	
		xbaropen[p] = malloc(maxj[p]*sizeof(double));//*
		twopidmudxbaropen[p] = malloc(maxj[p]*sizeof(double));
		twopidmudyopen[p] = malloc(maxj[p]*sizeof(double));

		xbaropen[p][0] = xbarcrit[p]; 
		chiMopen[p][0] = chiMcrit[p]; 
		muopen[p][0] = 0.0; 
		twopidmudxbaropen[p][0] = 0.0; 
		twopidmudyopen[p][0] = 0.0; 
		//printf("sizexbar = %d\n", sizexbar);
		//printf("p=%d/%d\n", p, size_ygrid);
		for (j=0;j<maxj[p]-1;j++) {
			//printf("lowerlimit[%d] = %d and upper[%d][0] = %d, upperlimit = %d\n", j, lowerlimit[j],j, upper[j][0], upperlimit[j]);
			mu[p][j][upperlimit[p][j]] = 0.0;
			chiMopen[p][j+1] = chiMax[p][j];
			muopen[p][j+1] = mu[p][j][0]; 
			xbaropen[p][j+1] = xbar[p][j];
			if (DEBUG == 1)
				printf("jmopen[%d] = %d\n", j, jmopen[p][j]); 
			if (j==0) 
				twopidmudxbar[p][j] = (2.0*M_PI) * ( ((xbar[p][j] - xbarcrit[p])/(xbar[p][j+1] - xbarcrit[p])) *(mu[p][j+1][0] - mu[p][j][0])/(xbar[p][j+1] - xbar[p][j]) + ((xbar[p][j+1] - xbar[p][j])/(xbar[p][j+1] - xbarcrit[p])) * (mu[p][j][0] - 0.0)/(xbar[p][j] - xbarcrit[p]) );	
			else if (j==maxj[p]-2) 
				twopidmudxbar[p][j] = twopidmudxbar[p][j-1];	
			else 
				twopidmudxbar[p][j] = (2.0*M_PI) * (((xbar[p][j] - xbar[p][j-1])/(xbar[p][j+1] - xbar[p][j-1])) *(mu[p][j+1][0] - mu[p][j][0])/(xbar[p][j+1] - xbar[p][j]) + ((xbar[p][j+1] - xbar[p][j])/(xbar[p][j+1] - xbar[p][j-1])) * (mu[p][j][0] - mu[p][j-1][0])/(xbar[p][j] - xbar[p][j-1]) );	
			if ( (p!=0) && (p!=size_ygrid-1) ) 
				twopidmudy[p][j] = (2.0*M_PI) * (mu[p+1][j][0] - mu[p-1][j][0])/(y_grid[p+1] - y_grid[p-1]);
			twopidmudxbaropen[p][j+1] = twopidmudxbar[p][j];
			twopidmudxbarantycal = 4.0*M_PI*xbar[p][j];
			if (DEBUG == 1)
				printf("%f %f %f %f %f\n", xbar[p][j], mu[p][j][0], Uperp[p][j][0], twopidmudxbar[p][j], twopidmudxbarantycal); 
		}
		muopen[p][maxj[p]-1] = 5000.0;
		chiMopen[p][maxj[p]-1] = 5000.0; 
		xbaropen[p][maxj[p]-1] = 100.0;
		twopidmudxbaropen[p][maxj[p]-1] = 100.0;
		twopidmudyopen[p][maxj[p]-1] = 0.0;
	}
	for (j=0;j<maxj[size_ygrid-1]-1;j++) 
		twopidmudy[size_ygrid-1][j] = 0.0; //(2.0*M_PI) * (mu[0][j][0] - mu[size_ygrid-2][j][0])/deltay;
	for (j=0;j<maxj[0]-1;j++) 
		twopidmudy[0][j] = 0.0 ; //(2.0*M_PI) * (mu[size_ygrid-1][j][0] - mu[1][j][0])/deltay;

	for (p=0; p<size_ygrid; p++) {
		if (DEBUG == 1) {
			printf("~~~~~The second element of FF is %f~~~~~\n",FF[p][0][1]);
			printf("~~~~~The second element of UU is %f~~~~~\n",UU[1]);
			printf("~~~~~The fourth element of mu is %f~~~~~\n",mumu[3]);
		}
		i=0;
		clock_t int1 = clock(); // Finds the time of the computation so far
		double inttime  = (double)(int1 - begin) / CLOCKS_PER_SEC;
		if (DEBUG == 1) 
			printf("in iondens2d: Array filling DONE: time is %f\n", inttime);

		/* DENSITY INTEGRALS 
		This part calculates the density integrals and outputs the result of the integration to a file fout and also the yz distribution function to three files one containing the distribution function the other two containing the velocity grid */
		////////////////////////////////////////////
		// calculate normalization at infinity
		deltax_inf = (xx[size_finegrid-1] - xx[size_finegrid-2]);
		j_inf = (int) sqrt(Ucap)/deltax_inf ;
		//printf("j_inf = %d, deltax_inf = %f\n", j_inf, deltax_inf);
		muinf[p] = malloc(j_inf*sizeof(double)); 
		vxinf[p] = malloc(j_inf*sizeof(double)); 
		chiinf[p] = malloc(j_inf*sizeof(double));
		for (j=j_inf-1; j>=0; j--) {	
			vxinf[p][j] = (double*)calloc(j_inf,sizeof(double*)); 
			chiinf[p][j] =  0.5*deltax_inf*j*deltax_inf*j;
			for (k=j; k<j_inf; k++) {
				vxinf[p][j][k] = sqrt(2.0*(chiinf[p][k] - chiinf[p][j]));
			}
			muinf[p][j] = chiinf[p][j];
		}

		intdxbar = 0.0;
		intdxbarflow = 0.0;
		intdvx = 0.0;
		intdvxflow = 0.0;
		for (j=0; j<j_inf; j++) {
			vxnew = 0.0;
			intdvxold = intdvx;
			intdvxflowold = intdvxflow;
			intdvx = 0.0;
			intdvxflow = 0.0;
			intdU = 0.0;
			intdUflow = 0.0;
			for (k=j; k<j_inf; k++) {
				vxold = vxnew;
				intdUold = intdU;
				intdUflowold = intdUflow;
				intdU = 0.0;
				intdUflow = 0.0;
				munew = muinf[p][k]; 
				vxnew = vxinf[p][j][k];
				Uperpnew = chiinf[p][k];
				//if ((j==3) && (k==4)) printf("munew = %f\tvxnew = %f\tUperpnew = %f\n", munew, vxnew, Uperpnew);
				sizeU = (int) 2.0*sqrt(2.0*(Ucap - Uperpnew))/dvz;
				Fold = F = 0.0;
				for (l=0; l < sizeU; l++) {	
					Fold = F;
					vzmin = -sqrt(2.0*(Ucap - Uperpnew)); 
					vzmin = 0.0;
					vz = vzmin + dvz*l;
					U = Uperpnew + 0.5*pow(vz, 2.0);
					//ystar = y_grid[p] - vz/alpha; 
					//if (U<=munew) vzinf = 0.0;
					//else vzinf = sqrt(2.0*(U - munew)); 
					//phiinfval = lin_interp(y_grid, phiinf, ystar + vzinf/alpha, size_ygrid, 1);
					//dphiinfdyval = lin_interp(y_grid, dphiinfdy, ystar + vzinf/alpha, size_ygrid, 1);
					//if (U<=munew + phiinfval) vzinf = 0.0;
					//else vzinf = sqrt(2.0*(U - munew - phiinfval)); 
					//vzinf = vzinf + (1.0/alpha)*dphiinfdyval;
					//do {
					//	vzinfold = vzinf;
					//}  while (fabs(vzinf - vzinfold) > TINY);
					vzinf = vz;// + (1.0/alpha)*dphiinfdyval;
					//if (vzinf <= 0.0) vzinf = 0.0;
					//yinf = ystar + vzinf/alpha;
					yinf = y_grid[p];
					//if (yinf >= y_grid[size_ygrid-1]) yinf = y_grid[size_ygrid-1] - TINY;
					//if (yinf <= y_grid[0]) yinf = y_grid[0] + TINY;
					//printf("vz = %f\t vzinf=%f\n", vz, vzinf);
					//if ( (U > munew) && (U - 0.5*vz*vz + 0.5*(vz-dvz)*(vz-dvz) < munew) ) {
					//	frac = (vz - sqrt(2.0*(munew - Uperpnew)))/dvz;
					//	Fold = trilin_interp(ystar, munew, 0.0, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
					//	F = trilin_interp(ystar, munew, U, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
					//}
					//else if (U > munew){
					//frac = 1.0;
					//F = trilin_interp(ystar, munew, U, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
					//}
					//else {
					//	frac = 1.0;
					//	F = 0.0;
					//}
					//intdU += frac*dvz*(F+Fold);
					F = trilin_interp(yinf, munew, vzinf, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
					//if (k==1 && j == 0) printf("F = %f\tFF=%f\t\tyinf=%f\tmunew=%f\tvzinf=%f\ty_grid = %f\tmumu = %f\tUU=%f\n", F, FF[p][1][1], yinf, munew, vzinf, y_grid[p], mumu[1], UU[1]);
					//if ((j==3) && (k==4))
					//	printf("F[%d] = %f\n", l, F, vzinf, munew);
					if (l!=0) {
						intdU += dvz*(F+Fold);
						intdUflow += alpha*(0.5*vz*vz - 0.5*(vz - dvz)*(vz-dvz))*(F+Fold);
					}
				}
				intdUantycal = exp(-chiinf[p][k])*(1.0/(2.0*M_PI));// result with phi =0
				if (DEBUG == 1) printf("Analytical intdU is %f, numerical one is %f\n", intdUantycal, intdU);
				if (k!=j) {	
					dvx = vxnew - vxold; 
					intdvx += 4.0*0.5*dvx*(intdU+intdUold); // CHECK NORMALIZATION
					intdvxflow += 4.0*0.5*dvx*(intdUflow+intdUflowold); // CHECK NORMALIZATION
					//printf("k = %d j = %d, intdvx = %f\n", k, j, intdvx);
				}
				intdvxantycal = (2.0/(2.0*sqrt(M_PI)))*exp(-deltax_inf*j*deltax_inf*j);//*erf(sqrt(xbar[j]*xbar[j]-(pos-xbar[j])*(pos-xbar[j])));
			}
			if (DEBUG == 1)	printf("intdvx is %f, analytical one is %f\n", intdvx, intdvxantycal); 
			//printf("intdvx is %f, while the analytical one is %f\n", intdvx, intdvxantycal);
			//intdvx = intdvxantycal;
			
			dxbar = deltax_inf;
			if (j!=0) {
				intdxbar += (intdvx+intdvxold)*dxbar; // multiply by two because only half the xbars are contemplated here
				intdxbarflow += (intdvxflow+intdvxflowold)*dxbar; // multiply by two because only half the xbars are contemplated here
			}
		}
		n_inf[p] = intdxbar;
		fluxinf2[p] = intdxbarflow;
		if (p==0) n_norm = n_inf[0];	
		if (DEBUG == 0) {	
			printf("n_inf = %f\n", n_inf[p]);
			if ( (n_inf[p]  != n_inf[p]) || (n_inf[p] < TINY) ) { 	
				printf("n_inf = %f\n", n_inf[p]);
				exit(-1);
			}
		}

		// density profile
		i=0; // start from x[0] = 0 i.e. the wall
		stop = 0; // set stop index to zero; it turns to 1 if density exceeds threshold in the input (expressed as fraction of density at infinity)
		ic = 0;
		while (stop == 0) {
		//for (ic=0; ic<size_ngrid; ic++) 
			i = ic*zoomfactor;
			intdxbar = 0.0;
			intdxbaropen = 0.0;
			intdxbaropenflow = 0.0;
			for (j=0; j<sizexbar; j++) { //*
				//printf("ic=%d/%d\tj=%d/%d\n", ic, size_phigrid, j, sizexbar);
				vxnew = 0.0;
				intdUopenold = intdUopen;
				intdUopenflowold = intdUopenflow; intdUopenBohmold = intdUopenBohm;
				intdUopenxbarold = intdUopenxbar; intdUopensquareold = intdUopensquare;
				intdUopen = 0.0;
				intdUopenflow = 0.0; intdUopenBohm = 0.0;
				intdUopenxbar = 0.0; intdUopensquare = 0.0;
				if (j == jmopen[p][i]) {
					oorbintgrd = oorbintgrdold = 0.0;
					oorbintgrdflow = oorbintgrdflowold = 0.0; oorbintgrdBohm = oorbintgrdBohmold = 0.0;
					oorbintgrdxbar = oorbintgrdxbarold = 0.0; oorbintgrdsquare = oorbintgrdsquareold = 0.0;
					sizeU = (int) sqrt(2.0*(Ucap - chiMax[p][j]))/dvzopen;
					if (j==0) { //????????????
						munew = 0.0;
						openorbitnew = 0.0;
						yopenorbitnew = 0.0;
					}
					else {
						munew = mu[p][j+1][0] + ( (chiMax[p][j+1] - chi[p][j+1][i]) / (chiMax[p][j+1] - chi[p][j+1][i] + chi[p][j][i] - chiMax[p][j]) ) * (mu[p][j][0] - mu[p][j+1][0]);
						openorbitnew = twopidmudxbar[p][j+1] + ( (chiMax[p][j+1] - chi[p][j+1][i]) / (chiMax[p][j+1] - chi[p][j+1][i] + chi[p][j][i] - chiMax[p][j]) ) * (twopidmudxbar[p][j] - twopidmudxbar[p][j+1]);
						yopenorbitnew = twopidmudy[p][j+1] + ( (chiMax[p][j+1] - chi[p][j+1][i]) / (chiMax[p][j+1] - chi[p][j+1][i] + chi[p][j][i] - chiMax[p][j]) ) * (twopidmudy[p][j] - twopidmudy[p][j+1]);
					}
					for (l=0; l < sizeU; l++) {
						oorbintgrdold = oorbintgrd;
						oorbintgrdflowold = oorbintgrdflow; oorbintgrdBohmold = oorbintgrdBohm;
						oorbintgrdxbarold = oorbintgrdxbar; oorbintgrdsquareold = oorbintgrdsquare;
						vzmin = 0.001;
						square = (xbar[p][j] - xx[imax[p][j]])*yopenorbitnew/(alpha*openorbitnew);
						if (vzmin < square) {
							printf("WARNING: vzmin = %f\tsquare = %f\n", vzmin, square);
							vzmin = square;
						}
						vz = vzmin + dvzopen*l;
						if (j!=0) { ///???????
							chinew = chi[p][j+1][i] + ( (chiMax[p][j+1] - chi[p][j+1][i]) / (chiMax[p][j+1] - chi[p][j+1][i] + chi[p][j][i] - chiMax[p][j]) ) * (chi[p][j][i] - chi[p][j+1][i]) ;
							vx0open = 0.0;
							vx0open = sqrt(TINY);
						}
						else {
							chinew = chi[p][j][i];
							if (chiMax[p][j] > chinew) 
								vx0open = sqrt(2.0*(chiMax[p][j] - chinew) + TINY);
							else vx0open = 0.0;
						}
						U = chinew + 0.5*pow(vz, 2.0); 
						ystar = y_grid[p] - vz/alpha; 
						if (U<=munew) vzinf = 0.0;
						else vzinf = sqrt(2.0*(U - munew)); 
						phiinfval = lin_interp(y_grid, phiinf, ystar + vzinf/alpha, size_ygrid, 1);
						dphiinfdyval = lin_interp(y_grid, dphiinfdy, ystar + vzinf/alpha, size_ygrid, 1);
						if (U<=munew + phiinfval) vzinf = 0.0;
						else vzinf = sqrt(2.0*(U - munew - phiinfval)); 
						vzinf = vzinf + (1.0/alpha)*dphiinfdyval;
						yinf = ystar + vzinf/alpha;
						if (yinf > y_grid[size_ygrid-1]) yinf = y_grid[size_ygrid-1];
						if (yinf < y_grid[0]) yinf = y_grid[0];
						//printf("yopenorbitnew = %f\n", yopenorbitnew);
						if ( (U > munew) && (U - 0.5*vz*vz + 0.5*(vz-dvzopen)*(vz-dvzopen) < munew) ) {
							//frac = (vz - sqrt(2.0*(munew - chinew)))/dvzopen;
							Fopen = trilin_interp(yinf, munew, 0.0, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							oorbintgrdold = ( sqrt(vx0open*vx0open + 2.0*alpha*sqrt(2.0*(munew - chinew))*openorbitnew - 2.0*(xbar[p][j] - xx[imax[p][j]])*yopenorbitnew) - vx0open )*Fopen;
							oorbintgrdflowold = (alpha*sqrt(2.0*(munew - chinew))*openorbitnew - (xbar[p][j] - xx[imax[p][j]])*yopenorbitnew)*Fopen;
							Fopen = trilin_interp(yinf, munew, U-munew, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							oorbintgrd = ( sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbitnew - 2.0*(xbar[p][j] - xx[imax[p][j]])*yopenorbitnew) - vx0open )*Fopen;
							oorbintgrdflow = (alpha*vz*openorbitnew - (xbar[p][j] - xx[imax[p][j]])*yopenorbitnew)*Fopen;
							if (fabs(Fopen) < TINY) oorbintgrd = 0.0;
						}
						else {
							//frac = 1.0;
							Fopen = trilin_interp(yinf, munew, vzinf, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							oorbintgrd = ( sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbitnew - 2.0*(xbar[p][j] - xx[imax[p][j]])*yopenorbitnew) - vx0open )*Fopen;
							oorbintgrdflow = (alpha*vz*openorbitnew - (xbar[p][j] - xx[imax[p][j]])*yopenorbitnew)*Fopen;
							if (fabs(Fopen) < TINY) oorbintgrd = 0.0;
						}
						if (oorbintgrd != oorbintgrd) {	
							printf("vx0open = %f, twopidmudxbar[%d] = %f, imaginary oorbintgrd in initial piece of integral due to negative value of twopidmudxbar?\n", vx0open, j, twopidmudxbar[p][j]); 
							printf("Fopen = %f\n", Fopen); 
							exit(-1);
						}
						if (i==0) 
							oorbintgrdantycal = sqrt(2.0*alpha*vz*M_PI*xbar[p][j])*(U-chiMax[p][j])*exp(-U)/pow(M_PI, 1.5);  // CHECK NORMALIZATION
							// for a flat potential, oorbintgrd can be calculated analytically
							//oorbintgrd = oorbintgrdantycal;
							//printf("oorbintgrd is %f, analytical is %F\n", oorbintgrd, oorbintgrdantycal);
						if (l!=0) {	
							intdUopen += 2.0*dvzopen*(oorbintgrd+oorbintgrdold);
							intdUopenflow += 2.0*dvzopen*(oorbintgrdflow+oorbintgrdflowold);
							intdUopenBohm += 2.0*dvzopen*(oorbintgrdBohm+oorbintgrdBohmold);
							intdUopenxbar += 2.0*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
							intdUopensquare += 2.0*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
						}
						else {	
							intdUopen += 0.0;
							intdUopenflow += 0.0;
							intdUopenBohm += 0.0;
							intdUopenxbar += 2.0*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
							intdUopensquare += 2.0*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
						} 
					}
					if (j == 0) { //*
						intdUopenold = 0.0;
						//printf("intdUopenold = %f\n", intdUopenold);
						dxbar = xbar[0] - xbarcrit;
						dxbar = 0.0;
						intdxbaropen += 0.5*(intdUopen+intdUopenold)*dxbar;
						intdxbaropenflow += 0.5*(intdUopenflow+intdUopenflowold)*dxbar;
						intdxbaropenBohm += 0.5*(intdUopenBohm+intdUopenBohmold)*dxbar;
						intdxbarxbar += 0.5*(intdUopenxbar+intdUopenxbarold)*dxbar;
						intdxbarsquare += 0.5*(intdUopensquare+intdUopensquareold)*dxbar; 
					}
				}
				else if (j > jmopen[p][i]) {
					oorbintgrd = oorbintgrdold = 0.0;
					oorbintgrdflow = oorbintgrdflowold = 0.0;
					oorbintgrdBohm = oorbintgrdBohmold = 0.0;
					oorbintgrdxbar = oorbintgrdxbarold = 0.0;
					oorbintgrdsquare = oorbintgrdsquareold = 0.0;
					sizeU = (int) sqrt(2.0*(Ucap - chiMax[p][j]))/dvzopen;
					for (l=0; l < sizeU; l++) {
						oorbintgrdold = oorbintgrd;
						oorbintgrdflowold = oorbintgrdflow;
						oorbintgrdBohmold = oorbintgrdBohm;
						oorbintgrdxbarold = oorbintgrdxbar;
						oorbintgrdsquareold = oorbintgrdsquare;
						vzmin = - sqrt(2.0*(Ucap - chiMax[p][j]));
						vzmin = 0.001;
						square = (xbar[p][j] - xx[imax[p][j]])*twopidmudy[p][j]/(alpha*twopidmudxbar[p][j]);
						if (vzmin < square) {
							printf("WARNING: vzmin = %f\tsquare = %f\n", vzmin, square);
							vzmin = square;
						}
						vz = vzmin + dvzopen*l;
						U = chiMax[p][j] + 0.5*pow(vz, 2.0);
						ystar = y_grid[p] - vz/alpha; 
						if (U<=munew) vzinf = 0.0;
						else vzinf = sqrt(2.0*(U - munew)); 
						phiinfval = lin_interp(y_grid, phiinf, ystar + vzinf/alpha, size_ygrid, 1);
						dphiinfdyval = lin_interp(y_grid, dphiinfdy, ystar + vzinf/alpha, size_ygrid, 1);
						if (U<=munew + phiinfval) vzinf = 0.0;
						else vzinf = sqrt(2.0*(U - munew - phiinfval)); 
						//if ( ((p==1) && (i==0)) && ( (j == 16) && (k==1) )) 
						//	printf("vzinf = %f, U= %f, munew = %f, phiinfval = %f\n", vzinf, U, munew, phiinfval);
						vzinf = vzinf + (1.0/alpha)*dphiinfdyval;
						//if ( ((p==1) && (i==0)) && ( (j == 16) && (k==1) )) 
						//	printf("vzinf = %f, alpha= %f, dphiinfdyval = %f\n", vzinf, alpha, dphiinfdyval);
						//vzinf = sqrt(2.0*(U - munew)); 
						//}  while (fabs(vzinf - vzinfold) > TINY);
						yinf = ystar + vzinf/alpha;
						if (yinf > y_grid[size_ygrid-1]) yinf = y_grid[size_ygrid-1];
						if (yinf < y_grid[0]) yinf = y_grid[0];
						//vx0open = sqrt(TINY + chiMax[j] - chi[j][i]);
						vx0open = sqrt(2.0*(chiMax[p][j] - chi[p][j][i]) + TINY); 
						if (vx0open != vx0open) {
							printf("HERE imaginary vx0open, j = %d, i is %d, chi[j][i] = %f, chiMax[j] = %f\n", j, i, chi[p][j][i], chiMax[p][j]); 
							exit(-1);
						}
						if ( (U >= mu[p][j][0]) && (U - 0.5*vz*vz + 0.5*(vz-dvzopen)*(vz-dvzopen) < mu[p][j][0]) ) {
							Fopen = trilin_interp(yinf, mu[p][j][0], 0.0, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							oorbintgrdold = sqrt(2.0*alpha*sqrt(2.0*(mu[p][j][0] - chiMax[p][j]))*twopidmudxbar[p][j] - 2.0*(xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j])*Fopen;
							oorbintgrdflowold = (alpha*sqrt(2.0*(mu[p][j][0] - chiMax[p][j]))*twopidmudxbar[p][j] - (xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j])*Fopen;
							oorbintgrdBohmold = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + 2.0*alpha*sqrt(2.0*(mu[p][j][0] - chiMax[p][j]))*twopidmudxbar[p][j] - 2.0*(xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j]))*Fopen;
							oorbintgrdxbarold = xbar[p][j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + 2.0*alpha*sqrt(2.0*(mu[p][j][0] - chiMax[p][j]))*twopidmudxbar[p][j] - 2.0*(xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j]))*Fopen;
							oorbintgrdsquareold = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + 2.0*alpha*sqrt(2.0*(mu[p][j][0] - chiMax[p][j]))*twopidmudxbar[p][j] - 2.0*(xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j]), 1.5))*Fopen;
							Fopen = trilin_interp(yinf, mu[p][j][0], vzinf, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							oorbintgrd = sqrt(2.0*alpha*vz*twopidmudxbar[p][j] - 2.0*(xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j])*Fopen;
							oorbintgrdflow = (alpha*vz*twopidmudxbar[p][j] - (xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j])*Fopen;
							oorbintgrdBohm = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + 2.0*alpha*vz*twopidmudxbar[p][j] - 2.0*(xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j]))*Fopen;
							oorbintgrdxbar = xbar[p][j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + 2.0*alpha*vz*twopidmudxbar[p][j] - 2.0*(xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j]))*Fopen;
							oorbintgrdsquare = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + 2.0*alpha*vz*twopidmudxbar[p][j] - 2.0*(xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j]), 1.5))*Fopen;
							//if (fabs(Fopen) < TINY) oorbintgrd = 0.0;

						}
						else {
							Fopen = trilin_interp(yinf, mu[p][j][0], vzinf, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							oorbintgrd = (sqrt(vx0open*vx0open + 2.0*alpha*vz*twopidmudxbar[p][j] - 2.0*(xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j]) - vx0open)*Fopen;
							oorbintgrdflow = alpha*vz*twopidmudxbar[p][j]*Fopen;
							oorbintgrdBohm = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + 2.0*alpha*vz*twopidmudxbar[p][j] - 2.0*(xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j]))*Fopen;
							oorbintgrdxbar = xbar[p][j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + 2.0*alpha*vz*twopidmudxbar[p][j] - 2.0*(xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j]))*Fopen;
							oorbintgrdsquare = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + 2.0*alpha*vz*twopidmudxbar[p][j] - 2.0*(xbar[p][j]-xx[imax[p][j]])*twopidmudy[p][j]), 1.5))*Fopen;
							//if (fabs(Fopen) < TINY) oorbintgrd = 0.0;
						}
						if (oorbintgrd != oorbintgrd) {
							printf("vx0open = %f, twopidmudxbar[%d] = %f, twopidmudy = %.15f, vz = %.15f, vy = %f, vzinf = %f, alpha = %f, dphidyinfval = %f, phiinfval = %f, U = %f, munew = %f\nHERE imaginary oorbintgrd in initial piece of integral due to negative value of twopidmudxbar?\n", vx0open, j, twopidmudxbar[p][j], twopidmudy[p][j], vz, xbar[p][j]-xx[imax[p][j]], vzinf, alpha, dphiinfdyval, phiinfval, U, munew); 
							printf("p=%d, i=%d, j=%d, k=%d\n", p, i, j, k);
							printf("Fopen = %f, dphiinfdyval = %f\n", Fopen, dphiinfdyval); 
							//printf("vzmin=%f\tsquare = %f\n", vzmin, square);
							exit(-1);
						}
						if (i==0) {
							oorbintgrdantycal = sqrt(2.0*alpha*vz*M_PI*xbar[p][j])*(U-chiMax[p][j])*exp(-U)/pow(M_PI, 1.5); 
							if (DEBUG == 1) 
								printf("twopidmudxbar = %f (should be %f with flat potential profile)\n", oorbintgrd, oorbintgrdantycal);
							// check not working
						}
		// for a flat potential, oorbintgrd can be calculated analytically
						//oorbintgrd = oorbintgrdantycal;
		//printf("oorbintgrd is %f, analytical is %F\n", oorbintgrd, oorbintgrdantycal);
						if (l!=0) {
							intdUopen       += 2.0*dvzopen*(oorbintgrd+oorbintgrdold);
							intdUopenflow   += 2.0*dvzopen*(oorbintgrdflow+oorbintgrdflowold);
							intdUopenBohm   += 2.0*dvzopen*(oorbintgrdBohm+oorbintgrdBohmold);
							intdUopenxbar   += 2.0*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
							intdUopensquare += 2.0*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
						}
						else {
							intdUopen       += 0.0;
							intdUopenflow   += 0.0;
							intdUopenBohm   += 0.0;
							intdUopenxbar   += 2.0*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
							intdUopensquare += 2.0*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
						} 
					}
					if (intdUopen != intdUopen)  {
						printf("HERE, j is %d\n", j); 
						exit(-1);
					}
					if (i==0) {	
						intdUopenantycal = 2.0*0.919*sqrt(2.0*alpha*xbar[p][j])*exp(-xbar[p][j]*xbar[p][j])/M_PI; 
						if (DEBUG == 1) { 
							printf("xbar[p][%d] = %f, analytical is %f, numerical is %f\n", j, xbar[p][j], intdUopenantycal, intdUopen); printf("vx0open is %f (should be zero)\n", vx0open); 
						}
					} 
					dxbar = xbar[p][j] - xbar[p][j-1];
					if ( (j == jmopen[p][i]+1) &&  (j != 1) ) { // in this case dxbar is different //*
						dxbar = (xbar[p][j] - xbar[p][j-1])*(chiMax[p][j] - chi[p][j][i])/ (chiMax[p][j] - chi[p][j][i] + chi[p][j-1][i] - chiMax[p][j-1]);// open orbit density does not need to be so accurate at this point
						if (DEBUG == 1)
							printf("dxbar = %f\n", dxbar);
					}
					intdxbaropen += 0.5*(intdUopen+intdUopenold)*dxbar;
					intdxbaropenflow += 0.5*(intdUopenflow+intdUopenflowold)*dxbar;
					intdxbaropenBohm += 0.5*(intdUopenBohm+intdUopenBohmold)*dxbar;
					intdxbarxbar += 0.5*(intdUopenxbar+intdUopenxbarold)*dxbar;
					intdxbarsquare += 0.5*(intdUopensquare+intdUopensquareold)*dxbar; 
					if (intdxbaropen != intdxbaropen) {
						printf("PROBLEM HERE, j is %d\n", j); 
						exit(-1);
					} 
				}
				if ( j>=jmclosed[p][i] ) {	
					/* We have entered the closed orbit integral */
					intdvxold = intdvx;
					//intdvxflowold = intdvxflow;
					//intdvxflow = 0.0;
					intdvx = 0.0;
					if (j == jmclosed[p][i]) {
						intdvx = 0.0;
						//intdxbar += 0.0; 
					}
					else if (j > jmclosed[p][i]) {	
						if (DEBUG ==1) {
							printf("lowerlimit[%d] = %d\tupper[%d][%d] = %d\n", j, lowerlimit[p][j], j, i, upper[p][j][i]);
							printf("mu = %f\n", mu[p][j][lowerlimit[p][j]]);
							printf("imax[p][%d] = %d\n", j, imax[p][j]);
						}
						for (k=lowerlimit[p][j]; k<upper[p][j][i]+1; k++) {	
							vxold = vxnew;
							intdUold = intdU;
							intdU = 0.0;
							if (k == upper[p][j][i]) {
								Uperpnew = chi[p][j][i];
								if (lowerlimit[p][j] == upper[p][j][i] ) {
									//munew = 0.0;
									munew = mu[p][j][k]; 
									//printf("why would I ever be here? i = %d mu = %f\n\n\n", i, munew);
								}
								else if (k == upperlimit[p][j]) // correct but unnecessary as taken ino account below
									munew = 0.0;
								else {
									munew = ((chi[p][j][i] - Uperp[p][j][k])*mu[p][j][k-1] + (Uperp[p][j][k-1] - chi[p][j][i])*mu[p][j][k])/(Uperp[p][j][k-1] - Uperp[p][j][k]);
								}
								vxnew = 0.0;
								dvx = vxold - vxnew; 
							}
							else {	
								Uperpnew = Uperp[p][j][k];
								vxnew = vx[p][j][i][k];
								dvx = vxold - vxnew;
								munew = mu[p][j][k]; 
							}
							sizeU = (int) ( 2.0*sqrt(2.0*(Ucap - Uperpnew))/dvz );
							for (l=0; l < sizeU; l++) {	
								Fold = F;
								vzmin = -sqrt(2.0*(Ucap - Uperpnew)); 
								vzmin = 0.0;
								vz = vzmin + dvz*l;
								U = Uperpnew + 0.5*vz*vz;
								ystar = y_grid[p] - vz/alpha; 
								if (U<=munew) vzinf = 0.0;
								else vzinf = sqrt(2.0*(U - munew)); 
								phiinfval = lin_interp(y_grid, phiinf, ystar + vzinf/alpha, size_ygrid, 1);
								dphiinfdyval = lin_interp(y_grid, dphiinfdy, ystar + vzinf/alpha, size_ygrid, 1);
								if (U<=munew + phiinfval) vzinf = 0.0;
								else vzinf = sqrt(2.0*(U - munew - phiinfval)); 
								vzinf = vzinf + (1.0/alpha)*dphiinfdyval;
								yinf = ystar + vzinf/alpha;
								if (yinf > y_grid[size_ygrid-1]) yinf = y_grid[size_ygrid-1];
								if (yinf < y_grid[0]) yinf = y_grid[0];
								F = trilin_interp(yinf, munew, vzinf, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
								if (l!=0) 
									intdU += dvz*(F+Fold);
							} 
							intdUantycal = exp(-Uperpnew)*(1.0/(2.0*M_PI));// result with phi =0
							//intdU = intdUantycal;
							if (DEBUG == 1) 	
								printf("Analytical intdU is %f, numerical one is %f\n", intdUantycal, intdU);
							if (k==lowerlimit[p][j]) {	
								intdvx += 0.0; 
							}
							else {	
								//printf("k = %d\t intdU = %f \t dvx = %f\n", k, intdU, dvx);
								intdvx += 4.0*0.5*dvx*(intdU+intdUold);
							}
							if (intdvx != intdvx) {	
								printf("intdvx is NAN, j=%d, i=%d\n", j, i); 
								exit(-1);
							} 
						}
						//intdvxantycal = (2.0/(2.0*sqrt(M_PI)))*exp(-(xx[i]-xbar[j])*(xx[i]-xbar[j]))*erf(sqrt(xbar[j]*xbar[j]-(xx[i]-xbar[j])*(xx[i]-xbar[j])));
						//if (DEBUG == 1) 
						//	printf("pos=%f, i=%d, j=%d, Uperp is %f, chi is %f\nvxnew and vxold are %f and %f and and vx is %f, dvx is %f\nintdvx is %f, analytical one is %f, upper is %d, upperlimit is %d\n", xx[i], i, j, Uperpnew, chi[j][i], vxnew, vxold, vx[j][i][k], dvx, intdvx, intdvxantycal, upper[j][i], upperlimit[j]); 
						//printf("intdvx is %f, while the analytical one is %f\n", intdvx, intdvxantycal);
						//intdvx = intdvxantycal;
						if (j==jmclosed[p][i]+1 && i!=0) {
							//dxbar = (xx[imax[j]] - xx[imax[j]])*(xbar[j] - xbar[j-1])/(xx[imax[j]] - xx[i]);
							//dxbar = (chiMax[j]-chi[j][i])/(2.0*(xx[i]-xx[imax[j]]));
							if (i==imax[p][j]) {	
								dxbar = xbar[p][j] - xbar[p][j-1]; 
							}
							else {	
								//dxbar = xbar[j] - (xbar[j-1]*(chiMax[j] - chi[j][i]) - xbar[j]*(chiMax[j-1] - chi[j-1][i]))/(-chiMax[j-1] + chi[j-1][i] + chiMax[j] - chi[j][i]); 
								dxbar = (xbar[p][j] - xbar[p][j-1])*(chiMax[p][j] - chi[p][j][i])/ (chiMax[p][j] - chi[p][j][i] + chi[p][j-1][i] - chiMax[p][j-1]);// open orbit density does not need to be so accurate at this point
							}
							intdxbar += 0.5*(intdvx+intdvxold)*dxbar;
							if (intdxbar != intdxbar) {	
								printf("intdxbar is NAN, j=%d, i=%d\n", j, i);  
								exit(-1);
							} 
						}
						else {
							dxbar = xbar[p][j] - xbar[p][j-1];
							intdxbar += 0.5*(intdvx+intdvxold)*dxbar;
							if (intdxbar != intdxbar) {
								printf("intdxbar is NAN, j=%d, i=%d\n", j, i); 
								exit(-1);
							} 
						}
					} 
				} 
			}
			n_grid[p][ic] = intdxbar + intdxbaropen;
			if (ic == 0) {
				Bohm = intdxbaropenBohm/n_grid[p][0];
				printf("flow velocity at x=0 = %f\n", (fluxinf[p])/n_grid[p][0]);
				//printf("flux evaluated at x=0 is %f\n", intdxbaropenflow/(n_inf[p]*alpha));
				printf("flux evaluated at x=0 is %f\n", intdxbaropenflow);
				//printf("flux evaluated at x=0 is %f\n", *flux);
				//printf("intdxbarsquare = %f\tintdxbarxbar = %f\n", intdxbarsquare, intdxbarxbar);
				//kBohmsquared = intdxbarxbar/(intdxbarsquare - 0.5*n_grid[0]); 
				fluxzero[p] = intdxbaropenflow;
			}
			if ( n_grid[p][ic] >= stopdens*n_inf[p]) {	
				printf("ic = %d/%d\n", ic, size_phigrid);
				stop = 1;
				size_ngrid[p] = ic; 
				printf("stopping because density larger than %f times the density at infinity\n", stopdens);
			}
			else if (x_grid[ic] > x_grid[size_phigrid-1] - limit_rho) {
				printf("ic = %d/%d\n", ic, size_phigrid);
				stop = 1;
				size_ngrid[p] = ic; 
				printf("WARNING in iondens2d.c: stopping for positive density\n");

			}
			if ( (DEBUG == 1) ) {	
				printf("%f, %f, %f is TOTAL, CLOSED and OPEN orbit density at position index %d, position %f, potential %f\n", n_grid[p][ic], intdxbar, intdxbaropen, ic, xx[i], phi_grid[p][ic]);
				if ( (n_grid[p][ic] != n_grid[p][ic]) || (n_grid[p][ic] < TINY) ) { 	
					printf("n_finorb[ic] = %f\n", n_grid[p][ic]);
					//exit(-1);
				} 
			} 
			//fprintf(fout, "%f %f %f %f\n", xx[i], n_grid[p][ic]/n_inf[p], intdxbar/n_inf[p], intdxbaropen/n_inf[p]);
			ic += 1;
		} 
		//fclose(fout);
		//if (stop == 0) { 
		//	printf("ERROR: the density never reached stopdens*n_inf\n"); 
		//	exit(-1); 
		//}
		for (ic = 0; ic < size_phigrid; ic++) {
			if (DEBUG == 1) 
				printf("before renormalizing n_finorb[%d] = %f\tn_inf = %f\tphi_grid[%d] = %f\n", ic, n_grid[p][ic], n_inf[p], ic, phi_grid[p][ic]);
			if (ic < size_ngrid[p]) {
				n_grid[p][ic] /= n_norm;
				if (DEBUG ==1) {
					printf("x_grid[%d] = %f\tn_finorb[%d] = %f\tphi_grid[%d] = %f\n", ic, x_grid[ic], ic, n_grid[p][ic], ic, phi_grid[p][ic]);
				}
			}
			else n_grid[p][ic] = 0.0;
		}
		fluxzero[p] /= n_norm;
		fluxinf2[p] /= n_norm;
		intfluxinf += fluxinf2[p];
		intfluxzero += fluxzero[p];


		if (DEBUG == 1) { // change this to 0 to debug
			printf("mu vy  chiM dmudvy\n");
		}
		for (j=0; j<sizemumu; j++) {
			vy_op[p][j] = lin_interp(muopen[p], xbaropen[p], mumu[j], maxj[p], 2700);
			chiMax_op[p][j] = lin_interp(muopen[p], chiMopen[p], mumu[j], maxj[p], 2700);
			dmudvy_op[p][j] = lin_interp(muopen[p], twopidmudxbaropen[p], mumu[j], maxj[p], 2700);
			if (DEBUG == 1) 
				printf("%d/%d %f %f %f %f\n", j, sizemumu, mumu[j], vy_op[p][j], chiMax_op[p][j], dmudvy_op[p][j]);
		}

		// If you love your variables (and your memory) set them free // wise words Robbie
		free(chiMopen[p]);//
		free(twopidmudxbaropen[p]);//
		free(twopidmudyopen[p]);//
		free(xbaropen[p]);//
		free(muopen[p]);//
		free(jmclosed[p]);//
		free(jmopen[p]);//
		free(xifunction[p]);//
		free(xbar[p]);//
		free(chimin[p]);//
		free(chiMax[p]);//
		free(chimpp[p]);//
		free(crossed_min[p]);//
		free(crossed_max[p]);//
		free(twopidmudxbar[p]);//
		free(twopidmudy[p]);//
		free(upperlimit[p]);//
		free(lowerlimit[p]);//
		free(xtop[p]);//
		free(imax[p]);//
		free(imin[p]);//

		for (w = 0; w < sizexbar; w++) {
			free(chi[p][w]);//
			free(Uperp[p][w]);//
			free(mu[p][w]);//
			free(upper[p][w]);//
			for (s = 0; s < size_finegrid; s++)
				free(vx[p][w][s]);//
			free(vx[p][w]);//
		}
		free(chi[p]);//
		free(Uperp[p]);//
		free(mu[p]);//
		free(upper[p]);//
		free(vx[p]);//

		free(phi[p]);//
		free(phi_x[p]);//
		free(phi_y[p]);//

		for (j=0; j<j_inf; j++) {
			free(vxinf[p][j]);//
		}
		free(vxinf[p]); //
		free(chiinf[p]); //
		free(muinf[p]);//

	}
	free(chiMopen);//
	free(twopidmudxbaropen);//
	free(twopidmudyopen);//
	free(xbaropen);//
	free(muopen);//
	free(jmclosed);//
	free(jmopen);//
	free(xifunction);//
	free(xbar);//
	free(chimin);//
	free(chiMax);//
	free(chimpp);//
	free(crossed_min);//
	free(crossed_max);//
	free(twopidmudxbar);//
	free(twopidmudy);//
	free(upperlimit);//
	free(lowerlimit);//
	free(xtop);//
	free(imax);//
	free(imin);//
	free(xx);//
	free(gg);//
	free(ff);//
	free(chi);//
	free(Uperp);//
	free(mu);//
	free(upper);//
	free(vx);//
	free(phi);//
	free(phi_x);//
	free(phi_y);//
	free(vxinf); //
	free(chiinf); //
	free(muinf);//
	free(dphiinfdy);//
	free(n_inf); //
	free(Chodura); //
	free(fluxinf); //
	free(fluxzero); //
	free(icrit); //
	clock_t end = clock(); // finds the end time of the computation
	double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
	intfluxinf /= size_ygrid;
	intfluxzero /= size_ygrid;
	printf("intfluxzero = %f\tintfluxinf = %f\n", intfluxzero, intfluxinf);
	printf("in iondens2d: module ran in %f seconds\n", jobtime);
	return;
}

// the set of functions below are necessary to model electron reflection from the infinitely thin Debye sheath obtain parallel velocity cutoff as a function of magnetic moment with assumption rho_e >> lambda_D
//double FF_e(double beta) {
//	double integrand, integral=0.0, beta_s = 0.005, betaval;
//	int ind, num_beta = (int) (beta/beta_s);
//	for (ind = 0; ind < num_beta+1; ind++) {
//		betaval = beta*ind/num_beta;
//		integrand = 0.5*(1.0-cos(2.0*betaval)) / ( M_PI - betaval + 0.5*sin(2.0*betaval) );
//		integral += integrand*beta/num_beta;
//	}
//        return integral;
//}
//double dFF_edbeta(double beta) {
//	double integrand = 0.5*(1.0-cos(2.0*beta)) / ( M_PI - beta + 0.5*sin(2.0*beta) );
//        return integrand;
//}
//double vparcut(double beta, double vcut) {
//	double vparcutval = sqrt((1.0-exp(-2.0*FF_e(beta))))*vcut/sin(beta);
//	return vparcutval;
//}
//double mucut(double beta, double vcut) {
//	double mucutval = 0.5*exp(-2.0*FF_e(beta))*pow(vcut/sin(beta), 2.0);
//	return mucutval;
//}
//double vparcut_mu(double mu, double vcut) {
//	double beta, musearch=-10000.0, vparcutn;
//	double betalow = 0.0, betaup = M_PI;
//	do {
//		beta = (betaup + betalow)/2.0;
//		musearch = 0.5*exp(-2.0*FF_e(beta))*pow(vcut/sin(beta), 2.0);
//		if (musearch < mu) {
//			betaup = beta;
//		}
//		else {
//			betalow = beta;
//		}
//		printf("%f %f\n", musearch, mu);
//	} while (fabs(musearch-mu) > TINY) ;
//	vparcutn = sqrt((1.0-exp(-2.0*FF_e(beta))))*vcut/sin(beta);
//	return vparcutn;
//}

//void newvcut(double *v_cut, double v_cutDS, double mioverme, double TioverTe, double u_i, double u_e, double current, double error_current, double weight) {
//	double old_v_cut = *v_cut;// u_etilde;
//	//u_etilde = u_e - sqrt(mioverme/M_PI)*0.5*exp(-0.5*old_v_cut*old_v_cut);
//	//*v_cut = (1.0-weight)*(old_v_cut) + weight*sqrt(-2.0*log(2.0*(-current + u_i - u_etilde)) - log(M_PI/mioverme) ); //OLD
//	printf("v_cut before = %f\n", v_cutDS);
//	*v_cut = sqrt( 2.0*( (1.0-weight)*0.5*old_v_cut*old_v_cut + weight* ( 0.5*old_v_cut*old_v_cut + (current - u_i + u_e)/( sqrt(mioverme/(TioverTe))*0.5*exp(-0.5*old_v_cut*old_v_cut) ) ) ) );
//	if ( (*v_cut*(*v_cut)*0.5 < 0.5*v_cutDS*v_cutDS ) )  {
//		printf("WARNING: v_cutDS > v_cut\n");
//		*v_cut = v_cutDS;
//	}
//	printf("v_cut after = %f\n", v_cutDS);
//}

void error_quasi(double *error, double *x_grid, double **ne_grid, double **ni_grid, double **phi_grid, int size_ygrid, int size_phigrid, int *size_ngrid) {
	int i, p, count=0;
	double res = 0.0, dev, devbig;
	// Calculate the residual of quasineutrality equation
	res = 0.0;
	devbig = dev = 0.0;
	for (p=0; p<size_ygrid; p++)  {
		for (i=size_ngrid[p]-1; i>=0; i--)  {
			//if (i==size_ngrid[p]-1)	printf("x phi ne ni rel_err \n");
			//printf("%f %f %f %f %f\n", x_grid[i], phi_grid[p][i], ne_grid[p][i], ni_grid[p][i], -ne_grid[p][i]/ni_grid[p][i] + 1.0);
			dev = fabs(-ne_grid[p][i]/ni_grid[p][i] + 1.0);
			res += fabs(-ne_grid[p][i]/ni_grid[p][i] + 1.0);
			if (dev > devbig)  devbig = dev; 
			count += 1;
		}
	}
	res /= (count);
	printf("res = %f\n", res);
	error[0] = res;
	error[1] = devbig;
	return;
}

void newguess(double *x_grid, double** ne_grid, double **ni_grid, double** phi_grid, double *phiinf_grid, int size_phigrid, int *size_ngrid, int size_ygrid, double pfac, double weight) {
	clock_t begin = clock(); 
	int i, j, s;
	double phiW_impose, **newphi, phi0, phip0, CC, pdec;

	pdec = 2.0/(1.0-pfac);
	printf("φ'' ~ φ^%f gives φ ~ x^%f\n", pfac, pdec);

	newphi = malloc(size_ygrid*sizeof(double)); // same as above 

	for (j=0; j < size_ygrid; j++) {
		newphi[j] = malloc(size_phigrid*sizeof(double)); // same as above 
		for (i=0; i< size_phigrid; i++) {
			if (i< size_ngrid[j]-1)
				//newphi[j][i] = ( ni_grid[j][i] - ne_grid[j][i] )*exp(-phi_grid[j][i]) + phi_grid[j][i];
				//newphi[j][i] = ( ni_grid[j][i] - ne_grid[j][i] ) + phi_grid[j][i];
				newphi[j][i] = log(ni_grid[j][i]);
			else if (i== size_ngrid[j]-1) {
				//newphi[j][i] = ( ni_grid[j][i] - ne_grid[j][i] )*exp(-phi_grid[j][i]) + phi_grid[j][i];
				//newphi[j][i] = ( ni_grid[j][i] - ne_grid[j][i] ) + phi_grid[j][i];
				newphi[j][i] = log(ni_grid[j][i]);
				phip0 = (newphi[j][i] - newphi[j][i-2])/(x_grid[i] - x_grid[i-2]);
				CC = pdec*(newphi[j][i-1] - phiinf_grid[j])/phip0 - x_grid[i-1];
				//CC = 0.0;
				phi0 =  (newphi[j][i-1] - phiinf_grid[j])/pow(x_grid[i-1]+CC, pdec) ;
				//phi0 =  phip0/(pdec*pow(x_grid[i-1]+CC, pdec-1)) ;
				printf("CC = %f and phi0 = %f\n", CC, phi0);
			}
			else 
				newphi[j][i] = phiinf_grid[j] + phi0*pow(x_grid[i]+CC, pdec);
		}
	}
	for (j=0; j<size_ygrid; j++) {
		for (i=size_phigrid-1; i>=0; i--) {
			//printf("j=%d/%d\ti=%d/%d\n", j, size_ygrid, i, size_phigrid);
			newphi[j][i] = weight*newphi[j][i] + (1.0-weight)*phi_grid[j][i] ; 
			phi_grid[j][i] = newphi[j][i];
		}
		free(newphi[j]);
	}
	free(newphi);

	clock_t end = clock(); // finds the end time of the computation
	double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Module to obtain the new guess for the electrostatic potential ran in %f\n", jobtime);
	return;
}
//void newguess(double *x_grid, double** ne_grid, double **ni_grid, double** phi_grid, int size_phigrid, int *size_ngrid, int size_ygrid, double pfac, double weight) {
//	clock_t begin = clock(); 
//	int i, j;
//	double phiW_impose;
//	double *phi_xp_red, **newphi;
//	//double res = 0.0, reslimit, dev, dev_0, devbig;
//	double phi0, phi_x0, CC, pdec;
//	pdec = 2.0/(1.0-pfac);
//	printf("φ'' ~ φ^%f gives φ ~ x^%f\n", pfac, pdec);
//	/* Below we initialize all arrays that contain functions of position x with the correct size n */
//	newphi = (double**)malloc(size_ygrid*sizeof(double)); // same as above 
//	for (j=0; j<size_ygrid; j++) 
//		newphi[j] = (double*)malloc(size_phigrid*sizeof(double)); // same as above 
//	for (j=0; j < size_ygrid; j++) {
//		for (i=0; i< size_phigrid; i++) {
//			if (i< size_ngrid[j]-1)
//				newphi[j][i] = ( ni_grid[j][i] - ne_grid[j][i] )*exp(-phi_grid[j][i]) + phi_grid[j][i];
//			else if (i== size_ngrid[j]-1) {
//				newphi[j][i] = ( ni_grid[j][i] - ne_grid[j][i] )*exp(-phi_grid[j][i]) + phi_grid[j][i];
//				phi_x0 = (newphi[j][i] - newphi[j][i-2])/(x_grid[i] - x_grid[i-2]);
//				CC = pdec*newphi[j][i-1]/phi_x0 - x_grid[i-1];
//				//CC = 0.0;
//				phi0 =  newphi[j][i-1] /pow(x_grid[i-1]+CC, pdec) ;
//				//phi0 =  phi_x0/(pdec*pow(x_grid[i-1]+CC, pdec-1)) ;
//				printf("CC = %f and phi0 = %f\n", CC, phi0);
//			}
//			else 
//				newphi[j][i] = phi0*pow(x_grid[i]+CC, pdec);
//		}
//	}
//
//	for (j=0; j < size_ygrid; j++) {
//		for (i=size_phigrid-1; i>=0; i--) {
//			newphi[j][i] = weight*newphi[j][i] + (1.0-weight)*phi_grid[j][i] ; 
//			phi_grid[j][i] = newphi[j][i];
//			//printf("newphi[%d][%d] = %f\n", j, i, phi_grid[j][i]);
//		}
//	}
//	for (j=0; j<size_ygrid; j++) free(newphi[j]);
//	free(newphi);
//
//	clock_t end = clock(); // finds the end time of the computation
//	double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
//	printf("Module to obtain the new guess for the electrostatic potential ran in %f\n", jobtime);
//	return;
//}


// The main function of MAGSHEATH
int main() {
	//double power = 1.0;
	clock_t begin_it = clock(); // Finds the start time of the computation
	double **vy_i_wall, **chiM_i, **twopidmudvy_i;
	double phiamp = 0.01;
	int size_mu_i, size_U_i, size_mu_e, size_vpar_e;
	int i, j, k, p, ncols=0, ndirname, fix_current=0;
	double error_MP[2], olderror_MP, Bohmshouldbe;
	double weight_j=WEIGHT_j, weight_MP=WEIGHT_MP;
	double TioverTe, phi0_init_MP = -0.0;
	double ***dist_i_GK, *mu_i, *U_i;
	double **dist_e_DK, *vpar_e, *mu_e;
	double dvpar = DV;
	char dirname[200];
	int convergence_MP = 0, convergence_j = 0;
	int N=0;
	double tot_time, v_cut , current, flux_e=0.0, *flux_i, mioverme;
	double *vpar_e_cut;
	double **ne_grid, **ni_grid, v_cutDS = 0.5;
	double target_current;
	double *storevals;
	// Pointers used to read the files
	double *x_grid, *y_grid, **phi_grid, *phiinf_grid;
	double Te, alpha, alpha_deg, deltaovalpha, delta, system_size;
	int size_phigrid, size_ygrid=SIZEY, *size_ngrid;
	char line_hundred[100];
	double grid_parameter, deltax, deltay;

	flux_i = malloc(size_ygrid*sizeof(double));
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
	if ((input = fopen("inputfile2d.txt", "r")) == NULL) { 
		printf("Cannot open %s\n", "inputfile2d.txt");
		fprintf(fout, "Cannot open %s\n", "inputfile2d.txt");
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
		if (i!= 4) {
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
		if (i==0) {
			// input in degrees only for convenience
			// should be below 5 degrees (0.1 rad) for asymptotic theory to be valid
			alpha_deg = *storevals; 
		}
		if (i==1) {
			// input in degrees only for convenience
			// should be below 5 degrees (0.1 rad) for asymptotic theory to be valid
			deltaovalpha = *storevals; 
		}
		else if (i==2) {
			TioverTe = *storevals; //linetodata(line_hundred, strlen(line_hundred), &ncols); 
		}
		else if (i==3) {
			mioverme = *storevals; 
		}
		else if (i==4) {
			fix_current = (int) (*storevals); 
		}
		else if (i==5) {
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
		if (i<4) {
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
	if (TioverTe < 0.6) weight_MP = WEIGHT_MP/3.0;
	printf("directory where output will be stored is %s\n", dirname);
	mkdir(dirname, S_IRWXU);
	Te = 1/TioverTe; // Te/Ti is used and labelled Te in code
	printf("Te = %f\n", Te);
	alpha = alpha_deg*M_PI/180; // alpha used in radians from now on
	delta = deltaovalpha*alpha;
	// initial iteration assumes flat potential profile in magnetic presheath
	// therefore, the parameter v_cutDS is equal to v_cut
	v_cutDS = v_cut;
	if (Te < 1.0) system_size = SYS_SIZ;
	else system_size = SYS_SIZ*sqrt(Te);

	i=0;
	printf("INPUT PARAMETERS:\n\tmagnetic field angle = α = %f deg (%f rad)\n", alpha_deg, alpha);
	fprintf(fout, "INPUT PARAMETERS:\n\tmagnetic field angle = α = %f deg (%f rad)\n", alpha_deg, alpha);
	printf("\tion temperature = τ = ZT_i/T_e for species %d = %f\n\tmass ratio = m_i/m_e for species ? = %f\n", i+1, TioverTe, mioverme);
	fprintf(fout, "\tion temperature = τ = ZT_i/T_e for species %d = %f\n\tmass ratio = m_i/m_e for species ? = %f\n", i+1, TioverTe, mioverme);
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
	size_ngrid = malloc(size_ygrid*sizeof(int));;

	printf("size of coarse potential grid in magnetic presheath = %d\n", size_phigrid);
	fprintf(fout, "size of coarse potential grid in magnetic presheath = %d\n", size_phigrid);
	//printf("size of coarse density grid in magnetic presheath = %d\n", size_ngrid);
	//fprintf(fout, "size of coarse potential grid in magnetic presheath = %d\n", size_ngrid);

	deltay = 2.0*M_PI/(delta*size_ygrid);
	printf("deltay = %f\n", deltay);
	//delta*size_ygrid*deltay = 2.0*M_PI

	// form x, phi, ne and ni grids for magnetic presheath
	x_grid = malloc(size_phigrid*sizeof(double));
	y_grid = malloc(size_ygrid*sizeof(double));
	for (p=0; p<size_ygrid; p++) y_grid[p] = deltay*p;

	phi_grid = malloc(size_ygrid*sizeof(double));
	phiinf_grid = malloc(size_ygrid*sizeof(double));
	ne_grid  = malloc(size_ygrid*sizeof(double));
	ni_grid  = malloc(size_ygrid*sizeof(double));
	for (p=0; p<size_ygrid; p++) {	
		phi_grid[p] = malloc(size_phigrid*sizeof(double));
		ne_grid[p] =  malloc(size_phigrid*sizeof(double));
		ni_grid[p] =  malloc(size_phigrid*sizeof(double));
	}

	make_phiinfgrid(0, phiamp, delta, y_grid, phiinf_grid, size_ygrid);
	make_phigrid(x_grid, y_grid, phi_grid, phiinf_grid, size_phigrid, deltax, size_ygrid, grid_parameter);

	for (p=0; p<size_ygrid; p++)
		printf("phiinfgrid[%d] = %f\n", p, phiinf_grid[p]);

	size_mu_i = (int) (MAXV/DV);
	size_U_i = (int)  (MAXV/DV);
	size_mu_e = (int) (MAXV/DV);
	size_vpar_e = (int) (MAXV/DV);
	dist_i_GK = malloc(size_ygrid*sizeof(double));
	mu_i = malloc(size_mu_i*sizeof(double));
	mu_e = malloc(size_mu_e*sizeof(double));
	U_i = malloc(size_U_i*sizeof(double));
	vpar_e = malloc(size_mu_e*sizeof(double));
	for (k=0; k < size_ygrid; k++) {
		dist_i_GK[k] = malloc(size_mu_i*sizeof(double));
		for (i=0; i < size_mu_i; i++) 
			dist_i_GK[k][i] = malloc(size_U_i*sizeof(double));
	}
	dist_e_DK = malloc(size_mu_e*sizeof(double));
	for (i=0; i < size_mu_e; i++) 
		dist_e_DK[i] = malloc(size_vpar_e*sizeof(double));
	Figen(mioverme, TioverTe, dist_i_GK, phiinf_grid, y_grid, mu_i, U_i, size_ygrid, size_U_i, size_mu_i, DV, DV, deltay);
	Fegen(dist_e_DK, mu_e, vpar_e, size_vpar_e, size_mu_e, DV, DV);
	vpar_e_cut  = malloc(size_mu_e*sizeof(double));
	vy_i_wall  = malloc(size_ygrid*sizeof(double));
	chiM_i  = malloc(size_ygrid*sizeof(double));
	twopidmudvy_i  = malloc(size_ygrid*sizeof(double));
	for (p=0; p<size_ygrid; p++) {
		vy_i_wall[p]  = malloc(size_mu_i*sizeof(double));
		chiM_i[p]  = malloc(size_mu_i*sizeof(double));
		twopidmudvy_i[p]  = malloc(size_mu_i*sizeof(double));
	} 

	// CARRY OUT A SINGLE ION DENSITY CALCULATION
	if (MAX_IT == 0) {	
		//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 5.0), 2.0) - grid_parameter ) / deltax );
		iondens2d(Te, alpha, size_ygrid, size_phigrid, size_ngrid, ni_grid, y_grid, x_grid, phi_grid, phiinf_grid, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, flux_i, ZOOM_MP, STOP_MP, vy_i_wall, chiM_i, twopidmudvy_i);
	}

	if (TEST==1) {
		v_cut = 2.0;
		//v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]) ;
		//for (i=0; i< size_mu_e; i++) 
		//	vpar_e_cut[i] = v_cutDS;
		printf("now evaluate electron density in MPS\n");
		//denszeroorb(-1.0, 1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, vpar_e_cut, 0.0, x_grid);
		for (p=0; p<size_ygrid; p++) 
			for (i=0; i<size_phigrid; i++) 
				ne_grid[p][i] = exp(phi_grid[p][i]);
		//for the moment assume Boltzmann electrons


		printf("now evaluate ion density in MPS\n");
		//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 5.0), 2.0) - grid_parameter ) / deltax );
		iondens2d(Te, alpha, size_ygrid, size_phigrid, size_ngrid, ni_grid, y_grid, x_grid, phi_grid, phiinf_grid, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, flux_i, ZOOM_MP, STOP_MP, vy_i_wall, chiM_i, twopidmudvy_i); 
		exit(0);
	}
	N=0;
	//ITERATE MAGNETIC PRESHEATH AND DEBYE SHEATH POTENTIAL TO FIND SOLUTION
	// SIMPLIFIED ELECTRON MODEL USED TO CALCULATE DEBYE SHEATH POTENTIAL DROP
	weight_j=WEIGHT_j;
	error_MP[0] = 10000.0;
	while ( ( (convergence_MP <= 1) || ( convergence_j <= 1 && fix_current == 1 ) ) && (N<MAX_IT) ) {
		olderror_MP = error_MP[0];
		printf("weight_MP = %f\n", weight_MP);
		printf("ITERATION # = %d\n", N);
		fprintf(fout, "ITERATION # = %d\n", N);
		current = target_current;
		printf("grid parameter = %f\n", grid_parameter);
		fprintf(fout, "grid parameter = %f\n", grid_parameter);
		printf("\t(phi_mp0, phi_ds0, phi_wall) = (%f, %f, %f)\n", phi_grid[0][0], -0.5*v_cut*v_cut - phi_grid[0][0], -0.5*v_cut*v_cut);
		fprintf(fout, "\t(phi_mp0, phi_ds0, phi_wall) = (%f, %f, %f)\n", phi_grid[0][0], -0.5*v_cut*v_cut - phi_grid[0][0], -0.5*v_cut*v_cut);

		printf("evaluate ion density in MPS\n");
		iondens2d(Te, alpha, size_ygrid, size_phigrid, size_ngrid, ni_grid, y_grid, x_grid, phi_grid, phiinf_grid, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, flux_i, ZOOM_MP, STOP_MP, vy_i_wall, chiM_i, twopidmudvy_i); 
		printf("iondens2d module ran\n");
		//printf("size of ion density grid = %d\n", size_ngrid);
		//fprintf(fout, "size of ion density grid = %d\n", size_ngrid);
		//if (v_cut*v_cut < - 2.0*phi_grid[0]) v_cut = sqrt(-2.0*phi_grid[0]) + TINY;
		//v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
		//if (v_cut*v_cut < - 2.0*phi_grid[0]) v_cutDS = 0.0;
		//printf("v_cutDS = %f\n", v_cutDS);
		//if (EGYRO == 1) {
		//	for (i=0; i< size_mu_e; i++) 
		//		vpar_e_cut[i] = vparcut_mu(mu_e[i], v_cutDS);
		//}
		//else {
		//	for (i=0; i< size_mu_e; i++) 
		//		vpar_e_cut[i] = v_cutDS;
		//}
		//denszeroorb(-1.0, 1.0, phi_grid, ne_grid, size_phigrid, &flux_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, vpar_e_cut, 0.0, x_grid);
		for (p=0; p<size_ygrid; p++) {
			//printf("p=%d/%d\n", p, size_ygrid);
			for (i=0; i<size_phigrid; i++) {
				//printf("i=%d/%d\n", i, size_phigrid);
				ne_grid[p][i] = exp(phi_grid[p][i]);
			}
		}
		//for the moment assume Boltzmann electrons
		Bohmshouldbe = (1.0/Te);
		printf("Bohm integral should be %f\n", Bohmshouldbe);
		fprintf(fout, "Bohm integral of converged MP solution should be %f\n", Bohmshouldbe);
		error_quasi(error_MP, x_grid, ne_grid, ni_grid, phi_grid, size_ygrid, size_phigrid, size_ngrid);
		printf("error_av = %f\terror_max = %f\n", error_MP[0], error_MP[1]);
		if ( (error_MP[0] < tol_MP[0]) && (error_MP[1] < tol_MP[1]) ) convergence_MP += 1 ;
		else convergence_MP = 0;
		if (convergence_MP == 0) {
			printf("MP not converged --> calculate new MP potential guess\n");
			fprintf(fout, "MP not converged --> calculate new MP potential guess\n");
			newguess(x_grid, ne_grid, ni_grid, phi_grid, phiinf_grid, size_phigrid, size_ngrid, size_ygrid, 1.5, weight_MP);// p, m);
		}
		else {
			printf("MP converged --> no iteration needed\n");
			fprintf(fout, "MP converged --> no iteration needed\n");
		}
		//current = flux_i - flux_e*sqrt(Te*mioverme);
		//if (fix_current == 1) {
		//	if (fabs((target_current - current)/flux_i) > tol_current) convergence_j = 0;
		//	else convergence_j += 1;
		//	printf("target current = %f +/- %f\n", target_current, tol_current*flux_i);
		//	fprintf(fout, "target current = %f +/- %f\n", target_current, tol_current*flux_i);
		//	printf("current = %f = ion_current (%f) - electron_current (%f) = %f x ion_current\n", current, flux_i, flux_e*sqrt(Te*mioverme), current/flux_i);
		//	fprintf(fout, "current = %f = ion_current (%f) - electron_current (%f) = %f x ion_current\n", current, flux_i, flux_e*sqrt(Te*mioverme), current/flux_i);
		//	if (convergence_j == 0) {
		//		fprintf(fout, "new WALL potential guess\n");
		//		printf("new WALL potential guess\n");
		//		newvcut(&v_cut, v_cutDS, mioverme, TioverTe, flux_i, flux_e*sqrt(Te*mioverme), target_current, tol_current, weight_j);
		//	}
		//	else {
		//		printf("current converged, NO wall potential iteration\n");
		//		fprintf(fout, "current converged, NO wall potential iteration\n");
		//	}
		//}
		//else {
		//	printf("current = %f = ion_current (%f) - electron_current (%f) = %f x ion_current\n", current, flux_i, flux_e*sqrt(Te*mioverme), current/flux_i);
		//	fprintf(fout, "current = %f = ion_current (%f) - electron_current (%f) = %f x ion_current\n", current, flux_i, flux_e*sqrt(Te*mioverme), current/flux_i);
		//}
		//printf("\tcurrent = %f (target = %f)\n", current, target_current);
		//fprintf(fout, "\tcurrent = %f (target = %f)\n", current, target_current);
		printf("convergence (MP, j) = (%d, %d)\n", convergence_MP, convergence_j);
		if ( (convergence_MP > 1) && ( convergence_j > 1 || fix_current == 0) ) {	
			//printf("ENTER ION DENSITY EVALUATION IN MPS\n");
			//The argument of iondens2d set to one makes the module compute the ion distribution function at x=0
			iondens2d(Te, alpha, size_ygrid, size_phigrid, size_ngrid, ni_grid, y_grid, x_grid, phi_grid, phiinf_grid, dist_i_GK, mu_i, U_i, size_mu_i, size_U_i, grid_parameter, flux_i, ZOOM_MP, 1.000, vy_i_wall, chiM_i, twopidmudvy_i); 
			error_quasi(error_MP, x_grid, ne_grid, ni_grid, phi_grid, size_ygrid, size_phigrid, size_ngrid);
			printf("error_av = %f\terror_max = %f\n", error_MP[0], error_MP[1]);
			//if ( error_MP[1] > 2.0*tol_MP[1] ) {
			//	convergence_MP = 0;
			//	N++;
			//}
		} 
		else N++;
	}
	clock_t end_it = clock(); // finds end time of last iteration
	tot_time = (double) (end_it - begin_it) / CLOCKS_PER_SEC;
	printf("At %dth iteration MP iteration with simplified DS model converged successfully in %f seconds\n", N, tot_time);
	fprintf(fout, "At %dth iteration MP iteration with simplified DS model converged successfully in %f seconds\n", N, tot_time);

	if (N== MAX_IT) {
		printf("No convergence after %d iterations :(\n", MAX_IT);
		fprintf(fout, "No convergence after %d iterations :(\n", MAX_IT);
	}
	fclose(fout);

	// Print potential and density profiles to files
	FILE *fp; 
	char fpstr[150];
	printf("ndirname = %d\n", ndirname);
	snprintf(fpstr, 150, "%s/phi_n_MP.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) printf("error when opening file\n");
	for (p=0; p<size_ygrid; p++) {
		for (i=0; i<size_phigrid; i++) {
			fprintf(fp, "%f %f %f %f\n", x_grid[i], phi_grid[p][i], ni_grid[p][i], ne_grid[p][i]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	snprintf(fpstr, 150, "%s/phi_MP2d.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) printf("error when opening file\n");
	for (i=0; i<size_phigrid; i++) {
		for (p=0; p<size_ygrid; p++) {
			fprintf(fp, "%f ", phi_grid[p][i]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	snprintf(fpstr, 150, "%s/x_MP2d.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) printf("error when opening file\n");
	for (i=0; i<size_phigrid; i++) {
		fprintf(fp, "%f ", x_grid[i]);
	}
	fclose(fp);
	snprintf(fpstr, 150, "%s/y_MP2d.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) printf("error when opening file\n");
	for (p=0; p<size_ygrid; p++) {
		fprintf(fp, "%f ", y_grid[p]);
	}
	fclose(fp);
	snprintf(fpstr, 150, "%s/inputfile2d.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file\n");
	}
	fprintf(fp, "%f\n", alpha_deg);
	fprintf(fp, "%f\n", deltaovalpha);
	fprintf(fp, "%f\n", TioverTe);
	fprintf(fp, "%f\n", mioverme);
	fprintf(fp, "%d\n", fix_current);
	if (fix_current == 1) fprintf(fp, "%f\n", target_current);
	else fprintf(fp, "%f\n", v_cut);
	fclose(fp);

	for (p=0; p<size_ygrid; p++) {
		free(phi_grid[p]);
		free(ne_grid[p]);
		free(ni_grid[p]);
	}
	free(phi_grid);
	free(vpar_e_cut);
	free(x_grid);
	free(ne_grid);
	free(ni_grid);
	free(size_ngrid);
	//return 0;
	exit(0);
}
//Authors: Alessandro Geraldini and Robbie Ewart 

//#include <stdlib.h>
