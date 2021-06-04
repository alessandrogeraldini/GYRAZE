// Author: Alessandro Geraldini
// MODIFIED 24 JUN 2019
/* This code calculates the lowest order ion density profile in the magnetic presheath for a given potential profile and entrance distribution function. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "mps.h"
#define TINY 1e-12
 //TINY is just a small number used to make some inequalities work numerically in case of exact equality. */

double tophat(double x1, double x2, double x) {
	double y;
	if ((x1 <= x ) && (x <= x2 ))
	{
		y = 1.0;
	}	
	else
	{
		y = 0.0;
	}
	return y;
}

int densfinorb(int sizevxopen, double Te, double alpha, int size_phigrid, int *size_ngrid, double* n_grid, double *x_grid, double* phi_grid, double charge, double **FF, double *mumu, double *UU, int sizemumu, int sizeUU, double grid_parameter, double *vx_i_DS, double *fx_i_DS, double *flux) {
// DECLARATIONS
//inpu: double rho_over_unitgrid = for ions in the magnetic presheath this is 1.0, for electrons in the Debye sheath this can be chosen as min(1.0, rho_e/lambda_D)
//double temp_over_unitgrid = 1.0;
clock_t begin = clock(); // Finds the start time of the computation
int zoomfactor = 2;
double phi_cut = -99999.0, n_inf=0.0; // keep phi_cut as large and negative as possible
int reflected = 0;
double Ucrit = 0.0, frac_reflected=0.0, frac = 0.0, vzcrit;
double deltax;
double *xx, *dxdf, pos, *phi, *phip, *ne, *niclosed, *niopen, *nitotal, **chi, **chiop, Telarge, stoplimit = 0.98;
/* pos is x, the position (distance from the wall); phi contains the values of phi(x), extracted from a file. phip is phi prime, first derivative of phi. ne is the electron density; niclosed is the closed orbit ion density, niopen is the open orbit ion density, nitotal is the sum of the two; newphi is the new electrostatic potential guess; chi is the effective potential;*/ 
int debug=0, stop = 0, sizexbar, problem = 0, noback = 0, maxj;
int icrit, ind;
/* n is the domain of the position x (the largest value of the index i, so e.g. x_grid[n] = L_2 in the paper); sizexbar is the domain of the position xbar (the largest value of the index j); */ 
double *openorbit, openorbitantycal, **mu, *muopen, *Ucritf, *xbar, *twopimuprime, *openintegral, *xifunction;
/* openorbit is the integral Delta_M in the paper (an integral over the last closed orbit, representing how much U_perp - chi_M change in one full last orbit); openorbitantycal is the analytical value of openorbit for a flat potential (the first potential guess usually); mu is the array containing values of mu(xbar, Uperp). The index j is for values of xbar, k is for values of Uperp; xbar is the grid of values used in the closed orbit integral; FF contains the distribution function, read from the file distfile.txt. UU and mumu contain the values of U and mu corresponding to the function FF (which is F(mu, U)); FFprime is the numerical first derivative of F with respect to U; Many of these arrays are pointers because we don't initially know the array size, until we read the files (e.g. the file with the potential and the one with the distribution function); Note: first index of FF and FFprime is mu, second one is U; Interesting output variables; xifunction is the function of x which defines the grid of values of xbar by finding a chi whose minimum lies exactly at each grid point x */
int *jmclosed, *jmopen,  i=0, ic=0, j=0, k=0, l=0, sizeU, size_finegrid;
/* jmclosed represent minimum values of xbar above which we integrate open and closed orbit density integrals respectively (xbar_m,o and xbar_m in the paper); i is an index usually representing the positin x  j is an index usually representing the orbit position xbar; k is an index usually representing the energy Uperp (or velocity vx). It's always used in conjunction with j (and sometimes i); l is an index (used in for loops) usually representing the total energy (or velocity vz). It's used only in the DENSITY INTEGRALS part of the code; sizeU is the size of the integration range over U (or velocity vz). It is set later on in the code */
int *crossed_max, *crossed_min, *kdrop;
/* crossed_max and crossed_min is non-zero when a minimum (or maximum) of chi is found for some xbar[j]; kdrop is an integer which is non-zero only if there is an additional potential barrier for the particle at x << rho e.g. for an ion if the sheath reverses Uperp is allowed to be above chiMax and therefore the k = 0 of Uperp[j][k] starts for Uperp[j][0] = chiMax[j] + phi_barrier instead of Uperp[j][0] = chiMax[j] which is the conventional way */
	
int temporary = 0;
/* the temporary flag is a fix that ensures the minimum value of chi is stored (otherwise the minimum flag which just turned on prevents it); */
int *lowerlimit, *upperlimit, **upper, *itop, *imax, *imin;
/* lowerlimit represents the lower limit of k in the integrals over Uperp (or vx). It's needed because some of the earlies energies; (which are the largest because thy are values of chi stored after the maximum is found); may be so large that they are associated with very small values of the distribution function. This avoids integrating in an empty portion of phase space; upperlimit[j] represents the largest value of k (the smallest stored energy Uperp = chi_minimum) associated with some value of j; upper[j][i] represents the value of k associated with the smallest value of vx when integrating over Uperp. Going above upperlimit[j][i] makes Uperp < chi so velocities imaginary; imax/imin[j] stores the position of the maximum/minimum of the effective potential chi (It's x_M/x_m in the paper, which depends on xbar). */
double **Uperp, ***vx, *chiMax, *chiMpp, *chimin, oorbintgrd, oorbintgrdantycal, oorbintgrdBohm, xbarmax, intdvxopen, intdvxopenBohm=0.0, intdvxopenfluidBohm=0.0, intdxbaropenold;
/* Uperp stores the possible values of Uperp associated with closed orbits, and so does vx; chiMax and chimin store the local maxima and minima of the effective potential maximum, oorbintgrd is the value of the integrand in the first open orbit integral (oorbintgrdantycal is the analytical result for flat potential); oorbintgrdBohm is the value of the integrand in the `Bohm' integral; intdvxopen and similars are where the integral of f_{0x} (v_x) over v_x, and its two important moments <v_x> (fluidBohm) and <1/v_x^2> (Bohm) are stored; They are a check that the extracted distribution function f_{0x} has the same moments it had when we carried out the density integral earlier; xbarmax is the maximum value of xbar in our grid*/
double vz, U, dvz = 0.2, dvzopen = 0.1, dvx, dvxopen, Deltavx, vxopen, dxbar, intdU=0.0, intdUopen=0.0, intdUopenBohm = 0.0;
	/* vz used in the density integral; U is the total energy, used in the density integral; dvz is the thickness of the vz grid used to take the integral over U (which is taken over vz in practice), dvzopen is the same for the open orbit piece; dvx is the thickness of the vx grid used to take the integral over Uperp ( which is taken over vx in practice). It must be evaluated because it depends on stored values of vx[j][i][k]; dvxopen is the thickness of the vx grid on which f_{0x} (v_x) is defined (and integrated to check consistency of its moments); dxbar is the thickness of the xbar grid; intdU is the value of the integral over U in the closed orbit density integration process; intdUopen is the same as above, for the open orbit integral; intdUopenBohm same, for Bohm integral */
double intdUold=0.0, intdvx=0.0, intdvxold = 0.0, intdxbar=0.0, intdxbaropen=0.0, intdxbaropenBohm = 0.0, Bohm = 0.0, F, Fold=0.0, Ucap;
	/* intdUold is a variable which stores the old intdU, so that the trapezium rule of integration can be applied (intdUold + intdU)*dvz; intdvx stores the integral over Uperp (hence over vx) in the closed orbit integral; intdxbar stores the value of the integral over xbar (which is the final result!), intdxbaropen does the same in the open orbit density integral; intdxbaropenBohm does the same for the Bohm integral; idealBohm is what the Bohm integral shoult be if Bohm condition is marginally satisfied; Bohm is the Bohm integral at the Debye sheath entrance x=0; F is the value of the distribution function evaluated in the density integrals by interpolating FF, and Fold is the `old' needed to apply the trapezium rule; Fprime is the bilinearly interpolated value of FFprime, and Fprimeold is the same at the previous grid point (needed for trapezium rule); used in INTEGRALS OF DISTRIBUTION FUNCTION AT INFINITY; Ucap is the topmost total energy integrated to */
double intdUopenflow = 0.0, intdUopenflowold = 0.0, intdxbaropenflow = 0.0, oorbintgrdflow = 0.0, oorbintgrdflowold = 0.0;
// values of various integrals
double oorbintgrdold=0.0, oorbintgrdBohmold=0.0, Fopen=0.0, intdUopenold=0.0, intdUopenBohmold;
double *chiprimetop, vx0open, aa, aaold; 
double intdUantycal=0.0, intdvxantycal=0.0, vxnew=0.0, vxold = 0.0, Uperpnew = 0.0, munew = 0.0, *xtop, intdUopenantycal=0.0;
/* intdUantycal is the integral over U (or v_z) for a flat potential profile (phi =0) for some value of xbar and Uperp; intdvxantycal  is the integral over Uperp (or vx) for a flat potential profile for some value of xbar; vxnew is the value of vx at the 'new' grid point, used in the vx integral (taken using the trapezium rule); vxold is the value of vx at the 'old' grid point, used in the vx integral; Uperpnew is the value of Uperp (used in the closed orbit density integral); munew is the valye of mu (used in the closed orbit density integral); xtop is the top bounce point x_t of the last closed orbit; intdUopenantycal is the analytical value of the integral over U  in the open orbit density integral */
double oorbintgrdxbar, oorbintgrdsquare, intdUopenxbar = 0.0, intdUopensquare = 0.0, intdxbarxbar=0.0, intdxbarsquare=0.0, kBohmsquared;
/* oorbintgrdxbar is an integral over the open orbit distribution function at x=0 which is needed to evaluate a coefficient that appears when; expanding quasineutrality near x=0. It was just for playing around and at the moment plays no role in the code; similarly with all other integrals here */
double oorbintgrdxbarold, oorbintgrdsquareold, intdUopenxbarold, intdUopensquareold;
// all used to take integrals above
double xi, *gg, *ff;
//double du, fluxinf1old, fluxinfintgrdold, fluxinfintgrd, u, Fprimeold, Fprime;
//double **distfunc_iprime, fluxinf, fluxinf1, densinf1, densinf, densinf1old;

printf("charge = %f\n", charge);

if (debug == 1) {
	for (i=0; i<size_phigrid;i++) {
		printf("index %d\tx = %f\tphi = %f\n", i, x_grid[i], phi_grid[i]);
	}
}

  gg = malloc(size_phigrid*sizeof(double));
  ff = malloc(size_phigrid*sizeof(double));
  for (i=0; i<size_phigrid; i++) {
    gg[i] = sqrt(x_grid[i]);
    ff[i] = pow(sqrt(grid_parameter)+gg[i], 2.0) - grid_parameter;
  }
  //size_finegrid = size_phigrid*zoomfactor;
  deltax = ff[1];
  size_finegrid = (int) (zoomfactor*ff[size_phigrid-1]/deltax);
  printf("size_finegrid = %d\n", size_finegrid);
  xx = (double*)calloc(size_finegrid,sizeof(double)); // xx = x has correct size
  dxdf = (double*)calloc(size_finegrid,sizeof(double)); // xx = x has correct size
  phi = (double*)calloc(size_finegrid,sizeof(double)); // phi now has correct size
  FILE *fp;
  i=0;
  //printf("size_finegrid = %d\n", size_finegrid);
  //printf("size_phigrid = %d\n", size_phigrid);
  //for (i = 0; i < size_phigrid; i += 1) printf("ff = %f\tphi = %f\n", ff[i], phi_grid[i]);
  //{
    gsl_interp_accel *acc
      = gsl_interp_accel_alloc ();
    gsl_spline *spline
      = gsl_spline_alloc (gsl_interp_cspline, size_phigrid);

    //gsl_spline_init (spline, gg, phi_grid, size_phigrid);
    gsl_spline_init (spline, ff, phi_grid, size_phigrid);
    //gsl_spline_init (spline, x_grid, phi_grid, size_phigrid);

    fp = fopen("phispline.txt", "w");
    if (fp == NULL) printf("Error: phispline not created\n");
    printf("phispline opened\n");
    //printf ("%f\n", gg[size_phigrid-1]);
    for (i = 0; i < size_finegrid; i += 1)
      {
	//printf("i = %d/%d\tdeltax=%f\tzoomfactor=%d\n", i, size_finegrid, deltax, zoomfactor);
	xi = TINY + i*deltax/zoomfactor;
	xx[i] = pow( pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), 2.0);
	dxdf[i] = sqrt(xx[i]/(grid_parameter+xi));
	//phi[i] = gsl_spline_eval (spline, xx[i]*rho_over_unitgrid, acc);
	//phi[i] = gsl_spline_eval (spline, sqrt(xx[i]*rho_over_unitgrid), acc);
	phi[i] = gsl_spline_eval (spline, xi, acc);
        fprintf(fp, "%f %f\n", pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), phi[i]);
	phi[i] *= charge;
	//printf("xx = %f\n", xx[i]);
	//xx[i] = xi*xi/rho_over_unitgrid;
        //fprintf(fp, "%f %f\n", xi, phi[i]);
      }
  fclose(fp);
  //}


// Introduce a number that equals Te when Te > 1 and 1 when Te<1. When Te is large, this number increases the number of grid points in tot. energy U to account for the fact that vz ~ vB > v_t,i when the electron temperature is large. Moreover, |v| ~ vB > v_t,i at the Debye sheath entrance, and so the number of grid points in xbar must also be increased.
if (Te>1.0) {	
	Telarge = Te;
}
else {
	Telarge = 1.0;
}
printf("Telarge = %f\n", Telarge);
xbarmax = 10.0*pow(Telarge, 0.5);
Ucap = 18.0 + 4.0*Te;


//
////if (phi[0] < 0.0) {
//distfunc_iprime = calloc(sizemumu,sizeof(double*));
//
///* The loop below evaluates the derivative of F with respect to total energy U at U = mu
// */
//for (i=0; i<sizemumu; i++)
//{
//	printf("i=%d/%d\n", i, sizemumu);
//	distfunc_iprime[i] = calloc(sizeUU,sizeof(double));
//	for (j=0; j < sizeUU; j++)
//	{
//		if ((j==sizeUU-1) || (i==sizemumu-1))
//		{
//			distfunc_iprime[i][j] = 0.0;	
//		}
//		else
//		{
//			distfunc_iprime[i][j] = (FF[i][j+1] - FF[i][j])/(UU[j+1] - UU[j]); 
//		}
//	}
//}
//
//fluxinf1old = fluxinf1 = fluxinf = 0.0;
//densinf1old = densinf1 = densinf = 0.0;
//du = 0.01;
//for (i=0; i<sizemumu; i++)
//{
//	if (i!=0)
//	{
//		fluxinf1old = fluxinf1;
//		fluxinf1 = 0.0;
//		densinf1old = densinf1;
//		densinf1 = 0.0;
//	}
//	Fprimeold = Fprime = 0.0;
//	Fold = F = 0.0;
//	fluxinfintgrdold = fluxinfintgrd = 0.0;
//	for (j=0; j< (int) 500*sqrt(1+Te); j++)
//	{
//		u = du*j;
//		fluxinfintgrdold = fluxinfintgrd;
//		Fprimeold = Fprime;
//		Fold = F;
//		U = pow(u, 2.0) + mumu[i];
//		Fprime = bilin_interp(mumu[i], U-mumu[i], distfunc_iprime, mumu, UU, sizemumu, sizeUU, -1, -1);
//		F = bilin_interp(mumu[i], U-mumu[i], FF, mumu, UU, sizemumu, sizeUU, -1, -1);
//		fluxinfintgrd = F*u*alpha;
//		fluxinf1 += 0.5*(fluxinfintgrd + fluxinfintgrdold)*du;
//		densinf1 += 0.5*(F + Fold)*du;
//	}
//	if (i!=0)
//	{
//		fluxinf += 0.5*(mumu[i]-mumu[i-1])*(fluxinf1 + fluxinf1old);
//		densinf += 0.5*(mumu[i]-mumu[i-1])*(densinf1 + densinf1old);
//	}
//}
//densinf *= (4.0*M_PI);
//fluxinf *= (4.0*M_PI);
//printf("densinf = %f\tfluxinf = %f\n", densinf, fluxinf);
////}


//printf("Ucap =%f\n", Ucap);
i=0;
/* We initialize all arrays that contain functions of position x with the correct size n */
ne = (double*)calloc(size_phigrid,sizeof(double)); // phi now has correct size
niclosed = (double*)calloc(size_phigrid,sizeof(double)); // array containing the density of closed orbits
niopen = (double*)calloc(size_phigrid,sizeof(double)); // array containing the density of open orbits
nitotal = (double*)calloc(size_phigrid,sizeof(double)); // array containing the total density of ions 
phip = (double*)calloc(size_finegrid,sizeof(double)); // phi now has correct size
jmclosed = (int*)calloc(size_finegrid,sizeof(int));
jmopen = (int*)calloc(size_finegrid,sizeof(int));
xifunction = (double*)calloc(size_finegrid,sizeof(double));

/* FORM XBAR GRIDS 
Take derivatives of phi and use them to obtain two grids for xbar, one to be used for closed orbits and one to be used for open orbits. */
xbar = (double*)calloc(size_finegrid,sizeof(double));	
//xbarop = (double*)calloc(n,sizeof(double)); // Open orbit grid defined such that every position x has a corresponding chi for which xM = x. This allows better resolution of the places where the open orbit integrand, oorbintgrd, is large (~alpha^(1/2))	
for (i=0; i<size_finegrid; i++)
{	
	//if (i!=0)	xiprime[i-1] = (xiprime[i]-xiprime[i-1])/(xx[i]-xx[i-1]); 
	// Evaluate derivative of phi
	jmopen[i] = jmclosed[i] = 0;
	if (i == 0)
	{	
		//phip[0] = (phi[1] - phi[0])/(xx[1]-xx[0]);
		phip[0] = (phi[1] - phi[0])/(xx[1]-xx[0]) + (xx[0] - xx[1])*((phi[2] - phi[1])/(xx[2]-xx[1]) - (phi[1] - phi[0])/(xx[1]-xx[0]))/(xx[2]-xx[1]) ; 
		//printf("phip[0] = %f\n", phip[0]);
		//phip[0] = (phi[1] - phi[0])/(deltax*dxdf[0]/zoomfactor);
		//printf("phip[0] = %f\n", phip[0]);
	}
	else if (i == size_finegrid-1)
	{	
		phip[i] = 0.5*(phi[i] - phi[i-1])/(xx[i] - xx[i-1]); 
		//printf("phip[%d] = %f\n", i, phip[i]);
		//phip[i] = (0.0 - phi[i-1])/(2.0*deltax*dxdf[i]/zoomfactor);
		//printf("phip[%d] = %f\n", i, phip[i]);
	}
		//phip[0] = phip[1] + (xx[0] - xx[1])*(phip[2] - phip[1])/(xx[2]-xx[1]); 
	else {	
		//phip[i] = (phi[i+1] - phi[i])/(xx[i+1] - xx[i]);
		phip[i] = ((xx[i] - xx[i-1])/(xx[i+1]- xx[i-1]))*(phi[i+1] - phi[i])/(xx[i+1] - xx[i]) + ((xx[i+1] - xx[i])/(xx[i+1]- xx[i-1]))*(phi[i] - phi[i-1])/(xx[i] - xx[i-1]); 
		//printf("phip[%d] = %f\n", i, phip[i]);
		//phip[i] = (phi[i+1] - phi[i-1])/(2.0*deltax*dxdf[i]/zoomfactor);
		//printf("phip[%d] = %f\n", i, phip[i]);
	}
	// xifunction is xbar corresponding to given position x being a stationary point
	xifunction[i] = xx[i] + phip[i]/2.0;
	if (debug == 1)
	 	printf("xifunction[%d] = %f\n", i, xifunction[i]); 
	if (i == 1) {	
		if  (xifunction[i] > xifunction[i-1]) { // immediately found that xi is increasing at x=0 telling us x_c = 0
			icrit = i-1;//=0;
		}
	}
	else if (i > 1) {	
		if ( (xifunction[i] > xifunction[i-1]) && (xifunction[i-1] < xifunction[i-2] ) ) {
		//if ( (xifunction[i] + TINY > xifunction[i-1]) && (xifunction[i-1] < xifunction[i-2] + TINY) ) 
			icrit = i-1;
		}
		//else if ( (xifunction[i] - TINY < xifunction[i-1]) && (xifunction[i-1] + TINY > xifunction[i-2]) ) 
		else if ( (xifunction[i] < xifunction[i-1]) && (xifunction[i-1]  > xifunction[i-2]) ) {
		// found maximum of xi: this only happens when phi has noise in second derivative
			//phip[i] = 0.5*(phi[i] - phi[i-1])/(xx[i] - xx[i-1]); 

			printf("ERROR in densfinorb: too much noise in second derivative\n");
			printf("i = %d\n", i);
			if (i<size_finegrid/2) {
				exit(-1);
			}
		}
	} 
}
printf("icrit = %d\n", icrit);
// Whole section below perhaps was overkill in xbar resolution
////////
//j=0;
//k = icrit+1; 
// the smallest value of xbar is not calculated at the critical point, but one step ahead
//for (i=icrit-1;i>=0;i--) {
//	while ( (xifunction[k] < xifunction[i]) && (k<size_finegrid) ) {	
//		xbar[j] = xifunction[k]; 
//		j++; k++; 
//	}	
//	xbar[j] = xifunction[i]; 
//	j++; 
//}
//while (k<size_finegrid) { 	
//	xbar[j] = xifunction[k]; 
//	k++; j++; 
//}	
////////
j=0;
for (k=icrit+1; k < size_finegrid; k++) {
	xbar[j] = xifunction[k]; 
	j++;
}
sizexbar = j;
if (debug == 1) {
	for (j = 0; j < sizexbar; j++)
		printf("xbar[%d] = %f\n", j, xbar[j]); 
}	
//printf("sizexbar = %d, size_finegrid (x) =%d\n", sizexbar, size_finegrid);
chi = (double**)calloc(sizexbar,sizeof(double*)); // chi(xbar, x) indices j and i 
chiop = (double**)calloc(sizexbar,sizeof(double*)); // chi(xbar, x) indices j and i 
Uperp = (double**)calloc(sizexbar,sizeof(double*)); 
mu = (double**)calloc(sizexbar,sizeof(double*)); 
muopen = (double*)calloc(sizexbar+1,sizeof(double));	
Ucritf = (double*)calloc(sizexbar+1,sizeof(double));	
vx = (double***)calloc(sizexbar,sizeof(double**)); 
upper = (int**)calloc(sizexbar,sizeof(int*)); 
chimin = (double*)calloc(sizexbar,sizeof(double)); 
chiMax = (double*)calloc(sizexbar,sizeof(double));
chiMpp = (double*)calloc(sizexbar,sizeof(double));
crossed_min = (int*)calloc(sizexbar,sizeof(int)); 
crossed_max = (int*)calloc(sizexbar,sizeof(int));
kdrop = (int*)calloc(sizexbar,sizeof(int));
twopimuprime = (double*)calloc(sizexbar,sizeof(double));
openintegral = (double*)calloc(sizexbar,sizeof(double));
openorbit = (double*)calloc(sizexbar,sizeof(double));
upperlimit = (int*)calloc(sizexbar,sizeof(int));
lowerlimit = (int*)calloc(sizexbar,sizeof(int));
chiprimetop = (double*)calloc(sizexbar,sizeof(double));
xtop = (double*)calloc(sizexbar,sizeof(double));
itop = (int*)calloc(sizexbar,sizeof(int));
imax = (int*)calloc(sizexbar,sizeof(int));
imin = (int*)calloc(sizexbar,sizeof(int));
/////////////////////////////////////////////////////
/* CLOSED ORBIT ARRAY FILLING
Set up the grid in xbar and also initialize all arrays that contain a different number at different values of xbar, indexed j */
for (j=0;j<sizexbar;j++)  {
	//printf("j=%d/%d\n", j, sizexbar);
	itop[j] = 0;
	imax[j] = imin[j] = -1;
	openorbit[j] = 0.0;
	crossed_min[j] = 0;
	crossed_max[j] = 0;
	xtop[j] = 0.0;
	chiprimetop[j] = 0.0;
	chimin[j] = 0.0;
	chiMax[j] = 0.0;
	upperlimit[j] = -1;
	lowerlimit[j] = 0; 
}
/* Loop below initializes all 2d arrays which are functions of xbar and x. It allocates the right amount of memory to arrays of pointers of size xbar. The result is a 2D array indexed j (size sizexbar) and i or k (size n, see above) */
for (j=0; j < sizexbar; j++) {
	chi[j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing chi(xbar, x) indexes j and i 
	chiop[j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing chi(xbar, x) indexes j and i 
	Uperp[j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing Uperp for closed orbits indexes j and k
	mu[j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing adiab invariant mu(xbar, Uperp), indexes j and k
	upper[j] = (int*)calloc(size_finegrid,sizeof(int)); // array containing the index of energy Uperp corresponding to chi(x)
	vx[j] = (double**)calloc(size_finegrid,sizeof(double*)); // see below
	/* The loop below initializes the 3d array containing the velocity of a particle at a given orbit position xbar, with particle position x and with energy Uperp, indexed i, j and k. */
	for (i = 0; i < size_finegrid; i++) {
		vx[j][i] = (double*)calloc(size_finegrid,sizeof(double)); 
	} 
} // 3D array containing value of vx at different xbar, x and Uperp
i=0;  // Re-set i to zero.
/* This for loop fills in the arrays, calculating the integrals where necessary */
for (j=0; j<sizexbar; j++) {	
	chi[j][i] =  pow((xx[i] - xbar[j]), 2.0) + phi[i];
}
	// chiop[j][i] =  pow((xx[i] - xbarop[j]), 2.0) + phi[i]; 
for (i=1; i<size_finegrid; i++) {
	for (j=0; j<sizexbar; j++)
	{	
		chi[j][i] =  pow((xx[i] - xbar[j]), 2.0) + phi[i];
//chiop[j][i] =  pow((xx[i] - xbarop[j]), 2.0) + phi[i];
/* 
FINDING MAXIMA/MINIMA
*/
/* Below, we use the array elements to define arrays for the effective potential maxima and minima that exist for every xbar (index j). We search through the chi curve from x=0 to the largest value of x but we search through every chi curve first (so first scan in xbar at fixed x, then move to the next x). We create arrays of mu, Uperp and vx.
Because I am comparing neighbouring values of x to find a maximum, I need to consider two distinct cases.The first one is treated in the if loop below, which the program should enter only if chi is decreasing at x=0. Implying that chi(x=0) is an effective potential maximum. The second one is the else if loop after that, which finds maxima of chi that are stationary points by comparing the value of the function before and after each point. Note that because we compare the function at a point with the function at points before and after, the point we consider at every iteration step in the index i is indexed (i-1), compared with (i-2) and i.  */
		if ( (i==1) && (chi[j][i] < chi[j][i-1]) )
		{	
			//if (phi_cut > chiMax[j]) {
		//		kdrop[j] = (int)((phi_cut - chiMax[j])/0.2); 
		//		printf("ever here?\n");
		//	}
		//	else kdrop[j] = 0;
			crossed_max[j] += 1;		
			imax[j] = 0;
			chiMax[j] = chi[j][i-1];
			Uperp[j][kdrop[j]] = chi[j][0];
			mu[j][0] = 0.0;
			vx[j][0][0] = 0.0; 
			//printf("charge = %f\n", charge);
			//printf("imax[%d]  = %d\n", j, imax[j]);
		} 
		else if ( (i > 1) && ((chi[j][i] < chi[j][i-1]) && (chi[j][i-1] > chi[j][i-2])) )
		{
			//printf("i = %d/%d\tj = %d/%d\n", i, size_finegrid, j, sizexbar);
			crossed_max[j] += 1;
// crossed_min and crossed_max let the code know whether we go up or down the effective potential well. We are starting to go down!
			if (crossed_max[j] > 1)
			{	
				printf("***WARNING*** There is more than one maximum!\n");
				problem = 1;
				printf("j is %d and maxima at %d and %d for first and second\n",  j, imax[j], i-1);
				for (k=imax[j]-1;k<i+1;k++) {
				printf("chi[%d][%d] = %f\n", j, k, chi[j][k]);}
				//if (i-1 - imin[j] > 2) 
				//{ 	exit(1); } 
				crossed_max[j] = 1; 
			}
			imax[j] = i-1; // Store the index of the maximum 
			chiMax[j] = chi[j][i-1]; //Store chi Maximum itself, a function of xbar (index j)
			//if (phi_cut > chiMax[j]) {
			//	kdrop[j] = (int)((phi_cut - chiMax[j])/0.2); 
			//	printf("ever here?\n");
			//}
			//else kdrop[j] = 0;
		}
/* We store the index corresponding to the position of a minimum for a given value of xbar.*/
		//watch out HERE
		else if ( (i > 1) && ((chi[j][i] > chi[j][i-1]) && (chi[j][i-1] < chi[j][i-2])) )
		{	
			crossed_min[j] += 1;
			//printf("i-1 = %d, imax[%d] = %d, crossed_min[%d] = %d, crossed_max[%d] = %d\n", i-1, j,imax[j], j, crossed_min[j], j, crossed_max[j]);
			if (crossed_min[j] > 1)
			{	
				printf("***WARNING*** There is more than one minimum!\n");
				problem = 1;
				printf("j is %d and minima at %d and %d for first and second\n",  j, imin[j], i-1);
				for (k=imin[j]-1;k<i+1;k++) {
				printf("chi[%d][%d] = %f\n", j, k, chi[j][k]);}
				//if (i-1 - imin[j] > 2) 	
				//{	exit(1); } 
				crossed_min[j] = 1; 
			}
			imin[j] = i-1; // Store the index of the maximum 
			upperlimit[j] = imin[j] - imax[j];
			chimin[j] = chi[j][i-1];  //Store chi Maximum itself, a function of xbar (index j)
			crossed_max[j] = 1;
/* We will now start going up the effective potential well! The temporary flag below is used because we still have to store the effective potential minimum as a possible value of Uperp. If we don't have this flag we miss the region near the minimum of chi. */
			temporary = 1; 
		}
/* FILLING IN ARRAYS */
/* Once we cross a maximum, we start filling in arrays for mu, Uperp and vx related to the orbit we are in. As we go down the maximum we store the values of Uperp = chi(x) we encounter which will form a grid of (unevenly spaced) allowed values of Uperp. We also store the value of the small deltamu associated with every value of Uperp above the current value of chi(x), and add it to the previous values */                                                                               	
		if ( (crossed_max[j] == 1  && crossed_min[j] == 0) || (temporary == 1) )
		{
			temporary = 0;
			Uperp[j][i-1-imax[j]+kdrop[j]] = chi[j][i-1];
			if (Uperp[j][i-1-imax[j]+kdrop[j]] < 18.0*Telarge && Uperp[j][i-2-imax[j]+kdrop[j]] > 18.0*Telarge)
				lowerlimit[j] = i-1-imax[j]+kdrop[j]; 
			mu[j][i-1-imax[j]+kdrop[j]] = 0.0;
			//printf("kdrop[j] = %d\n", kdrop[j]);
			upper[j][i-1] = i-1-imax[j]+kdrop[j];
/* Note that the size of the dimension of the array with values of Uperp is set to n, which is larger than the size it will turn out to be. Unfortunately in C I have no way to append elements to arrays as I go along, but need to give enough memory to the array from the start. n is the largest possible size the array could have. */
			for (k=0;k<=upper[j][i-1]; k++) 
			{	
				// replaced k with upper below
				//if (upper[j][i-1] == 0)
				if (k == upper[j][i-1])
				{	
					vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));// equiv 0
					mu[j][k] += 0.0; 
				}
				else if (k== upper[j][i-1] -1) {
					vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));
					//chiprime = ;
					mu[j][k] += (1.0/M_PI)*sqrt(chi[j][i-2]-chi[j][i-1])*(2.0/3.0)*(xx[i-1] - xx[i-2]);
					//printf("%f \n", sqrt(chi[j][i-2]-chi[j][i-1]));
//(1.0/M_PI)*(vx[j][i-1][k] + vx[j][i-2][k])*(xx[i-1] - xx[i-2]); 
				}
				else
				{	
					vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));
					if (debug==1) 
					{ 
						printf("vx[%d][%d][%d] = %f, Uperp = %f, chi = %f\n", j,i-1,k,vx[j][i-1][k], Uperp[j][k], chi[j][i-1]); 
					}
					//if (k==0)  {
						//printf("vx[%d][%d][0] = %f\n", j,i-1, vx[j][i-1][0]);
						//printf("vx[%d][%d][0] = %f\n", j,i-2, vx[j][i-2][0]);
					//}

					mu[j][k] += (1.0/M_PI)*(vx[j][i-1][k] + vx[j][i-2][k])*(xx[i-1] - xx[i-2]); 
				}
			}
		}
/* Once we cross the minimum, we stop creating array elements with values of Uperp. However, we keep storing the value of vx associated with any given point x on an effective potential curve with xbar, with energy Uperp and using this value to finish performing the mu integral. This should happen as long the effective potential at the point under consideration is smaller than the effective potential maximum. */
		//else if (crossed_min[j] == 1 && crossed_max[j] == 1 && chi[j][i-1] < chiMax[j] - TINY)
		else if (crossed_min[j] == 1 && crossed_max[j] == 1 && chi[j][i-1] < chiMax[j] )
		{	
			for (k=0;k <= upperlimit[j] ;k++)
			{	
				if (chi[j][i-1] < Uperp[j][k])
				{	
					vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));
					mu[j][k] += (1.0/M_PI)*(vx[j][i-1][k] + vx[j][i-2][k])*(xx[i-1] - xx[i-2]); 
				}
				// watch out here
				//else if (Uperp[j][k-1] > chi[j][i-1] + TINY && Uperp[j][k] < chi[j][i-1] - TINY)
				else if (Uperp[j][k-1] > chi[j][i-1]  && Uperp[j][k] <= chi[j][i-1] )
				{	
					upper[j][i-1] = k;
					//mu[j][k] += (2.0/M_PI)*(vx[j][i-2][k])*(xx[i-1] - xx[i-2])*(Uperp[j][k] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]); 
					ind = 0;
					while (Uperp[j][k] < chi[j][i-2-ind]) {
						ind++;
					}
					mu[j][k] += (1.0/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[j][k] - chi[j][i-2-ind], 1.5)/(chi[j][i-1-ind] - chi[j][i-2-ind]); 
					//printf("YO %f: i = %d, j = %d, k = %d \nchi = %f, %f\n", (xx[i-1] - xx[i-2])*pow(Uperp[j][k] - chi[j][i-2], 1.5)/(chi[j][i-1] - chi[j][i-2]), i, j, k, chi[j][i-1], chi[j][i-2]);
				}
				if (mu[j][k] != mu[j][k])
				{	if (debug == 1)
					{	
						printf("mu is NAN\n"); 
					}
					problem = 1;  
				}  
			}
		}
/* When the effective potential at the iteration point (i-1) under consideration becomes larger than the effective potential maximum, we finish performing the mu integral. We also store the position of the top of the orbit which has chi = chiMax, in order to perform the open orbit integral. If the loop below is accessed, a switch it turned off to signify the no more closed orbits can be present */
		//else if ( ( crossed_min[j] == 1 ) && ( crossed_max[j] == 1) && ( chi[j][i-1] > chiMax[j] - TINY) ) {	
		else if ( ( crossed_min[j] == 1 ) && ( crossed_max[j] == 1) && ( chi[j][i-1] > chiMax[j] ) ) {	
			itop[j] = i-2;
			xtop[j] = xx[i-2] + ((chiMax[j] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]))*(xx[i-1] - xx[i-2]);
			chiprimetop[j] = (chi[j][i-1] - chi[j][i-2])/(xx[i-1] - xx[i-2]);
			for (k=0; k<upper[j][i-2]; k++)
			{
				mu[j][k] += (1.0/M_PI)*(vx[j][i-2][k])*(xx[i-1] - xx[i-2])*(Uperp[j][k] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]);
				if (mu[j][k] != mu[j][k])  printf("mu is NAN\n");
			}
			crossed_max[j] = 0; }
		if (j!=0)
		{	
			//if ( ((chiMax[j-1] < chi[j-1][i-1] + TINY) && (chiMax[j] + TINY > chi[j][i-1])) ) { // || (imax[j-1] == -1 && imax[j] != -1) )
			if ( ((chiMax[j-1] < chi[j-1][i-1] ) && (chiMax[j] > chi[j][i-1])) ) { // || (imax[j-1] == -1 && imax[j] != -1) )
				jmclosed[i-1] = j-1;
				if (i - 1 > icrit) { 	
					jmopen[i-1] = j-1; 
					if (debug == 1) {
						printf("i = %d, j = %d, icrit = %d\n", i-1, j-1, icrit); 
						printf("jmopen[%d] = %d\n", i-1, jmopen[i-1]); 
					}
				} 
			}
			if (i==size_finegrid-1 && crossed_max[j] == 1 && crossed_min[j] == 1) {	
				itop[j] = size_finegrid-1;
				if (noback == 0) {
					maxj = j;	
					noback = 1; 
				} 
			} 
		} 
	} 
}	
printf("xbar[maxj=%d] = %f\n", maxj, xbar[maxj]);	
// OPEN ORBIT INTEGRAL
/* Now we perform the open orbit integral. We use a change of variables which makes the integrand smooth at the top bounce point. The change of variables is to some var = sqrt(x_t - x) */
muopen[0] = 0.0;
for (j=0;j<maxj;j++) {
	muopen[j+1] = mu[j][0];
	Ucritf[j+1] = chiMax[j] - mu[j][0] ;
	if (debug == 1)
	{	
		printf("jmopen[%d] = %d\n", j, jmopen[j]); 
	}
	//METHOD 1
	if (j != 0)
	{	
		//twopimuprime[j-1] = 2.0*M_PI*(mu[j][0] - mu[j-1][0])/(xbar[j] - xbar[j-1]);	
		twopimuprime[j] = M_PI*(mu[j+1][0] - mu[j][0])/(xbar[j+1] - xbar[j]) + M_PI*(mu[j][0] - mu[j-1][0])/(xbar[j] - xbar[j-1]);	
	
	}
	else
	{
		//twopimuprime[j] = 2.0*M_PI*(mu[j][0] - 0.0)/(xbar[j+1] - xbar[j]);	
		twopimuprime[j] = M_PI*(mu[j+1][0] - mu[j][0])/(xbar[j+1] - xbar[j]);	
	}
// METHOD 2
	//if (imax[j] != 0) {
	//	openintegral[j] += 4.0*(xx[imax[j]+1] - xx[imax[j]])*(xx[imax[j]+1] - xx[imax[j]])/vx[j][i][0];
	//}
	//else if (imax[j] == 0) {
	//	openintegral[j] += 2.0*(xx[imax[j]+1] - xx[imax[j]])*(xx[imax[j]+1] - xx[imax[j]])/vx[j][i][0];
	//}
	//aaold = 8.0*sqrt(xtop[j] - xx[imax[j]])/sqrt(-chiMpp[j]);
	//aaold = 0.0;
	for (i=imax[j]+2;i<itop[j];i++)
	{
		openintegral[j] += 2.0*( (xx[i] - xx[imax[j]])/vx[j][i][0] + (xx[i-1] - xx[imax[j]])/vx[j][i-1][0] )*(xx[i] - xx[i-1]);
		//aa = 8.0*(xx[i] - xx[imax[j]])*sqrt(xtop[j] - xx[i])/vx[j][i][0];
		//if (fabs(vx[j][i][0]) < 0.0001)
		//	printf("HERE INF, j=%d, i=%d\n", j, i); 
		//if (aa != aa)
		//{
		//	printf("aa is %f\ni is %d, j is %d\nitop is %d\n", aa, i, j, itop[j]);
		//	printf("xtop[%d] = %f, xx[%d] = %f, xx[imax=%d]=%f\n", j, xtop[j], i, xx[i], imax[j], xx[imax[j]]);
		//	printf("vx = %f\n", vx[j][i][0]);
		//	openintegral[j] += 0.0;
		//}
		//else
		//{
		//	openintegral[j] += 0.5*(aa + aaold)*(sqrt(xtop[j] - xx[i-1]) - sqrt(xtop[j] - xx[i]));
		//	aaold = aa;
		//}
	}
	//aa = 8.0*(xtop[j] - xx[imax[j]])/sqrt(chiprimetop[j]);
	//if (aa != aa)
	//{
	//		printf("aa= %f, chiprimetop = %f, j = %d\n", aa, chiprimetop[j], j);
	//		aa = 0.0;
	//}
	//openintegral[j] += 0.5*(aa + aaold)*sqrt(xtop[j] - xx[itop[j]-1]);
	//openorbit[j] = openintegral[j];
	openorbit[j] = twopimuprime[j];  
	openorbitantycal = 4.0*M_PI*xbar[j];
	if (debug ==1)
	{	
		printf("twopimuprime[%d] = %f, openintegral[%d] = %f\n", j, twopimuprime[j], j, openintegral[j]); 
		printf("At xbar[%d] = %f, mu is %f (should be %f if phi is flat), openorbit is %f (should be %f if phi is flat)\n", j, xbar[j], mu[j][0], Uperp[j][0], openorbit[j], openorbitantycal); 
	} 
}
// temporary
muopen[maxj] = 1000.0;
//
Ucritf[0] = Ucritf[1] - muopen[1]*(Ucritf[2] - Ucritf[1])/(muopen[2] - muopen[1]);
//Ucritf[0] = phi[0];
if (debug == 1) {
	printf("~~~~~The second element of FF is %f~~~~~\n",FF[0][1]);
	printf("~~~~~The second element of UU is %f~~~~~\n",UU[1]);
	printf("~~~~~The fourth element of mu is %f~~~~~\n",mumu[3]);
}

i=0;
clock_t int1 = clock(); // Finds the time of the computation so far
double inttime  = (double)(int1 - begin) / CLOCKS_PER_SEC;
if (debug == 1) 
	printf("in densfinorb: Array filling DONE: time is %f\n", inttime);
/* DENSITY INTEGRALS 
This part calculates the density integrals and outputs the result of the integration to a file fout and also the yz distribution function to three files one containing the distribution function the other two containing the velocity grid */
FILE *fout; 
if ((fout = fopen("PostProcessing/densfinorb_out.txt", "w")) == NULL)
{	
	printf("Cannot open densfinorb_out.txt");
	exit(EXIT_FAILURE);
}
i=0;
stop = 0;
ic = 0;
while (stop == 0) 
{	
	i = ic*zoomfactor;
	pos = xx[i];
	intdxbar = 0.0;
	intdxbaropen = 0.0;
	intdxbaropenflow = 0.0;
	for (j=0; j<sizexbar; j++) {
		vxnew = 0.0;
		intdUopenold = intdUopen;
		intdUopenflowold = intdUopenflow; intdUopenBohmold = intdUopenBohm;
		intdUopenxbarold = intdUopenxbar; intdUopensquareold = intdUopensquare;
		intdUopen = 0.0;
		intdUopenflow = 0.0; intdUopenBohm = 0.0;
		intdUopenxbar = 0.0; intdUopensquare = 0.0;
		if (j == jmopen[i]) {	
			oorbintgrd = oorbintgrdold = 0.0;
			oorbintgrdflow = oorbintgrdflowold = 0.0; oorbintgrdBohm = oorbintgrdBohmold = 0.0;
			oorbintgrdxbar = oorbintgrdxbarold = 0.0; oorbintgrdsquare = oorbintgrdsquareold = 0.0;
			sizeU = (int) sqrt(Ucap - chiMax[j])/dvzopen;
			for (l=0; l < sizeU; l++) {	
				oorbintgrdold = oorbintgrd;
				oorbintgrdflowold = oorbintgrdflow; oorbintgrdBohmold = oorbintgrdBohm;
				oorbintgrdxbarold = oorbintgrdxbar; oorbintgrdsquareold = oorbintgrdsquare;
				vz = dvzopen*l;
				U = chi[j][k] + pow(vz, 2.0); 
				if ( (U > mu[j][0]) && (U - vz*vz + (vz-dvzopen)*(vz-dvzopen) < mu[j][0]) ) {
					frac = (vz - sqrt(mu[j][0] - chi[j][k]))/dvzopen;
					Fopen = bilin_interp(mu[j][0], 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					//printf("frac = %f\n", frac);
					oorbintgrdold = sqrt(alpha*sqrt(mu[j][0] - chi[j][k])*openorbit[j])*Fopen;
					oorbintgrdflowold = 0.5*alpha*sqrt(mu[j][0] - chi[j][k])*openorbit[j]*Fopen;
					Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					oorbintgrd = sqrt(alpha*vz*openorbit[j])*Fopen;
					oorbintgrdflow = 0.5*alpha*vz*openorbit[j]*Fopen;
					//if (oorbintgrd != oorbintgrd) 
					//	printf("AHA\n\n\n\n\n");
				}
				else {
					frac = 1.0;
					Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					oorbintgrd = sqrt(alpha*vz*openorbit[j])*Fopen;
					oorbintgrdflow = 0.5*alpha*vz*openorbit[j]*Fopen;
					//if (oorbintgrd != oorbintgrd) 
					//	printf("AHA\n\n\n\n\n");
				}

					
				// is it worth increasing accuracy?
				//if (i > icrit)
				//	Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
				//else
				//	Fopen = 0.0; 	

				if (oorbintgrd != oorbintgrd)
				{	
					problem = 1; 
					if (debug == 1)
					{
						printf("vx0open = %f, openorbit[%d] = %f, imaginary oorbintgrd in initial piece of integral due to negative value of openorbit?\n", vx0open, j, openorbit[j]); 
					}
				}
				if (i==0) 
					oorbintgrdantycal = sqrt(2.0*alpha*vz*M_PI*xbar[j])*(U-chiMax[j])*exp(-U)/pow(M_PI, 1.5); 
// for a flat potential, oorbintgrd can be calculated analytically
					//oorbintgrd = oorbintgrdantycal;
					//printf("oorbintgrd is %f, analytical is %F\n", oorbintgrd, oorbintgrdantycal);
				if (l!=0) {	
					intdUopen += 2.0*frac*dvzopen*(oorbintgrd+oorbintgrdold);
					intdUopenflow += 2.0*frac*dvzopen*(oorbintgrdflow+oorbintgrdflowold);
					intdUopenBohm += 2.0*frac*dvzopen*(oorbintgrdBohm+oorbintgrdBohmold);
					intdUopenxbar += 2.0*frac*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
					intdUopensquare += 2.0*frac*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
				}
				else {	
					intdUopen += 0.0;
					intdUopenflow += 0.0;
					intdUopenBohm += 0.0;
					intdUopenxbar += 2.0*frac*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
					intdUopensquare += 2.0*frac*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
				} 
			} 
			intdxbaropen += 0.5*(intdUopen+intdUopenold)*dxbar;
			intdxbaropenflow += 0.5*(intdUopenflow+intdUopenflowold)*dxbar;
			intdxbaropenBohm += 0.5*(intdUopenBohm+intdUopenBohmold)*dxbar;
			intdxbarxbar += 0.5*(intdUopenxbar+intdUopenxbarold)*dxbar;
			intdxbarsquare += 0.5*(intdUopensquare+intdUopensquareold)*dxbar; 
		}
		if (j > jmopen[i]) {	
			oorbintgrd = oorbintgrdold = 0.0;
			oorbintgrdflow = oorbintgrdflowold = 0.0;
			oorbintgrdBohm = oorbintgrdBohmold = 0.0;
			oorbintgrdxbar = oorbintgrdxbarold = 0.0;
			oorbintgrdsquare = oorbintgrdsquareold = 0.0;
			sizeU = (int) sqrt(Ucap - chiMax[j])/dvzopen;
			for (l=0; l < sizeU; l++) {
				oorbintgrdold = oorbintgrd;
				oorbintgrdflowold = oorbintgrdflow;
				oorbintgrdBohmold = oorbintgrdBohm;
				oorbintgrdxbarold = oorbintgrdxbar;
				oorbintgrdsquareold = oorbintgrdsquare;
				vz = dvzopen*l;
				U = chiMax[j] + pow(vz, 2.0);
				//vx0open = sqrt(TINY + chiMax[j] - chi[j][i]);
				vx0open = sqrt(chiMax[j] - chi[j][i]);
				if (vx0open != vx0open) {		
					problem = 1; 
					if (debug == 1)
						printf("HERE imaginary vx0open, j = %d, i is %d, chi[j][i] = %f, chiMax[j] = %f\n", j, i, chi[j][i], chiMax[j]); 
				}
				if ( (U >= mu[j][0]) && (U - vz*vz + (vz-dvzopen)*(vz-dvzopen) < mu[j][0]) ) {
					frac = (vz - sqrt(mu[j][0] - chiMax[j]))/dvzopen;
					Fopen = bilin_interp(mu[j][0], 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					oorbintgrdold = sqrt(alpha*sqrt(mu[j][0] - chiMax[j])*openorbit[j])*Fopen;
					oorbintgrdflowold = 0.5*alpha*sqrt(mu[j][0] - chiMax[j])*openorbit[j]*Fopen;
					oorbintgrdBohmold = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*sqrt(mu[j][0] - chiMax[j])*openorbit[j]))*Fopen;
					oorbintgrdxbarold = xbar[j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*sqrt(mu[j][0] - chiMax[j])*openorbit[j]))*Fopen;
					oorbintgrdsquareold = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + alpha*sqrt(mu[j][0] - chiMax[j])*openorbit[j]), 1.5))*Fopen;
					Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					oorbintgrd = sqrt(alpha*vz*openorbit[j])*Fopen;
					oorbintgrdflow = 0.5*alpha*vz*openorbit[j]*Fopen;
					oorbintgrdBohm = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*vz*openorbit[j]))*Fopen;
					oorbintgrdxbar = xbar[j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*vz*openorbit[j]))*Fopen;
					oorbintgrdsquare = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + alpha*vz*openorbit[j]), 1.5))*Fopen;
					//if (oorbintgrd != oorbintgrd) 
					if (frac != frac) 
						printf("AHA\n\n\n\n\n");

				}
				else {
					frac = 1.0;
					Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					//if (i ==0) printf("l = %d, Fopen = %f\n", l, Fopen);
					oorbintgrd = (sqrt(vx0open*vx0open + alpha*vz*openorbit[j]) - vx0open)*Fopen;
					oorbintgrdflow = 0.5*alpha*vz*openorbit[j]*Fopen;
					oorbintgrdBohm = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*vz*openorbit[j]))*Fopen;
					oorbintgrdxbar = xbar[j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*vz*openorbit[j]))*Fopen;
					oorbintgrdsquare = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + alpha*vz*openorbit[j]), 1.5))*Fopen;
					//if (oorbintgrd != oorbintgrd) 
					if (frac != frac) 
						printf("AHA\n\n\n\n\n");
				}
				//Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1);
				if (oorbintgrd != oorbintgrd) {
					problem = 1; 	
					if (debug == 1) {
						printf("vx0open = %f, openorbit[%d] = %f, HERE imaginary oorbintgrd in initial piece of integral due to negative value of openorbit?\n", vx0open, j, openorbit[j]); 
					} 
				}
				if (i==0) {
					oorbintgrdantycal = sqrt(2.0*alpha*vz*M_PI*xbar[j])*(U-chiMax[j])*exp(-U)/pow(M_PI, 1.5); 
					if (debug == 1) 
					printf("openorbit = %f (should be %f with flat potential profile)\n", oorbintgrd, oorbintgrdantycal);
					// check not working
				}
// for a flat potential, oorbintgrd can be calculated analytically
				//oorbintgrd = oorbintgrdantycal;
//printf("oorbintgrd is %f, analytical is %F\n", oorbintgrd, oorbintgrdantycal);
				if (l!=0) {
					//if (oorbintgrdold < TINY ) {
					//	intdUopen += 4.0*(2.0/3.0) *pow(dvzopen, 1.5) * sqrt(alpha*openorbit[j])*bilin_interp(mu[j][0], 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
;
					//}
					
					intdUopen += 2.0*frac*dvzopen*(oorbintgrd+oorbintgrdold);
					intdUopenflow += 2.0*frac*dvzopen*(oorbintgrdflow+oorbintgrdflowold);
					intdUopenBohm += 2.0*frac*dvzopen*(oorbintgrdBohm+oorbintgrdBohmold);
					intdUopenxbar += 2.0*frac*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
					intdUopensquare += 2.0*frac*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
				}
				else {
					intdUopen += 0.0;
					intdUopenflow += 0.0;
					intdUopenBohm += 0.0;
					intdUopenxbar += 2.0*frac*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
					intdUopensquare += 2.0*frac*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
				} 
			}
			if (intdUopen != intdUopen)  {
				problem = 1; 
				if (debug == 1)
					printf("HERE, j is %d\n", j); 
			}
			if (i==0) {	
				intdUopenantycal = 2.0*0.919*sqrt(2.0*alpha*xbar[j])*exp(-xbar[j]*xbar[j])/M_PI; 
			if (debug == 1) { 
				printf("xbar[%d] = %f, analytical is %f, numerical is %f\n", j, xbar[j], intdUopenantycal, intdUopen); printf("vx0open is %f (should be zero)\n", vx0open); 
			}
			} 
			dxbar = xbar[j] - xbar[j-1];
			if ( (j == jmopen[i]+1) && (i > icrit) ) {
				if ( fabs(pos - xx[imax[j]]) < TINY )
					intdxbaropen += 0.0; 
				else { 
					dxbar = xbar[j] - (xbar[j-1]*(chiMax[j] - chi[j][i]) - xbar[j]*(chiMax[j-1] - chi[j-1][i]))/(-chiMax[j-1] + chi[j-1][i] + chiMax[j] - chi[j][i]); 
					if (debug == 1)
						printf("dxbar = %f\n", dxbar);
				}
			}
//(chi[j][k]-chi[j][i])/(2.0*(pos-xx[k])); 
			intdxbaropen += 0.5*(intdUopen+intdUopenold)*dxbar;
			intdxbaropenflow += 0.5*(intdUopenflow+intdUopenflowold)*dxbar;
			intdxbaropenBohm += 0.5*(intdUopenBohm+intdUopenBohmold)*dxbar;
			intdxbarxbar += 0.5*(intdUopenxbar+intdUopenxbarold)*dxbar;
			intdxbarsquare += 0.5*(intdUopensquare+intdUopensquareold)*dxbar; 
			if (intdxbaropen != intdxbaropen) {
				problem = 1; 
				if (debug == 1)
					printf("PROBLEM HERE, j is %d\n", j); 
			} 
		}
		if ( j>=jmclosed[i] )  /* We have entered the closed orbit integral */
		{	
			intdvxold = intdvx;
			//intdvxflowold = intdvxflow;
			//intdvxflow = 0.0;
			intdvx = 0.0;
			if (j == jmclosed[i]) {
				intdvx = 0.0;
				//intdxbar += 0.0; 
			}
			else if (j > jmclosed[i]) {	
				if (debug ==1) {
					printf("lowerlimit[%d] = %d\tupper[%d][%d] = %d\n", j, lowerlimit[j], j, i, upper[j][i]);
					printf("imax[%d] = %d\n", j, imax[j]);
				}
				for (k=lowerlimit[j]; k<upper[j][i]+1; k++) {	
					vxold = vxnew;
					intdUold = intdU;
					intdU = 0.0;
					if (k == upper[j][i]) {
						Uperpnew = chi[j][i];
						if (lowerlimit[j] == upper[j][i] )
							munew = 0.0; 
						else
							munew = ((chi[j][i] - Uperp[j][k])*mu[j][k-1] + (Uperp[j][k-1] - chi[j][i])*mu[j][k])/(Uperp[j][k-1] - Uperp[j][k]);
						if (munew != munew)
						{ 	
							problem = 1; 
							if (debug == 1)
							{	
								printf("munew is NAN, j = %d, i = %d, k = %d, lowerlimit = %d, upper = %d\n", j, i, k, lowerlimit[j], upper[j][i]); 
								printf("mu[%d][%d] = %f, mu[%d][%d] = %f\n", j, k-1, mu[j][k-1], j, k, mu[j][k]); 
							} 
						}
						vxnew = 0.0;
						dvx = vxold - vxnew; 
					}
					else
					{	
						Uperpnew = Uperp[j][k];
						vxnew = vx[j][i][k];
						dvx = vxold - vxnew;
						munew = mu[j][k]; 
					}
					Ucrit = lin_interp(muopen, Ucritf, munew, sizexbar, 791);
					sizeU = (int) sqrt(Ucap - Uperp[j][k])/dvz;
					reflected = 0;
					for (l=0; l < sizeU; l++)
					{	
						if (l!=0)
						{	
							Fold = F;
							vz = dvz*l;
							U = Uperpnew + pow(vz, 2.0);
							if ( (U > munew) && (U - vz*vz + (vz-dvz)*(vz-dvz) < munew) ) {
								frac = (vz - sqrt(munew - Uperpnew))/dvz;
								//printf("frac = %f\n", frac);
								Fold = bilin_interp(munew, 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
								F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
							}
							else {
								frac = 1.0;
								F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
							}
							//U = Uperpnew + pow(vz, 2.0);
							if ( (U-munew < Ucrit) && (reflected == 1) ) {
								frac_reflected = 1.0;
							}
							else if ( (U-munew > Ucrit) && (reflected == 1) ) {
								reflected = 0;
								vzcrit = sqrt(Ucrit + munew - Uperpnew);
								frac_reflected = (vzcrit - (vz - dvz))/dvz;
							}
							else frac_reflected = 0.0;
									
							//if (reflected == 1) printf("particles are being reflected\n");
							//F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
							//printf("frac_reflected = %f\n", frac_reflected);
							//printf("F = %f\n", F);
							intdU += (1.0 + frac_reflected)*frac*dvz*(F+Fold);
						}
						else {	
							//vz = dvz*l;
							U = Uperpnew ;//+ pow(vz, 2.0);
							if (U-munew < Ucrit) {
								reflected = 1;
							}
							F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
							intdU += 0.0;
						} 
					}
					intdUantycal = exp(-Uperp[j][k])*(1.0/(2.0*M_PI));// result with phi =0
					//intdU = intdUantycal;
					if (debug == 1) 	
						printf("Analytical intdU is %f, numerical one is %f\n", intdUantycal, intdU);
					if (k==lowerlimit[j])
					{	intdvx += 0.0; }
					else
					{	
						intdvx += 4.0*0.5*dvx*(intdU+intdUold);
					}
					if (intdvx != intdvx)
					{	problem = 1; 
						if (debug ==1)
							printf("intdvx is NAN, j=%d, i=%d\n", j, i); 
					} 
				}
				intdvxantycal = (2.0/(2.0*sqrt(M_PI)))*exp(-(pos-xbar[j])*(pos-xbar[j]))*erf(sqrt(xbar[j]*xbar[j]-(pos-xbar[j])*(pos-xbar[j])));
				if (debug == 1) {
					printf("pos=%f, i=%d, j=%d, Uperp is %f, chi is %f\nvxnew and vxold are %f and %f and and vx is %f, dvx is %f\nintdvx is %f, analytical one is %f, upper is %d, upperlimit is %d\n", pos, i, j, Uperpnew, chi[j][i], vxnew, vxold, vx[j][i][k], dvx, intdvx, intdvxantycal, upper[j][i], upperlimit[j]); 
				}
				//printf("intdvx is %f, while the analytical one is %f\n", intdvx, intdvxantycal);
				//intdvx = intdvxantycal;
				if (j==jmclosed[i]+1 && i!=0)
				{
					//dxbar = (xx[imax[j]] - xx[imax[j]])*(xbar[j] - xbar[j-1])/(xx[imax[j]] - pos);
				//	dxbar = (chiMax[j]-chi[j][i])/(2.0*(pos-xx[imax[j]]));
					if (i==imax[j])
					{	dxbar = xbar[j] - xbar[j-1]; }
					else
					{	dxbar = xbar[j] - (xbar[j-1]*(chiMax[j] - chi[j][i]) - xbar[j]*(chiMax[j-1] - chi[j-1][i]))/(-chiMax[j-1] + chi[j-1][i] + chiMax[j] - chi[j][i]); }
					intdxbar += 0.5*(intdvx+intdvxold)*dxbar;
					if (intdxbar != intdxbar) {	
						problem = 1; 
						if (debug == 1)
						{	
							printf("intdxbar is NAN, j=%d, i=%d\n", j, i);  
						} 
					} 
				}
				else
				{
					dxbar = xbar[j] - xbar[j-1];
					intdxbar += 0.5*(intdvx+intdvxold)*dxbar;
					if (intdxbar != intdxbar)
					{
						problem = 1; 
						if (debug == 1)
						{
							printf("intdxbar is NAN, j=%d, i=%d\n", j, i); 
						}
					} 
				}
			} 
		} 
	}

	ne[ic] = exp(phi[i]/Te); 
	niclosed[ic] = intdxbar;
	niopen[ic] = intdxbaropen;
	nitotal[ic] = niopen[ic] + niclosed[ic];
	if (ic == 0) {
		Bohm = intdxbaropenBohm/nitotal[ic];
		*flux = intdxbaropenflow/alpha; 
		printf("in densfinorb.c: flow velocity at x=0 = %f\n", (*flux)/nitotal[ic]);
		printf("intdxbarsquare = %f\tintdxbarxbar = %f\n", intdxbarsquare, intdxbarxbar);
		kBohmsquared = intdxbarxbar/(intdxbarsquare - 0.5*nitotal[0]); 
	}
	if ( (charge > 0.0) && ((nitotal[ic] >= stoplimit) || (ic == size_phigrid-1)) ) {	
		stop = 1;
		*size_ngrid = ic; 
		printf("in densfinorb.c: stopping for positive ion density");
	}
	if ( (ic != 0) && (charge < 0.0) ) {
		if ( (nitotal[ic] > n_inf ) ) {// || (ic == size_phigrid-1) ) 
			//stop = 1;
			n_inf = nitotal[ic];
			*size_ngrid = ic+1; 
			//stop = 1;
			//printf("in densfinorb.c: stopping for electron density with n_inf = %f\n", n_inf);
		}
		if (ic == size_phigrid-1 - 20) {
			stop = 1;
			//n_inf = nitotal[ic-1];
		}
	}
	if (debug == 1)
	{	
		printf("%f, %f, %f is TOTAL, CLOSED and OPEN orbit density at position index %d, position %f\n", nitotal[ic], niclosed[ic], niopen[ic], ic, pos);
		if ( (nitotal[ic] != nitotal[ic]) || (nitotal[ic] < TINY) )
		{ 	
			problem = 1; 
			printf("nitotal[ic] = %f\n", nitotal[ic]);
		} 
	} 
	fprintf(fout, "%f %f %f %f\n", pos, nitotal[ic], niclosed[ic], niopen[ic]);
	//i += 1; 
	ic += 1;
} 
if (charge < 0.0) {
	printf("ninf = %f\n", n_inf);
	//for (ic = 0; ic < *size_ngrid; ic++) {
	//	nitotal[ic] /= n_inf;
	//}
}
printf("size_ngrid = %d\n", *size_ngrid);
fclose(fout);
clock_t int2 = clock(); // Finds the start time of the computation
double inttime1  = (double)(int2 - int1) / CLOCKS_PER_SEC;
if (debug == 1) 
	printf("in densfinorb: After density evaluation time is %f\n", inttime1);

if (problem == 1)
{	
	printf("***WARNING*** there is a problem with the iteration\n"); 
	// If a problem is found, re-start iteration from a simple first guess
}

if (sizevxopen > 0)
{
	dvxopen = (5.0/sizevxopen)*sqrt(Telarge);
	Ucap = 18.0 + 2.0*Te;
	/* The part below performs the integration over the open orbit distribution function at x=0 in a way that allows to extract the distribution function in v_x (because v_x is the only velocity component that matters in the Debye sheath).*/
	FILE *output, *outputyz, *outputyzvy, *outputyzvz;
	if ((output = fopen("PostProcessing/f0x.txt", "w")) == NULL)
	{        // Check for presence of file
		printf("Cannot open %s\n", "outputfile.txt");
		exit(EXIT_FAILURE);
	}
	if ((outputyz = fopen("PostProcessing/f0yz.txt", "w")) == NULL)
	{
		// Check for presence of file
		printf("Cannot open %s\n", "outputyz.txt");
		exit(EXIT_FAILURE);
	}
	if ((outputyzvy = fopen("PostProcessing/vy0.txt", "w")) == NULL)
	{
		// Check for presence of file
		printf("Cannot open %s\n", "outputyzvy.txt");
		exit(EXIT_FAILURE);
	}
	if ((outputyzvz = fopen("PostProcessing/vz0.txt", "w")) == NULL)
	{
		// Check for presence of file
		printf("Cannot open %s\n", "outputyzvz.txt");
		exit(EXIT_FAILURE);
	}
	intdxbaropen = 0.0;
	intdvxopen = 0.0;
	intdUopen = 0.0;
	F = 0.0;
	dvzopen = 0.1;
	sizeU = sqrt(Ucap)/dvzopen;
	if (debug == 1) printf("maxj = %d, sizeU = %d\n", maxj, sizeU);
	for (i=0; i<sizevxopen; i++) {	
		vxopen = i*dvxopen;
		intdxbaropenold = intdxbaropen;
		intdxbaropen = 0.0;
		for (j=0; j < maxj; j++)
		{	
			if (i==0)
			{
				fprintf(outputyzvy, "%f\n", xbar[j]);
			}
			intdUopenold = intdUopen;	
			intdUopen = 0.0;
			F = 0.0;
			//sizeU = (int) ( sqrt(Ucap - pow(vxopen, 2.0) - chi[j][0])/dvzopen );
			for (l=0; l<sizeU; l++)
			{	
				vz = dvzopen*l;
				if (j == 0 && i == 0)
				{	fprintf(outputyzvz, "%f\n", vz); }
				Fold = F;
				vx0open = sqrt(chiMax[j] - chi[j][0]);
				//printf("vx0open = %f\n", vx0open);
				U = pow(vz, 2.0) + pow(vx0open, 2.0) + chi[j][0];
				Deltavx = sqrt(pow(vx0open, 2.0) + alpha*vz*openorbit[j]) -  vx0open;
				aa = tophat(vx0open, vx0open + Deltavx, vxopen);
				//Fopen = bilin_interp(mu[j][0], U, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
				Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1);
				F = Fopen*aa;
				if (i==0)
				{
					if (l==sizeU - 1 && i == 0)
					{	
						fprintf(outputyz, "%f\n", Fopen*Deltavx); 
					}
					else 
					{	
						fprintf(outputyz, "%f ", Fopen*Deltavx); 
					}
				}
				if (debug == 1)
					printf("tophat = %f when vx0open = %f, Deltavx = %f and vxopen = %f\n, F = %f\n", aa, vx0open, Deltavx, vxopen, F);
				if (l!=0) {	
					intdUopen += 2.0*(F+Fold)*dvzopen; 
				}
			}
			if (j != 0) {	
				intdxbaropen += 0.5*(intdUopen + intdUopenold)*(xbar[j] - xbar[j-1]); 
			} 
		}
		if (debug == 1)
			printf("FBohm = %f\n", intdxbaropen);
		if (i!=0)
		{
			intdvxopen += 0.5*(intdxbaropen + intdxbaropenold)*dvxopen;
			intdvxopenfluidBohm += 0.5*(intdxbaropen*vxopen + intdxbaropenold*(vxopen-dvxopen))*dvxopen;
			if (i>1)
			{
			intdvxopenBohm += 0.5*(intdxbaropen/pow(vxopen, 2.0) + intdxbaropenold/pow(vxopen-dvxopen, 2.0))*dvxopen;
			//intdvxopenBohm += 0.5*(intdxbaropen + intdxbaropenold)*(1/(vxopen-dvxopen) - 1/vxopen );
			}
			else if (i==1)
			{
			intdvxopenBohm += 0.5*intdxbaropen/pow(vxopen, 2.0)*dvxopen;
			}
		}
		vx_i_DS[i] = vxopen;
		fx_i_DS[i] = intdxbaropen/nitotal[0];
		fprintf(output, "%f %f\n", vxopen, intdxbaropen/nitotal[0]);
	}
	fclose(output);
	fclose(outputyz);
	fclose(outputyzvy);
	fclose(outputyzvz);

	// Calculating correction due to non boltzmann electrons
	// AG: to the Bohm condition
	/*the gradient of the function */
	//double ne_p0=0.0;
	//ne_p0 = (ne_grid[1] - ne_grid[0]) / (phi[1] - phi[0]);
	//printf("ne_p0 = %f (%f)\n", ne_p0, (ne_grid[1] - ne_grid[0]) / (phi_grid[1] - phi_grid[0]));

	printf("in densfinorb: The density at x=0 obtained from the extracted distribution function is %f\nThe density at x=0 obtained from the ion density integral (more direct) is %f\n", intdvxopen, nitotal[0]);
	//printf("Bohm integral = %f (obtained from the extracted distribution function at x=0)\nBohm integral = %f (obtained directly from the distribution function at infinity)\nThe Bohm condition is: Bohm integral = %f\n", intdvxopenBohm/intdvxopen, Bohm, (2.0*ne_p0)/(Te*nitotal[0])); 
	printf("in densfinorb: The flow at x=0 obtained from the extracted distribution function is %f\nThe flow at x=0 obtained from the distribution function at infinity is %f\nThe Bohm speed is %f\n", intdvxopenfluidBohm/intdvxopen, *flux, 1.0/sqrt(2.0)); 
	printf("in densfinorb: kBohmsquared= %f (should be >0)\n", kBohmsquared);

}
else if (sizevxopen < 0) {
	FILE *output;
	if ((output = fopen("PostProcessing/f0x.txt", "w")) == NULL)
	{       
		printf("Cannot open %s\n", "outputfile.txt");
		exit(EXIT_FAILURE);
	}
	printf("muopen Ucritf\n");
	for (j=0;j<(-sizevxopen);j++) {
		if (debug == 0) {
			printf("%f %f\n", muopen[j], Ucritf[j]); 
		}
		//fprintf(output, "%f %f\n", vxopen, intdxbaropen/nitotal[0]);
		if (muopen[maxj] < vx_i_DS[j]) {
			fx_i_DS[j] = sqrt(2.0*Ucritf[maxj]);
			//printf("Here? vparcut = %f\n", fx_i_DS[j]);
		}
		else fx_i_DS[j] = sqrt(2.0*lin_interp(muopen, Ucritf, vx_i_DS[j], maxj, 986));
		fprintf(output, "%f %f\n", vx_i_DS[j], fx_i_DS[j]);
		printf("j=%d\tmu_crit = %f\tvpar_crit = %f\n", j, vx_i_DS[j], fx_i_DS[j]);
	}
	fclose(output);
}
// If you love your variables (and your memory) set them free
free(ne);//
free(xx);//
free(dxdf);//
free(niclosed);//
free(niopen);//
free(jmclosed);//
free(jmopen);//
free(xifunction);//
free(xbar);//
free(chimin);//
free(chiMax);//
free(chiMpp);//
free(crossed_min);//
free(crossed_max);//
free(twopimuprime);
free(openintegral);
free(openorbit);
free(upperlimit);
free(lowerlimit);
free(chiprimetop);
free(xtop);
free(itop);
free(imax);
free(imin);


for (int w = 0; w < sizexbar; w++)
{
	free(chi[w]);
	free(chiop[w]);
	free(Uperp[w]);
	free(mu[w]);
	free(upper[w]);//
	for (int s = 0; s < size_finegrid; s++)
	{
		free(vx[w][s]);//
	}
	free(vx[w]);//
}
free(chi);//
free(chiop);//
free(Uperp);//
free(mu);
free(upper);
free(vx);

for (i=0; i<*size_ngrid; i++) {
	//printf("iondens.c: phi = %f\tnitotal = %f\n", phi[i], nitotal[i]);
	n_grid[i] = nitotal[i];
}
free(phi);//
free(gg);//
free(ff);//
free(phip);//
free(nitotal);//

free(Ucritf);//
free(muopen);//
free(kdrop);//
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

clock_t end = clock(); // finds the end time of the computation
double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
printf("in densfinorb: module ran in %f seconds\n", jobtime);
return problem;
} // closes densfinorb function 
