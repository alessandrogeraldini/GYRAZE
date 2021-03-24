// Author: Alessandro Geraldini
// MODIFIED 24 JUN 2019
/* This code calculates the lowest order ion density profile in the magnetic presheath for a given potential profile and entrance distribution function. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mps.h"

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

int densfinorb(int dodistfunc,double* ne_grid, double* phi_grid, int *p_sizepoint)
{// DECLARATIONS
double rho_over_unitgrid = 1.0; // this needs to be an input; for ions in the magnetic presheath this is 1.0, for electrons in the Debye sheath this can be chosen as min(1.0, rho_e/lambda_D)
double temp_over_unitgrid = 1.0;
clock_t begin = clock(); // Finds the start time of the computation
double phi_cut = -999.0;
int reflected = 0;
int upgraded = 1, p_size = *p_sizepoint;
double Ucrit = 0.0;
double alpha, *gg, *xx, pos, *phi, *phip, *ne, *niclosed, *niopen, *nitotal, **chi, **chiop, Te, Telarge, Tesmall, small = 1e-6, stoplimit = 0.96;
/* alpha is the angle the magnetic field makes with the wall, read from inputfile.txt; gg = sqrt(x). I often refer to it as g but call it gg to make searching for it easier; pos is x, the position (distance from the wall); phi contains the values of phi(x), extracted from a file. phip is phi prime, first derivative of phi. phipp is the second derivative; phipg and phippg are derivatives of phi wrt g = sqrt(x). ne is the electron density; niclosed is the closed orbit ion density, niopen is the open orbit ion density, nitotal is the sum of the two; newphi is the new electrostatic potential guess; chi is the effective potential; small is just a small number used to make some inequalities work numerically in case of exact equality. */
int debug=0, n, stop = 0, sizexbar, icritlow = 0, icritup = 0, ilocalmax = 0, L1, problem = 0, noback = 0, maxj;
/* n is the domain of the position x (the largest value of the index i, so e.g. gg[n]*gg[n] = L_2 in the paper); sizexbar is the domain of the position xbar (the largest value of the index j); L1 is the largest value of position x ( = gg^2) at which we evaluate the ion density */
double *openorbit, openorbitantycal, **mu, *muopen, *Ucritf, *xbar, *xbarop, **FF, *UU, *mumu, *twopimuprime, *openintegral, *xifunction;
/* openorbit is the integral Delta_M in the paper (an integral over the last closed orbit, representing how much U_perp - chi_M change in one full last orbit); openorbitantycal is the analytical value of openorbit for a flat potential (the first potential guess usually); mu is the array containing values of mu(xbar, Uperp). The index j is for values of xbar, k is for values of Uperp; xbar is the grid of values used in the closed orbit integral; xbarop is an alternative grid used for open orbit integral; FF contains the distribution function, read from the file distfile.txt. UU and mumu contain the values of U and mu corresponding to the function FF (which is F(mu, U)); FFprime is the numerical first derivative of F with respect to U; Many of these arrays are pointers because we don't initially know the array size, until we read the files (e.g. the file with the potential and the one with the distribution function); Note: first index of FF and FFprime is mu, second one is U; Interesting output variables; xifunction is the function of x which defines the grid of values of xbar by finding a chi whose minimum lies exactly at each grid point x */
int *jmclosed, *jmopen,  i=0, j=0, k=0, l=0, sizeU;
/* jmclosed represent minimum values of xbar above which we integrate open and closed orbit density integrals respectively (xbar_m,o and xbar_m in the paper); i is an index usually representing the positin x (= gg^2); j is an index usually representing the orbit position xbar; k is an index usually representing the energy Uperp (or velocity vx). It's always used in conjunction with j (and sometimes i); l is an index (used in for loops) usually representing the total energy (or velocity vz). It's used only in the DENSITY INTEGRALS part of the code; sizeU is the size of the integration range over U (or velocity vz). It is set later on in the code */
int *crossed_max, *crossed_min, *kdrop;
/* crossed_max and crossed_min is non-zero when a minimum (or maximum) of chi is found for some xbar[j]; kdrop is an integer which is non-zero only if there is an additional potential barrier for the particle at x << rho e.g. for an ion if the sheath reverses Uperp is allowed to be above chiMax and therefore the k = 0 of Uperp[j][k] starts for Uperp[j][0] = chiMax[j] + phi_barrier instead of Uperp[j][0] = chiMax[j] which is the conventional way */
	
int temporary = 0, nrows_distfile = 0, ncols_distfile = 0, ncols_Umufile = 0, sizeUU, sizemumu;
/* the temporary flag is a fix that ensures the minimum value of chi is stored (otherwise the minimum flag which just turned on prevents it); nrows_distfile counts the number of rows in the input file with the distribution function; ncols_distfile does the same for columns; ncols_Umufile  finds the number of columns for both rows or Umufile.txt. The first rows contains values of U, the second of mu. The columns are the possible values; sizeUU and sizemumu stores the final number of columns */
int *lowerlimit, *upperlimit, **upper, cols, rows, *itop, *imax, *imin;
/* lowerlimit represents the lower limit of k in the integrals over Uperp (or vx). It's needed because some of the earlies energies; (which are the largest because thy are values of chi stored after the maximum is found); may be so large that they are associated with very small values of the distribution function. This avoids integrating in an empty portion of phase space; upperlimit[j] represents the largest value of k (the smallest stored energy Uperp = chi_minimum) associated with some value of j; upper[j][i] represents the value of k associated with the smallest value of vx when integrating over Uperp. Going above upperlimit[j][i] makes Uperp < chi so velocities imaginary; imax/imin[j] stores the position of the maximum/minimum of the effective potential chi (It's x_M/x_m in the paper, which depends on xbar). */
double **Uperp, ***vx, *chiMax, *chiMpp, *chimin, oorbintgrd, oorbintgrdantycal, oorbintgrdBohm, xbarmax, intdvxopen, intdvxopenBohm=0.0, intdvxopenfluidBohm=0.0, intdxbaropenold, sizevxopen;
/* Uperp stores the possible values of Uperp associated with closed orbits, and so does vx; chiMax and chimin store the local maxima and minima of the effective potential maximum, oorbintgrd is the value of the integrand in the first open orbit integral (oorbintgrdantycal is the analytical result for flat potential); oorbintgrdBohm is the value of the integrand in the `Bohm' integral; intdvxopen and similars are where the integral of f_{0x} (v_x) over v_x, and its two important moments <v_x> (fluidBohm) and <1/v_x^2> (Bohm) are stored; They are a check that the extracted distribution function f_{0x} has the same moments it had when we carried out the density integral earlier; xbarmax is the maximum value of xbar in our grid*/
double vz, U, dvz = 0.2, dvzopen = 0.2, dvx, dvxopen, Deltavx, vxopen, dxbar, intdU=0.0, intdUopen=0.0, intdUopenBohm = 0.0;
	/* vz used in the density integral; U is the total energy, used in the density integral; dvz is the thickness of the vz grid used to take the integral over U (which is taken over vz in practice), dvzopen is the same for the open orbit piece; dvx is the thickness of the vx grid used to take the integral over Uperp ( which is taken over vx in practice). It must be evaluated because it depends on stored values of vx[j][i][k]; dvxopen is the thickness of the vx grid on which f_{0x} (v_x) is defined (and integrated to check consistency of its moments); dxbar is the thickness of the xbar grid; intdU is the value of the integral over U in the closed orbit density integration process; intdUopen is the same as above, for the open orbit integral; intdUopenBohm same, for Bohm integral */
double intdUold=0.0, intdvx=0.0, intdvxold = 0.0, intdxbar=0.0, intdxbaropen=0.0, intdxbaropenBohm = 0.0, Bohm = 0.0, F, Fold=0.0, Ucap;
	/* intdUold is a variable which stores the old intdU, so that the trapezium rule of integration can be applied (intdUold + intdU)*dvz; intdvx stores the integral over Uperp (hence over vx) in the closed orbit integral; intdxbar stores the value of the integral over xbar (which is the final result!), intdxbaropen does the same in the open orbit density integral; intdxbaropenBohm does the same for the Bohm integral; idealBohm is what the Bohm integral shoult be if Bohm condition is marginally satisfied; Bohm is the Bohm integral at the Debye sheath entrance x=0; F is the value of the distribution function evaluated in the density integrals by interpolating FF, and Fold is the `old' needed to apply the trapezium rule; Fprime is the bilinearly interpolated value of FFprime, and Fprimeold is the same at the previous grid point (needed for trapezium rule); used in INTEGRALS OF DISTRIBUTION FUNCTION AT INFINITY; Ucap is the topmost total energy integrated to */
double *storevals, *storevalspot, *storevalsinput;
// Pointers used to read the files
double intdUflow = 0.0, intdUflowold = 0.0, intdUopenflow = 0.0, intdUopenflowold = 0.0, intdxbaropenflow = 0.0, oorbintgrdflow = 0.0, oorbintgrdflowold = 0.0;
// values of various integrals
double oorbintgrdold=0.0, oorbintgrdBohmold=0.0, Fopen=0.0, intdUopenold=0.0, intdUopenBohmold;
double *chiprimetop, vx0open, aa, aaold; 
double intdUantycal=0.0, intdvxantycal=0.0, vxnew=0.0, vxold = 0.0, Uperpnew = 0.0, munew = 0.0, *xtop, intdUopenantycal=0.0;
/* intdUantycal is the integral over U (or v_z) for a flat potential profile (phi =0) for some value of xbar and Uperp; intdvxantycal  is the integral over Uperp (or vx) for a flat potential profile for some value of xbar; vxnew is the value of vx at the 'new' grid point, used in the vx integral (taken using the trapezium rule); vxold is the value of vx at the 'old' grid point, used in the vx integral; Uperpnew is the value of Uperp (used in the closed orbit density integral); munew is the valye of mu (used in the closed orbit density integral); xtop is the top bounce point x_t of the last closed orbit; intdUopenantycal is the analytical value of the integral over U  in the open orbit density integral */
double flux0; 
/* flux0 stores the flux at x=0 */
double oorbintgrdxbar, oorbintgrdsquare, intdUopenxbar = 0.0, intdUopensquare = 0.0, intdxbarxbar=0.0, intdxbarsquare=0.0, kBohmsquared;
/* oorbintgrdxbar is an integral over the open orbit distribution function at x=0 which is needed to evaluate a coefficient that appears when; expanding quasineutrality near x=0. It was just for playing around and at the moment plays no role in the code; similarly with all other integrals here */
double oorbintgrdxbarold, oorbintgrdsquareold, intdUopenxbarold, intdUopensquareold;
// all used to take integrals above
char line_distfile[200000], line_Umufile[200000], line_pot[200], line_input[20]; 
// The chars above are used to read the files with distribution function, potential input file and file containing grid of U and mu space; I can probably use one same variable instead of these four
/////////////////////////////////////////////////////////////
/* READ INPUT FILE */
FILE *input;
if ((input = fopen("inputfile.txt", "r")) == NULL)
{ printf("Cannot open %s\n", "inputfile.txt"); exit(EXIT_FAILURE); }
i=0;
while (fgets(line_input, 20, input) != NULL)
{	
	storevalsinput = linetodata(line_input, strlen(line_input), &ncols_Umufile);
	if (i==0)
	{	alpha = *storevalsinput;}
	else if (i==1)
	{	Te = *storevalsinput;}    
	i += 1;	
}
fclose(input);
// Introduce a number that equals Te when Te > 1 and 1 when Te<1. When Te is large, this number increases the number of grid points in tot. energy U to account for the fact that vz ~ vB > v_t,i when the electron temperature is large. Moreover, |v| ~ vB > v_t,i at the Debye sheath entrance, and so the number of grid points in xbar must also be increased.
if (Te>1.0)
{	Telarge = Te;
	Tesmall = 1.0;}
else
{	Telarge = 1.0;
	Tesmall = Te;}
xbarmax = 13.0*pow(Telarge, 0.5);
Ucap = 18.0 + 4.0*Te;
printf("Ucap =%f\n", Ucap);
i=0;
printf("alpha = %f\n", alpha);
/* READ POTENTIAL FILE */
// Now we read the file containing the potential profile and the grid of values of g = sqrt(x) on which is is defined
FILE *fp;
if ((fp = fopen("phidata.txt", "r")) == NULL)
{printf("Cannot open %s\n", "fp.txt"); exit(EXIT_FAILURE);}
/* The while loop counts the lines in the file to determine the size of the arrays to be created */
i=0;
while (fgets(line_pot, 200, fp) != NULL)
{i += 1;}
n=i; // n is the array size
rewind(fp); // go back to beginning of file
/* We initialize all arrays that contain functions of position x with the correct size n */
ne = (double*)calloc(n,sizeof(double)); // phi now has correct size
phi = (double*)calloc(n,sizeof(double)); // phi now has correct size
phip = (double*)calloc(n,sizeof(double)); // phi now has correct size
gg = (double*)calloc(n,sizeof(double)); // gg = sqrt(x) has correct size
xx = (double*)calloc(n,sizeof(double)); // xx = x has correct size
niclosed = (double*)calloc(n,sizeof(double)); // array containing the density of closed orbits
niopen = (double*)calloc(n,sizeof(double)); // array containing the density of open orbits
nitotal = (double*)calloc(n,sizeof(double)); // array containing the total density of ions 
jmclosed = (int*)calloc(n,sizeof(int));
jmopen = (int*)calloc(n,sizeof(int));
xifunction = (double*)calloc(n,sizeof(double));
// Begin reading potential file here
i=0;
while (fgets(line_pot, 200, fp) != NULL)
{
	storevalspot = linetodata(line_pot, strlen(line_pot), &ncols_Umufile);
	gg[i] = *storevalspot;
	xx[i] = gg[i]*gg[i];
	phi[i] = *(storevalspot+1);
	i += 1;
}
fclose(fp); // Close file

if (upgraded == 0) 
	Fgenerator(phi_grid[0]);
else 
	Fgenerator(phi[0]);

/* FORM XBAR GRIDS 
Take derivatives of phi and use them to obtain two grids for xbar, one to be used for closed orbits and one to be used for open orbits. */
xbar = (double*)calloc(n,sizeof(double));	
xbarop = (double*)calloc(n,sizeof(double)); // Open orbit grid defined such that every position x has a corresponding chi for which xM = x. This allows better resolution of the places where the open orbit integrand, oorbintgrd, is large (~alpha^(1/2))	
for (i=0; i<n; i++)
{	
	//if (i!=0)	xiprime[i-1] = (xiprime[i]-xiprime[i-1])/(xx[i]-xx[i-1]); 
	// Evaluate derivative of phi
	jmopen[i] = jmclosed[i] = 0;
	if (i == 0)
	{	
		//phip[i] = 0.5*(phi[i+2] - phi[i+1])/(xx[i+2] - xx[i+1]) + 0.5*(phi[i+1] - phi[i])/(xx[i+1] - xx[i]); 
		//phip[0] = 0.5*(phi[2] - phi[1])/(xx[2] - xx[1]) + 0.5*(phi[1] - phi[0])/(xx[1]-xx[0]) + (xx[0] - xx[1])*(0.5*(phi[3] - phi[2])/(xx[3]-xx[2]) - 0.5*(phi[1] - phi[0])/(xx[1]-xx[0]))/(xx[2]-xx[1]) ; 
		phip[0] = (phi[1] - phi[0])/(xx[1]-xx[0]) + (xx[0] - xx[1])*((phi[2] - phi[1])/(xx[2]-xx[1]) - (phi[1] - phi[0])/(xx[1]-xx[0]))/(xx[2]-xx[1]) ; }
		
	else if (i == n-1)
	{	
		phip[i] = 0.5*(phi[i] - phi[i-1])/(xx[i] - xx[i-1]); 
	}
		//phip[0] = phip[1] + (gg[0]*gg[0] - gg[1]*gg[1])*(phip[2] - phip[1])/(gg[2]*gg[2]-gg[1]*gg[1]); 
	else 
	{	
		//phip[i] = (phi[i+1] - phi[i])/(xx[i+1] - xx[i]); 
		phip[i] = 0.5*(phi[i+1] - phi[i])/(xx[i+1] - xx[i]) + 0.5*(phi[i] - phi[i-1])/(xx[i] - xx[i-1]); 
	}
	// xifunction is xbar corresponding to given position x being a stationary point
	xifunction[i] = gg[i]*gg[i] + phip[i]/2.0;
	if (debug == 1)
	{ 	printf("xifunction[%d] = %f\n", i, xifunction[i]); }
	if (i == 1)
	{	if  (xifunction[i] > xifunction[i-1]) // immediately found that xi is increasing at x=0 telling us x_c = 0
		{	icritup = i-1; } }
	else if (i > 1)
	{	if ( (xifunction[i] + small > xifunction[i-1]) && (xifunction[i-1] < xifunction[i-2] + small) ) 
		{	icritup = icritlow = i-1 ; }
		else if ( (xifunction[i] - small < xifunction[i-1]) && (xifunction[i-1] + small > xifunction[i-2]) ) 
		// found maximum of xi: this only happens when phi has noise in second derivative
		{	if ( xifunction[i-1] > xifunction[ilocalmax] || ilocalmax == 0 )  
			{	ilocalmax = i-1;} }
		else if ( (ilocalmax !=0) && (xifunction[ilocalmax] - small < xifunction[i-1]) && (xifunction[ilocalmax] + small > xifunction[i-2]) )
		{	icritup = i-1 ; } } }
//	if (i - 1 >= icritup)
//	{	xbar[j] = xifunction[i-1];
//		j += 1; } 
printf("icritup = %d\n", icritup);
if (ilocalmax != 0)
{	i=0;
	printf("ilocalmax = %d\n", ilocalmax);
	while (xifunction[i] > xifunction[ilocalmax])
	{	i+=1; }
	icritlow = i-1; }
j=0;
//if (icritlow == icritup)
//{ 	xbar[0] = xifunction[icritlow]; j = 1; }
k = icritup+1;	
for (i=icritlow-1;i>=0;i--)
{	while ( (xifunction[k] < xifunction[i]) && (k<n) )
	{	//if (xifunction[k] > xifunction[ilocalmax] + 0.01)
		//{	
		xbar[j] = xifunction[k]; j++; 
		//} 
		k++; }	
	//if (xifunction[i] > xifunction[ilocalmax] + 0.01)
	//{	
	xbar[j] = xifunction[i]; j++; 
	//} 
	}
while (k<n) { 	
	xbar[j] = xifunction[k]; 
	k++; j++; 
}	
sizexbar = j;
if (debug == 1)
{	for (j = 0; j < sizexbar; j++)
	{	printf("xbar[%d] = %f\n", j, xbar[j]); } }
printf("icritlow = %d, icritup = %d\n", icritlow, icritup);
printf("sizexbar = %d, n =%d\n", sizexbar, n);
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
for (j=0;j<sizexbar;j++) 
{	itop[j] = 0;
	imax[j] = imin[j] = -1;
	openorbit[j] = 0.0;
	crossed_min[j] = 0;
	crossed_max[j] = 0;
	xtop[j] = 0.0;
	chiprimetop[j] = 0.0;
	chimin[j] = 0.0;
	chiMax[j] = 0.0;
	upperlimit[j] = -1;
	lowerlimit[j] = 0; }
/* phi will be read from the file phidata.txt which contains values of phi and the corresponding value of sqrt(x) */
if ((fp = fopen("phidata.txt", "r")) == NULL)
{	printf("Cannot open %s\n", "fp.txt");
	exit(EXIT_FAILURE); }
/* Loop below initializes all 2d arrays which are functions of xbar and x. It allocates the right amount of memory to arrays of pointers of size xbar. The result is a 2D array indexed j (size sizexbar) and i or k (size n, see above) */
for (j=0; j < sizexbar; j++)
{	chi[j] = (double*)calloc(n,sizeof(double)); // array containing chi(xbar, x) indexes j and i 
	chiop[j] = (double*)calloc(n,sizeof(double)); // array containing chi(xbar, x) indexes j and i 
	Uperp[j] = (double*)calloc(n,sizeof(double)); // array containing Uperp for closed orbits indexes j and k
	mu[j] = (double*)calloc(n,sizeof(double)); // array containing adiab invariant mu(xbar, Uperp), indexes j and k
	upper[j] = (int*)calloc(n,sizeof(int)); // array containing the index of energy Uperp corresponding to chi(x)
	vx[j] = (double**)calloc(n,sizeof(double*)); // see below
	/* The loop below initializes the 3d array containing the velocity of a particle at a given orbit position xbar, with particle position x and with energy Uperp, indexed i, j and k. */
	for (i = 0; i < n; i++)
	{	vx[j][i] = (double*)calloc(n,sizeof(double)); } } // 3D array containing value of vx at different xbar, x and Uperp
i=0;  // Re-set i to zero.
/* This for loop fills in the arrays, calculating the integrals where necessary */
//printf("n=%d\n", n);
for (j=0; j<sizexbar; j++)
{	
	chi[j][i] =  pow((pow(gg[i],2.0) - xbar[j]), 2.0) + phi[i];
}
	// chiop[j][i] =  pow((pow(gg[i],2.0) - xbarop[j]), 2.0) + phi[i]; }
for (i=1; i<n; i++)
{	
	for (j=0; j<sizexbar; j++)
	{	
		chi[j][i] =  pow((pow(gg[i],2.0) - xbar[j]), 2.0) + phi[i];
//chiop[j][i] =  pow((pow(gg[i],2.0) - xbarop[j]), 2.0) + phi[i];
/* 
FINDING MAXIMA/MINIMA
*/
/* Below, we use the array elements to define arrays for the effective potential maxima and minima that exist for every xbar (index j). We search through the chi curve from x=0 to the largest value of x but we search through every chi curve first (so first scan in xbar at fixed x, then move to the next x). We create arrays of mu, Uperp and vx.
Because I am comparing neighbouring values of x to find a maximum, I need to consider two distinct cases.The first one is treated in the if loop below, which the program should enter only if chi is decreasing at x=0. Implying that chi(x=0) is an effective potential maximum. The second one is the else if loop after that, which finds maxima of chi that are stationary points by comparing the value of the function before and after each point. Note that because we compare the function at a point with the function at points before and after, the point we consider at every iteration step in the index i is indexed (i-1), compared with (i-2) and i.  */
		if ( (i==1) && (chi[j][i] < chi[j][i-1]) )
		{	
			if (phi_cut > chiMax[j]) 
				kdrop[j] = (int)((phi_cut - chiMax[j])/0.2); 
			else kdrop[j] = 0;
			crossed_max[j] += 1;		
			imax[j] = 0;
			chiMax[j] = chi[j][i-1];
			Uperp[j][kdrop[j]] = chi[j][0];
			mu[j][0] = 0.0;
			vx[j][0][0] = 0.0; } 
		else if ( (i > 1) && ((chi[j][i] < chi[j][i-1]) && (chi[j][i-1] > chi[j][i-2])) )
		{
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
			if (phi_cut > chiMax[j]) 
				kdrop[j] = (int)((phi_cut - chiMax[j])/0.2); 
			else kdrop[j] = 0;
		}
/* We store the index corresponding to the position of a minimum for a given value of xbar.*/
		else if ( (i > 1) && ((chi[j][i] > chi[j][i-1] - small ) && (chi[j][i-1] < chi[j][i-2] + small)) )
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
			if (Uperp[j][i-1-imax[j]+kdrop[j]] < 10.0*Telarge && Uperp[j][i-2-imax[j]+kdrop[j]] > 10.0*Telarge)
				lowerlimit[j] = i-1-imax[j]+kdrop[j]; 
			mu[j][i-1-imax[j]+kdrop[j]] = 0.0;
			upper[j][i-1] = i-1-imax[j]+kdrop[j];
/* Note that the size of the dimension of the array with values of Uperp is set to n, which is larger than the size it will turn out to be. Unfortunately in C I have no way to append elements to arrays as I go along, but need to give enough memory to the array from the start. n is the largest possible size the array could have. */
			for (k=0;k<=upper[j][i-1]; k++) 
			{	
				// replaced k with upper below
				if (upper[j][i-1] == 0)
				{	
					vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));
					mu[j][k] += 0.0; }
				else
				{	
					vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));
					if (debug==1) 
					{ 
						printf("vx[%d][%d][%d] = %f, Uperp = %f, chi = %f\n", j,i-1,k,vx[j][i-1][k], Uperp[j][k], chi[j][i-1]); 
					}
					mu[j][k] += (1.0/M_PI)*(vx[j][i-1][k] + vx[j][i-2][k])*(pow(gg[i-1], 2.0) - pow(gg[i-2], 2.0)); 
				}
			}
		}
/* Once we cross the minimum, we stop creating array elements with values of Uperp. However, we keep storing the value of vx associated with any given point x on an effective potential curve with xbar, with energy Uperp and using this value to finish performing the mu integral. This should happen as long the effective potential at the point under consideration is smaller than the effective potential maximum. */
		else if (crossed_min[j] == 1 && crossed_max[j] == 1 && chi[j][i-1] < chiMax[j] - small)
		{	
			for (k=0;k <= upperlimit[j] ;k++)
			{	
				if (chi[j][i-1] < Uperp[j][k])
				{	
					vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));
					mu[j][k] += (1.0/M_PI)*(vx[j][i-1][k] + vx[j][i-2][k])*(pow(gg[i-1], 2.0) - pow(gg[i-2], 2.0)); 
				}
				else if (Uperp[j][k-1] > chi[j][i-1] - small && Uperp[j][k] < chi[j][i-1] + small)
				{	
					upper[j][i-1] = k;
					mu[j][k] += (2.0/M_PI)*(vx[j][i-2][k])*(pow(gg[i-1], 2.0) - pow(gg[i-2], 2.0))*(Uperp[j][k] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]); 
				}
				if (mu[j][k] != mu[j][k])
				{	if (debug == 1)
					{	
						printf("mu is NAN\n"); 
					}
					problem = 1;  
				}  
			if (mu[j][k] != mu[j][k])
				{	
					if (debug ==1)
					{ 
						printf("mu is NAN\n"); 
					}
					problem = 1; 
				}  
			}
		}
/* When the effective potential at the iteration point (i-1) under consideration becomes larger than the effective potential maximum, we finish performing the mu integral. We also store the position of the top of the orbit which has chi = chiMax, in order to perform the open orbit integral. If the loop below is accessed, a switch it turned off to signify the no more closed orbits can be present */
		else if ( ( crossed_min[j] == 1 ) && ( crossed_max[j] == 1) && ( chi[j][i-1] > chiMax[j] - small) )
		{	itop[j] = i-2;
			xtop[j] = gg[i-2]*gg[i-2] + ((chiMax[j] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]))*(gg[i-1]*gg[i-1] - gg[i-2]*gg[i-2]);
			chiprimetop[j] = (chi[j][i-1] - chi[j][i-2])/(gg[i-1]*gg[i-1] - gg[i-2]*gg[i-2]);
			for (k=0; k<upper[j][i-2]; k++)
			{
				mu[j][k] += (1.0/M_PI)*(vx[j][i-2][k])*(pow(gg[i-1], 2.0) - pow(gg[i-2], 2.0))*(Uperp[j][k] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]);
				if (mu[j][k] != mu[j][k])  printf("mu is NAN\n");
			}
			crossed_max[j] = 0; }
		if (j!=0)
		{	
			if ( ((chiMax[j-1] < chi[j-1][i-1] + small) && (chiMax[j] + small > chi[j][i-1])) ) // || (imax[j-1] == -1 && imax[j] != -1) )
			{	jmclosed[i-1] = j-1;
				if (i - 1 > icritup)
				{ 	if (debug == 1) { printf("i = %d, j = %d, icritup = %d, icritlow = %d\n", i-1, j-1, icritup, icritlow); }
					jmopen[i-1] = j-1; 
					//printf("jmopen[%d] = %d\n", i-1, jmopen[i-1]);
					if (debug == 1)
					{	printf("jmopen[%d] = %d\n", i-1, jmopen[i-1]); 
					} 
				} 
			}
			if (i==n-1 && crossed_max[j] == 1 && crossed_min[j] == 1)
			{	itop[j] = n-1;
				if (noback == 0)
				{	maxj = j;	
					noback = 1; 
				} 
			} 
		} 
	} 
}
	/* Below loop is just here to place a cap on the maximum value of j at which the last orbit integral needs to be computed. If not present, the code still works if the open orbit density evaluation piece of the code limits the upper limit on xbar to some value (4). However, it will spit out error messages to do with the last orbit integral failing at large values of xbar, where the top bounce point of the last orbit lies so far out that it is beyond the x domain that we are considering. */
		//if (sswitch[j] != 0)
		//{
		//	if ((jmopen[i-1] == 0) && (j == sizexbar -1))
		//	{	maxj = j;	
		//		noback = 1; } } } } // closes all loops
		//if (i==n-1 && sswitch[j] != 0)
		//{	itop[j] = n-1;
		//	if (noback == 0)
		//	{	maxj = j;	
		//		noback = 1; } } } } // closes all loops
fclose(fp); // closes file
// OPEN ORBIT INTEGRAL
/* Now we perform the open orbit integral. We use a change of variables which makes the integrand smooth at the top bounce point. The change of variables is to some var = sqrt(x_t - x) */
muopen[0] = 0.0;
for (j=0;j<maxj;j++)
{
	muopen[j+1] = mu[j][0];
	Ucritf[j+1] = chiMax[j];
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
	aaold = 8.0*sqrt(xtop[j] - gg[imax[j]]*gg[imax[j]])/sqrt(-chiMpp[j]);
	//aaold = 0.0;
	for (i=imax[j]+1;i<itop[j];i++)
	{
		aa = 8.0*(gg[i]*gg[i] - gg[imax[j]]*gg[imax[j]])*sqrt(xtop[j] - gg[i]*gg[i])/vx[j][i][0];
		if (fabs(vx[j][i][0]) < 0.0001)
		{ printf("HERE INF, j=%d, i=%d\n", j, i); }
		if (aa != aa)
		{
			printf("aa is %f\ni is %d, j is %d\nitop is %d\n", aa, i, j, itop[j]);
			printf("xtop[%d] = %f, x[%d] = %f\n", j, xtop[j], i, gg[i]*gg[i]);
			printf("vx = %f\n", vx[j][i][0]);
			openintegral[j] += 0.0;
		}
		else
		{
		openintegral[j] += 0.5*(aa + aaold)*(sqrt(xtop[j] - gg[i-1]*gg[i-1]) - sqrt(xtop[j] - gg[i]*gg[i]));
		aaold = aa;
		}
	}
	aa = 8.0*(xtop[j] - gg[imax[j]]*gg[imax[j]])/sqrt(chiprimetop[j]);
	if (aa != aa)
	{
			printf("aa= %f, chiprimetop = %f, j = %d\n", aa, chiprimetop[j], j);
			aa = 0.0;
	}
	openintegral[j] += 0.5*(aa + aaold)*sqrt(xtop[j] - gg[itop[j]-1]*gg[itop[j]-1]);
	//openorbit[j] = openintegral[j];
	openorbit[j] = twopimuprime[j];  
	openorbitantycal = 4.0*M_PI*xbar[j];
	if (debug ==1)
	{	
		printf("twopimuprime[%d] = %f, openintegral[%d] = %f\n", j, twopimuprime[j], j, openintegral[j]); 
		printf("At xbar[%d] = %f, mu is %f, openorbit is %f while analytical one for a flat electrostatic potential is %f\n", j, xbar[j],mu[j][0], openorbit[j], openorbitantycal); 
	} 
}
Ucritf[0] = Ucritf[1] - muopen[1]*(Ucritf[2] - Ucritf[1])/(muopen[2] - muopen[1]);
//for (j=0;j<sizexbar;j++)
//{	
//	if (j != 0)
//	{	if (j!=1)
//		{
//		twopimuprime[j-1] = 2.0*M_PI*( 0.5*(mu[j][0] - mu[j-1][0])/(xbar[j] - xbar[j-1]) + 0.5*(mu[j-1][0] - mu[j-2][0])/(xbar[j-1] - xbar[j-2]) );	
//		}
//		else
//		{
//		twopimuprime[0] = 2.0*M_PI*( 0.5*(mu[1][0] - mu[0][0])/(xbar[j] - xbar[j-1]) ); 
//		}
//		openorbit[j-1] = twopimuprime[j-1];  
//		openorbitantycal = 4.0*M_PI*xbar[j-1];
//		if (debug ==1)
//		{	printf("openorbit[%d] = %f\n", j-1, openorbit[j-1]); 
//			printf("At xbar[%d] = %f, mu is %f, openorbit is %f while analytical one is %f\n", j-1, xbar[j-1],mu[j-1][0], openorbit[j-1], openorbitantycal); } } }

/* This part of the code extracts the distribution function from a file. We import the distribution function into a 2 dimensional array, F(mu,U). */
FILE *distfile, *Umufile;
if ((distfile = fopen("distfuncin.txt", "r")) == NULL)
{	
	printf("cannot open file %s\n", "distfuncin.txt");
	exit(-1); 
}
/* Count the number of rows in the distfuncin.txt file */
while(fgets(line_distfile, 20000, distfile) != NULL)
{	
	nrows_distfile += 1; 
}
/* Allocate the right amount of memory to the distribution function pointer */
FF = (double**)calloc(nrows_distfile,sizeof(double*));
/* The number of rows is also the the size of the array in UU */
//UU = malloc(nrows_distfile*sizeof(double));
cols = nrows_distfile;
nrows_distfile = 0; // Set number of rows counter to zero again
fclose(distfile); // Close file
/* Check file exists and can be opened
 */
if ((distfile = fopen("distfuncin.txt", "r")) == NULL)
{
	printf("cannot open file %s\n", "distfuncin.txt");
	exit(-1);
}
/* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and counting the columns. Once the columns are counted we allocate memory to FF and assign each number to a different column of FF (for a fixed row). The while loop then jumps to a new line and repeats the process */
while (fgets(line_distfile, 200000, distfile) != NULL)
{	//string = malloc(strlen(line_distfile)*sizeof(char));
	//string = line_distfile;
	//printf("The string is %s\n", string);
	//storevals = linetodata(string, strlen(string), &ncols_distfile);
	storevals = linetodata(line_distfile, strlen(line_distfile), &ncols_distfile);
	FF[nrows_distfile] = (double*)calloc(ncols_distfile,sizeof(double));
	mumu = (double*)calloc(ncols_distfile,sizeof(double));
	rows = ncols_distfile;
	FF[nrows_distfile] = storevals;
	nrows_distfile +=1; 
}
fclose(distfile);
printf("~~~~~The second element of FF is %f~~~~~\n",FF[0][1]);
/* Now we extract the two 1D arrays representing the grid points on which FF (the distribution function array) is defined. */
if ((Umufile = fopen("Umufile.txt", "r")) == NULL)
{	printf("cannot open file %s\n", "Umufile.txt");
	exit(-1); }
/* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and assigning the values of the first line to mumu, second to UU */
i=0;
while (fgets(line_Umufile, 200000, Umufile) != NULL)
{	//string = malloc(strlen(line_Umufile)*sizeof(char));
	//string = line_Umufile;
	//printf("The string is %s\n", string);
	//storevals = linetodata(string, strlen(string), &ncols_Umufile);
	storevals = linetodata(line_Umufile, strlen(line_Umufile), &ncols_Umufile);
	//index = 1;
	if (i == 0)
	{
	sizemumu = ncols_Umufile;
	mumu = storevals;
	//printf("mumu[%d] is %f and number of elements is %d\n", index, mumu[index], sizemumu);
	}
	else
	{
	sizeUU = ncols_Umufile;
	//UU = malloc(ncols_Umufile*sizeof(double));
	UU = storevals;
	//printf("UU[%d] is %f and number of elements is %d\n", index, UU[index], sizeUU);
	}
	i += 1; }
fclose(Umufile);
i=0;
printf("Array filling DONE\n");
clock_t int1 = clock(); // Finds the time of the computation so far
double inttime  = (double)(int1 - begin) / CLOCKS_PER_SEC;
printf("After array filling time is %f\n", inttime);
/* DENSITY INTEGRALS 
This part calculates the density integrals and outputs the result of the integration to a file fout and also the yz distribution function to three files one containing the distribution function the other two containing the velocity grid */
FILE *fout; 
if ((fout = fopen("PostProcessing/niout.txt", "w")) == NULL)
{	
	printf("Cannot open niout.txt");
	exit(EXIT_FAILURE);
}
i=0;
stop = 0;
printf("sizexbar = %d\n", sizexbar);
printf("Starting density evaluation\n");
while (stop == 0) // i<L1; i++)
{	pos = gg[i]*gg[i];
	intdxbar = 0.0;
	intdxbaropen = 0.0;
	intdxbaropenflow = 0.0;
	for (j=1; j<sizexbar; j++)
	{	vxnew = 0.0;
		intdUopenold = intdUopen;
		intdUopenflowold = intdUopenflow;
		intdUopenBohmold = intdUopenBohm;
		intdUopenxbarold = intdUopenxbar;
		intdUopensquareold = intdUopensquare;
		intdUopen = 0.0;
		intdUopenflow = 0.0;
		intdUopenBohm = 0.0;
		intdUopenxbar = 0.0;
		intdUopensquare = 0.0;
		if (j == jmopen[i])
		{	oorbintgrd = oorbintgrdold = 0.0;
			oorbintgrdflow = oorbintgrdflowold = 0.0;
			oorbintgrdBohm = oorbintgrdBohmold = 0.0;
			oorbintgrdxbar = oorbintgrdxbarold = 0.0;
			oorbintgrdsquare = oorbintgrdsquareold = 0.0;
			//if (icritlow - j < 0)
			//{	k = 0; }
			//else
			//{	k = icritlow - j; }
			//sizeU = (int) sqrt(Ucap - chi[j][k])/dvzopen;
			sizeU = (int) sqrt(Ucap - chiMax[j])/dvzopen;
			for (l=0; l < sizeU; l++)
			{	oorbintgrdold = oorbintgrd;
				oorbintgrdflowold = oorbintgrdflow;
				oorbintgrdBohmold = oorbintgrdBohm;
				oorbintgrdxbarold = oorbintgrdxbar;
				oorbintgrdsquareold = oorbintgrdsquare;
				vz = dvzopen*l;
				U = chi[j][k] + pow(vz, 2.0); // is it worth increasing accuracy?
				if (i > icritup)
				{	Fopen = bilin_interp(mu[j][0], U, FF, mumu, UU, sizemumu, sizeUU, -1, -1); }
				else
				{	Fopen = 0; }	
				oorbintgrd = sqrt(alpha*vz*openorbit[j])*Fopen;
				oorbintgrdflow = 0.5*alpha*vz*openorbit[j]*Fopen;
				//oorbintgrdBohm = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*vz*openorbit[j]))*Fopen;
				//oorbintgrdBohm = (alpha*vz*openorbit[j]/(pow(vx0open, 3.0)))*Fopen;
				//oorbintgrdBohm = oorbintgrd/pow(vx0open, 2.0);
				//oorbintgrdxbar = xbar[j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*vz*openorbit[j]))*Fopen;
				//oorbintgrdsquare = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + alpha*vz*openorbit[j]), 1.5))*Fopen;
				//if (j == jmopen[i])
				//{	oorbintgrdold = sqrt(alpha*vz*(openorbit[j-1]*(chiMax[j] - chi[j][i]) + openorbit[j]*(chiMax[j-1] - chi[j-1][i]))/(chiMax[j-1] - chi[j-1][i] + chiMax[j] - chi[j][i]))*Fopen; }
				if (oorbintgrd != oorbintgrd)
				{	
					problem = 1; 
					if (debug == 1)
					{
						printf("vx0open = %f, openorbit[%d] = %f, imaginary oorbintgrd in initial piece of integral due to negative value of openorbit?\n", vx0open, j, openorbit[j]); 
					}
				}
				if (i==0) // for a flat potential, oorbintgrd can be calculated analytically
				{	oorbintgrdantycal = sqrt(2.0*alpha*vz*M_PI*xbar[j])*(U-chiMax[j])*exp(-U)/pow(M_PI, 1.5); }
					//oorbintgrd = oorbintgrdantycal;
					//printf("oorbintgrd is %f, analytical is %F\n", oorbintgrd, oorbintgrdantycal);
				if (l!=0)
				{	intdUopen += 2.0*dvzopen*(oorbintgrd+oorbintgrdold);
					intdUopenflow += 2.0*dvzopen*(oorbintgrdflow+oorbintgrdflowold);
					intdUopenBohm += 2.0*dvzopen*(oorbintgrdBohm+oorbintgrdBohmold);
					intdUopenxbar += 2.0*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
					intdUopensquare += 2.0*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); }
				else
				{	intdUopen += 0.0;
					intdUopenflow += 0.0;
					intdUopenBohm += 0.0;
					intdUopenxbar += 2.0*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
					intdUopensquare += 2.0*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); } } }
		if (j > jmopen[i]) 
		{	oorbintgrd = oorbintgrdold = 0.0;
			oorbintgrdflow = oorbintgrdflowold = 0.0;
			oorbintgrdBohm = oorbintgrdBohmold = 0.0;
			oorbintgrdxbar = oorbintgrdxbarold = 0.0;
			oorbintgrdsquare = oorbintgrdsquareold = 0.0;
			sizeU = (int) sqrt(Ucap - chiMax[j])/dvzopen;
			for (l=0; l < sizeU; l++)
			{	oorbintgrdold = oorbintgrd;
				oorbintgrdflowold = oorbintgrdflow;
				oorbintgrdBohmold = oorbintgrdBohm;
				oorbintgrdxbarold = oorbintgrdxbar;
				oorbintgrdsquareold = oorbintgrdsquare;
				vz = dvzopen*l;
				U = chiMax[j] + pow(vz, 2.0);
				vx0open = sqrt(small + chiMax[j] - chi[j][i]);
				if (vx0open != vx0open)
				{		
					problem = 1; 
					if (debug == 1)
					{
						printf("HERE imaginary vx0open, j = %d, i is %d, chi[j][i] = %f, chiMax[j] = %f\n", j, i, chi[j][i], chiMax[j]); 
					}
				}
				Fopen = bilin_interp(mu[j][0], U, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
				oorbintgrd = (sqrt(vx0open*vx0open + alpha*vz*openorbit[j]) - vx0open)*Fopen;
				oorbintgrdflow = 0.5*alpha*vz*openorbit[j]*Fopen;
				oorbintgrdBohm = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*vz*openorbit[j]))*Fopen;
				//oorbintgrdBohm = (alpha*vz*openorbit[j]/(pow(vx0open, 3.0)))*Fopen;
				//oorbintgrdBohm = oorbintgrd/pow(vx0open, 2.0);
				oorbintgrdxbar = xbar[j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*vz*openorbit[j]))*Fopen;
				oorbintgrdsquare = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + alpha*vz*openorbit[j]), 1.5))*Fopen;
				//if (j == jmopen[i])
				//{	oorbintgrdold = sqrt(alpha*vz*(openorbit[j-1]*(chiMax[j] - chi[j][i]) + openorbit[j]*(chiMax[j-1] - chi[j-1][i]))/(chiMax[j-1] - chi[j-1][i] + chiMax[j] - chi[j][i]))*Fopen; }
				if (oorbintgrd != oorbintgrd)
				{	problem = 1; 	
					if (debug == 1)
					{	printf("vx0open = %f, openorbit[%d] = %f, HERE imaginary oorbintgrd in initial piece of integral due to negative value of openorbit?\n", vx0open, j, openorbit[j]); } }
				if (i==0) // for a flat potential, oorbintgrd can be calculated analytically
				{	oorbintgrdantycal = sqrt(2.0*alpha*vz*M_PI*xbar[j])*(U-chiMax[j])*exp(-U)/pow(M_PI, 1.5); }
					//oorbintgrd = oorbintgrdantycal;
					//printf("oorbintgrd is %f, analytical is %F\n", oorbintgrd, oorbintgrdantycal);
				if (l!=0)
				{	intdUopen += 2.0*dvzopen*(oorbintgrd+oorbintgrdold);
					intdUopenflow += 2.0*dvzopen*(oorbintgrdflow+oorbintgrdflowold);
					intdUopenBohm += 2.0*dvzopen*(oorbintgrdBohm+oorbintgrdBohmold);
					intdUopenxbar += 2.0*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
					intdUopensquare += 2.0*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); }
				else
				{	intdUopen += 0.0;
					intdUopenflow += 0.0;
					intdUopenBohm += 0.0;
					intdUopenxbar += 2.0*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
					intdUopensquare += 2.0*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); } }
			if (intdUopen != intdUopen) 
			{	problem = 1; 
				if (debug == 1)
				{	printf("HERE, j is %d\n", j); } }
			if (i==0) // at x=0 this can be evaluated analytically
			{	intdUopenantycal = 2.0*0.919*sqrt(2.0*alpha*xbar[j])*exp(-xbar[j]*xbar[j])/M_PI; } 
				//printf("xbar[%d] = %f, analytical is %f, numerical is %f\n", j, xbar[j], intdUopenantycal, intdUopen); printf("vx0open is %f (should be zero)\n", vx0open); 
			dxbar = xbar[j] - xbar[j-1];
			if ( (j == jmopen[i]+1) && (i > icritup) )
			{	if ( fabs(pos - gg[imax[j]]*gg[imax[j]]) < small )
				{	intdxbaropen += 0.0; }
				else
				{	dxbar = xbar[j] - (xbar[j-1]*(chiMax[j] - chi[j][i]) - xbar[j]*(chiMax[j-1] - chi[j-1][i]))/(-chiMax[j-1] + chi[j-1][i] + chiMax[j] - chi[j][i]); } }
//(chi[j][k]-chi[j][i])/(2.0*(pos-gg[k]*gg[k])); 
			intdxbaropen += 0.5*(intdUopen+intdUopenold)*dxbar;
			intdxbaropenflow += 0.5*(intdUopenflow+intdUopenflowold)*dxbar;
			intdxbaropenBohm += 0.5*(intdUopenBohm+intdUopenBohmold)*dxbar;
			intdxbarxbar += 0.5*(intdUopenxbar+intdUopenxbarold)*dxbar;
			intdxbarsquare += 0.5*(intdUopensquare+intdUopensquareold)*dxbar; 
			if (intdxbaropen != intdxbaropen) 
			{	problem = 1; 
				if (debug == 1)
				{	printf("HERE, j is %d\n", j); } } }
		if ( j>=jmclosed[i] )  /* We have entered the closed orbit integral */
		{	
			intdvxold = intdvx;
			//intdvxflowold = intdvxflow;
			//intdvxflow = 0.0;
			intdvx = 0.0;
			if (j == jmclosed[i])
			{	intdvx = 0.0;
				intdxbar += 0.0; }
			else if (j > jmclosed[i])
			{	for (k=lowerlimit[j]; k<upper[j][i]+1; k++)
				{	vxold = vxnew;
					intdUold = intdU;
					intdUflowold = intdUflow;
					intdUflow = 0.0;
					intdU = 0.0;
					if (k == upper[j][i])
					{	
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
					sizeU = (int) sqrt(Ucap - Uperp[j][k])/dvz;
					for (l=0; l < sizeU; l++)
					{	if (l!=0)
						{	
							Fold = F;
							vz = dvz*l;
							U = Uperpnew + pow(vz, 2.0);
							Ucrit = lin_interp(muopen, Ucritf, munew, sizexbar, 791);
							if (U < Ucrit) reflected = 1;
							if (reflected == 1) printf("Houston, we have a problem\n");
							reflected = 0;
							F = bilin_interp(munew, U, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
							intdU += (1.0+ (double)reflected)*dvz*(F+Fold);
							intdUflow += (1.0- (double)reflected)*alpha*dvz*vz*(F+Fold); }
						else
						{	vz = dvz*l;
							U = Uperpnew + pow(vz, 2.0);
							F = bilin_interp(munew, U, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
							intdU += 0.0;
							intdUflow += 0.0; } }
					intdUantycal = exp(-Uperp[j][k])*(1.0/(2.0*M_PI));// result with phi =0
					//intdU = intdUantycal;
					//printf("Analytical intdU is %f, numerical one is %f\n", intdUantycal, intdU);
					if (k==lowerlimit[j])
					{	intdvx += 0.0; }
					else
					{	
						intdvx += 4.0*0.5*dvx*(intdU+intdUold);
						//intdvxflow += 4.0*0.5*dvx*(intdUflow + intdUflowold); 
					}
					if (intdvx != intdvx)
					{	problem = 1; 
						if (debug ==1)
						{	printf("intdvx is NAN, j=%d, i=%d\n", j, i); } } }
					intdvxantycal = (2.0/(2.0*sqrt(M_PI)))*exp(-(pos-xbar[j])*(pos-xbar[j]))*erf(sqrt(xbar[j]*xbar[j]-(pos-xbar[j])*(pos-xbar[j])));
					if (debug == 1) {
					printf("pos=%f, i=%d, j=%d, Uperp is %f, chi is %f, vxnew and vxold are %f and %f and and vx is %f, dvx is %f, intdvx is %f, analytical one is %f, upper is %d, upperlimit is %d\n", pos, i, j, Uperpnew, chi[j][i], vxnew, vxold, vx[j][i][k], dvx, intdvx, intdvxantycal, upper[j][i], upperlimit[j]); }
				//printf("intdvx is %f, while the analytical one is %f\n", intdvx, intdvxantycal);
				//intdvx = intdvxantycal;
				if (j==jmclosed[i]+1 && i!=0)
				{
					//dxbar = (gg[imax[j]]*gg[imax[j]] - gg[imax[j]]*gg[imax[j]])*(xbar[j] - xbar[j-1])/(gg[imax[j]]*gg[imax[j]] - pos);
				//	dxbar = (chiMax[j]-chi[j][i])/(2.0*(pos-gg[imax[j]]*gg[imax[j]]));
					if (i==imax[j])
					{	dxbar = xbar[j] - xbar[j-1]; }
					else
					{	dxbar = xbar[j] - (xbar[j-1]*(chiMax[j] - chi[j][i]) - xbar[j]*(chiMax[j-1] - chi[j-1][i]))/(-chiMax[j-1] + chi[j-1][i] + chiMax[j] - chi[j][i]); }
					intdxbar += 0.5*(intdvx+intdvxold)*dxbar;
					if (intdxbar != intdxbar)
					{	problem = 1; 
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
	ne[i] = exp(phi[i]/Te); 
	niclosed[i] = intdxbar;
	niopen[i] = intdxbaropen;
	nitotal[i] = niopen[i] + niclosed[i];
	if (i == 0)
	{	Bohm = intdxbaropenBohm/nitotal[i];
		flux0 = (intdxbaropenflow)/nitotal[i]; 
		printf("intdxbarsquare = %f\n", intdxbarsquare);
		kBohmsquared = intdxbarxbar/(intdxbarsquare - 0.5*nitotal[0]); }
	if ((nitotal[i] >= stoplimit) || (i == n-1))
	{	
		stop = 1;
		L1 = i+1; 
	}
	if (debug == 1)
	{	
		printf("%f is CLOSED orbit density at position index %d, position %f\n", niclosed[i], i, pos);
		printf("%f is OPEN orbit density at position index %d, position %f\n", niopen[i], i, pos);
		printf("%f is TOTAL orbit density at position index %d, position %f\n", nitotal[i], i, pos); 
		if ( (nitotal[i] != nitotal[i]) || (nitotal[i] < small) )
		{ 	
			problem = 1; 
		} 
	} 
	fprintf(fout, "%f %f %f %f\n", gg[i], nitotal[i], niclosed[i], niopen[i]);
	i += 1; 
} 
fclose(fout);
clock_t int2 = clock(); // Finds the start time of the computation
double inttime1  = (double)(int2 - int1) / CLOCKS_PER_SEC;
printf("After density evaluation time is %f\n", inttime1);
clock_t end = clock(); // finds the end time of the computation
double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;

if (problem == 1)
{	
	printf("***WARNING*** there is a problem with the iteration\n"); 
	// If a problem is found, re-start iteration from a simple first guess
}

if (dodistfunc != 0)
{
	sizevxopen = (int) 100*sqrt(Telarge)/sqrt(Tesmall);
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
	printf("maxj = %d, sizeU = %d\n", maxj, sizeU);
	for (i=0; i<sizevxopen; i++)
	{	vxopen = i*dvxopen;
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
				Fopen = bilin_interp(mu[j][0], U, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
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
				//printf("tophat = %f when vx0open = %f, Deltavx = %f and vxopen = %f\n, F = %f\n", aa, vx0open, Deltavx, vxopen, F);
				if (l!=0)
				{	intdUopen += 2.0*(F+Fold)*dvzopen; } }
			if (j != 0)
			{	intdxbaropen += 0.5*(intdUopen + intdUopenold)*(xbar[j] - xbar[j-1]); } }
		//printf("FBohm = %f\n", intdxbaropen);
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
	fprintf(output, "%f %f\n", vxopen, intdxbaropen/nitotal[0]);
	}
	fclose(output);
	fclose(outputyz);
	fclose(outputyzvy);
	fclose(outputyzvz);

	// Calculating correction due to non boltzmann electrons
	// AG: to the Bohm condition

	/*the gradient of the function */
	int cont = 1;
	int i = 0;
	double ne_p0=0.0;

	if (upgraded == 0) {
		while (cont == 1)
		{
			if (i == p_size - 1)
			{
				printf("The phi value exceeds the given phi grid\n");
				printf("The phi value was: %.10f", phi[0]);
				exit(EXIT_FAILURE);
			}
			if (phi[0]/Te >= phi_grid[i])
			{
				if (phi[0]/Te < phi_grid[i + 1])
				{
					if (i != 0)
					{
						ne_p0 = (ne_grid[i + 1] - ne_grid[i - 1]) / (phi_grid[i + 1] - phi_grid[i - 1]);//+(phi[0] - phi_grid[i])* ((ne_grid[i + 1] + ne_grid[i - 1] - (2.0 * ne_grid[i])) / (pow((phi_grid[i + 1] - phi_grid[i]), 2.0)));
					}
					else
					{
						ne_p0 = (ne_grid[i + 1] - ne_grid[i]) / (phi_grid[i + 1] - phi_grid[i]);
					}
					cont = 0;
				}
			}
			i++;
		}
	}
	else {
		ne_p0 = (ne_grid[1] - ne_grid[0]) / (phi[1] - phi[0]);
		printf("ne_p0 = %f (%f)\n", ne_p0, (ne_grid[1] - ne_grid[0]) / (phi_grid[1] - phi_grid[0]));
	}

	printf("The density at x=0 obtained from the extracted distribution function is %f\nThe density at x=0 obtained from the ion density integral (more direct) is %f\n", intdvxopen, nitotal[0]);
	printf("Bohm integral = %f (obtained from the extracted distribution function at x=0)\nBohm integral = %f (obtained directly from the distribution function at infinity)\nThe Bohm condition is: Bohm integral = %f\n", intdvxopenBohm/intdvxopen, Bohm, (2.0*ne_p0)/(Te*nitotal[0])); 
	printf("The flow at x=0 obtained from the extracted distribution function is %f\nThe flow at x=0 obtained from the distribution function at infinity is %f\nThe Bohm speed is %f\n", intdvxopenfluidBohm/intdvxopen, flux0, 1.0/sqrt(2.0)); 
	printf("kBohmsquared= %f (should be >0)\n", kBohmsquared);

	clock_t end = clock(); // finds the end time of the computation
	double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Evaluation of ion distribution function at Debye sheath entrance ran in %f seconds\n", jobtime);
}
// If you love your variables (and your memory) set them free
free(ne);
//free(phi);
//free(phip);
free(gg);
free(xx);
free(niclosed);
free(niopen);
free(nitotal);
free(jmclosed);
free(jmopen);
free(xifunction);
free(xbar);
free(xbarop);

free(chimin);
free(chiMax);
free(chiMpp);
free(crossed_min);
free(crossed_max);
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
	free(upper[w]);
	for (int s = 0; s < n; s++)
	{
		free(vx[w][s]);
	}
	free(vx[w]);
}
for (int w = 0; w < nrows_distfile; w++)
{
	free(FF[w]);
}
free(FF);

if (upgraded == 1) {
	for (i=0; i<n; i++) {
		phi_grid[i] = phi[i];
		//printf("iondens.c: phi = %f\n", phi[i]);
	}
	*p_sizepoint = n;
}
free(phi);
free(phip);

printf("Iondens module ran in %f seconds\n", jobtime);
return problem;
} // closes main function 
