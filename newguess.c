
// LAST SUBSTANTIAL MODIFICATION MADE 13 NOV 2018
/* This code calculates the next electrostatic potential guess in the iteration to obtain the self-consistent magnetic presheath electrostatic potential profile */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mps.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#define TINY 1e-8

void newvcut(double *v_cut, double v_cutDS, double mioverme, double u_i, double u_e, double current, int *convergence) {
	double weight = 0.3;
	double old_v_cut = *v_cut, u_etilde;
	u_etilde = u_e - sqrt(mioverme/M_PI)*0.5*exp(-0.5*old_v_cut*old_v_cut);
	//printf("Old v_cut = %f\n", *v_cut); //printf("u_etilde = %f\n", u_etilde);
	printf("target current = %f\n", current);
	printf("current = %f = u_i (%f) - u_e (%f)\n", u_i-u_e, u_i, u_e);
	if ( (-current + u_i - u_etilde) < 0.0 ) {
		printf("WARNING: The log below would be negative, iteration avoided\n");
		*v_cut = old_v_cut;
	}
	else { 
		*v_cut = (1.0-weight)*(*v_cut) + weight*sqrt(-2.0*log(2.0*(-current + u_i - u_etilde)) - log(M_PI/mioverme) );
	}
	if ( (*v_cut*(*v_cut)*0.5 < 0.5*v_cutDS*v_cutDS ) )  *v_cut = v_cutDS;
	if ( (fabs(current - u_i + u_e)/u_i > 0.01) && (*convergence == 2) ) *convergence = 1;

}

void newguess(int *convergence, int problem, double Te, double alpha, double *x_grid, double* ne_grid, double *ni, double* phi_grid,int size_phigrid, int *psize_ngrid, double gammasq, double **distfunc_i, double *UU, double *mumu, int sizeUU, int sizemumu, double v_cutDS) { // gsl_permutation *p, gsl_matrix *m) {
//DECLARATIONS
	int size_ngrid = *psize_ngrid, new_size_ngrid = 0;
	double phiinf, phipinf, phippinf, decay, amp;
	double temp = 0.0, dirichlet = 1;
	clock_t begin = clock(); // Finds the start time of the computation
	double *phipp_red, deltaxsq;
	double fac=1.0;
	int debug = 0;
	int clip_index = 0;
	int row, col=0, i, j;
	char line[20000];
	double *storevals, k32, k32denom, k32num, k1DK, *Tevect, *k1GKvect; 
	double *phip, *phipp;;
	double k1DKnum1old, densinf, fluxinf, k1GK, **distfunc_iprime, k32denom1, k1DKnum1, k1DKnum, *newphi, *newphi_dir, k32denom1old, fluxinf1old, fluxinf1;
	double *deltanewphi;
	double densinf1old, densinf1, *gg;
	double CC, *flux, *ne, k1DKthres = 0.3, Fprime, Fprimeold, F, Fold, Teor1, fluxinfintgrdold, fluxinfintgrd;
	double phi1, weight, U, u, du;
	double sig=0.0, res = 0.0, av_dev=0.0, siglimit, siglimit1, reslimit, reslimit1, dev, devbig;
	double random_value;
	double negradinf, dev_0;
	//int size_phippit_DS = -1;
	//double *ni_DS_model;
	//double u_etilde, u_i, old_v_cut;
	FILE *fout2, *fpotold;

//if (*convergence == 1) method = 2;

printf("gammasqe = %f\n", gammasq);
printf("size_phigrid = %d\nsize_ngrid = %d\n", size_phigrid, size_ngrid);
	if (problem == 1) 
	{ 	
		*convergence = 0; 
	}

/* Below we initialize all arrays that contain functions of position x with the correct size n */
ne = (double*)calloc(    size_phigrid,sizeof(double)); // phi now has correct size
gg = (double*)calloc(    size_phigrid,sizeof(double)); // gg now has correct size
phip = (double*)calloc(  size_phigrid,sizeof(double)); // phi now has correct size
phipp = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size
newphi = (double*)calloc(size_phigrid,sizeof(double)); // same as above 
newphi_dir = (double*)calloc(size_phigrid,sizeof(double)); // same as above 
deltanewphi = (double*)calloc(size_phigrid,sizeof(double)); // same as above 
flux = (double*)calloc(  size_phigrid,sizeof(double)); // same as above 

row = 0;
if (Te>1.001)
{	
	Teor1 = Te;
}
else 
{	
	Teor1 = 1.0; 
}
siglimit1 = 0.2;//0.03
siglimit = 0.002;//0.007
reslimit1 = 0.01;//0.03
reslimit = 0.001;//0.03
if (gammasq > TINY) {
	siglimit = 0.005;
	siglimit1 = 0.1;
}
if (Te>2.1)
{
	siglimit1 = 0.012;
	siglimit = 0.012;
}
// Assign value to weight depending on whether convergence flag 

random_value = (double) rand()/RAND_MAX;//float in range 0 to 1
//printf("The random number between 0 and 1 is %f\n", random_value);
if (alpha < 3*M_PI/180.0)
	weight = 0.2;// 0.55
else if (gammasq < 1.0) weight = 2.0*(gammasq/(1.0+gammasq))*0.4;
else if (gammasq < TINY) weight = 0.4;
else weight = 0.2;
weight = 0.1;
printf("WEIGHT = %f\n", weight);

//if (alpha > 0.02) weight = 0.4;
//else weight = 0.2;


row = 0;
if ((fpotold = fopen("PostProcessing/phidataold.txt", "w")) == NULL)
{
	printf("Cannot open phidataold.txt");
	exit(EXIT_FAILURE);
}
for (i=0; i<size_phigrid; i++)
{
	gg[i] = sqrt(x_grid[i]);
	fprintf(fpotold, "%f %f\n", gg[i], phi_grid[i]);
}
fclose(fpotold); 
// INTEGRALS OF DISTRIBUTION FUNCTION AT INFINITY
/* In the section below I will evaluate the integrals that are necessary at the Chodura sheath entrance. They will be used in the potential iteration to extrapolate the electrostatic potential at large values of x which inevitably the code cannot solve for (boundary condition at infinity). */

//FILE *fpk1;
//if ((fpk1 = fopen("k1gyrokinetic.txt", "r")) == NULL) {	
//// Check for presence of file
//	printf("Cannot open %s\n", "k1gyrokinetic.txt");
//	exit(EXIT_FAILURE); 
//}
///* The while loop counts the lines in the file to determine the size of the arrays to be created */
//row=0;
//// Count the number of rows in the distfuncin.txt file
//while(fgets(line, 20000, fpk1) != NULL)
//	row += 1; 
//k1GKvect = (double*)calloc(row,sizeof(double));
//Tevect = (double*)calloc(row,sizeof(double));
//rewind(fpk1);
//row=0;
//col = 2;
//while (fgets(line, 20000, fpk1) != NULL) {
//	printf("i=%d\n", row);
//	printf("%s\n", line);
//	//storevals = linetodata(line, strlen(line), &col);
//	storevals = linetodatanew(line, &col);
//	printf("i=%d\n", row);
//
//	Tevect[row] = *storevals;
//	k1GKvect[row] = *(storevals+1);
//	printf("i=%d, Tevect=%f and k1GKvect=%f\n", row, Tevect[row], k1GKvect[row]);
//	row += 1;
//}
//fclose(fpk1); // close file


//printf("k1GK = %f\n", k1GK);

// Note: first index of distfunc_i and distfunc_iprime is mu, second one is U

distfunc_iprime = calloc(sizemumu,sizeof(double*));


/* The loop below evaluates the derivative of F with respect to total energy U at U = mu
 */
for (i=0; i<sizemumu; i++)
{
	distfunc_iprime[i] = calloc(sizeUU,sizeof(double));
	for (j=0; j < sizeUU; j++)
	{
		if ((j==sizeUU-1) || (i==sizemumu-1))
		{
			distfunc_iprime[i][j] = 0.0;	
		}
		else
		{
			distfunc_iprime[i][j] = (distfunc_i[i][j+1] - distfunc_i[i][j])/(UU[j+1] - UU[j]); 
		}
	}
}

k32num = 0.0;
for (j=0; j<sizemumu; j++)
{
	if (j!=0)
	{
		//k32num += 0.5*0.5*(16.0*M_PI/3)*(distfunc_iprime[j][j]+distfunc_iprime[j-1][j-1])*(mumu[j] - mumu[j-1]);
		k32num += 0.5*0.5*(16.0*M_PI/3)*(distfunc_iprime[j][0]+distfunc_iprime[j-1][0])*(mumu[j] - mumu[j-1]);
		//printf("k32num = %f\n", k32num);
	}
}
// Numerator ok

k1DKnum1old = 0.0;
k32denom1 = 0.0;
k1DKnum1 = 0.0;
k32denom1old = 0.0;
k32denom = k1DKnum = 0.0;
fluxinf1old = fluxinf1 = fluxinf = 0.0;
densinf1old = densinf1 = densinf = 0.0;
du = 0.01;
for (i=0; i<sizemumu; i++)
{
	if (i!=0)
	{
	k32denom1old = k32denom1;
	k32denom1 = 0.0;
	k1DKnum1old = k1DKnum1;
	k1DKnum1 = 0.0;
	fluxinf1old = fluxinf1;
	fluxinf1 = 0.0;
	densinf1old = densinf1;
	densinf1 = 0.0;
	}
	Fprimeold = Fprime = 0.0;
	Fold = F = 0.0;
	fluxinfintgrdold = fluxinfintgrd = 0.0;
	for (j=0; j< (int) 500*sqrt(Teor1); j++)
	{
		u = du*j;
		fluxinfintgrdold = fluxinfintgrd;
		Fprimeold = Fprime;
		Fold = F;
		U = pow(u, 2.0) + mumu[i];
		//Fprime = bilin_interp(mumu[i], U, distfunc_iprime, mumu, UU, sizemumu, sizeUU, -1, -1);
		//F = bilin_interp(mumu[i], U, distfunc_i, mumu, UU, sizemumu, sizeUU, -1, -1);
		Fprime = bilin_interp(mumu[i], U-mumu[i], distfunc_iprime, mumu, UU, sizemumu, sizeUU, -1, -1);
		F = bilin_interp(mumu[i], U-mumu[i], distfunc_i, mumu, UU, sizemumu, sizeUU, -1, -1);
		fluxinfintgrd = F*u*alpha;
		k32denom1 += 0.5*(Fprimeold + Fprime)*du;
		k1DKnum1 += 0.5*(Fprimeold + Fprime)*du;
		fluxinf1 += 0.5*(fluxinfintgrd + fluxinfintgrdold)*du;
		densinf1 += 0.5*(F + Fold)*du;
		//printf("k32denom1 = %f\n", k32denom1);
	}
	k32denom1 *= mumu[i];
	if (i!=0)
	{
	k32denom += 0.5*(mumu[i]-mumu[i-1])*(k32denom1 + k32denom1old);
	k1DKnum += 0.5*(mumu[i]-mumu[i-1])*(k1DKnum1 + k1DKnum1old);
	fluxinf += 0.5*(mumu[i]-mumu[i-1])*(fluxinf1 + fluxinf1old);
	densinf += 0.5*(mumu[i]-mumu[i-1])*(densinf1 + densinf1old);
	}
}

negradinf = (ne_grid[size_phigrid - 3] - ne_grid[size_phigrid - 1]) / (phi_grid[size_phigrid - 3] - phi_grid[size_phigrid - 1]);
if (fabs(phi_grid[size_ngrid]) < 0.05) negradinf = 1.0;
else negradinf = (ne[size_ngrid-2] - 1.0)/phi_grid[size_ngrid-2];
//negradinf = (ne_grid[size_ngrid-1] - 1.0)/phi_grid[size_ngrid-1];
negradinf = 1.0;
printf("negradinf = %.5f\n", negradinf);

//printf("Tevect = %f\n", Tevect[0]);
//i = 0;
//if (Te <= 1)
//{
//	while (Tevect[i] < Te)
//	{
//		//	printf("Tevect = %f\n", Tevect[i]);
//		i += 1;
//	}
//	k1GK = (k1GKvect[i] * (-Tevect[i - 1] + Te) + k1GKvect[i - 1] * (Tevect[i] - Te)) / (Tevect[i] - Tevect[i - 1]);
//
//}
//else if (Te > 1)
//{
//	k1GK = 0; //~~need to change this just didn't want the value looking correct in case someone thought it was~~
//}

k32denom = 1.0 + 4.0*M_PI*k32denom;
k1DKnum = negradinf - 4.0*M_PI*k1DKnum*Te;

k32 = k32num/k32denom;
k1DK = 2.0*k1DKnum/(Te*k32denom);
k1DK *= Te;
densinf *= (4*M_PI);
fluxinf *= (4*M_PI); 

//printf("k32num = %f\n", k32num);
printf("k1DK = %f\nk32 = %f\n", k1DK, k32);
printf("fluxinf = %f\n", fluxinf);
printf("densinf = %f\n", densinf);
for (i=0; i<size_phigrid; i++)
{	
	if ( (i != size_phigrid-1) && (i!=0) )
	{	
		//phipp[i] = ((x_grid[i] - x_grid[i-1])/(x_grid[i+1] - x_grid[i-1]))*(phip[i+1] - phip[i])/(x_grid[i+1] - x_grid[i]);
		//+ ((x_grid[i+1] - x_grid[i])/(x_grid[i+1] - x_grid[i-1]))*(phip[i] - phip[i-1])/(x_grid[i] - x_grid[i-1]);
		phip[i] = (phi_grid[i+1] - phi_grid[i]) / (x_grid[i+1] - x_grid[i]);
		phipp[i] = (phi_grid[i+1] - 2.0*phi_grid[i] + phi_grid[i-1]) / pow(x_grid[i+1] - x_grid[i], 2.0);
	}
	else if (i==size_phigrid-1) phipp[size_phigrid-1] = 0.0; //(phip[size_phigrid-1] - phip[size_phigrid-2])/(x_grid[size_phigrid-1] - x_grid[size_phigrid-2]);

	if (i==2) phipp[0] = phipp[1] - x_grid[1]*(phipp[2] - phipp[1])/(x_grid[2] - x_grid[1]);
	if (i==0) phip[i] = (phi_grid[i+1] - phi_grid[i]) / (x_grid[i+1] - x_grid[i]);
	if (i==size_phigrid-1) phip[i] = 0.0;
}




sig = 0.0;
res = 0.0;
devbig = dev = 0.0;
dev_0 = 0.0;
// The part that calculates the new electrostatic potential
for (i=0; i<size_ngrid; i++)  {
	if (gammasq > TINY) {
		//ni[i] /= ni[new_size_ngrid-1];
		//ne_grid[i] /= ne_grid[size_ngrid-1];
	}

	//if (ne_grid[i] > 1) {
	//	printf("ERROR: ne_grid[%d] = %f > 1.0\n", i, ne_grid[i]);
	//}
	if (ne_grid[i] < 0) {
		printf("ERROR: ne_grid[i] is less than zero\n");
		exit(-1);
	}
	//if ( (ni[i] > 0.95) && (new_size_ngrid == 0) ) {
	//	new_size_ngrid = i;
	//}
	if (debug == 0) {
		if (i==0) 
		printf("phi phipp ne ni=\n");
		printf("%d %f %f %f %f %f\n", i, phi_grid[i], gammasq*phipp[i], ne_grid[i], ni[i], (-ne_grid[i] + gammasq*phipp[i]) + ni[i]);
	}
	if ( (gammasq < TINY) || ( (i!=0) && (i!=size_ngrid-1) ) ) {
		dev = fabs((-ne_grid[i] + gammasq*phipp[i]) + ni[i]);
		sig += pow((-ne_grid[i] + gammasq*phipp[i]) + ni[i], 2.0);
		res += fabs((-ne_grid[i] + gammasq*phipp[i]) + ni[i]);
		av_dev += (-ne_grid[i] + gammasq*phipp[i] + ni[i]);
		if (dev > devbig)  devbig = dev; 
		if (i==0) dev_0 = fabs((-ne_grid[i] + gammasq*phipp[i]) + ni[i]);
	}
	if (i==size_ngrid-1)
	{
		if (problem == 0)
		{
			//newphi_dir[i] = Te * log( ni[i] - ne_grid[i] + exp(phi_grid[i]/Te) - gammasq*phipp[i]); 
			newphi_dir[i] = Te * log( ni[i] - ne_grid[i] + exp(phi_grid[i]/Te) ); 
			//newphi_dir[i] = Te * ( ni[i] - ne_grid[i] + phi_grid[i]/Te ); 
			//if (newphi_dir[i]/phi_grid[0] < 0.0) newphi_dir[i] = 0.0;
			if (newphi_dir[i] > 0.0) newphi_dir[i] = 0.0;
			printf("%f\n", Te * log(ni[i] - ne_grid[i] + exp(-fabs(phi_grid[i])/Te)) ) ; 
			phipp[i] = (ne_grid[i] - ni[i])/gammasq;
			if (phipp[i] > 0.0) phipp[i] = 0.0;
			//printf("%f\n", phipp[i]) ; 
			//printf("phip[%d] = %f\n", i, phip[i]);
			printf("newphi_dir[%d] = %f\n", i, newphi_dir[i]);
			
			phip[i] = 0.0;
		}
		else
		{
			newphi_dir[i] = 0.0;
		}
		if (fabs(k1DK) < k1DKthres)
		{
			CC = sqrt(20.0/k32)/pow(-0.5*newphi_dir[i], 0.25) - x_grid[i];
			printf("k32 = %f\t newphi_dir[%d] = %f\tx_grid[%d] = %f\n", k32, i, newphi_dir[i], i, x_grid[i]);
			printf("C_3/2 = %f\n", CC);
			if (fabs(newphi_dir[i]) < 1.0e-10) CC = 0.0;
		}
		else if (fabs(k1DK) > k1DKthres)
		{
			printf("~~~Chodura Oversatistied~~~\n");
			//phi1 = newphi[i]*exp(-k1GK*x_grid[i]);
			//phi1 = newphi[i] * exp(k1DK * x_grid[i]);
			phi1 = newphi_dir[i] * exp(x_grid[i]);
		}
	}
	else  {
		// This is where the iteration happens
		if (problem == 0) {
			//newphi_dir[i] = Te * log( ni[i] - ne_grid[i] + exp(phi_grid[i]/Te) - gammasq*phipp[i] ) ; 
			newphi_dir[i] = Te * log( ni[i] - ne_grid[i] + exp(phi_grid[i]/Te) ) ; 
			//newphi_dir[i] = Te * ( ni[i] - ne_grid[i] + phi_grid[i]/Te ); 
			if (newphi_dir[i] > 0.0) newphi_dir[i] = 0.0;
			//printf("%f\n", Te * log(ni[i] - ne_grid[i] + exp(-fabs(phi_grid[i])/Te)) ) ; 
			printf("newphi_dir[%d] = %f\n", i, newphi_dir[i]);
			phipp[i] = (ne_grid[i] - ni[i])/gammasq;
			if (phipp[i] > 0.0) phipp[i] = 0.0;
			//printf("%f\n", phipp[i]) ; 
		}
		else {
			newphi[i] = 0.0; 
		}
	}
}
if (new_size_ngrid == 0) new_size_ngrid = size_ngrid;
sig /= new_size_ngrid;
sig = sqrt(sig);
res /= new_size_ngrid;
av_dev /= new_size_ngrid;

if (gammasq > TINY)
	printf("test boundary value of phi: %f\n\n\n\n", newphi_dir[new_size_ngrid-1] );

phiinf = -0.00000;
deltaxsq = x_grid[1]*x_grid[1];

if ( (gammasq > TINY) ) {
	if (dirichlet == 1) {

	 gsl_vector *newphi_gsl = gsl_vector_alloc (new_size_ngrid-2);
	  gsl_matrix * m = gsl_matrix_alloc (new_size_ngrid-2, new_size_ngrid-2);

	  for (i = 0; i < new_size_ngrid-2; i++)
	    for (j = 0; j < new_size_ngrid-2; j++) {
		if (i == j) 
	      gsl_matrix_set (m, i, j, -2.0);
	      //gsl_matrix_set (m, i, j, -2.0 - deltaxsq*phi_grid[i]/gammasq);
		else if ( (i==j+1) || (i==j-1) )
	      gsl_matrix_set (m, i, j, 1.0);
		else 
	      gsl_matrix_set (m, i, j, 0.0);
		}

	phipp_red = malloc((new_size_ngrid-2)*sizeof(double));
	for (i=1;i<new_size_ngrid-1; i++) {
		phipp_red[i-1] = deltaxsq*phipp[i];// - deltaxsq*phi_grid[i]/gammasq;
		if (i==1) phipp_red[i-1] += 0.5*v_cutDS*v_cutDS;
		//if (i==new_size_ngrid-2) phipp_red[i-1] -= newphi_dir[new_size_ngrid-1];
		if (i==new_size_ngrid-2) phipp_red[i-1] -= phiinf;
	}

	  gsl_vector_view phipp_gsl
	    = gsl_vector_view_array (phipp_red, new_size_ngrid-2);


	  int s;
	  gsl_permutation * p = gsl_permutation_alloc (new_size_ngrid-2);
	  gsl_linalg_LU_decomp (m, p, &s);

	  gsl_linalg_LU_solve (m, p, &phipp_gsl.vector, newphi_gsl);

	  //printf ("newphi_gsl = \n");
	  //printf("%f\n", -0.5*v_cutDS*v_cutDS);
	  //gsl_vector_fprintf (stdout, newphi_gsl, "%g");

	  newphi[0] = - 0.5*v_cutDS*v_cutDS;
	  printf("%f\n", newphi[0]);
	  for (i=0; i<new_size_ngrid-2; i++) {
		temp = gsl_vector_get(newphi_gsl, i);
		//newphi[i+1] = temp + phi_grid[i+1];
		newphi[i+1] = temp;
		//if (newphi[i+1] > Te*log(gammasq*phipp[i] - ne_grid[i] + exp(phi_grid[i]) + ni[i])) newphi[i] = Te*log(gammasq*phipp[i] - ne_grid[i] + exp(phi_grid[i]) + ni[i]);
		//if ( (newphi[i+1] > 0.0) && (clip_index==0) ) 
		//if ( (newphi[i+1] > newphi_dir[new_size_ngrid-1]) && (clip_index==0) ) 
		if ( (newphi[i+1] > phiinf) && (clip_index==0) ) 
		{
			//newphi[i+1] = 0.0;
			printf("newphi at edge = %f\n", newphi_dir[new_size_ngrid-1]);
			printf("clipped\n\n\n\n\n\n\n");
			clip_index = i+1;
		}
		printf("%f\n", newphi[i+1]);
		//printf("i = %d/%d\n", i, new_size_ngrid-1);
	  }
	  //newphi[new_size_ngrid-1] = Te*log(gammasq*phipp[i] - ne_grid[i] + exp(phi_grid[i]) + ni[i]);
	  //newphi[new_size_ngrid-1] = newphi_dir[new_size_ngrid-1];
	  //newphi[new_size_ngrid-1] = 0.0;
	  printf("is this equal to clip_index %d\n", new_size_ngrid);
	  printf("%f\n", newphi[new_size_ngrid-1]);


	  gsl_permutation_free (p);
	  gsl_matrix_free (m);
	  gsl_vector_free (newphi_gsl);
	  free(phipp_red);

	}
	else {
		newphi[size_ngrid-1] = 0.0;
		newphi[size_ngrid-2] = 0.0;
		for (i=new_size_ngrid-3; i>=0; i--) {
			newphi[i] = phipp[i+1]*deltaxsq + 2.0*newphi[i+1] - newphi[i+2];
		}
		for (i=new_size_ngrid-3; i>=0; i--) {
			newphi[i] *= (-0.5*v_cutDS*v_cutDS/newphi[0]);
			printf("newphi in neumann bc gives %f\n", newphi[i]);
		}
	}
}

if (clip_index == 0) clip_index = new_size_ngrid;


printf("sig = %f, res = %f, devbig = %f, dev_0 = %f\n", sig, res, devbig, dev_0);

printf("conditions %d %d %d\n", (res < reslimit), devbig < 6.0*reslimit, dev_0 < 3.0*reslimit);
if (res < reslimit && devbig < 6.0*reslimit && dev_0 < 3.0*reslimit)  {
	*convergence = 2;
}
else {
	if ((fout2 = fopen("phidata.txt", "w")) == NULL)
	{
		printf("Cannot open phidata.txt");
		exit(EXIT_FAILURE);
	}
	for (i=0; i<size_phigrid; i++)
	{	
		fprintf(fout2, "%f %f\n", gg[i], newphi[i]);
	}
	fclose(fout2);
	if (*convergence == 0 && ( res < 2.0*reslimit1 ) && (dev_0 < 5.0*siglimit1) ) //|| devbig < 0.2) )
	{
		*convergence = 1; 
		printf("Now don't smooth = %d, sig = %f\n", *convergence, sig);
	}
	else if ((*convergence == 1) && (res > 2.0*reslimit1))
	{
		*convergence = 0;
	}
}

if (gammasq < TINY) 
	printf("Bohm integral should be %f\n", 2.0*Te*(ne_grid[1] - ne_grid[0])/((phi_grid[1] - phi_grid[0])*ni[0]));
//phipinf = (newphi[clip_index-1] - newphi[clip_index-2])/x_grid[1];
phipinf = (newphi[clip_index-1] - newphi[clip_index-2])/(x_grid[clip_index-1] - x_grid[clip_index-2]);
//phippinf = (newphi[clip_index-1] - 2.0*newphi[clip_index-2] + newphi[clip_index-3])/(pow(x_grid[1], 2.0));
//phiinf = 3.0*phipinf*phipinf/(2.0*phippinf);
//for (i=clip_index-1; i>=0; i--) newphi[i] += phiinf;
//phipinf += phipp[clip_index-1]*(x_grid[clip_index] - x_grid[clip_index-1]);
//decay = -(phipinf/newphi[clip_index-1]);
//amp = newphi[clip_index-1]*pow(x_grid[clip_index-1] + decay, 2.0);
printf("newphi  = %f, previous newphi = %f\n", newphi[clip_index-1], newphi[clip_index-2]);
//amp = - phipinf/newphi[clip_index-1];
//decay = -log(-newphi[clip_index-1]/amp)/x_grid[clip_index-1];

decay = -(newphi[clip_index-1]/phipinf)*2.0 - x_grid[clip_index-1];
//decay = sqrt((newphi[clip_index-1]/fabs(phipp[clip_index-1]))*6.0) - x_grid[clip_index-1];
if (x_grid[clip_index-1] + decay < 0.0) {
	printf("warning\n"); 
	decay = 0.0;
}
amp = newphi[clip_index-1]*pow(x_grid[clip_index-1] + decay, 2.0);
//amp = newphi[clip_index-1]*pow(x_grid[clip_index-1], 2.0);
//decay = 0.0;
if (dirichlet == 0) { amp = 0.0; decay = 0.0;}
printf("amp = %f, decay = %f, newphi = %f %f, phipinf = %f\n", amp, decay, newphi[clip_index-1], newphi[clip_index-2], phipinf);
printf("clip_index = %d\n", clip_index);
for (i=clip_index; i<size_phigrid; i++) {
	if (gammasq < TINY) {
		if (problem == 0) {
			if (fabs(k1DK) < k1DKthres) {
				newphi_dir[i] = -2.0*(400.0/pow(k32, 2.0))/pow(x_grid[i] + CC, 4.0);
			}
			else if (fabs(k1DK) > k1DKthres)
			{
				//newphi[i] = phi1*exp(k1GK*x_grid[i]); // k1GK should be negative
				//newphi[i] = phi1 * exp(-k1DK * x_grid[i]);
				newphi_dir[i] = phi1 * exp(-1.0 * x_grid[i]);
			}
			//printf("phi = %f, newphi[i] = %f\n", phi_grid[i], newphi[i]);
		}
		else {
			newphi_dir[i] = 0.0;
		}
		//newphi_dir[i] = 0.0;
	}
	else {
		//newphi_dir[i] = 0.0;
		//newphi[i] = -amp * exp(-decay * x_grid[i]);
		newphi[i] = amp/pow(decay+x_grid[i], 2.0);
		//newphi[i] = newphi[clip_index-1]*(/pow(decay+x_grid[i], 2.0);
		printf("newphi[%d] = %f\n", i, newphi[i]);
	}
}
for (i=size_phigrid-1; i>=0; i--) {
	if ( (gammasq < TINY) ) // || (i > icut-1) )
		newphi[i] = newphi_dir[i] ;
	//else phi_grid[i] *= (-0.5*v_cutDS*v_cutDS/phi_grid[0]);

	newphi[i] = weight*newphi[i] + (1.0-weight)*phi_grid[i] ; 
	//printf("%f %f %f\n", phi_grid[i], newphi[i], newphi_dir[i]);
	phi_grid[i] = newphi[i];
}
free(phip);
free(phipp);
free(flux);
free(newphi);
free(deltanewphi);
free(ne);
free(gg);

for (i=0; i<sizemumu; i++) {
	 free(distfunc_iprime[i]);
}
free(distfunc_iprime);

free(newphi_dir);

printf("av_dev = %f\n\n\n\n", av_dev);

//*psize_ngrid = clip_index;
clock_t end = clock(); // finds the end time of the computation
double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
printf("Module to obtain the new guess for the electrostatic potential ran in %f\n", jobtime);
return;
}


