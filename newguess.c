
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

void newguess(int *convergence, int problem, double Te, double alpha, double *x_grid, double* ne_grid, double *ni, double* phi_grid,int size_phigrid, int size_ngrid, double gammasq, double **distfunc_i, double *UU, double *mumu, int sizeUU, int sizemumu, double v_cutDS) {
//DECLARATIONS
	double phipinf, decay, amp, fict_charge = 0.0;
	double temp = 0.0;
	clock_t begin = clock(); // Finds the start time of the computation
	double *phipp_red, deltaxsq;
	double fac=1.0;
	int method = 2;
	int debug = 0, gone_negative = 0;
	int row, col=0, i, ii, j;
	char line[20000];
	int newmethod = 1;
	double deltax, deltaphi;
	double newweight = 0.2, newphiprime;
	double *storevals, k32, k32denom, k32num, k1DK, *Tevect, *k1GKvect; 
	double *phipg, *phippg, *phip, *phipp;;
	double *altphipp, *altphipsqovtwo, *altphip;
	double k1DKnum1old, densinf, fluxinf, k1GK, **distfunc_iprime, k32denom1, k1DKnum1, k1DKnum, *newphi, k32denom1old, fluxinf1old, fluxinf1;
	double *deltanewphi;
	double densinf1old, densinf1, *gg;
	double CC, *flux, *ne, k1DKthres = 0.3, Fprime, Fprimeold, F, Fold, Teor1, fluxinfintgrdold, fluxinfintgrd;
	double phi1, weight, U, u, du;
	double smoothweight = 2.0/3.0;
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

//ni_DS_model = malloc(size_ngrid*sizeof(double));
/* Below we initialize all arrays that contain functions of position x with the correct size n */
ne = (double*)calloc(    size_phigrid,sizeof(double)); // phi now has correct size
gg = (double*)calloc(    size_phigrid,sizeof(double)); // gg now has correct size
phip = (double*)calloc(  size_phigrid,sizeof(double)); // phi now has correct size
phipg = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size
phipp = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size
phippg = (double*)calloc(size_phigrid,sizeof(double)); // phi now has correct size
newphi = (double*)calloc(size_phigrid,sizeof(double)); // same as above 
deltanewphi = (double*)calloc(size_phigrid,sizeof(double)); // same as above 
flux = (double*)calloc(  size_phigrid,sizeof(double)); // same as above 

altphipp = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size
altphipsqovtwo = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size
altphip = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size

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
reslimit1 = 0.02;//0.03
reslimit = 0.003;//0.03
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
if (*convergence == 0)
{	
	//weight = 0.2 + 0.2*random_value/pow(Teor1, 1.0); 
	//weight = 0.1 + 0.2*random_value/pow(Teor1, 1.0); 
	weight = 0.5;
}
else
{	
	//weight = 0.1 + 0.2*random_value/pow(Teor1, 1.0); 
	weight = 0.2;
}
//if (gammasq > TINY) 
if (alpha < 2*M_PI/180.0)
	weight = 0.2;// 0.55
else weight = 0.4;

//if (alpha<0.05) weight = 0.1; 

printf("WEIGHT = %f\n\n\n\n\n", weight);

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

k32denom = 1 + 4.0*M_PI*k32denom;
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

//if (fabs(gammasq) > TINY) {

	//k1DK = 1.0; // temporary	
	//weight = 0.05;
	//for (i=0; i<size_phigrid; i++)
	//{	
	//	if ( (i != size_phigrid-1) && (i!=0) )
	//	{	
	//		//phip[i] = ((x_grid[i] - x_grid[i-1])/(x_grid[i+1] - x_grid[i-1]))*(phi_grid[i+1] - phi_grid[i])/(x_grid[i+1] - x_grid[i]);
	//		//+ ((x_grid[i+1] - x_grid[i])/(x_grid[i+1] - x_grid[i-1]))*(phi_grid[i] - phi_grid[i-1])/(x_grid[i] - x_grid[i-1]);
	//		//phipg[i] = (phi_grid[i+1] - phi_grid[i])/(gg[i+1] - gg[i]);
	//		//if (i!=0) {	
	//		//	phipp[i] = 2.0*(phip[i] - phip[i-1])/(x_grid[i+1] - x_grid[i-1]);
	//		//	phippg[i] = 2.0*(phipg[i] - phipg[i-1])/(gg[i+1] - gg[i-1]); 
	//		//}	
	//		//if (i==2)  {	
	//		//	phipp[0] = phipp[1] - x_grid[1]*(phipp[2] - phipp[1])/(x_grid[2] - x_grid[1]);
	//		//	phippg[0] = phippg[1] - sqrt(x_grid[1])*(phippg[2] - phippg[1])/(sqrt(x_grid[2]) - sqrt(x_grid[1]))  ; 
	//		//} 
	//	}
	//	else if (i==0) phip[0] = (phi_grid[1] - phi_grid[0])/(x_grid[1] - x_grid[0]);
	//	else if (i==size_phigrid-1) phip[size_phigrid-1] = (phi_grid[size_phigrid-1] - phi_grid[size_phigrid-2])/(x_grid[size_phigrid-1] - x_grid[size_phigrid-2]);
	//	//else {	
	//	//	phip[i] = phip[i-1];
	//	//	phipg[i] = phipg[i-1];
	//	//	phipp[i] = phipp[i-1];
	//	//	phippg[i] = phippg[i-1];
	//	//} 
	//}
	for (i=0; i<size_phigrid; i++)
	{	
		if ( (i != size_phigrid-1) && (i!=0) )
		{	
			//phipp[i] = ((x_grid[i] - x_grid[i-1])/(x_grid[i+1] - x_grid[i-1]))*(phip[i+1] - phip[i])/(x_grid[i+1] - x_grid[i]);
			//+ ((x_grid[i+1] - x_grid[i])/(x_grid[i+1] - x_grid[i-1]))*(phip[i] - phip[i-1])/(x_grid[i] - x_grid[i-1]);
			phip[i] = (phi_grid[i+1] - phi_grid[i]) / (x_grid[i+1] - x_grid[i]);
			phipp[i] = (phi_grid[i+1] - 2.0*phi_grid[i] + phi_grid[i-1]) / pow(x_grid[i+1] - x_grid[i], 2.0);
			//phipg[i] = (phi_grid[i+1] - phi_grid[i])/(gg[i+1] - gg[i]);
			//if (i!=0) {	
			//	phipp[i] = 2.0*(phip[i] - phip[i-1])/(x_grid[i+1] - x_grid[i-1]);
			//	phippg[i] = 2.0*(phipg[i] - phipg[i-1])/(gg[i+1] - gg[i-1]); 
			//}	
			//if (i==2)  {	
			//	phipp[0] = phipp[1] - x_grid[1]*(phipp[2] - phipp[1])/(x_grid[2] - x_grid[1]);
			//	phippg[0] = phippg[1] - sqrt(x_grid[1])*(phippg[2] - phippg[1])/(sqrt(x_grid[2]) - sqrt(x_grid[1]))  ; 
			//} 
		}
		else if (i==size_phigrid-1) phipp[size_phigrid-1] = 0.0; //(phip[size_phigrid-1] - phip[size_phigrid-2])/(x_grid[size_phigrid-1] - x_grid[size_phigrid-2]);

		if (i==2) phipp[0] = phipp[1] - x_grid[1]*(phipp[2] - phipp[1])/(x_grid[2] - x_grid[1]);
		if (i==0) phip[i] = (phi_grid[i+1] - phi_grid[i]) / (x_grid[i+1] - x_grid[i]);
		if (i==size_phigrid-1) phip[i] = 0.0;
		
		//if (i==2)  {	
		//	phipp[0] = phipp[1] - x_grid[1]*(phipp[2] - phipp[1])/(x_grid[2] - x_grid[1]);
		//}
		//else {	
		//	phip[i] = phip[i-1];
		//	phipg[i] = phipg[i-1];
		//	phipp[i] = phipp[i-1];
		//	phippg[i] = phippg[i-1];
		//} 
	}
	//for (i=1; i<size_phigrid-1; i++) 
	//	phipp[i] = 0.4*phipp[i] + 0.3*phipp[i-1] + 0.3*phipp[i+1];
	//i=0;
	//while (weight * (ni[i] - gammasq*phipp[i] - ne_grid[i] + exp(phi_grid[i]/Te)) + (1.0-weight)*exp(phi_grid[i]/Te) > 1.0) {
	//	size_ngrid = i;
	//	i++;
	//}
//}




sig = 0.0;
res = 0.0;
devbig = dev = 0.0;
dev_0 = 0.0;
// The part that calculates the new electrostatic potential
for (i=0; i<size_ngrid; i++)  {
	if (gammasq > TINY) {
		//ne_grid[i] /= ne_grid[size_ngrid-1];
		ni[i] /= ni[size_ngrid-1];
		// this is necessary when solving Poisson's equation with Dirichlet boundary conditions
		// at x=infinity
		//printf("hellow\n\n\n\n\n");
	}

	if (ne_grid[i] > 1) {
		printf("ERROR: ne_grid[%d] = %f > 1.0\n", i, ne_grid[i]);
	}
	if (ne_grid[i] < 0) {
		printf("ERROR: ne_grid[i] is less than zero\n");
	}
	if (debug == 0) {
		if (i==0) 
		printf("phi phipp ne ni=\n");
		printf("%f %f %f %f %f\n", phi_grid[i], gammasq*phipp[i], ne_grid[i], ni[i], (-ne_grid[i] + gammasq*phipp[i]) + ni[i]);
	}
	//dev = fabs((ne_grid[i] + gammasq*phipp[i])/ni[i] - 1.0);
	//sig += pow((ne_grid[i] + gammasq*phipp[i])/ni[i] - 1.0, 2.0);
	//res += fabs((ne_grid[i] + gammasq*phipp[i])/ni[i] - 1.0);
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
			if (fabs(gammasq) < TINY) {
				//newphi[i] = - ((phi_grid[0]-TINY)/fabs(phi_grid[0]-TINY)) * Te * log(weight * (ni[i] - gammasq*phipp[i] - ne_grid[i] + exp(phi_grid[i]/Te)) + (1.0-weight)*exp(-fabs(phi_grid[i])/Te)); 
				newphi[i] = - ((phi_grid[0]-TINY)/fabs(phi_grid[0]-TINY)) * Te * log(ni[i] + gammasq*phipp[i] - ne_grid[i] + exp(-fabs(phi_grid[i])/Te)) ; 
				//newphi[i] = Te * log( (ni[i] + gammasq*phipp[i] - ne_grid[i] + exp(phi_grid[i]/Te)) ); 
				altphipp[i] = - (-ne_grid[i] + ni[i] + fict_charge + fac*exp(-fabs(phi_grid[i])/Te));
				if (newphi[i]/phi_grid[0] < 0.0) newphi[i] = 0.0;
			}
			else {
				//phipp[i] = (-ne_grid[i] + ni[i])/gammasq;
				phipp[i] = (ne_grid[i] - ni[i])/gammasq - phipp[i];
				//if (phipp[i] > 0.0) phipp[i] = 0.0;
				//altphipp[i] = -phip[i]*(ne_grid[i] - ni[i])/gammasq;
				altphipp[i] = (ne_grid[i] - ni[i] + fac*exp(-fabs(phi_grid[i])/Te));
				printf("altphipp[%d] = %f\n", i, altphipp[i]);
				printf("phip[%d] = %f\n", i, phip[i]);
				
				phip[i] = 0.0;
			}
		}
		else
		{
			newphi[i] = 0.0;
		}
		if (fabs(k1DK) < k1DKthres)
		{
			CC = sqrt(20.0/k32)/pow(-0.5*newphi[i], 0.25) - x_grid[i];
			printf("k32 = %f\t newphi[%d] = %f\tx_grid[%d] = %f\n", k32, i, newphi[i], i, x_grid[i]);
			printf("C_3/2 = %f\n", CC);
			if (fabs(newphi[i]) < 1.0e-10) CC = 0.0;
		}
		else if (fabs(k1DK) > k1DKthres)
		{
			printf("~~~Chodura Oversatistied~~~\n");
			//phi1 = newphi[i]*exp(-k1GK*x_grid[i]);
			//phi1 = newphi[i] * exp(k1DK * x_grid[i]);
			phi1 = newphi[i] * exp(x_grid[i]);
		}
	}
	else  {
		// This is where the iteration happens
		if (problem == 0) {
			if (fabs(gammasq) < TINY) {
			//if (weight * (ni[i] - gammasq*phipp[i] - ne_grid[i] + exp(phi_grid[i]/Te)) + (1.0-weight)*exp(-fabs(phi_grid[i])/Te) < 1.0) {
				altphipp[i] = (ne_grid[i] - ni[i] + fict_charge + fac*exp(-fabs(phi_grid[i])/Te));
				//newphi[i] = - ((phi_grid[0]-TINY)/fabs(phi_grid[0]-TINY)) * Te * log(weight * (ni[i] - gammasq*phipp[i] - ne_grid[i] + exp(phi_grid[i]/Te)) + (1.0-weight)*exp(-fabs(phi_grid[i])/Te)); 
				newphi[i] = - ((phi_grid[0]-TINY)/fabs(phi_grid[0]-TINY)) * Te * log(ni[i] + gammasq*phipp[i] - ne_grid[i] + exp(-fabs(phi_grid[i])/Te)) ; 
				  //newphi[i] = Te * log( (ni[i] + gammasq*phipp[i] - ne_grid[i] + exp(phi_grid[i]/Te)) ); 
			  //newphi[i] = - ((phi_grid[0]-TINY)/fabs(phi_grid[0]-TINY)) * Te * log(weight * (ni[i] + gammasq*phipp[i] - ne_grid[i] + exp(phi_grid[i]/Te)) + (1.0-weight)*exp(-fabs(phi_grid[i])/Te)); 
			}
			//else {
			//  newphi[i] = 0.0;
			//}
			else {
				//phipp[i] = (ne_grid[i] - ni[i])/gammasq;
				phipp[i] = (ne_grid[i] - ni[i])/gammasq - phipp[i];
				altphipp[i] = (ne_grid[i] - ni[i] + fac*exp(phi_grid[i]/Te) + fict_charge);
				//if (phipp[i] > 0.0) phipp[i] = 0.0;
			}
		}
		else {
			newphi[i] = 0.0; 
		}
	}
}
sig /= size_ngrid;
sig = sqrt(sig);
res /= size_ngrid;
av_dev /= size_ngrid;

deltaxsq = x_grid[1]*x_grid[1];
 gsl_vector *newphi_gsl = gsl_vector_alloc (size_ngrid-2);
if ( (gammasq > TINY) ){
  gsl_matrix * m = gsl_matrix_alloc (size_ngrid-2, size_ngrid-2);

  for (i = 0; i < size_ngrid-2; i++)
    for (j = 0; j < size_ngrid-2; j++) {
	if (i == j) 
      gsl_matrix_set (m, i, j, -2.0);
	else if ( (i==j+1) || (i==j-1) )
      gsl_matrix_set (m, i, j, 1.0);
	else 
      gsl_matrix_set (m, i, j, 0.0);
	}

phipp_red = malloc((size_ngrid-2)*sizeof(double));
//phipp_red[0] = deltaxsq*phipp[0];
for (i=1;i<size_ngrid-1; i++) {
	//phipp_red[i] = deltaxsq*phipp[i];
	phipp_red[i-1] = deltaxsq*(ne_grid[i] - ni[i])/gammasq ;
	if (i==1) phipp_red[i-1] += 0.5*v_cutDS*v_cutDS;
}

  gsl_vector_view phipp_gsl
    = gsl_vector_view_array (phipp_red, size_ngrid-2);


  int s;

  gsl_permutation * p = gsl_permutation_alloc (size_ngrid-2);

  //gsl_linalg_LU_decomp (&m.matrix, p, &s);

  gsl_linalg_LU_decomp (m, p, &s);

  //gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
  gsl_linalg_LU_solve (m, p, &phipp_gsl.vector, newphi_gsl);
  //gsl_linalg_LU_solve (m, p, b, x);

  //printf ("newphi_gsl = \n");
  //printf("%f\n", -0.5*v_cutDS*v_cutDS);
  //gsl_vector_fprintf (stdout, newphi_gsl, "%g");

  newphi[0] = - 0.5*v_cutDS*v_cutDS;
  for (i=0; i<size_ngrid-2; i++) {
	temp = gsl_vector_get(newphi_gsl, i);
	if (temp < 0.0) {
		newphi[i+1] = temp;
	}
	else newphi[i] = 0.0;
  	//printf("%f\n", newphi[i]);
  	//printf("i = %d/%d\n", i, size_ngrid-1);
  }
  newphi[size_ngrid-1] = 0.0;

  gsl_permutation_free (p);

//printf("2nd element is %f\n", newphi_gsl->data[1]);
 //exit(0);
}


printf("sig = %f, res = %f, devbig = %f, dev_0 = %f\n", sig, res, devbig, dev_0);


altphipsqovtwo[size_ngrid] = 0.0;
//devbig = 0.0;
//printf("energy comparison\n");
if (gammasq > TINY) {
	for (i=size_ngrid-1; i>=0; i--) {
		deltaphi = fabs(phi_grid[i+1] - phi_grid[i]);
		altphipsqovtwo[i] = altphipsqovtwo[i+1] + deltaphi*altphipp[i]; 
		altphip[i] = altphipsqovtwo[i] - 0.5*gammasq*phip[i]*phip[i] - fac* ( 1.0 - exp(-fabs(phi_grid[i]/Te)) ) ;
		//res += fabs( altphip[i]/altphipsqovtwo[i] ); //
		//printf("negative? %f\n", altphipsqovtwo[i] - gammasq*0.5*phip[i]*phip[i] - 1.0);
		//printf("altphipp[%d] = %f\n", i, altphipp[i]);
		//if (i!= size_ngrid-1){
		//dev = ( altphipsqovtwo[i] - 0.5*gammasq*phip[i]*phip[i] - 1.0 + exp(-fabs(phi_grid[i])/Te) ) ;/// altphipsqovtwo[i];
	}
	//else if (method == 2) {
	//	phip[size_ngrid] = 0.0;
	//	for (i=size_ngrid-1; i>=0; i--) {
	//		deltax = x_grid[i+1] - x_grid[i];
	//		phip[i] = phip[i+1] - phipp[i]*deltax;
	//		printf("phip[%d] = %f\n", i, phip[i]);
	//		printf("phipp[%d] = %f\n", i, phipp[i]);
	//		printf("deltax = %f\n", deltax);
	//	}

	//	//for (ii=i; ii>=0; ii--)
	//	//{	
	//	//	//phip[ii] -= deltax*(phipp[i+1] + phipp[i])*0.5; 
	//	//	altphipsqovtwo[ii] -= deltax*altphipp[i]; 
	//	//} 
	//} 
}

printf("conditions %d %d %d\n", (res < reslimit), devbig < 6.0*reslimit, dev_0 < 3.0*reslimit);
if (res < reslimit && devbig < 6.0*reslimit && dev_0 < 3.0*reslimit)  {
	*convergence = 2;
}
else {
	// SMOOTHING OF PHI
	// Now we take derivatives of phi_grid
	//for (i=0; i<size_phigrid; i++)
	//{	
	//	if (i != size_phigrid-1)
	//	{	
	//		//phip[i] = (newphi[i+1] - newphi[i])/(x_grid[i+1] - x_grid[i]);
	//		phipg[i] = (newphi[i+1] - newphi[i])/(gg[i+1] - gg[i]);
	//		if (i!=0) {	
	//			//phipp[i] = 2.0*(phip[i] - phip[i-1])/(x_grid[i+1] - x_grid[i-1]);
	//			phippg[i] = 2.0*(phipg[i] - phipg[i-1])/(gg[i+1] - gg[i-1]); 
	//		}	
	//		if (i==1)  {	
	//			//phipp[0] = phipp[1];
	//			phippg[0] = phippg[1]; 
	//		} 
	//	}
	//	else {	
	//		//phip[i] = phip[i-1];
	//		phipg[i] = phipg[i-1];
	//		//phipp[i] = phipp[i-1];
	//		phippg[i] = phippg[i-1];
	//	} 
	//}
	/* Smooth second derivative of phi wrt sqrt(x) */
	//if (*convergence == 0)
	//{ 	
		//for (i=0; i<size_ngrid; i++) 
			//newphi[i] = phi_grid[size_ngrid];
		//for (i=0; i<size_ngrid; i++) {	
		//	newphi[i] = phi_grid[size_ngrid];
		//	phipg[i] = phipg[size_ngrid];
		//	if (i!=0)
		//	{	
		//		phippg[i] = ((1.0-smoothweight)*phippg[i-1]*(gg[i+1] - gg[i]) + smoothweight*phippg[i]*(gg[i+1] - gg[i-1]) + (1.0-smoothweight)*phippg[i+1]*(gg[i] - gg[i-1]))/((gg[i+1] - gg[i-1]));
		//		//phipp[i] = ((1.0-smoothweight)*phipp[i-1]*(x_grid[i+1] - x_grid[i]) + smoothweight*phipp[i]*(x_grid[i+1] - x_grid[i-1]) + (1.0-smoothweight)*phipp[i+1]*(x_grid[i] - x_grid[i-1]))/((x_grid[i+1] - x_grid[i-1])); 
		//	} 
		//}
		//// Re-integrate smoothed second derivative to first
		//for (i=size_ngrid-1; i>=0; i--)
		//{	for (ii=i; ii>=0; ii--)
		//	{	
		//		phipg[ii] -= (gg[i+1] - gg[i])*(phippg[i+1] + phippg[i])*0.5; 
		//	} 
		//}  
		//		//phip[ii] -= (x_grid[i+1] - x_grid[i])*(phipp[i+1] + phipp[i])*0.5; 
		//// Re-integrate first derivative to phi
		//for (i=size_ngrid-1; i>=0; i--)
		//{	for (ii=i; ii>=0; ii--)
		//	{	
		//		newphi[ii] -= (gg[i+1] - gg[i])*(phipg[i+1] + phipg[i])*0.5; 
		//	} 
		//} 
	//}
	
	//if (size_phippit_DS < 0) size_phippit_DS = size_ngrid;
	// can get rid of loop below now
	if (gammasq > 100000.0+TINY)  {
		deltax = x_grid[0];
		////printf("newphi[%d] = %f\n", size_ngrid-1, newphi[size_ngrid-1]);
		//for (i=1; i<size_ngrid; i++) {
		//	deltax = (x_grid[i] - x_grid[i-1]);
		//	newphi[i] = newphi[i-1]*(1.0-deltax*exp(-fabs(phi_grid[i-1]))/(gammasq*phip[i-1])) + deltax*altphip[i-1]/(phip[i-1]*gammasq); // - fac*deltax*exp(-fabs(phi_grid[i]/Te))); 
		//	printf("newphi = %f\n", newphi[i]);
		//	//newphi[i] = newphi[i-1] + sqrt(2.0* (altphipsqovtwo[i-1] + exp(-fabs(newphi[i-1])) - 1.0)/gammasq)*deltax;
		//}
		if (method == 1) {
			deltanewphi[size_ngrid] = 0.0;
			for (i=size_ngrid-1; i>=0; i--)
			{	
				deltax = (x_grid[i+1] - x_grid[i]);
				//printf("altphip[%d] = %f\n", i, altphip[i]);
				//printf("phip[%d] = %f\n", i, phip[i]);
				deltanewphi[i] = (phip[i]*gammasq*newphi[i+1] - altphip[i]*deltax)/(phip[i]*gammasq - fac*deltax*exp(phi_grid[i]/Te)); 
				//printf("deltanewphi[%d] = %f\n", i, deltanewphi[i]);
			} 
		}
		else if (method == 2) {
			//deltanewphi[0] = 0.0;
			//for (i=1; i<=size_ngrid; i++)
			//{	
			//	deltax = (x_grid[i+1] - x_grid[i]);
			//	//printf("altphip[%d] = %f\n", i, altphip[i]);
			//	//printf("phip[%d] = %f\n", i, phip[i]);
			//	//deltanewphi[i] = newphi[i-1]*(1.0 - deltax*exp(-phi_grid[i-1])/(gammasq*phip[i-1])) +  altphip[i-1]*deltax/(phip[i-1]*gammasq);// - fac*deltax*exp(-fabs(phi_grid[i]/Te))); 
			//	deltanewphi[i] = deltanewphi[i] 
			//	//deltanewphi[i-1] = deltanewphi[i] - phip[i-1]*deltax;
			//	
			//	printf("deltanewphi[%d] = %f\n", i, deltanewphi[i]);
			//} 

			//deltanewphi[size_ngrid-1] = 0.0;
			//deltanewphi[size_ngrid] = 0.0;
			//for (i=size_ngrid-2; i>=0; i--)
			//{	
			//	deltax = (x_grid[i+1] - x_grid[i]);
			//	deltanewphi[i] = 2.0*deltanewphi[i+1] - deltanewphi[i+2] + deltax*deltax*phipp[i+1];
			//	//deltanewphi[i-1] = deltanewphi[i] - phip[i-1]*deltax;
			//	
			//	printf("deltanewphi[%d] = %f\n", i, deltanewphi[i]);
			//} 

			altphip[size_ngrid] = 0.0;
			for (i=size_ngrid-1; i>=0; i--)
			{	
				deltax = (x_grid[i+1] - x_grid[i]);
				altphip[i] = altphip[i+1] - deltax*phipp[i];
				//deltanewphi[i-1] = deltanewphi[i] - phip[i-1]*deltax;
				
				//printf("deltanewphi[%d] = %f\n", i, deltanewphi[i]);
			} 
			deltanewphi[0] = 0.0;
			for (i=1; i<size_ngrid; i++)
			{	
				deltax = (x_grid[i+1] - x_grid[i]);
				deltanewphi[i] = deltanewphi[i-1] + deltax*altphip[i];
				//deltanewphi[i-1] = deltanewphi[i] - phip[i-1]*deltax;
				
				printf("deltanewphi[%d] = %f\n", i, deltanewphi[i]);
			} 
		}
		for (i=0; i<=size_ngrid; i++) {
			newphi[i] = phi_grid[i] + deltanewphi[i]; // = newphi[i-1]*(1.0 - deltax*exp(-phi_grid[i-1])/(gammasq*phip[i-1])) +  altphip[i-1]*deltax/(phip[i-1]*gammasq);// - fac*deltax*exp(-fabs(phi_grid[i]/Te))); 
			printf("newphi[%d] = %f\n", i, newphi[i]);
			if (newphi[i] > 0.0) newphi[i] = 0.0;
		}
	}

	//if (newmethod == 1)
	//{ 	

	//	for (i=0; i<size_ngrid; i++) {	
	//		newphi[i] = phi_grid[size_ngrid];
	//		phipg[i] = phipg[size_ngrid];
	//		//phippg[i] = (ni[i] - ne_grid[i] + newweight*phippg[i])/newweight;
	//		//if (i!=0)
	//		//{	
	//		//	phippg[i] = ((1.0-smoothweight)*phippg[i-1]*(gg[i+1] - gg[i]) + smoothweight*phippg[i]*(gg[i+1] - gg[i-1]) + (1.0-smoothweight)*phippg[i+1]*(gg[i] - gg[i-1]))/((gg[i+1] - gg[i-1]));
	//		//	//phipp[i] = ((1.0-smoothweight)*phipp[i-1]*(x_grid[i+1] - x_grid[i]) + smoothweight*phipp[i]*(x_grid[i+1] - x_grid[i-1]) + (1.0-smoothweight)*phipp[i+1]*(x_grid[i] - x_grid[i-1]))/((x_grid[i+1] - x_grid[i-1])); 
	//		//} 
	//	}
	//	// Re-integrate smoothed second derivative to first
	//	for (i=size_ngrid-1; i>=0; i--)
	//	{	for (ii=i; ii>=0; ii--)
	//		{	
	//			phipg[ii] -= (gg[i+1] - gg[i])*(phippg[i+1] + phippg[i])*0.5; 
	//		} 
	//	}  
	//	//phip[ii] -= (x_grid[i+1] - x_grid[i])*(phipp[i+1] + phipp[i])*0.5; 
	//	//for (i=size_phigrid-1; i>=size_ngrid; i--)
	//	//	printf("i = %d, phi = %f, newphi = %f\n", i, phi_grid[i], newphi[i]);
	//	// Re-integrate first derivative to phi
	//	for (i=size_ngrid-1; i>=0; i--)
	//	{	for (ii=i; ii>=0; ii--)
	//		{	
	//			newphi[ii] -= (gg[i+1] - gg[i])*(phipg[i+1] + phipg[i])*0.5; 
	//		} 
	//		//printf("i = %d, phi = %f, newphi = %f\n", i, phi_grid[i], newphi[i]);
	//	} 
	//}

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

//printf("old new\n");
for (i=0; i<size_ngrid; i++) {
	//printf("(old phi, unchanged new phi) = %f, %f\n", phi_grid[i], newphi[i]);
	//if ( gammasq > TINY ) {
	//	newphi[i] = weight*newphi[i] + phi_grid[i];
	//}
	//else {
	newphi[i] = weight*newphi[i] + (1.0-weight)*phi_grid[i];
	//}
	//printf("%f %f\n", phi_grid[i], newphi[i]);
	phi_grid[i] = newphi[i];
}
phipinf = (newphi[size_ngrid-1] - newphi[size_ngrid-2])/x_grid[1];
//decay = -(phipinf/newphi[size_ngrid-1]);
decay = -(newphi[size_ngrid-1]/phipinf)*2.0 - x_grid[size_ngrid-1];
amp = newphi[size_ngrid-1]*pow(x_grid[size_ngrid-1] + decay, 2.0);
for (i=size_ngrid; i<size_phigrid; i++) {
	if (fabs(gammasq) < TINY) {
		if (problem == 0) {
			if (fabs(k1DK) < k1DKthres)
			{
				newphi[i] = -2.0*(400.0/pow(k32, 2.0))/pow(x_grid[i] + CC, 4.0);
			}
			else if (fabs(k1DK) > k1DKthres)
			{
				//newphi[i] = phi1*exp(k1GK*x_grid[i]); // k1GK should be negative
				//newphi[i] = phi1 * exp(-k1DK * x_grid[i]);
				newphi[i] = phi1 * exp(-1.0 * x_grid[i]);
			}
			//printf("phi = %f, newphi[i] = %f\n", phi_grid[i], newphi[i]);
		}
		else {
			newphi[i] = 0.0;
		}
	}
	else {
		newphi[i] = 0.0;
		newphi[i] = newphi[size_ngrid-1] * exp(-1.0 * x_grid[i]);
		newphi[i] = amp/pow(decay+x_grid[i], 2.0);
		//printf("%f %f\n", phi_grid[i], newphi[i]);
		phi_grid[i] = newphi[i];
	}
}
//for (i=0; i<size_phigrid; i++) {
//	printf("(old phi, unchanged new phi) = %f, %f\n", phi_grid[i], newphi[i]);
//	if ( gammasq > TINY ) {
//		newphi[i] = weight*newphi[i] + phi_grid[i];
//	}
//	//if ( (gammasq > TINY) && (*convergence == 1) ) {
//	//	if (i==0) newphi[i] = 0.5*v_cut*v_cut;
//	//	else if (i>= size_ngrid-2) newphi[i] = 0.0;
//	//	else newphi[i] = weight * (newphi_gsl->data[i]) + (1.0-weight)*phi_grid[i];
//	//	if (newphi[i] < 0.0) gone_negative = 1;
//	//	if (gone_negative == 1) newphi[i] = 0.0;
//	//}
//	//if ( (i < size_ngrid) && (*convergence != 2) )
//		printf("(old phi, new phi) = %f, %f\n", phi_grid[i], newphi[i]);
//	phi_grid[i] = newphi[i];
//	
//}
free(phip);
free(phipp);
free(phipg);
free(phippg);
free(flux);
free(newphi);
free(deltanewphi);
free(ne);
free(gg);
//free(ni_DS_model);

for (i=0; i<sizemumu; i++) {
	 free(distfunc_iprime[i]);
}
free(distfunc_iprime);


free(altphipp);
free(altphipsqovtwo);
free(altphip);
//gsl_vector_free (newphi_gsl);

printf("av_dev = %f\n\n\n\n", av_dev);

clock_t end = clock(); // finds the end time of the computation
double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
printf("Module to obtain the new guess for the electrostatic potential ran in %f\n", jobtime);
return;
}


