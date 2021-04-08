
// LAST SUBSTANTIAL MODIFICATION MADE 13 NOV 2018
/* This code calculates the next electrostatic potential guess in the iteration to obtain the self-consistent magnetic presheath electrostatic potential profile */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mps.h"

int newguess(int convergence, int problem, double Te, double alpha, double *x_grid, double* ne_grid, double *ni, double* phi_grid,int size_phigrid, int size_ngrid, double u_e, double *v_cut, int quasineutral, double **distfunc_i, double *UU, double *mumu, int sizeUU, int sizemumu) {
//DECLARATIONS
	clock_t begin = clock(); // Finds the start time of the computation
	int row, col=0, i, ii, j;
	char line[20000];
	int newmethod = 1;
	double newweight = 0.2, newphiprime;
	double *storevals, k32, k32denom, k32num, k1DK, *Tevect, *k1GKvect; 
	double *phipg, *phippg, *phip, *phipp;;
	double k1DKnum1old, densinf, fluxinf, k1GK, **distfunc_iprime, k32denom1, k1DKnum1, k1DKnum, *newphi, k32denom1old, fluxinf1old, fluxinf1;
	double densinf1old, densinf1, *gg;
	double CC, *flux, *ne, k1DKthres = 0.3, Fprime, Fprimeold, F, Fold, Teor1, fluxinfintgrdold, fluxinfintgrd;
	double phi1, weight, U, u, du;
	double smoothweight = 2.0/3.0;
	double sig=0.0, siglimit, siglimit1, dev, devbig;
	double random_value;
	double negradinf, dev_0;
	double u_etilde, u_i, current = 0.0, mioverme = 1836.0, set_current = 1, old_v_cut;
	FILE *fout2, *fpotold;

printf("size_phigrid = %d\nsize_ngrid = %d\n", size_phigrid, size_ngrid);
	if (problem == 1) 
	{ 	
		convergence = 0; 
	}
/* Below we initialize all arrays that contain functions of position x with the correct size n */
ne = (double*)calloc(    size_phigrid,sizeof(double)); // phi now has correct size
gg = (double*)calloc(    size_phigrid,sizeof(double)); // gg now has correct size
phip = (double*)calloc(  size_phigrid,sizeof(double)); // phi now has correct size
phipg = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size
phipp = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size
phippg = (double*)calloc(size_phigrid,sizeof(double)); // phi now has correct size
newphi = (double*)calloc(size_phigrid,sizeof(double)); // same as above 
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
siglimit1 = 0.03;//0.03
siglimit = 0.002;//0.007
if (Te>2.1)
{
	siglimit1 = 0.012;
	siglimit = 0.012;
}
// Assign value to weight depending on whether convergence flag 

random_value = (double) rand()/RAND_MAX;//float in range 0 to 1
printf("The random number between 0 and 1 is %f\n", random_value);
if (convergence == 0)
{	
	weight = 0.2 + 0.4*random_value/pow(Teor1, 1.0); 
	weight = 0.4;
}
else
{	
	weight = 0.1 + 0.2*random_value/pow(Teor1, 1.0); 
	weight = 0.4;
}


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
//printf("YOLO\n");
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
//printf("YOLO\n");


//printf("k1GK = %f\n", k1GK);

// Note: first index of distfunc_i and distfunc_iprime is mu, second one is U

distfunc_iprime = (double**)calloc(sizemumu,sizeof(double*));
//printf("distfunc_i[30][29] = %f and distfunc_i[29][30] =%f \n", distfunc_i[30][29], distfunc_i[29][30]);


/* The loop below evaluates the derivative of F with respect to total energy U at U = mu
 */
for (i=0; i<sizemumu; i++)
{
	distfunc_iprime[i] = (double*)calloc(sizeUU,sizeof(double));
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
	//distfunc_iprime[i] = (1/pow(M_PI, 1.5))*exp(-mumu[i]);
		//printf("distfunc_iprime = %f\n", distfunc_iprime[i][j]);
	}
}

k32num = 0.0;
for (j=0; j<sizemumu; j++)
{
	if (j!=0)
	{
		k32num += 0.5*0.5*(16.0*M_PI/3)*(distfunc_iprime[j][j]+distfunc_iprime[j-1][j-1])*(mumu[j] - mumu[j-1]);
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
		Fprime = bilin_interp(mumu[i], U, distfunc_iprime, mumu, UU, sizemumu, sizeUU, -1, -1);
		F = bilin_interp(mumu[i], U, distfunc_i, mumu, UU, sizemumu, sizeUU, -1, -1);
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

printf("k32num = %f\n", k32num);
printf("k32 = %f\n", k32);
printf("k1DK = %f, k1DK threshold is %f\n", k1DK, k1DKthres);
printf("Te = %f\n", Te);
printf("fluxinf = %f\n", fluxinf);

printf("Integrals of distribution function at infinity DONE\n");

if (quasineutral == 0) {
	for (i=0; i<size_phigrid; i++)
	{	
		if (i != size_phigrid-1)
		{	
			phip[i] = (newphi[i+1] - newphi[i])/(x_grid[i+1] - x_grid[i]);
			phipg[i] = (newphi[i+1] - newphi[i])/(gg[i+1] - gg[i]);
			if (i!=0) {	
				phipp[i] = 2.0*(phip[i] - phip[i-1])/(x_grid[i+1] - x_grid[i-1]);
				phippg[i] = 2.0*(phipg[i] - phipg[i-1])/(gg[i+1] - gg[i-1]); 
			}	
			if (i==1)  {	
				phipp[0] = phipp[1];
				phippg[0] = phippg[1]; 
			} 
		}
		else {	
			phip[i] = phip[i-1];
			phipg[i] = phipg[i-1];
			phipp[i] = phipp[i-1];
			phippg[i] = phippg[i-1];
		} 
	}
}

sig = 0.0;
devbig = dev = 0.0;
dev_0 = 0.0;
// The part that calculates the new electrostatic potential
for (i=0; i<size_ngrid; i++)
{
	if (ne_grid[i] > 1) {
		printf("ERROR: ne_grid[i] is greater than one\n");
	}
	if (ne_grid[i] < 0) {
		printf("ERROR: ne_grid[i] is less than zero\n");
	}
	//printf("phi=%f\tne=%f\tni=%f\n", phi_grid[i], ne_grid[i], ni[i]);
	dev = fabs(ni[i]/ne_grid[i] - 1.0);
	if (dev > devbig)  devbig = dev; 
	if (quasineutral == 1) sig += pow((ni[i]/ne_grid[i] - 1.0), 2.0);
	else sig += pow((ne_grid[i] - phipp[i])/ni[i] - 1.0, 2.0);
	if (i==size_ngrid-1)
	{
		if (problem == 0)
		{
			printf("quasineutral = %d\n", quasineutral);
			if (quasineutral == 1) 
				newphi[i] = Te * log(weight * (ni[i] - ne_grid[i] + exp(phi_grid[i]/Te)) + (1.0-weight)*exp(phi_grid[i]/Te)); 
			else
				//newphi[i] = Te * log(weight * (ni[i] - ne_grid[i] + exp(phi_grid[i]/Te) - phipp[i]) + (1.0-weight)*exp(phi_grid[i]/Te)); 
				phipp[i] = ne_grid[i] - ni[i];
		}
		else
		{
			newphi[i] = 0.0;
		}
		if (fabs(k1DK) < k1DKthres)
		{
			if (quasineutral == 1) {
				CC = sqrt(20.0/k32)/pow(-0.5*newphi[i], 0.25) - x_grid[i];
				printf("C_3/2 = %f, newphi[i] = %f\n", CC, newphi[i]);
			}
		}
		else if (fabs(k1DK) > k1DKthres)
		{
			printf("k1DK = %f\n, k1DKthres=%f\n", k1DK, k1DKthres);
			printf("~~~Chodura Oversatistied~~~\n");
			//phi1 = newphi[i]*exp(-k1GK*x_grid[i]);
			phi1 = newphi[i] * exp(k1DK * x_grid[i]);
		}
	}
	else  {
		// This is where the iteration happens
		if (problem == 0) {
			if (quasineutral == 1) 
				newphi[i] = Te * log(weight * (ni[i] - ne_grid[i] + exp(phi_grid[i]/Te)  ) + (1.0-weight)*exp(phi_grid[i]/Te)); 
			else
				//newphi[i] = Te * log(weight * (ni[i] - ne_grid[i] + exp(phi_grid[i]/Te) - phipp[i]  ) + (1.0-weight)*exp(phi_grid[i]/Te)); 
				phipp[i] = ne_grid[i] - ni[i];
		}
		else {
			newphi[i] = 0.0; 
		}
	}
	//printf("phi = %f, newphi[i] = %f\n", phi_grid[i], newphi[i]);
}
for (i=size_ngrid; i<size_phigrid; i++)
{
	if (quasineutral == 1) {
		if (fabs(k1DK) < k1DKthres)
		{
			newphi[i] = -2.0*(400.0/pow(k32, 2.0))/pow(x_grid[i] + CC, 4.0);
		}
		else if (fabs(k1DK) > k1DKthres)
		{
			//newphi[i] = phi1*exp(k1GK*x_grid[i]); // k1GK should be negative
			newphi[i] = phi1 * exp(-k1DK * x_grid[i]);
		}
		//printf("phi = %f, newphi[i] = %f\n", phi_grid[i], newphi[i]);
	}
}
sig /= size_ngrid;
sig = sqrt(sig);
dev_0 = fabs((ne_grid[0] / ni[0]) - 1.0);
printf("phi_grid[0] is %f\n", phi_grid[0]);

if (sig < siglimit && devbig < 5.0*siglimit && dev_0 < siglimit1)  {
	convergence = 2;
	printf("convergence = %d, sig = %f, devbig = %f, dev_0 = %f\n", convergence, sig, devbig, dev_0);
}
else {
	// SMOOTHING OF PHI
	// Now we take derivatives of phi_grid
	for (i=0; i<size_phigrid; i++)
	{	
		if (i != size_phigrid-1)
		{	
			//phip[i] = (newphi[i+1] - newphi[i])/(x_grid[i+1] - x_grid[i]);
			phipg[i] = (newphi[i+1] - newphi[i])/(gg[i+1] - gg[i]);
			if (i!=0) {	
				//phipp[i] = 2.0*(phip[i] - phip[i-1])/(x_grid[i+1] - x_grid[i-1]);
				phippg[i] = 2.0*(phipg[i] - phipg[i-1])/(gg[i+1] - gg[i-1]); 
			}	
			if (i==1)  {	
				//phipp[0] = phipp[1];
				phippg[0] = phippg[1]; 
			} 
		}
		else {	
			//phip[i] = phip[i-1];
			phipg[i] = phipg[i-1];
			//phipp[i] = phipp[i-1];
			phippg[i] = phippg[i-1];
		} 
	}
	/* Smooth second derivative of phi wrt sqrt(x) */
	if (convergence == 0)
	{ 	
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
	}

	if (quasineutral == 0) {
		for (i=size_ngrid-1; i>=0; i--)
		{	for (ii=i; ii>=0; ii--)
			{	
				phip[ii] -= (x_grid[i+1] - x_grid[i])*(phipp[i+1] + phipp[i])*0.5; 
			} 
		}  
		for (i=size_ngrid-1; i>=0; i--)
		{	for (ii=i; ii>=0; ii--)
			{	
				newphi[ii] -= (x_grid[i+1] - x_grid[i])*(phip[i+1] + phip[i])*0.5; 
			} 
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

				//newphi[ii] -= (x_grid[i+1] - x_grid[i])*(phip[i+1] + phip[i])*0.5; 
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
	printf("convergence = %d, sig = %f, devbig = %f\n", convergence, sig, devbig);
	if (convergence == 0 && ( sig < siglimit1 ) && (dev_0 < siglimit1) ) //|| devbig < 0.2) )
	{
		convergence = 1; 
		printf("now don't smooth = %d, sig = %f\n", convergence, sig);
	}
	else if ((convergence == 1) && (sig > siglimit1))
	{
		convergence = 0;
	}
}
printf("After smoothing phi[0] is: %f\n", newphi[0]);
for (i=0; i<size_phigrid; i++) {
	phi_grid[i] = newphi[i];
	//printf("phi = %f\n", phi_grid[i]);
}
free(phip);
free(phipp);
free(phipg);
free(phippg);
printf("YOLO\n");
free(flux);
free(newphi);
//for (int w = 0;w < sizemumu;w++)
//{
//	free(distfunc_iprime[w]);
//
//}
//free(distfunc_iprime);
//for (int m = 0;m < sizeUU;m++)
//{
//	free(distfunc_i[m]);
//
//}
//free(distfunc_i);

u_i = fluxinf/alpha;
if (set_current == 1) {
	printf("Old v_cut = %f\n", *v_cut);
	old_v_cut = *v_cut;
	u_etilde = u_e - sqrt(mioverme/M_PI)*0.5*exp(-0.5*old_v_cut*old_v_cut);
	printf("u_etilde = %f\n", u_etilde);
	printf("u_e = %f\n", u_e);
	printf("u_i = %f\n", u_i);
	printf("current = %f\n", current);
	if ( (-current +u_i - u_etilde) < 0.0 ) {
		printf("The log below will be negative\n");
		*v_cut = old_v_cut;
	}
	else { 
		*v_cut = 0.5*(*v_cut) + 0.5*sqrt(-2.0*log(2.0*(-current + u_i - u_etilde)) - log(M_PI/mioverme) );
	}
	//if ((u_i - u_e - current)/u_i > 0.01) {
	//	*v_cut -= 0.01;	
	//}
	//*v_cut = sqrt(0.8*(*v_cut)*(*v_cut) - 0.2*2.0*log(2.0*(-current + u_i - u_etilde)) - log(M_PI/mioverme) );
	//printf("New v_cut = %f\n", *v_cut);
	if ( (*v_cut*(*v_cut)*0.5 < - phi_grid[0]) )  *v_cut = sqrt(- 2.0*phi_grid[0]);
	printf("New v_cut = %f\n", *v_cut);
	printf("New phi_cut = %f\n", 0.5*(*v_cut)*(*v_cut));
	//if ( (fabs(old_v_cut - (*v_cut)) > 0.02) && (convergence == 2) ) convergence = 1;
	if ( (fabs(current - u_i + u_e)/u_i > 0.002) && (convergence == 2) ) convergence = 1;
}
clock_t end = clock(); // finds the end time of the computation
double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
printf("Module to obtain the new guess for the electrostatic potential ran in %f\n", jobtime);
return convergence;
}
