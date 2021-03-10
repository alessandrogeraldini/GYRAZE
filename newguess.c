// LAST SUBSTANTIAL MODIFICATION MADE 13 NOV 2018
/* This code calculates the next electrostatic potential guess in the iteration to obtain the self-consistent magnetic presheath electrostatic potential profile */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "Headerfiles/linetodata.h"
#include "Headerfiles/newguess.h"
#include "Headerfiles/bilin_int.h"
#include "Headerfiles/cspline.h"
#include "Headerfiles/makelookup.h"
#define M_PI acos(-1)

//struct spline
//{
//	double* a;
//	double* b;
//	double* c;
//	double* d;
//
//};

int newguess(int convergence, int problem, double* ne_grid, double* phi_grid,int p_size,double* ne_new_grid, double * phi_new_grid, double* ne_cut_grid, double* phi_cut_grid, int cut_points,int *use_new, double u_e, double *v_cut)
{//DECLARATIONS
	clock_t begin = clock(); // Finds the start time of the computation
	int debug = 0, lookup_inversion = 0, upgraded = 1;
	int row, col, i, ii, sizemumu, sizeUU, j, L1=0, n, cont = 1, h;
	char line[200000], *string;
	double *storevals, *gg, *xx, *ni, **FF, *mumu, *UU, k32, k32denom, k32num, k1DK, *Tevect, *k1GKvect; 
	double *phipg, *phippg, *phip, *phipp;;
	double k1DKnum1old, densinf, fluxinf, k1GK, **FFprime, k32denom1, k1DKnum1, k1DKnum, *newphi, k32denom1old, fluxinf1old, fluxinf1;
	double densinf1old, densinf1;
	double CC, *flux, *ne, k1DKthres = 0.3, Fprime, Fprimeold, F, Fold, Te, Teor1, fluxinfintgrdold, fluxinfintgrd;
	double phi1, alpha, weight, U, u, du, *phi;
	double smoothweight = 2.0/3.0;
	double sig=0.0, siglimit, siglimit1, dev, devbig;
	double random_value;
	double negradinf, dev_0;
	double u_etilde, u_i, current = 0.0, mioverme = 1836.0, set_current = 1, old_v_cut;
	*use_new = 0;


	FILE *distfile, *densfile, *Umufile, *fout2, *input, *fp, *fpotold;
	if (problem == 1) 
	{ 	
		convergence = 0; 
	}
if ((fp = fopen("phidata.txt", "r")) == NULL)
{// Check for presence of file
	printf("Cannot open %s\n", "phidata.txt");
	exit(EXIT_FAILURE); }
/* The while loop counts the lines in the file to determine the size of the arrays to be created */
row=0;
while(fgets(line, 20000, fp) != NULL)
{ // Count the number of rows in the distfuncin.txt file
	row += 1; }
n=row;
row=0;
/* Below we initialize all arrays that contain functions of position x with the correct size n */
ne = (double*)calloc(n,sizeof(double)); // phi now has correct size
ni = (double*)calloc(n,sizeof(double)); // phi now has correct size
phi = (double*)calloc(n,sizeof(double)); // phi now has correct size
phip = (double*)calloc(n,sizeof(double)); // phi now has correct size
phipg = (double*)calloc(n,sizeof(double)); // phi now has correct size
phipp = (double*)calloc(n,sizeof(double)); // phi now has correct size
phippg = (double*)calloc(n,sizeof(double)); // phi now has correct size
gg = (double*)calloc(n,sizeof(double)); // gg = sqrt(x) has correct size
xx = (double*)calloc(n,sizeof(double)); // xx = x has correct size
newphi = (double*)calloc(n,sizeof(double)); // same as above 
flux = (double*)calloc(n,sizeof(double)); // same as above 
row = 0;
rewind(fp);
while (fgets(line, 2000, fp) != NULL)
{	
	//Fill phi and gg arrays
	//string = malloc(strlen(line)*sizeof(char));
	//string = line;
	//storevals = linetodata(string, strlen(string), &col);
	storevals = linetodata(line, strlen(line), &col);
	gg[row] = *storevals;
	phi[row] = *(storevals+1);
	xx[row] = gg[row]*gg[row];
	if (debug == 1)
	{	printf("at row = %d, phi = %f, gg = %f, xx = %f\n", row, phi[row], gg[row], xx[row]); }
	row += 1; 
}
fclose(fp);
if ((input = fopen("inputfile.txt", "r")) == NULL)
{ // Check for presence of file
       	printf("Cannot open %s\n", "inputfile.txt");
        exit(1); }
i=0;
while (fgets(line, 20, input) != NULL)
{	
	//string = malloc(strlen(line)*sizeof(char));
    //string = line;
    //printf("The string is %s\n", string);
	//storevals = linetodata(string, strlen(string), &col);
	storevals = linetodata(line, strlen(line), &col);

	if (i==0)
	{	alpha = *storevals; }
	else if (i==1)
	{	Te = *storevals; }
	i += 1; 
}
fclose(input);
if (Te>1.001)
{	
	Teor1 = Te;
}
else 
{	
	Teor1 = 1.0; 
}
siglimit1 = 0.03;//0.03
siglimit = 0.009;//0.007
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
	weight = 0.8*random_value/pow(Teor1, 1.0); 
}
else
{	
	weight = 0.4*random_value/pow(Teor1, 1.0); 
}
if ((Umufile = fopen("Umufile.txt", "r")) == NULL)
{
	printf("cannot open file %s\n", "Umufile.txt");
	exit(-1);
}
/* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and assigning the values of the first line to mumu, second to UU  */
row=0;
while (fgets(line, 200000, Umufile) != NULL)
{	
	//string = malloc(strlen(line)*sizeof(char));
	//string = line;
	//storevals = linetodata(string, strlen(string), &col);
	storevals = linetodata(line, strlen(line), &col);
	if (row == 0)
	{	mumu = storevals;
		sizemumu = col; }
	else
	{	UU = storevals;
		sizeUU = col; }
	row += 1; 
}
row=0;
fclose(Umufile);
/* This part of the code extracts the distribution function from a file. We import the distribution function into a 2 dimensional array, F(mu,U). */
if ((distfile = fopen("distfuncin.txt", "r")) == NULL)
{	printf("cannot open file %s\n", "distfuncin.txt");
	exit(-1); }
while(fgets(line, 20000, distfile) != NULL)
{ // Count the number of rows in the distfuncin.txt file
	row += 1; }
// Allocate the right amount of memory to the distribution function pointer
FF = (double**)calloc(row,sizeof(double*));
// The number of rows is also the the size of the array in UU
row = 0; // Set number of rows counter to zero again
fclose(distfile); // Close file
if ((distfile = fopen("distfuncin.txt", "r")) == NULL)
{	// Check file exists and can be opened
	printf("cannot open file %s\n", "distfuncin.txt");
	exit(-1); }
while (fgets(line, 200000, distfile) != NULL)
{ 
	/* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and counting the columns. Once the columns are counted we allocate memory to FF and assign each number to a different column of FF (for a fixed row). The  while loop then jumps to a new line and repeats the process. */
	//string = malloc(strlen(line)*sizeof(char));
	//string = line;
	//printf("The string is %s\n", string);
	//storevals = linetodata(string, strlen(string), &col);
	storevals = linetodata(line, strlen(line), &col);

	FF[row] = (double*)calloc(col,sizeof(double));
	FF[row] = storevals;
	row +=1; 
}
fclose(distfile);
row = 0;
if ((densfile = fopen("PostProcessing/niout.txt", "r")) == NULL)
{	
	printf("Cannot open file %s\n", "niout.txt"); 
	exit(-1); 
}
while(fgets(line, 20000, densfile) != NULL)
{ // Count the number of rows in the distfuncin.txt file
	row += 1; }
ni = (double*)calloc(row,sizeof(double));
row = 0;
rewind(densfile);
while (fgets(line, 20000, densfile) != NULL)
{
	/* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and counting the columns. Once the columns are counted we allocate memory to FF and assign each number to a different column of FF (for a fixed row). The  while loop then jumps to a new line and repeats the process. */
	//string = malloc(strlen(line)*sizeof(char));
	//string = line;
	//storevals = linetodata(string, strlen(string), &col);
	storevals = linetodata(line, strlen(line), &col);

	ni[row] = *(storevals+1);
	//printf("ni[%d] = %f\n", row, ni[row]);
	row += 1; 
}
fclose(densfile);
L1 = row;
printf("L1 = %d\n", L1);
if ((fpotold = fopen("PostProcessing/phidataold.txt", "w")) == NULL)
{
	printf("Cannot open phidataold.txt");
	exit(EXIT_FAILURE);
}
for (i=0; i<n; i++)
{
	fprintf(fpotold, "%f %f\n", gg[i], phi[i]);
}
fclose(fpotold); 
// INTEGRALS OF DISTRIBUTION FUNCTION AT INFINITY
/* In the section below I will evaluate the integrals that are necessary at the Chodura sheath entrance. They will be used in the potential iteration to extrapolate the electrostatic potential at large values of x which inevitably the code cannot solve for (boundary condition at infinity). */
FILE *fpk1;
if ((fpk1 = fopen("k1gyrokinetic.txt", "r")) == NULL)
{	// Check for presence of file
	printf("Cannot open %s\n", "k1gyrokinetic.txt");
	exit(EXIT_FAILURE); 
}
/* The while loop counts the lines in the file to determine the size of the arrays to be created */
row=0;
while(fgets(line, 20000, fpk1) != NULL)
{ // Count the number of rows in the distfuncin.txt file
	row += 1; }
k1GKvect = (double*)calloc(row,sizeof(double));
Tevect = (double*)calloc(row,sizeof(double));
rewind(fpk1);
row=0;
while (fgets(line, 2000, fpk1) != NULL)
{
	//string = malloc(strlen(line)*sizeof(char));
	//string = line;
	//storevals = linetodata(string, strlen(string), &col);
	storevals = linetodata(line, strlen(line), &col);

	Tevect[row] = *storevals;
	k1GKvect[row] = *(storevals+1);
	//printf("i=%d, Tevect=%f and k1GKvect=%f\n", i, Tevect[i], k1GKvect[i]);
	row += 1;
}
fclose(fpk1); // close file


//printf("k1GK = %f\n", k1GK);

// Note: first index of FF and FFprime is mu, second one is U

FFprime = (double**)calloc(sizemumu,sizeof(double*));
//printf("FF[30][29] = %f and FF[29][30] =%f \n", FF[30][29], FF[29][30]);


/* The loop below evaluates the derivative of F with respect to total energy U at U = mu
 */
for (i=0; i<sizemumu; i++)
{
	FFprime[i] = (double*)calloc(sizeUU,sizeof(double));
	for (j=0; j < sizeUU; j++)
	{
		if ((j==sizeUU-1) || (i==sizemumu-1))
		{
			FFprime[i][j] = 0.0;	
		}
		else
		{
			FFprime[i][j] = (FF[i][j+1] - FF[i][j])/(UU[j+1] - UU[j]); 
		}
	//FFprime[i] = (1/pow(M_PI, 1.5))*exp(-mumu[i]);
		//printf("FFprime = %f\n", FFprime[i][j]);
	}
}

k32num = 0.0;
for (j=0; j<sizemumu; j++)
{
	if (j!=0)
	{
		k32num += 0.5*0.5*(16.0*M_PI/3)*(FFprime[j][j]+FFprime[j-1][j-1])*(mumu[j] - mumu[j-1]);
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
		Fprime = bilin_interp(mumu[i], U, FFprime, mumu, UU, sizemumu, sizeUU, -1, -1);
		F = bilin_interp(mumu[i], U, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
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
if (upgraded == 0) {
if (phi[0] / Te < phi_grid[0])
{
	printf("##The potential exceeds the cutoff potentia##\n");
	makelookup(phi_new_grid, ne_new_grid, p_size, sqrt(-2 * phi[0]/Te) + pow(10.0,-10.0), 0.0, &u_e);
	//generate_cspline(spl_new_netophi, ne_new_grid, phi_new_grid, p_size);
	//generate_cspline(spl_new_phitone, phi_new_grid, ne_new_grid, p_size);
	printf("New phi cut off = %.5f\nNew ne cut off = %.5f\n", Te*phi_new_grid[0], ne_new_grid[0]);
	*use_new = 1;
}
if (*use_new == 1)
{
	if (ni[0] < ne_new_grid[0])
	{
		printf("Ion density lower than new electron density\nhence due to weighting the cutoff must be lowered\n");
		cont = 1;
		h = 0;
		while (cont == 1)
		{
			if (h == cut_points - 1)
			{
				printf("We've got a problem");
			}
			if ((ne_new_grid[0] + weight * (ni[0] - ne_new_grid[0]) <= ne_cut_grid[h]) && (ne_new_grid[0] + weight * (ni[0] - ne_new_grid[0]) > ne_cut_grid[h + 1]))
			{
				makelookup(phi_new_grid, ne_new_grid, p_size, sqrt(fabs(2.0*phi_cut_grid[h + 1])), sqrt(fabs(2.0*phi_cut_grid[h + 1] - 2.0*Te*phi_new_grid[0])), &u_e);
				//generate_cspline(spl_new_netophi, ne_new_grid, phi_new_grid, p_size);
				//generate_cspline(spl_new_phitone, phi_new_grid, ne_new_grid, p_size);
				*use_new = 1;
				printf("New phi cut off = %.5f\nNew ne cut off = %.5f\n", Te*phi_new_grid[0], ne_new_grid[0]);
				cont = 0;
			}
			h++;
		}
	}
}
else if (*use_new == 0) // I think I need to remove this now
{
	if (ni[0] < ne_grid[0])
	{
		printf("Ion density lower than new electron density\nhence due to weighting the cutoff must be lowered\n");
		cont = 1;
		h = 0;
		//printf("The weighted density is %.5f\n", ne_grid[0] + weight * (ni[0] - ne_grid[0]));
		while (cont == 1)
		{
			if (h == cut_points - 1)
			{
				printf("We've got a problem");
			}
			if ((ne_grid[0] + weight * (ni[0] - ne_grid[0]) <= ne_cut_grid[h]) && (ne_grid[0] + weight * (ni[0] - ne_grid[0]) > ne_cut_grid[h + 1]))
			{
				//printf("The reduced weighted density is %.5f\n", ne_cut_grid[h + 1]);
				makelookup(phi_new_grid, ne_new_grid, p_size, sqrt(fabs(2.0*phi_cut_grid[h + 1])), sqrt(fabs(2.0*phi_cut_grid[h + 1] - 2.0*Te*phi_new_grid[0])), &u_e);
				//generate_cspline(spl_new_netophi, ne_new_grid, phi_new_grid, p_size);
				//generate_cspline(spl_new_phitone, phi_new_grid, ne_new_grid, p_size);
				*use_new = 1;
				printf("New phi cut off = %.5f\nNew ne cut off = %.5f\n", Te*phi_new_grid[0],ne_new_grid[0]);
				cont = 0;
			}
			h++;
		}
	}
}
}
if (*use_new == 0)
{
	negradinf = (ne_grid[p_size - 3] - ne_grid[p_size - 1]) / (phi_grid[p_size - 3] - phi_grid[p_size - 1]);
}
else if (*use_new == 1)
{
	negradinf = (ne_new_grid[p_size - 3] - ne_new_grid[p_size - 1]) / (phi_new_grid[p_size - 3] - phi_new_grid[p_size - 1]);
}
printf("negradinf = %.5f%\n", negradinf);

printf("Tevect = %f\n", Tevect[0]);
i = 0;
if (Te <= 1)
{
	while (Tevect[i] < Te)
	{
		//	printf("Tevect = %f\n", Tevect[i]);
		i += 1;
	}
	k1GK = (k1GKvect[i] * (-Tevect[i - 1] + Te) + k1GKvect[i - 1] * (Tevect[i] - Te)) / (Tevect[i] - Tevect[i - 1]);

}
else if (Te > 1)
{
	k1GK = 0; //~~need to change this just didn't want the value looking correct in case someone thought it was~~
}

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

sig = 0.0;
devbig = dev = 0.0;
dev_0 = 0.0;
// The part that calculates the new electrostatic potential
for (i=0; i<L1; i++)
{
	if (*use_new == 0)
	{
		//ne[i] = call_spline(spl_phitone, phi_grid, phi[i] / Te, p_size,421);
		//if (phi[i]/Te < phi_grid[0]) ne[i] = ne_grid[0];
		//else
		//printf("%f %f \n", phi_grid[i], phi[i]/Te);
		//ne[i] = lin_interp(phi_grid, ne_grid, phi[i] / Te, p_size, 462);
		ne[i] = ne_grid[i];
		//printf("%f %f \n", ne_grid[i], ne[i]/Te);
		//printf("i=%d\nne[i] = %f\nphi[i]/Te = %f\nphi_grid[0] = %f\nne_grid[0] = %f\n", i, ne[i], phi[i]/Te, phi_grid[0], ne_grid[0]);
		if (ne[i] > 1)
		{
			printf("ne[i] is greater than one\n");
		}
		if (ne[i] < 0)
		{
			printf("ne[i] is less than zero\n");
		}
	}
	else if (*use_new == 1)
	{
		//ne[i] = call_spline(spl_new_phitone, phi_new_grid, phi[i] / Te, p_size,425);
		ne[i] = lin_interp(phi_new_grid, ne_new_grid, phi[i] / Te, p_size, 475);
		if (ne[i] > 1)
		{
			printf("ne[i] is greater than one\n");
		}
		if (ne[i] < 0)
		{
			printf("ne[i] is less than zero\n");
		}
	}
	dev = fabs(ni[i]/ne[i] - 1.0);
	if (dev > devbig)
	{	devbig = dev; }
	sig += pow((ni[i]/ne[i] - 1.0), 2.0);
	if (i==L1-1)
	{
		//printf("ni= %f, weight = %f, phi[i] = %f\n", ni[i], weight, phi[i]);
		//newphi[i] = phi[i] - weight*(ne[i] - ni[i]);
		if (problem == 0)
		{
			if (*use_new == 0)
			{
				//newphi[i] = Te * call_spline(spl_netophi, ne_grid, ne[i] + weight * (ni[i] - ne[i]), p_size, 439);
				if (lookup_inversion == 1)
					newphi[i] = Te * lin_interp(ne_grid, phi_grid, ne[i] + weight * (ni[i] - ne[i]), p_size, 498);
				else 
					newphi[i] = Te * log(weight * (ni[i] - ne[i] + exp(phi[i]/Te)) + (1.0-weight)*exp(phi[i]/Te)); 
			}
			else if (*use_new == 1)
			{
				//newphi[i] = Te * call_spline(spl_new_netophi, ne_new_grid, ne[i] + weight * (ni[i] - ne[i]), p_size, 443);
				if (lookup_inversion == 1)
					newphi[i] = Te * lin_interp(ne_new_grid, phi_new_grid, ne[i] + weight * (ni[i] - ne[i]), p_size, 503);
				else 
					newphi[i] = Te * log(weight * (ni[i] - ne[i] + exp(phi[i]/Te)) + (1.0-weight)*exp(phi[i]/Te)); 
			}
		}
		else
		{
			newphi[i] = 0.0;
		}
		//newphi[i] = Te*log(weight*ni[i] + (1-weight)*ne[i]);
		if (fabs(k1DK) < k1DKthres)
		{
			CC = sqrt(20.0/k32)/pow(-0.5*newphi[i], 0.25) - gg[i]*gg[i];
			printf("C_3/2 = %f, newphi[i] = %f\n", CC, newphi[i]);
		}
		else if (fabs(k1DK) > k1DKthres)
		{
			printf("~~~Chodura Oversatistied~~~\n");
			//phi1 = newphi[i]*exp(-k1GK*xx[i]);
			phi1 = newphi[i] * exp(k1DK * xx[i]);
		}
	}
	else 
	{
		// This is where the iteration happens
		//newphi[i] = phi[i] - weight*(ne[i] - ni[i]);
		if (problem == 0)
		{
			//printf("%.5f\n", ne[i] + weight * (ni[i] - ne[i]));
			if (*use_new == 0)
			{
				//newphi[i] = Te * call_spline(spl_netophi, ne_grid, ne[i] + weight * (ni[i] - ne[i]), p_size, 469);
				if (lookup_inversion == 1)
					newphi[i] = Te * lin_interp(ne_grid, phi_grid, ne[i] + weight * (ni[i] - ne[i]), p_size, 533);
				else 
					newphi[i] = Te * log(weight * (ni[i] - ne[i] + exp(phi[i]/Te)) + (1.0-weight)*exp(phi[i]/Te)); 
			}
			if (*use_new == 1)
			{
				//printf("%d\n", i);
				//newphi[i] = Te * call_spline(spl_new_netophi, ne_new_grid, ne[i] + weight * (ni[i] - ne[i]), p_size, 473);
				if (lookup_inversion == 1)
					newphi[i] = Te * lin_interp(ne_new_grid, phi_new_grid, ne[i] + weight * (ni[i] - ne[i]), p_size, 539);
				else
					newphi[i] = Te * log(weight * (ni[i] - ne[i] + exp(phi[i]/Te)) + (1.0-weight)*exp(phi[i]/Te)); 
			}
			
		}
		else
		{
			newphi[i] = 0.0; 
		}
		//newphi[i] = Te*log(weight*ni[i] + (1-weight)*ne[i]); (exactly the same as the above)
		//printf("ni= %f\n", ni[i]);
		//printf("newphi[i] = %f\n", newphi[i]);
	}
}
for (i=L1; i<n; i++)
{
	if (fabs(k1DK) < k1DKthres)
	{
		newphi[i] = -2.0*(400.0/pow(k32, 2.0))/pow(gg[i]*gg[i] + CC, 4.0);
	}
	else if (fabs(k1DK) > k1DKthres)
	{
		//newphi[i] = phi1*exp(k1GK*xx[i]); // k1GK should be negative
		newphi[i] = phi1 * exp(-k1DK * xx[i]);
		//printf("phi = %f, newphi[i] = %f\n", phi[i], newphi[i]);
	}
	if (*use_new == 0)
	{
		//ne[i] = call_spline(spl_phitone, phi_grid, phi[i] / Te, p_size, 500);
		printf("%f\n", phi[i]/Te);
		//ne[i] = lin_interp(phi_grid, ne_grid, phi[i] / Te, p_size, 568);
		ne[i] = ne_grid[i];
	}
	else if (*use_new == 1)
	{
		//ne[i] = call_spline(spl_new_phitone, phi_new_grid, phi[i] / Te, p_size, 504);
		ne[i] = lin_interp(phi_new_grid, ne_new_grid, phi[i] / Te, p_size, 574);
	}
	
}
sig /= L1;
sig = sqrt(sig);
dev_0 = fabs((ne[0] / ni[0]) - 1.0);
printf("phi[0] is %f\n", phi[0]);
if (sig < siglimit && devbig < 5.0*siglimit && dev_0 < siglimit1) // || devbig < 0.04)
{
	convergence = 2;
	printf("convergence = %d, sig = %f\n", convergence, sig);
}
else
{
	// SMOOTHING OF PHI
	// Now we take derivatives of phi
	for (i=0; i<n; i++)
	{	
		if (i != n-1)
		{	phip[i] = (newphi[i+1] - newphi[i])/(gg[i+1]*gg[i+1] - gg[i]*gg[i]);
			phipg[i] = (newphi[i+1] - newphi[i])/(gg[i+1] - gg[i]);
			if (i!=0)
			{	phipp[i] = 2.0*(phip[i] - phip[i-1])/(gg[i+1]*gg[i+1] - gg[i-1]*gg[i-1]);
				phippg[i] = 2.0*(phipg[i] - phipg[i-1])/(gg[i+1] - gg[i-1]); }	
			if (i==1)
			{	phipp[0] = phipp[1];
				phippg[0] = phippg[1]; } }
		else
		{	
			phip[i] = phip[i-1];
			phipg[i] = phipg[i-1];
			phipp[i] = phipp[i-1];
			phippg[i] = phippg[i-1];
		} 
	}
	/* Smooth second derivative of phi wrt sqrt(x)
	*/
	if (convergence == 0)
	{ 	
		for (i=0; i<L1; i++)
		{	newphi[i] = phi[L1];
			phipg[i] = phipg[L1];
			if (i!=0)
			{	
				phippg[i] = ((1.0-smoothweight)*phippg[i-1]*(gg[i+1] - gg[i]) + smoothweight*phippg[i]*(gg[i+1] - gg[i-1]) + (1.0-smoothweight)*phippg[i+1]*(gg[i] - gg[i-1]))/((gg[i+1] - gg[i-1]));
				//phipp[i] = ((1.0-smoothweight)*phipp[i-1]*(xx[i+1] - xx[i]) + smoothweight*phipp[i]*(xx[i+1] - xx[i-1]) + (1.0-smoothweight)*phipp[i+1]*(xx[i] - xx[i-1]))/((xx[i+1] - xx[i-1])); 
			} 
		}
		// Re-integrate smoothed second derivative to first
		for (i=L1-1; i>=0; i--)
		{	for (ii=i; ii>=0; ii--)
			{	phipg[ii] -= (gg[i+1] - gg[i])*(phippg[i+1] + phippg[i])*0.5; 
			} 
		}  
				//phip[ii] -= (xx[i+1] - xx[i])*(phipp[i+1] + phipp[i])*0.5; 
		// Re-integrate first derivative to phi
		for (i=L1-1; i>=0; i--)
		{	for (ii=i; ii>=0; ii--)
			{	
				newphi[ii] -= (gg[i+1] - gg[i])*(phipg[i+1] + phipg[i])*0.5; 
			} 
		} 
	}
				//newphi[ii] -= (xx[i+1] - xx[i])*(phip[i+1] + phip[i])*0.5; 
	if ((fout2 = fopen("phidata.txt", "w")) == NULL)
	{
		printf("Cannot open phidata.txt");
		exit(EXIT_FAILURE);
	}
	for (i=0; i<n; i++)
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
free(ne);
free(ni);
free(phi);
free(phip);
free(phipp);
free(phipg);
free(phippg);
free(gg);
free(xx);
free(flux);
for (int w = 0;w < sizemumu;w++)
{
	free(FFprime[w]);

}
free(FFprime);
for (int m = 0;m < sizeUU;m++)
{
	free(FF[m]);

}
free(FF);

u_i = fluxinf/alpha;
if (set_current == 1) {
	printf("Old v_cut = %f\n", *v_cut);
	old_v_cut = *v_cut;
	u_etilde = u_e - sqrt(mioverme/M_PI)*0.5*exp(-0.5*old_v_cut*old_v_cut);
	printf("u_etilde = %f\n", u_etilde);
	printf("u_e = %f\n", u_e);
	printf("u_i = %f\n", u_i);
	printf("current = %f\n", current);
	*v_cut = 0.8*(*v_cut) + 0.2*sqrt(-2.0*log(2.0*(-current + u_i - u_etilde)) - log(M_PI/mioverme) );
	//printf("New v_cut = %f\n", *v_cut);
	if ( (*v_cut*(*v_cut)*0.5 < - phi_grid[0]) )  *v_cut = sqrt(- 2.0*phi_grid[0]);
	printf("New v_cut = %f\n", *v_cut);
	printf("New phi_cut = %f\n", 0.5*(*v_cut)*(*v_cut));
	//if ( (fabs(old_v_cut - (*v_cut)) > 0.02) && (convergence == 2) ) convergence = 1;
	if ( (fabs(current - u_i + u_e)/u_i > 0.008) && (convergence == 2) ) convergence = 1;
}

clock_t end = clock(); // finds the end time of the computation
double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
printf("Module to obtain the new guess for the electrostatic potential ran in %f\n", jobtime);
return convergence;
}
