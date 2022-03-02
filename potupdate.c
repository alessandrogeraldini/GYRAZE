// LAST SUBSTANTIAL MODIFICATION MADE 16 JAN 2022
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

void newvcut(double *v_cut, double v_cutDS, double mioverme, double u_i, double u_e, double current, double error_current, int *convergence, double weight) {
	double old_v_cut = *v_cut, u_etilde;
	u_etilde = u_e - sqrt(mioverme/M_PI)*0.5*exp(-0.5*old_v_cut*old_v_cut);
	//printf("Old v_cut = %f\n", *v_cut); //printf("u_etilde = %f\n", u_etilde);
	printf("target current = %f\n", current);
	printf("current = %f = u_i (%f) - u_e (%f)\n", u_i-u_e, u_i, u_e);
	if (*convergence == 0) {
		if ( (-current + u_i - u_etilde) < TINY ) {
			printf("WARNING: The log below would be negative, iteration avoided\n");
			*v_cut = old_v_cut;
		}
		else { 
			//*v_cut = (1.0-weight)*(old_v_cut) + weight*sqrt(-2.0*log(2.0*(-current + u_i - u_etilde)) - log(M_PI/mioverme) );
			*v_cut = sqrt(2.0*( (1.0-weight)*0.5*old_v_cut*old_v_cut + weight* ( 0.5*old_v_cut*old_v_cut + (current - u_i + u_e)/( sqrt(mioverme)*0.5*exp(-0.5*old_v_cut*old_v_cut) ) ) ) );
		}
		if ( (*v_cut*(*v_cut)*0.5 < 0.5*v_cutDS*v_cutDS ) )  {
			printf("WARNING: v_cutDS is larger than v_cut\n");
			*v_cut = v_cutDS;
		}
	}

	if ( (fabs(current - u_i + u_e)/u_i > error_current) ) *convergence = 0;
	else *convergence += 1;
}

void newguess(int *convergence, double Te, double *x_grid, double* ne_grid, double *ni_grid, double* phi_grid,int size_phigrid, int *psize_ngrid, double invgammasq, double v_cutDS, double pfac, double weight) {
	
clock_t begin = clock(); 
int size_ngrid = *psize_ngrid;
int i, j;
double Bohmshouldbe, phiW_impose, temp;
double *phip, *phipp, *phipp_red, *newphi, *gg;
double res = 0.0, reslimit, dev, dev_0, devbig;
double phi0, phip0, CC, pdec, deltaxsq;
//FILE *fout;

pdec = 2.0/(1.0-pfac);
printf("decay power is %f\n", pdec);

Bohmshouldbe = 2.0*Te*(ne_grid[1] - ne_grid[0])/((phi_grid[1] - phi_grid[0])*ni_grid[0]);
printf("Bohm integral should be %f\n", Bohmshouldbe);

/* Below we initialize all arrays that contain functions of position x with the correct size n */
gg = (double*)calloc(    size_phigrid,sizeof(double)); // gg now has correct size
phip = (double*)calloc(  size_phigrid,sizeof(double)); // phi now has correct size
phipp = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size
newphi = (double*)calloc(size_phigrid,sizeof(double)); // same as above 

//if (invgammasq < TINY) reslimit = 0.0005;
reslimit = 0.002;

//if (invgammasq > TINY) reslimit *= invgammasq;

printf("weight = %f\n", weight);

//if ((fout = fopen("phidataold.txt", "w")) == NULL)
//{
//	printf("Cannot open phidataold.txt");
//	exit(EXIT_FAILURE);
//}
//for (i=0; i<size_phigrid; i++)
//{
//	gg[i] = sqrt(x_grid[i]);
//	fprintf(fout, "%f %f\n", gg[i], phi_grid[i]);
//}
//fclose(fout); 

// Calculate first and second derivatives of potential
for (i=0; i<size_phigrid; i++) {	
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

// Calculate the residual of Poisson's/quasineutrality equation
res = 0.0;
devbig = dev = 0.0;
dev_0 = 0.0;
//if (invgammasq < TINY) {
//	fout = fopen("OUTPUT/iterationdata_MP.txt", "w");
//	fout = fopen("OUTPUT/iterationdata_DS.txt", "w");
//	if (fout == NULL) {
//		printf("Cannot open iterationdata.txt");
//		exit(EXIT_FAILURE);
//	}
//}
for (i=size_ngrid-1; i>=0; i--)  {
	//if (i==size_ngrid-1)	printf("x phi phipp phipp/gammasq ne ni err \n");
	//printf(fout, "%f %f %f %f %f %f %f\n", invgammasq, x_grid[i], phi_grid[i], phipp[i], ne_grid[i], ni_grid[i], (-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
	//fprintf(fout, "%f %f %f %f %f %f %f\n", invgammasq, x_grid[i], phi_grid[i], phipp[i], ne_grid[i], ni_grid[i], (-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
	if ( (i!=0) && (i!=size_ngrid-1) ) {
		dev = fabs((-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
		res += fabs((-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
		if (dev > devbig)  devbig = dev; 
		if (i==1) dev_0 = fabs((-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
	}
}
//fclose(fout);

res /= (size_ngrid-2);

deltaxsq = x_grid[1]*x_grid[1];

if ( invgammasq > TINY ) {
	phi0 = (ne_grid[size_ngrid-1] - ni_grid[size_ngrid-1])/(pow(phi_grid[size_ngrid-1], pfac));
	
	gsl_vector *newphi_gsl = gsl_vector_alloc (size_phigrid-2);
	gsl_matrix * m = gsl_matrix_alloc (size_phigrid-2, size_phigrid-2);

	for (i = 0; i < size_phigrid-2; i++) {
		for (j = 0; j < size_phigrid-2; j++) {
			if (i == j) 
				gsl_matrix_set (m, i, j, -2.0*invgammasq/deltaxsq - 1.0*exp(phi_grid[i+1]/Te)); 
			else if ( (i==j+1) || (i==j-1) )
				gsl_matrix_set (m, i, j, 1.0*invgammasq/deltaxsq);
			else 
				gsl_matrix_set (m, i, j, 0.0);
		}
	}

	phipp_red = malloc((size_phigrid-2)*sizeof(double));
	for (i=0;i<size_phigrid-2; i++) {
		if (i< size_ngrid-2)
			phipp_red[i] = (ne_grid[i+1] - ni_grid[i+1]) - phi_grid[i+1]* exp(phi_grid[i+1]/Te);
		else phipp_red[i] = phi0*pow(phi_grid[i+1], pfac) - phi_grid[i+1]* exp(phi_grid[i+1]/Te);
	}
	phiW_impose = (-0.5*v_cutDS*v_cutDS); 
	//printf("phi0 = %f\tCC=%f\n", phi0, CC);
	//printf("phiW_impose = %f\n", phiW_impose);
	phipp_red[0] -= (phiW_impose*invgammasq/deltaxsq);
	//phi0 = phip[size_ngrid-1]/(pdec*pow(x_grid[size_ngrid-1] + CC, pdec-1));
	//phipp_red[size_phigrid-1] -= ( phi0*pow(x_grid[size_phigrid-1] + CC, pdec)*invgammasq/deltaxsq );

	gsl_vector_view phipp_gsl
	    = gsl_vector_view_array (phipp_red, size_phigrid-2);

	int s;
	gsl_permutation * p = gsl_permutation_alloc (size_phigrid-2);
	clock_t t1 = clock(); // finds the end time of the computation
	gsl_linalg_LU_decomp (m, p, &s);
	clock_t t2 = clock(); // finds the end time of the computation
	double decomptime  = (double)(t2 - t1) / CLOCKS_PER_SEC;
	printf("decomptime = %f\n", decomptime);
	gsl_linalg_LU_solve (m, p, &phipp_gsl.vector, newphi_gsl);

	//printf ("newphi_gsl = \n");
	//printf("%f\n", -0.5*v_cutDS*v_cutDS);
	//gsl_vector_fprintf (stdout, newphi_gsl, "%g");

	newphi[0] = - 0.5*v_cutDS*v_cutDS;
	printf("%f\n", newphi[0]);
	for (i=0; i<size_ngrid; i++) {
		temp = gsl_vector_get(newphi_gsl, i);
		newphi[i+1] = temp;
		printf("%f\n", newphi[i+1]);
	}
	phip0 = (newphi[size_ngrid-1] - newphi[size_ngrid-2])/(x_grid[size_ngrid-1] - x_grid[size_ngrid-2]);
	CC = pdec*newphi[size_ngrid-1]/phip0 - x_grid[size_ngrid-1];
	if (CC < -x_grid[size_ngrid-1]) CC = 0.0;
	phi0 =  newphi[size_ngrid-1] /pow(x_grid[size_ngrid-1]+CC, pdec) ;
	printf("CC = %f and phi0 = %f\n", CC, phi0);
	for (i=size_ngrid; i<size_phigrid; i++) {
		newphi[i] = phi0*pow(x_grid[i]+CC, pdec);
		printf("%f\n", newphi[i]);
	}

	//phi0 = (ne_grid[size_ngrid-1] - ni_grid[size_ngrid-1])/pow(phi_grid[size_ngrid-1], pfac);
	//CC = 0.0;
	//newphi[size_phigrid-1] = phi0*pow(x_grid[size_phigrid-1] + CC, pdec);

	printf("%f\n", newphi[size_phigrid-1]);

	gsl_permutation_free (p);
	gsl_matrix_free (m);
	gsl_vector_free (newphi_gsl);
	free(phipp_red);
}
else {
	for (i=0; i< size_phigrid; i++) {
		if (i< size_ngrid-1)
			//newphi[i] = Te * log( ni_grid[i] - ne_grid[i] + exp(phi_grid[i]/Te)); // + phi_grid[i];
			//newphi[i] = phi_grid[i] + 1.0*Te* ( ni_grid[i] - ne_grid[i] ); // + exp(phi_grid[i]/Te)); // + phi_grid[i];
			newphi[i] = Te * ( ni_grid[i] - ne_grid[i] )*exp(-phi_grid[i]/Te) + phi_grid[i];
			//newphi[i] = (1.0/ne_grid[i])*Te * ( ni_grid[i] - ne_grid[i] ) + phi_grid[i];
		else if (i== size_ngrid-1) {
			//newphi[i] = Te * log( ni_grid[i] - ne_grid[i] + exp(phi_grid[i]/Te)); // + phi_grid[i];
			//newphi[i] = phi_grid[i] + 1.0*Te* ( ni_grid[i] - ne_grid[i] ); // + exp(phi_grid[i]/Te)); // + phi_grid[i];
			newphi[i] = Te * ( ni_grid[i] - ne_grid[i] )*exp(-phi_grid[i]/Te) + phi_grid[i];
			//newphi[i] = (1.0/ne_grid[i])*Te * ( ni_grid[i] - ne_grid[i] ) + phi_grid[i];
			phip0 = (newphi[i] - newphi[i-2])/(x_grid[i] - x_grid[i-2]);
			CC = pdec*newphi[i-1]/phip0 - x_grid[i-1];
			//CC = 0.0;
			phi0 =  newphi[i-1] /pow(x_grid[i-1]+CC, pdec) ;
			//phi0 =  phip0/(pdec*pow(x_grid[i-1]+CC, pdec-1)) ;
			printf("CC = %f and phi0 = %f\n", CC, phi0);
		}
		else 
			newphi[i] = phi0*pow(x_grid[i]+CC, pdec);
	// + newphi[size_ngrid-2] - phi0*pow(x_grid[size_ngrid-2], pdec);
		//printf("newphi[%d] = %f\n", i, newphi[i]);
	}
}

printf("res = %f, devbig = %f, dev_0 = %f\n", res, devbig, dev_0);

printf("conditions %d %d %d\n", res < reslimit, devbig < 5.0*reslimit, dev_0 < 5.0*reslimit);
if ( res < reslimit && devbig < 5.0*reslimit && dev_0 < 5.0*reslimit )  {
	*convergence += 1 ;
}
else {
	//if ((fout = fopen("phidata.txt", "w")) == NULL)
	//{
	//	printf("Cannot open phidata.txt");
	//	exit(EXIT_FAILURE);
	//}
	//for (i=0; i<size_phigrid; i++)
	//{	
	//	fprintf(fout, "%f %f\n", gg[i], newphi[i]);
	//}
	//fclose(fout);
	*convergence = 0;
}

//if (*convergence == 0) {
for (i=size_phigrid-1; i>=0; i--) {
	newphi[i] = weight*newphi[i] + (1.0-weight)*phi_grid[i] ; 
	phi_grid[i] = newphi[i];
}
//}
free(phip);
free(phipp);
free(newphi);
free(gg);

clock_t end = clock(); // finds the end time of the computation
double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
printf("Module to obtain the new guess for the electrostatic potential ran in %f\n", jobtime);
return;
}
