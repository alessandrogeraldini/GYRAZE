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

void newvcut(double *v_cut, double v_cutDS, double mioverme, double TioverTe, double u_i, double u_e, double current, double error_current, double weight) {
	double old_v_cut = *v_cut;// u_etilde;
	//u_etilde = u_e - sqrt(mioverme/M_PI)*0.5*exp(-0.5*old_v_cut*old_v_cut);
	//*v_cut = (1.0-weight)*(old_v_cut) + weight*sqrt(-2.0*log(2.0*(-current + u_i - u_etilde)) - log(M_PI/mioverme) ); //OLD
	printf("v_cut before = %f\n", v_cutDS);
	*v_cut = sqrt( 2.0*( (1.0-weight)*0.5*old_v_cut*old_v_cut + weight* ( 0.5*old_v_cut*old_v_cut + (current - u_i + u_e)/( sqrt(mioverme/(TioverTe))*0.5*exp(-0.5*old_v_cut*old_v_cut) ) ) ) );
	if ( (*v_cut*(*v_cut)*0.5 < 0.5*v_cutDS*v_cutDS ) )  {
		printf("WARNING: v_cutDS > v_cut\n");
		*v_cut = v_cutDS;
	}
	printf("v_cut after = %f\n", v_cutDS);
}

void error_quasi(double *error, double *x_grid, double *ne_grid, double *ni_grid, double *phi_grid, int size_phigrid, int size_ngrid) {
	int i;
	double res = 0.0, dev, devbig;
	// Calculate the residual of quasineutrality equation
	res = 0.0;
	devbig = dev = 0.0;
	for (i=size_ngrid-1; i>=0; i--)  {
		if (i==size_ngrid-1)	printf("gamma^-2 x phi phipp ne ni err \n");
		printf("%f %f %f %f %f\n", x_grid[i], phi_grid[i], ne_grid[i], ni_grid[i], -ne_grid[i]/ni_grid[i] + 1.0);
		dev = fabs(-ne_grid[i]/ni_grid[i] + 1.0);
		res += fabs(-ne_grid[i]/ni_grid[i] + 1.0);
		if (dev > devbig)  devbig = dev; 
	}
	res /= (size_ngrid);
	error[0] = res;
	error[1] = devbig;
	return;
}

void newguess(double *x_grid, double** ne_grid, double **ni_grid, double** phi_grid, double *phiinf_grid, int size_phigrid, int *psize_ngrid, int size_ygrid, double pfac, double weight) {
	clock_t begin = clock(); 
	int size_ngrid = *psize_ngrid;
	int i, j, s;
	double phiW_impose, temp;
	double *phipp_red, *newphi;
	//double res = 0.0, reslimit, dev, dev_0, devbig;
	double phi0, phip0, CC, pdec, deltaxsq;
	//FILE *fout;

	//if (size_ngrid == size

	pdec = 2.0/(1.0-pfac);
	printf("φ'' ~ φ^%f gives φ ~ x^%f\n", pfac, pdec);

	/* Below we initialize all arrays that contain functions of position x with the correct size n */
	newphi = (double*)calloc(size_phigrid,sizeof(double)); // same as above 

	deltaxsq = x_grid[1]*x_grid[1];

	for (j=0; j < size_ygrid; j++) {
		for (i=0; i< size_phigrid; i++) {
			if (i< size_ngrid-1)
				newphi[j][i] = ( ni_grid[j][i] - ne_grid[j][i] )*exp(-phi_grid[j][i]) + phi_grid[j][i];
			else if (i== size_ngrid-1) {
				newphi[i] = ( ni_grid[j][i] - ne_grid[j][i] )*exp(-phi_grid[j][i]) + phi_grid[j][i];
				phip0 = (newphi[j][i] - newphi[j][i-2])/(x_grid[i] - x_grid[i-2]);
				CC = pdec*(newphi[j][i-1] - phiinf_grid[i-1])/phip0 - x_grid[i-1];
				//CC = 0.0;
				phi0 =  (newphi[j][i-1] - phiinf_grid[i-1])/pow(x_grid[i-1]+CC, pdec) ;
				//phi0 =  phip0/(pdec*pow(x_grid[i-1]+CC, pdec-1)) ;
				printf("CC = %f and phi0 = %f\n", CC, phi0);
			}
			else 
				newphi[j][i] = phiinf_grid[i] + phi0*pow(x_grid[i]+CC, pdec);
		}
	}

	//printf("res = %f, devbig = %f, dev_0 = %f\n", res, devbig, dev_0);
	//printf("conditions %d %d %d\n", res < reslimit, devbig < 5.0*reslimit, dev_0 < 5.0*reslimit);

	for (i=size_phigrid-1; i>=0; i--) {
		newphi[j][i] = weight*newphi[j][i] + (1.0-weight)*phi_grid[j][i] ; 
		phi_grid[j][i] = newphi[j][i];
	}
	free(newphi);

	clock_t end = clock(); // finds the end time of the computation
	double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Module to obtain the new guess for the electrostatic potential ran in %f\n", jobtime);
	return;
}
