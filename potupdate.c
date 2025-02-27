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

void newvcut(double *v_cut, double v_cutDS, double u_i, double u_e, double current, double error_current, double weight) {
	double old_v_cut = *v_cut;// u_etilde;
	//u_etilde = u_e - sqrt(mioverme/M_PI)*0.5*exp(-0.5*old_v_cut*old_v_cut);
	//*v_cut = (1.0-weight)*(old_v_cut) + weight*sqrt(-2.0*log(2.0*(-current + u_i - u_etilde)) - log(M_PI/mioverme) ); //OLD
	printf("phi_wall before = %f --> \t", 0.5*(*v_cut)*(*v_cut));
	*v_cut = sqrt( 2.0*( (1.0-weight)*0.5*old_v_cut*old_v_cut + weight* ( 0.5*old_v_cut*old_v_cut + (current - u_i + u_e)/( 0.5*exp(-0.5*old_v_cut*old_v_cut) ) ) ) );
	if ( (*v_cut*(*v_cut)*0.5 < 0.5*v_cutDS*v_cutDS ) )  {
		printf("WARNING: v_cutDS > v_cut\n");
		*v_cut = v_cutDS;
	}
	printf("v_cut after = %f\n", *v_cut);
	printf("phi_wall after = %f --> \t", 0.5*(*v_cut)*(*v_cut));
}

void error_Poisson(double *error, double *x_grid, double *ne_grid, double *ni_grid, double *nioverne, double *phi_grid, int size_phigrid, int size_ngrid, double invgammasq) {
	int i;
	double *phip, *phipp;
	double res = 0.0, dev, devbig;
	/* Below we initialize all arrays that contain functions of position x with the correct size n */
	phip = (double*)calloc(  size_phigrid,sizeof(double)); // phi now has correct size
	phipp = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size
	// Calculate first and second derivatives of potential
	for (i=0; i<size_phigrid; i++) {	
		if ( (i != size_phigrid-1) && (i!=0) ) {	
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
	for (i=size_ngrid-1; i>=0; i--)  {
		//if (i==size_ngrid-1)	printf("gamma^-2 x phi phipp ne ni err \n");
		//printf("%f %f %f %f %f %f %f\n", invgammasq, x_grid[i], phi_grid[i], phipp[i], ne_grid[i], ni_grid[i], (-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
		if ( (i!=0) && (i!=size_ngrid-1) ) {
			dev = fabs((-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
			res += fabs((-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
			//dev = fabs((-ne_grid[i] + phipp[i]*invgammasq) + ni_grid[i]);
			//res += fabs((-ne_grid[i] + phipp[i]*invgammasq) + ni_grid[i]);
			if (dev > devbig)  {
				devbig = dev; 
				//printf("index corresponding to largest error = %d\tposition =%f\n", i, x_grid[i]);
			}
		}
	}
	res /= (size_ngrid-2);
	error[0] = res;
	error[1] = devbig;
	free(phip); free(phipp);
	return;
}

void newguess(double *x_grid, double* ne_grid, double *ni_grid, double* phi_grid, int size_phigrid, int size_ngridin, double invgammasq, double v_cutDS, double pfac, double weight) {
//attempt_at_adapting_grid failed. It was an attempt at refining grid near x=0 using a maximum Delta x maximum Delta phi gridding
int i, j, s, attempt_at_adapting_grid=0, size_ngrid=size_ngridin;
double phiW_impose, temp;
double *new_x;
double *phipp_red, *newphi, *phipp;
//double res = 0.0, reslimit, dev, dev_0, devbig;
double phi0, phip0, CC, pdec, deltaxsq, deltax = 0.3, deltaphi = 0.1;
//FILE *fout;

pdec = 2.0/(1.0-pfac);
printf("φ'' ~ φ^%f gives φ ~ x^%f\n", pfac, pdec);
printf("0.5*v_cutDS*v_cutDS = %f\n", 0.5*v_cutDS*v_cutDS);

/* Below we initialize all arrays that contain functions of position x with the correct size n */
newphi = (double*)calloc(size_phigrid,sizeof(double)); // same as above 
//phip = (double*)calloc(  size_phigrid,sizeof(double)); // phi now has correct size
phipp = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size

// Calculate first and second derivatives of potential
for (i=0; i<size_phigrid; i++) {	
	if ( (i != size_phigrid-1) && (i!=0) )
	{	
		//phipp[i] = ((x_grid[i] - x_grid[i-1])/(x_grid[i+1] - x_grid[i-1]))*(phip[i+1] - phip[i])/(x_grid[i+1] - x_grid[i]);
		//+ ((x_grid[i+1] - x_grid[i])/(x_grid[i+1] - x_grid[i-1]))*(phip[i] - phip[i-1])/(x_grid[i] - x_grid[i-1]);
		//phip[i] = (phi_grid[i+1] - phi_grid[i]) / (x_grid[i+1] - x_grid[i]);
		phipp[i] = (phi_grid[i+1] - 2.0*phi_grid[i] + phi_grid[i-1]) / pow(x_grid[i+1] - x_grid[i], 2.0);
	}
	else if (i==size_phigrid-1) phipp[size_phigrid-1] = 0.0; //(phip[size_phigrid-1] - phip[size_phigrid-2])/(x_grid[size_phigrid-1] - x_grid[size_phigrid-2]);

	if (i==2) phipp[0] = phipp[1] - x_grid[1]*(phipp[2] - phipp[1])/(x_grid[2] - x_grid[1]);
	//if (i==0) phip[i] = (phi_grid[i+1] - phi_grid[i]) / (x_grid[i+1] - x_grid[i]);
	//if (i==size_phigrid-1) phip[i] = 0.0;
}

deltaxsq = x_grid[1]*x_grid[1];

if ( invgammasq > TINY ) {
	phi0 = (ne_grid[size_ngrid-1] - ni_grid[size_ngrid-1])/(pow(phi_grid[size_ngrid-1], pfac));
	gsl_vector *newphi_gsl = gsl_vector_alloc (size_phigrid-2);
	gsl_matrix * m = gsl_matrix_alloc (size_phigrid-2, size_phigrid-2);

	for (i = 0; i < size_phigrid-2; i++) {
		for (j = 0; j < size_phigrid-2; j++) {
			if (i == j) {
				gsl_matrix_set (m, i, j, -2.0*invgammasq/deltaxsq - 1.0*exp(phi_grid[i+1])); 
			}
			else if ( (i==j+1) || (i==j-1) ) {
				gsl_matrix_set (m, i, j, 1.0*invgammasq/deltaxsq);
			}
			else { 
				gsl_matrix_set (m, i, j, 0.0);
			}
		}
	}
	//phibar = ( - (0.5*v_cutDS*v_cutDS) - (1.0-weight)*phi_grid[0] )/weight;
	phipp_red = malloc((size_phigrid-2)*sizeof(double));
	//printf("size_phigrid = %d\nsize_ngrid=%d\n\n", size_phigrid, size_ngrid);
	for (i=0;i<size_phigrid-2; i++) {
		//if ( (ne_grid[i+1] > ni_grid[i+1]) && (stop == 0) ) {
		//	size_ngrid = i + 1;
		//	stop = 1;
		//}
		if (i< size_ngrid-2)
			phipp_red[i] = (ne_grid[i+1] - ni_grid[i+1]) - 1.0*phi_grid[i+1]* exp(phi_grid[i+1]);
		else phipp_red[i] = phi0*pow(phi_grid[i+1], pfac) - 1.0*phi_grid[i+1]* exp(phi_grid[i+1]);
		//printf("phipp_red[%d/%d/%d] = %f\n", i, size_ngrid-3, size_phigrid-3, phipp_red[i]);
		//if (i< size_ngrid-2)
		//	printf("i = %d,ne, ni =  %f, %f\n", i, ne_grid[i+1], ni_grid[i+1]); 
	}
	printf("weight= %f, phi_grid[0] = %f\n", weight, phi_grid[0]);
	//phiW_impose = ( - (0.5*v_cutDS*v_cutDS) - (1.0-weight)*phi_grid[0] )/weight;
	phiW_impose = - (0.5*v_cutDS*v_cutDS);
	//printf("phi0 = %f\tCC=%f\n", phi0, CC);
	//printf("phiW_impose = %f\n", phiW_impose);
	phipp_red[0] -= (phiW_impose*invgammasq/deltaxsq);
	//phi0 = phip[size_ngrid-1]/(pdec*pow(x_grid[size_ngrid-1] + CC, pdec-1));
	//phipp_red[size_phigrid-1] -= ( phi0*pow(x_grid[size_phigrid-1] + CC, pdec)*invgammasq/deltaxsq );

	gsl_vector_view phipp_gsl
	    = gsl_vector_view_array (phipp_red, size_phigrid-2);

	gsl_permutation * p = gsl_permutation_alloc (size_phigrid-2);
	clock_t t1 = clock(); // finds the end time of the computation
	gsl_linalg_LU_decomp (m, p, &s);
	clock_t t2 = clock(); double decomptime  = (double)(t2 - t1) / CLOCKS_PER_SEC;
	printf("LU decomposition time = %f\n", decomptime);
	gsl_linalg_LU_solve (m, p, &phipp_gsl.vector, newphi_gsl);

	//printf ("newphi_gsl = \n");
	//printf("%f\n", -0.5*v_cutDS*v_cutDS);
	//gsl_vector_fprintf (stdout, newphi_gsl, "%g");

	newphi[0] = phiW_impose;
	//printf("newphi[0/%d] = %f (imposed)\n", size_ngrid, newphi[0]);
	for (i=0; i<size_ngrid-1; i++) {
		temp = gsl_vector_get(newphi_gsl, i);
		if (temp > 0.0) {
			printf("WARNING: φ > 0.0 and thus non-monotonic in the Debye sheath\n\tFor the moment, this code cannot handle non-monotonic profiles => Exiting code\n");
			temp = 0.0;
			exit(-1);
		}
		if (temp < newphi[i]) {
			printf("WARNING: φ is non-monotonic in the Debye sheath at x = %f\n\tFor the moment, this code cannot handle non-monotonic profiles\n", x_grid[i+1]);
			//exit(-1);
		}
		newphi[i+1] = temp;
		//printf("newphi[%d/%d] = %f\n", i+1, size_ngrid, newphi[i+1]);
	}
	phip0 = (newphi[size_ngrid-1] - newphi[size_ngrid-2])/(x_grid[size_ngrid-1] - x_grid[size_ngrid-2]);
	//printf("phip0 = %f\n\n", phip0);
	CC = pdec*newphi[size_ngrid-1]/phip0 - x_grid[size_ngrid-1];
	//CC = sqrt(pdec*(pdec-1)*newphi[size_ngrid-1]/(invgammasq*(ne_grid[size_ngrid-1] - ni_grid[size_ngrid-1]))) - x_grid[size_ngrid-1];
	printf("CC1 = %f\n\n", CC);
	//CC = (x_grid[size_ngrid-1]*pow(newphi[size_ngrid-1]/newphi[size_ngrid-2], -1.0/pdec) - x_grid[size_ngrid-2])/(1.0 - pow(newphi[size_ngrid-1]/newphi[size_ngrid-2], -1.0/pdec));
	//printf("CC2 = %f\n\n", CC);
	if (CC < -x_grid[size_ngrid-1]) {
		printf("Warning: CC = 0.0\n");
		CC=0.0;
	}
	//if (CC != CC)  CC = pdec*newphi[size_ngrid-1]/phip0 - x_grid[size_ngrid-1];
	//if (pdec*(pdec-1)*newphi[size_ngrid-1]/(ne_grid[size_ngrid-1] - ni_grid[size_ngrid-1]) < 0.0) CC = 0.0; 
	phi0 =  newphi[size_ngrid-1] /pow(x_grid[size_ngrid-1]+CC, pdec) ;
	printf("pdec = %f\n", pdec);
	//printf("C_ds = %f and a_ds = %f\n", CC, phi0);
	for (i=size_ngrid-1; i<size_phigrid; i++) {
		newphi[i] = phi0*pow(x_grid[i]+CC, pdec);
		//printf("newphi[%d/%d] = %f (analytical decay)\n", i, size_ngrid, newphi[i]);
	}
	//printf("%f last\n", newphi[size_phigrid-1]);
	gsl_permutation_free(p);
	gsl_matrix_free (m);
	gsl_vector_free (newphi_gsl);
	free(phipp_red);
}
else {
	for (i=0; i< size_phigrid; i++) {
		if (i< size_ngrid-1)
			//newphi[i] = log( ni_grid[i] - ne_grid[i] + exp(phi_grid[i])); // + phi_grid[i];
			//newphi[i] = phi_grid[i] + 1.0* ( ni_grid[i] - ne_grid[i] ); // + exp(phi_grid[i])); // + phi_grid[i];
			newphi[i] = ( ni_grid[i] - ne_grid[i] )*exp(-phi_grid[i]) + phi_grid[i];
			//newphi[i] = (1.0/ne_grid[i]) * ( ni_grid[i] - ne_grid[i] ) + phi_grid[i];
		else if (i== size_ngrid-1) {
			newphi[i] = log( ni_grid[i] - ne_grid[i] + exp(phi_grid[i])); // + phi_grid[i];
			//newphi[i] = phi_grid[i] + 1.0* ( ni_grid[i] - ne_grid[i] ); // + exp(phi_grid[i])); // + phi_grid[i];
			newphi[i] = ( ni_grid[i] - ne_grid[i] )*exp(-phi_grid[i]) + phi_grid[i];
			//newphi[i] = (1.0/ne_grid[i]) * ( ni_grid[i] - ne_grid[i] ) + phi_grid[i];
			phip0 = (newphi[i] - newphi[i-2])/(x_grid[i] - x_grid[i-2]);
			CC = pdec*newphi[i-1]/phip0 - x_grid[i-1];
			//CC = 0.0;
			phi0 =  newphi[i-1] /pow(x_grid[i-1]+CC, pdec) ;
			//phi0 =  phip0/(pdec*pow(x_grid[i-1]+CC, pdec-1)) ;
			//printf("CC = %f and phi0 = %f\n", CC, phi0);
		}
		else { 
			newphi[i] = phi0*pow(x_grid[i]+CC, pdec);
		}
	// + newphi[size_ngrid-2] - phi0*pow(x_grid[size_ngrid-2], pdec);
	}
}

//printf("res = %f, devbig = %f, dev_0 = %f\n", res, devbig, dev_0);
//printf("conditions %d %d %d\n", res < reslimit, devbig < 5.0*reslimit, dev_0 < 5.0*reslimit);

for (i=size_phigrid-1; i>=0; i--) {
	newphi[i] = weight*newphi[i] + (1.0-weight)*phi_grid[i] ; 
	//printf("newphi[%d] = %f\n", i, newphi[i]);
	phi_grid[i] = newphi[i];
}

//phifree(phip);
//free(phipp);

// remake grid in magnetic presheath
//FILE *fp;
if (attempt_at_adapting_grid==1) {
	new_x = malloc(size_phigrid*sizeof(double));
	printf("in remake_MPgrid: size_phigrid = %d\n", size_phigrid);
	newphi[0] = phi_grid[0];
	new_x[0] = x_grid[0];
	for (i=1; i<size_phigrid; i++) {
		if (phi_grid[i] - phi_grid[i-1] > deltaphi) {
			printf("equal phi interval");
			//gsl_interp_accel *acc
			//= gsl_interp_accel_alloc ();
			//gsl_spline *spline
			//= gsl_spline_alloc (gsl_interp_cspline, size_phigrid);
			//gsl_spline_init (spline, phi_grid, x_grid, size_phigrid);
			newphi[i] = newphi[i-1] + deltaphi;
			//new_x[i] = gsl_spline_eval (spline, newphi[i], acc);
			new_x[i] = lin_interp(phi_grid, x_grid, newphi[i], size_phigrid, 293);
			//gsl_spline_free (spline);
			//gsl_interp_accel_free (acc);
		}
		else {
			printf("equal x interval");
			new_x[i] = new_x[i-1] + deltax;
			if (new_x[i] < x_grid[size_phigrid-1]) {
				//gsl_interp_accel *acc
				//= gsl_interp_accel_alloc ();
				//gsl_spline *spline
				//= gsl_spline_alloc (gsl_interp_cspline, size_phigrid);
				//gsl_spline_init (spline, x_grid, phi_grid, size_phigrid);
				//newphi[i] = gsl_spline_eval (spline, new_x[i], acc);
				newphi[i] = lin_interp(x_grid, phi_grid, new_x[i], size_phigrid, 293);
				//gsl_spline_free (spline);
				//gsl_interp_accel_free (acc);
				printf("path A\n");
			}
			else {
				phip0 = (newphi[i] - newphi[i-1])/(new_x[i] - new_x[i-1]);
				CC = pdec*newphi[size_ngrid-1]/phip0 - new_x[i-1];
				phi0 =  newphi[i-1] /pow(new_x[i-1]+CC, pdec) ;
				newphi[i] = phi0*pow(new_x[i]+CC, pdec);
				printf("path B\n");
			}
			//printf(" i = %d\n", i);
		}
		printf("x phi = (%f %f)\n", new_x[i], newphi[i]);
		//fprintf(fp, "%f %f\n", g, new_phi[i]);
	}
	//*psize_phigrid = i;
	for (i=0; i < size_phigrid; i++) {
		phi_grid[i] = newphi[i];
		x_grid[i] = new_x[i];
	}
	free(new_x);
}
free(newphi);
return;

}
