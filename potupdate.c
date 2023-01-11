// LAST SUBSTANTIAL MODIFICATION MADE 16 JAN 2022
/* This code calculates the next electrostatic potential guess in the iteration to obtain the self-consistent magnetic presheath electrostatic potential profile */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mps_renorm.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>

void newvcut(double *v_cut, double v_cutDS, double u_i, double u_e, double current, double error_current, double weight) {
	double old_v_cut = *v_cut;// u_etilde;
	//u_etilde = u_e - sqrt(mioverme/M_PI)*0.5*exp(-0.5*old_v_cut*old_v_cut);
	//*v_cut = (1.0-weight)*(old_v_cut) + weight*sqrt(-2.0*log(2.0*(-current + u_i - u_etilde)) - log(M_PI/mioverme) ); //OLD
	printf("v_cut before = %f\n", v_cutDS);
	*v_cut = sqrt( 2.0*( (1.0-weight)*0.5*old_v_cut*old_v_cut + weight* ( 0.5*old_v_cut*old_v_cut + (current - u_i + u_e)/( 0.5*exp(-0.5*old_v_cut*old_v_cut) ) ) ) );
	if ( (*v_cut*(*v_cut)*0.5 < 0.5*v_cutDS*v_cutDS ) )  {
		printf("WARNING: v_cutDS > v_cut\n");
		*v_cut = v_cutDS;
	}
	printf("v_cut after = %f\n", v_cutDS);
}

void error_Poisson(double *error, double *x_grid, double *ne_grid, double *ni_grid, double *nioverne, double *phi_grid, int size_phigrid, int size_ngrid, double invgammasq) {
	clock_t begin = clock(); 
	int i;
	double *phip, *phipp;
	double res = 0.0, dev, devbig;
	/* Below we initialize all arrays that contain functions of position x with the correct size n */
	phip = (double*)calloc(  size_phigrid,sizeof(double)); // phi now has correct size
	phipp = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size
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
	for (i=size_ngrid-1; i>=0; i--)  {
		if (i==size_ngrid-1)	printf("gamma^-2 x phi phipp ne ni err \n");
		printf("%f %f %f %f %f %f %f\n", invgammasq, x_grid[i], phi_grid[i], phipp[i], ne_grid[i], ni_grid[i], (-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
		if ( (i!=0) && (i!=size_ngrid-1) ) {
			dev = fabs((-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
			res += fabs((-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
			if (dev > devbig)  devbig = dev; 
		}
	}
	res /= (size_ngrid-2);
	error[0] = res;
	error[1] = devbig;
	free(phip); free(phipp);
	return;
}

void newguess(double *x_grid, double* ne_grid, double *ni_grid, double* phi_grid, int size_phigrid, int size_ngrid, double invgammasq, double v_cutDS, double pfac, double weight) {
clock_t begin = clock(); 
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
//phip = (double*)calloc(  size_phigrid,sizeof(double)); // phi now has correct size
//phipp = (double*)calloc( size_phigrid,sizeof(double)); // phi now has correct size

//reslimit = 0.001;
//
////if (invgammasq > TINY) reslimit *= invgammasq;
//
//printf("weight = %f\n", weight);
//
////if ((fout = fopen("phidataold.txt", "w")) == NULL)
////{
////	printf("Cannot open phidataold.txt");
////	exit(EXIT_FAILURE);
////}
////for (i=0; i<size_phigrid; i++)
////{
////	gg[i] = sqrt(x_grid[i]);
////	fprintf(fout, "%f %f\n", gg[i], phi_grid[i]);
////}
////fclose(fout); 
//
//// Calculate first and second derivatives of potential
//for (i=0; i<size_phigrid; i++) {	
//	if ( (i != size_phigrid-1) && (i!=0) )
//	{	
//		//phipp[i] = ((x_grid[i] - x_grid[i-1])/(x_grid[i+1] - x_grid[i-1]))*(phip[i+1] - phip[i])/(x_grid[i+1] - x_grid[i]);
//		//+ ((x_grid[i+1] - x_grid[i])/(x_grid[i+1] - x_grid[i-1]))*(phip[i] - phip[i-1])/(x_grid[i] - x_grid[i-1]);
//		phip[i] = (phi_grid[i+1] - phi_grid[i]) / (x_grid[i+1] - x_grid[i]);
//		phipp[i] = (phi_grid[i+1] - 2.0*phi_grid[i] + phi_grid[i-1]) / pow(x_grid[i+1] - x_grid[i], 2.0);
//	}
//	else if (i==size_phigrid-1) phipp[size_phigrid-1] = 0.0; //(phip[size_phigrid-1] - phip[size_phigrid-2])/(x_grid[size_phigrid-1] - x_grid[size_phigrid-2]);
//
//	if (i==2) phipp[0] = phipp[1] - x_grid[1]*(phipp[2] - phipp[1])/(x_grid[2] - x_grid[1]);
//	if (i==0) phip[i] = (phi_grid[i+1] - phi_grid[i]) / (x_grid[i+1] - x_grid[i]);
//	if (i==size_phigrid-1) phip[i] = 0.0;
//}
//
//// Calculate the residual of Poisson's/quasineutrality equation
//res = 0.0;
//devbig = dev = 0.0;
//dev_0 = 0.0;
////if (invgammasq < TINY) {
////	fout = fopen("OUTPUT/iterationdata_MP.txt", "w");
////	fout = fopen("OUTPUT/iterationdata_DS.txt", "w");
////	if (fout == NULL) {
////		printf("Cannot open iterationdata.txt");
////		exit(EXIT_FAILURE);
////	}
////}
//for (i=size_ngrid-1; i>=0; i--)  {
//	//if (i==size_ngrid-1)	printf("x phi phipp phipp/gammasq ne ni err \n");
//	//printf(fout, "%f %f %f %f %f %f %f\n", invgammasq, x_grid[i], phi_grid[i], phipp[i], ne_grid[i], ni_grid[i], (-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
//	//fprintf(fout, "%f %f %f %f %f %f %f\n", invgammasq, x_grid[i], phi_grid[i], phipp[i], ne_grid[i], ni_grid[i], (-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
//	if ( (i!=0) && (i!=size_ngrid-1) ) {
//		dev = fabs((-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
//		res += fabs((-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
//		if (dev > devbig)  devbig = dev; 
//		if (i==1) dev_0 = fabs((-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
//	}
//}
////fclose(fout);
//
//res /= (size_ngrid-2);

deltaxsq = x_grid[1]*x_grid[1];

if ( invgammasq > TINY ) {
	phi0 = (ne_grid[size_ngrid-1] - ni_grid[size_ngrid-1])/(pow(phi_grid[size_ngrid-1], pfac));
	
	gsl_vector *newphi_gsl = gsl_vector_alloc (size_phigrid-2);
	gsl_matrix * m = gsl_matrix_alloc (size_phigrid-2, size_phigrid-2);

	for (i = 0; i < size_phigrid-2; i++) {
		for (j = 0; j < size_phigrid-2; j++) {
			if (i == j) {
				gsl_matrix_set (m, i, j, -2.0*invgammasq/deltaxsq - 1.0*exp(phi_grid[i+1])); 
				//if (j==size_phigrid-3) printf("%f\n", -2.0*invgammasq/deltaxsq - 1.0*exp(phi_grid[i+1]));
				//else printf("%f ", -2.0*invgammasq/deltaxsq - 1.0*exp(phi_grid[i+1]));
			}
			else if ( (i==j+1) || (i==j-1) ) {
				gsl_matrix_set (m, i, j, 1.0*invgammasq/deltaxsq);
				//if (j==size_phigrid-3) printf("%f\n", 1.0*invgammasq/deltaxsq);
				//else printf("%f ", 1.0*invgammasq/deltaxsq);
			}
			else { 
				gsl_matrix_set (m, i, j, 0.0);
				//if (j==size_phigrid-3) printf("0.0\n");
				//else printf("0.0 ");
			}
		}
	}

	phipp_red = malloc((size_phigrid-2)*sizeof(double));
	printf("size_phigrid = %d\nswize_ngrid=%d\n\n", size_phigrid, size_ngrid);
	for (i=0;i<size_phigrid-2; i++) {
		if (i< size_ngrid-2)
			phipp_red[i] = (ne_grid[i+1] - ni_grid[i+1]) - phi_grid[i+1]* exp(phi_grid[i+1]);
		else phipp_red[i] = phi0*pow(phi_grid[i+1], pfac) - phi_grid[i+1]* exp(phi_grid[i+1]);
		//else phipp_red[i] = phi0*pow(phi_grid[i+1], pfac) - phi_grid[i+1]* exp(phi_grid[i+1]);
		//printf("phipp_red[%d/%d/%d] = %f\n", i, size_ngrid-3, size_phigrid-3, phipp_red[i]);
		//printf("ne, ni =  %f, %f\n", ne_grid[i+1], ni_grid[i+1]); 
	}
	phiW_impose = (-0.5*v_cutDS*v_cutDS); 
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

	newphi[0] = - 0.5*v_cutDS*v_cutDS;
	printf("%f\n", newphi[0]);
	for (i=0; i<size_ngrid-1; i++) {
		temp = gsl_vector_get(newphi_gsl, i);
		if (temp > 0.0) {
			printf("WARNING: φ > 0.0 in the Debye sheath\n\tFor the moment, this code cannot handle non-monotonic profiles, although in the future perhaps...\n\tIn the meanwhile, exiting code.\n");
			temp = 0.0;
			//exit(-1);
		}
		newphi[i+1] = temp;
		printf("newphi[%d/%d] = %f\n", i+1, size_ngrid, newphi[i+1]);
	}
	phip0 = (newphi[size_ngrid-1] - newphi[size_ngrid-2])/(x_grid[size_ngrid-1] - x_grid[size_ngrid-2]);
	//CC = pdec*newphi[size_ngrid-1]/phip0 - x_grid[size_ngrid-1];
	CC = (x_grid[size_ngrid-1]*pow(newphi[size_ngrid-1]/newphi[size_ngrid-2], -1.0/pdec) - x_grid[size_ngrid-2])/(1.0 - pow(newphi[size_ngrid-1]/newphi[size_ngrid-2], -1.0/pdec));
	if (CC < -x_grid[size_ngrid-1]) CC = 0.0;
	phi0 =  newphi[size_ngrid-1] /pow(x_grid[size_ngrid-1]+CC, pdec) ;
	//printf("C_ds = %f and a_ds = %f\n", CC, phi0);
	for (i=size_ngrid-1; i<size_phigrid; i++) {
		newphi[i] = phi0*pow(x_grid[i]+CC, pdec);
		//printf("%f decay\n", newphi[i]);
	}
	printf("%f last\n", newphi[size_phigrid-1]);
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
			//newphi[i] = log( ni_grid[i] - ne_grid[i] + exp(phi_grid[i])); // + phi_grid[i];
			//newphi[i] = phi_grid[i] + 1.0* ( ni_grid[i] - ne_grid[i] ); // + exp(phi_grid[i])); // + phi_grid[i];
			newphi[i] = ( ni_grid[i] - ne_grid[i] )*exp(-phi_grid[i]) + phi_grid[i];
			//newphi[i] = (1.0/ne_grid[i]) * ( ni_grid[i] - ne_grid[i] ) + phi_grid[i];
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

//printf("res = %f, devbig = %f, dev_0 = %f\n", res, devbig, dev_0);
//printf("conditions %d %d %d\n", res < reslimit, devbig < 5.0*reslimit, dev_0 < 5.0*reslimit);

for (i=size_phigrid-1; i>=0; i--) {
	newphi[i] = weight*newphi[i] + (1.0-weight)*phi_grid[i] ; 
	phi_grid[i] = newphi[i];
}
//phifree(phip);
//free(phipp);
free(newphi);

clock_t end = clock(); // finds the end time of the computation
double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
printf("Module to obtain the new guess for the electrostatic potential ran in %f\n", jobtime);
return;
}
