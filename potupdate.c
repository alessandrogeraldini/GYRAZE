
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
#include <gsl/gsl_blas.h>
#if DS_SLP
#include <glpk.h>
#endif

/* 
	CALCULATE THE NEW TOTAL POTENTIAL DROP ACROSS THE COMBINED MAGNETIC PRESHEATH AND DEBYE SHEATH
*/
void newvcut(double *v_cut, double v_cutDS, double u_i, double u_e, double current, double error_current, double weight) {
	double old_v_cut = *v_cut;// u_etilde;
	//u_etilde = u_e - sqrt(mioverme/M_PI)*0.5*exp(-0.5*old_v_cut*old_v_cut);
	//*v_cut = (1.0-weight)*(old_v_cut) + weight*sqrt(-2.0*log(2.0*(-current + u_i - u_etilde)) - log(M_PI/mioverme) ); //OLD
	printf("phi_wall before = %f --> \t", 0.5*(*v_cut)*(*v_cut));
	*v_cut = sqrt( 2.0*( (1.0-weight)*0.5*old_v_cut*old_v_cut + weight* ( 0.5*old_v_cut*old_v_cut + (current - u_i + u_e)/( 0.5*exp(-0.5*old_v_cut*old_v_cut) ) ) ) );
	//*v_cut = sqrt( 2.0*( (1.0-weight)*0.5*old_v_cut*old_v_cut - weight* log(sqrt(2.0*M_PI)*(current + u_i - u_e) + exp(-0.5*old_v_cut*old_v_cut) ) ) );
	//*v_cut = sqrt( - 2.0*( log( 1.0*sqrt(2.0*M_PI)*(current + u_i - u_e) + 1.0*exp(-0.5*old_v_cut*old_v_cut) ) ) );
	//if ( (*v_cut*(*v_cut)*0.5 < 0.5*v_cutDS*v_cutDS ) )  {
	//	printf("WARNING: v_cutDS > v_cut\n");
	//	*v_cut = v_cutDS;
	//}
	if ( (*v_cut*(*v_cut)*0.5 < 0.5*old_v_cut*old_v_cut - 0.5*v_cutDS*v_cutDS ) )  {
		printf("WARNING: old_v_cut= %f\tnewv_cut = %f\t", old_v_cut, *v_cut);
		*v_cut = sqrt(old_v_cut*old_v_cut - v_cutDS*v_cutDS) + TINY;
		printf("new new_v_cut = %f\n", *v_cut);
	}

	printf("v_cut after = %f\n", *v_cut);
	printf("phi_wall after = %f --> \t", 0.5*(*v_cut)*(*v_cut));
}


/* 
	CALCULATE THE ERROR IN POISSON'S EQUATION OR QUASINEUTRALITY (if invgammasq = 0)
*/
int error_Poisson_imax = -1;

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
	error_Poisson_imax = -1;
	for (i=size_ngrid-1; i>=0; i--)  {
		//if (i==size_ngrid-1)	printf("gamma^-2 x phi phipp ne ni err \n");
		//printf("%f %f %f %f %f %f %f\n", invgammasq, x_grid[i], phi_grid[i], phipp[i], ne_grid[i], ni_grid[i], (-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
		if ( (i!=0) && (i!=size_ngrid-1) ) {
			dev = fabs((-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
			res += dev;//fabs((-ne_grid[i] + phipp[i]*invgammasq)/ni_grid[i] + 1.0);
			//dev = fabs((-ne_grid[i] + phipp[i]*invgammasq) + ni_grid[i]);
			//res += fabs((-ne_grid[i] + phipp[i]*invgammasq) + ni_grid[i]);
			if (dev > devbig)  {
				devbig = dev; error_Poisson_imax = i; 
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

/* 
	CALCULATE THE NEW POTENTIAL PROFILE ACROSS THE MAGNETIC PRESHEATH (if invgammasq=0) OR DEBYE SHEATH
*/
// The function below is momentarily only valid for Ti = Te, so TiovTe = 1 in the physical input file
void newguess(double *x_grid, double* ne_grid, double *ni_grid, double* phi_grid, int size_phigrid, int size_ngridin, double invgammasq, double v_cutDS, double pfac, double weight) {
//attempt_at_adapting_grid failed. It was an attempt at refining grid near x=0 using a maximum Delta x maximum Delta phi gridding
int i, j, s, attempt_at_adapting_grid=0, size_ngrid=size_ngridin;
double phiW_impose, temp, res = 1000.0;
double *new_x;
double *phipp_red, *newphi, *oldphi, *phipp;
double phi0, phip0, CC, pdec, deltaxsq, deltax = 0.3, deltaphi = 0.1;
//FILE *fout;

pdec = 2.0/(1.0-pfac);
printf("φ'' ~ φ^%f gives φ ~ x^%f\n", pfac, pdec);
printf("0.5*v_cutDS*v_cutDS = %f\n", 0.5*v_cutDS*v_cutDS);

/* Below we initialize all arrays that contain functions of position x with the correct size n */
newphi = (double*)calloc(size_phigrid,sizeof(double)); // same as above 
oldphi = (double*)calloc(size_phigrid,sizeof(double)); // same as above 
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

for (i=0;i<size_phigrid;i++)
	oldphi[i] = phi_grid[i];

if ( invgammasq > TINY ) { // DEBYE SHEATH ITERATION
	phi0 = (ne_grid[size_ngrid-1] - ni_grid[size_ngrid-1])/(pow(phi_grid[size_ngrid-1], pfac));
	printf("At beginning phi0 = %f\n", phi0);
	gsl_vector *newphi_gsl = gsl_vector_alloc (size_phigrid-2);
	gsl_matrix * m = gsl_matrix_alloc (size_phigrid-2, size_phigrid-2);


	while (res>0.005) { // This loop is taking advantage of Nicole's trick

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
		if (i< size_ngrid-2)
			phipp_red[i] = (ne_grid[i+1] - ni_grid[i+1]) - phi_grid[i+1]* exp(phi_grid[i+1]);
		else phipp_red[i] = phi0*pow(phi_grid[i+1], pfac) - phi_grid[i+1]* exp(phi_grid[i+1]);
	}
	printf("weight= %f, phi_grid[0] = %f\n", weight, phi_grid[0]);
	//phiW_impose = ( - (0.5*v_cutDS*v_cutDS) - (1.0-weight)*phi_grid[0] )/weight;
	phiW_impose = - (0.5*v_cutDS*v_cutDS);
	//printf("phi0 = %f\tCC=%f\n", phi0, CC);
	//printf("phiW_impose = %f\n", phiW_impose);
	phipp_red[0] -= (phiW_impose*invgammasq/deltaxsq);
	//phi0 = phip[size_ngrid-1]/(pdec*pow(x_grid[size_ngrid-1] + CC, pdec-1));
	//phipp_red[size_phigrid-1] -= ( phi0*pow(x_grid[size_phigrid-1] + CC, pdec)*invgammasq/deltaxsq );

	gsl_vector_view phipp_gsl = gsl_vector_view_array (phipp_red, size_phigrid-2);

	gsl_permutation * p = gsl_permutation_alloc (size_phigrid-2);
	clock_t t1 = clock(); // finds the end time of the computation
	gsl_linalg_LU_decomp (m, p, &s);
	clock_t t2 = clock(); double decomptime  = (double)(t2 - t1) / CLOCKS_PER_SEC;
	printf("LU decomposition time = %f\n", decomptime);
	gsl_linalg_LU_solve (m, p, &phipp_gsl.vector, newphi_gsl);

	newphi[0] = phiW_impose;
	//printf("newphi[0/%d] = %f (imposed)\n", size_ngrid, newphi[0]);
        printf("In Debye sheath, asymptotic result starts at x = %f\n\n", x_grid[size_ngrid-1]);
	for (i=0; i<size_ngrid-1; i++) {
		temp = gsl_vector_get(newphi_gsl, i);
		if (temp > 0.0) {
			printf("WARNING: φ > 0.0 at x = %f, and thus non-monotonic in the Debye sheath\n\tFor the moment, this code cannot handle non-monotonic profiles => Exiting code\n", x_grid[i+1]);
//			temp = 0.0;
//			exit(-1);
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
	printf("In Debye sheath CC1 = %f\n\n", CC);
	        printf("In Debye sheath, max x = %f\n\n", x_grid[size_phigrid-1]);
	printf("last phi is %f\n", newphi[size_ngrid-1]);
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
	printf("phi0 = %f\n", phi0);
	//printf("C_ds = %f and a_ds = %f\n", CC, phi0);
	for (i=size_ngrid-1; i<size_phigrid; i++) {
		newphi[i] = phi0*pow(x_grid[i]+CC, pdec);
		//printf("newphi[%d/%d] = %f (analytical decay)\n", i, size_ngrid, newphi[i]);
	}
	for (i=size_phigrid-1; i>= 0; i--) {
		res += fabs(phi_grid[i] - newphi[i]);
		phi_grid[i] = newphi[i];
	}
	res /= size_phigrid;
	gsl_permutation_free(p);
	printf("res = %f\n", res);
	res = 0.00001;
	}
	//printf("%f last\n", newphi[size_phigrid-1]);
	gsl_matrix_free (m);
	gsl_vector_free (newphi_gsl);
	free(phipp_red);
}
else { // MAGNETIC PRESHEATH ITERATION
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

int mixed = 0;
#if ANDERSON_M > 0
/* Anderson acceleration of the Debye sheath Picard map: with f = G(phi) - phi the residual of the
 * map, take phi + weight*f corrected by the last ANDERSON_M differences of phi and of f, which can
 * cancel the weakly unstable modes that plain damping cannot. The history is kept through a rising
 * error (that is when it is needed) and dropped only when the grid changes; the correction is
 * bounded by DS_DPHIMAX so a bad extrapolation cannot run away. */
if (invgammasq != 0.0) {
	static double *px = NULL, *pf = NULL, *dX[ANDERSON_M], *dF[ANDERSON_M];
	static int nh = 0, nprev = -1;
	int k, n = size_phigrid;
	double *f = malloc(n*sizeof(double)), gam[ANDERSON_M], E_act = 0.0;
	for (i=0; i<n; i++) f[i] = newphi[i] - oldphi[i];
	for (i=1; i<size_ngrid-1; i++) {
		double pp = (phi_grid[i+1] - 2.0*phi_grid[i] + phi_grid[i-1])/(x_grid[1]*x_grid[1]);
		E_act = fmax(E_act, fabs((-ne_grid[i] + pp*invgammasq)/ni_grid[i] + 1.0));
	}
	if (px == NULL || nprev != n) {
		free(px); free(pf);
		px = malloc(n*sizeof(double)); pf = malloc(n*sizeof(double));
		for (k=0; k<ANDERSON_M; k++) { if (nprev > 0) { free(dX[k]); free(dF[k]); } dX[k] = malloc(n*sizeof(double)); dF[k] = malloc(n*sizeof(double)); }
		nh = 0; nprev = n;
	}
	else {   /* newest difference first, oldest dropped */
		double *ox = dX[ANDERSON_M-1], *of = dF[ANDERSON_M-1];
		for (k=ANDERSON_M-1; k>0; k--) { dX[k] = dX[k-1]; dF[k] = dF[k-1]; }
		dX[0] = ox; dF[0] = of;
		for (i=0; i<n; i++) { dX[0][i] = oldphi[i] - px[i]; dF[0][i] = f[i] - pf[i]; }
		if (nh < ANDERSON_M) nh++;
	}
	for (i=0; i<n; i++) { px[i] = oldphi[i]; pf[i] = f[i]; }

	/* least squares min |f - dF*gam| via the (small) normal equations, lightly regularized */
	for (k=0; k<ANDERSON_M; k++) gam[k] = 0.0;
	if (nh > 0) {
		gsl_matrix *A = gsl_matrix_alloc(nh, nh);
		gsl_vector *b = gsl_vector_alloc(nh), *g = gsl_vector_alloc(nh);
		gsl_permutation *perm = gsl_permutation_alloc(nh);
		double tr = 0.0;
		for (k=0; k<nh; k++) {
			double bk = 0.0;
			for (i=0; i<n; i++) bk += dF[k][i]*f[i];
			gsl_vector_set(b, k, bk);
			for (j=0; j<nh; j++) {
				double a = 0.0;
				for (i=0; i<n; i++) a += dF[k][i]*dF[j][i];
				if (k == j) tr += a;
				gsl_matrix_set(A, k, j, a);
			}
		}
		for (k=0; k<nh; k++) gsl_matrix_set(A, k, k, gsl_matrix_get(A, k, k) + 1e-8*tr/nh);
		gsl_linalg_LU_decomp(A, perm, &s);
		gsl_linalg_LU_solve(A, perm, b, g);
		for (k=0; k<nh; k++) gam[k] = gsl_vector_get(g, k);
		gsl_permutation_free(perm); gsl_vector_free(b); gsl_vector_free(g); gsl_matrix_free(A);
	}
	for (i=0; i<n; i++) {
		newphi[i] = oldphi[i] + weight*f[i];
		for (k=0; k<nh; k++) newphi[i] -= gam[k]*(dX[k][i] + weight*dF[k][i]);
	}
	/* trust region: same limit on the change in phi as the Newton update */
	double dmax = 0.0;
	for (i=0; i<n; i++) dmax = fmax(dmax, fabs(newphi[i] - oldphi[i]));
	if (dmax > DS_DPHIMAX)
		for (i=0; i<n; i++) newphi[i] = oldphi[i] + (DS_DPHIMAX/dmax)*(newphi[i] - oldphi[i]);
	printf("Anderson: %d stored iterates, max|dphi| = %f, error %f\n", nh, fmin(dmax, DS_DPHIMAX), E_act);
	free(f);
	mixed = 1;
}
#endif
for (i=size_phigrid-1; i>=0; i--) {
	if (!mixed) newphi[i] = weight*newphi[i] + (1.0-weight)*oldphi[i] ;
	//printf("newphi[%d] = %f\n", i, newphi[i]);
	phi_grid[i] = newphi[i];
}

// remake grid in magnetic presheath // This attempt failed miserably, don't change flag value to 1
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
	free(oldphi);
	free(phipp);
	return;

}

/*
	NEWTON-RAPHSON UPDATE OF THE DEBYE SHEATH POTENTIAL

	Solves  phi'' - gamma^2*(ne - ni) = 0  for the Debye sheath region via one
	Newton-Raphson step:  J * dphi = -F,  phi += weight * dphi

	The Jacobian is tridiagonal:
	  off-diagonal:   1/dx^2
	  diagonal:      -2/dx^2  -  gamma^2 * (ne_corr_delta[i] + ne_corr_chiM[i] + ne_grid[i] - ni_corr[i])
	where i is the matrix row (0 = first interior point).
	ne_corr_chiM[0] is NaN; ne_corr_chiM[1] is used in its place.

	Boundary condition and asymptotic extension are kept identical to newguess.
*/
/* x = V diag(s/(s^2 + lam^2)) U^T b, for J = U diag(s) V^T: the Levenberg-Marquardt step, i.e. the solution
 * of (J^T J + lam^2 I) x = J^T b */
static void lm_solve(const gsl_matrix *U, const gsl_matrix *V, const gsl_vector *S, double lam,
                     const gsl_vector *b, gsl_vector *x)
{
	size_t n = S->size, i;
	gsl_vector *c = gsl_vector_alloc(n);
	gsl_blas_dgemv(CblasTrans, 1.0, U, b, 0.0, c);
	for (i = 0; i < n; i++) {
		double si = gsl_vector_get(S, i);
		gsl_vector_set(c, i, gsl_vector_get(c, i) * si / (si*si + lam*lam));
	}
	gsl_blas_dgemv(CblasNoTrans, 1.0, V, c, 0.0, x);
	gsl_vector_free(c);
}

#if DS_SLP
/* Endgame DS step by linear programming (DS_SLP). With r_i(d) = s_i (-F_i + (J d)_i) the linearised relative
 * Poisson residual after a step d (s_i = 1/(gamma^2 ni); F_vec holds -F, as newguess_NR builds it), solve
 *
 *     minimise t   s.t.  -t <= r_i(d) <= t            every row
 *                        a_i >= |r_i(d)|,  sum_{i<navg} a_i <= avgcap * navg     (error_Poisson's mean)
 *                        |d_j| <= bound_j                                         (trust region)
 *                        |d_{j-1} - 2 d_j + d_{j+1}| <= smooth
 *
 * That LP is infeasible when the current mean is above avgcap and the trust region is too small to bring
 * it under the cap in one step (then nothing was solved and the LM step was kept). So in that case the LP
 * is phase A instead: minimise the mean, with the max held at or below its current value (t <= max_i |r_i(0)|).
 * d = 0 satisfies every constraint of either phase, so the LP is always feasible. *phase is 1 (A) or 2.
 *
 * Returns 1 and fills d (and the predicted max t and mean) on an optimal solution, 0 otherwise. */
static int slp_solve(const gsl_matrix *Jc, const double *Fneg, const double *sc, int n, int navg,
                     const double *bound, double smooth, double avgcap, double *d, double *tpred, double *apred,
                     int *phase)
{
	int i, j, k, r = 0, ok;
	double m0 = 0.0, a0 = 0.0;                         /* max and mean of |r_i| at d = 0 */
	for (i = 0; i < n; i++) m0 = fmax(m0, fabs(sc[i] * Fneg[i]));
	for (k = 0; k < navg; k++) a0 += fabs(sc[k] * Fneg[k]);
	a0 /= navg;
	*phase = (a0 > avgcap) ? 1 : 2;
	const int ct = n + 1, ca = n + 2;                  /* 1-based columns: d_1..d_n, t, a_1..a_n */
	size_t cap = (size_t)4 * n * (n + 1) + (size_t)navg + (size_t)6 * n + 16, ne = 0;
	int *ia = malloc((cap + 1) * sizeof(int)), *ja = malloc((cap + 1) * sizeof(int));
	double *ar = malloc((cap + 1) * sizeof(double));
	glp_prob *lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MIN);
	glp_add_cols(lp, 2 * n + 1);
	for (j = 0; j < n; j++) {
		if (bound[j] > 0.0) glp_set_col_bnds(lp, j + 1, GLP_DB, -bound[j], bound[j]);
		else                glp_set_col_bnds(lp, j + 1, GLP_FX, 0.0, 0.0);
	}
	if (*phase == 1) {                                 /* minimise the mean, max held */
		glp_set_col_bnds(lp, ct, GLP_DB, 0.0, m0);
		for (k = 0; k < navg; k++) glp_set_obj_coef(lp, ca + k, 1.0 / navg);
	}
	else {                                             /* minimise the max, mean capped */
		glp_set_col_bnds(lp, ct, GLP_LO, 0.0, 0.0);
		glp_set_obj_coef(lp, ct, 1.0);
	}
	for (k = 0; k < n; k++) glp_set_col_bnds(lp, ca + k, GLP_LO, 0.0, 0.0);
	glp_add_rows(lp, 4 * n + 1 + 2 * (n - 2));
	for (i = 0; i < n; i++) {
		double c = sc[i] * Fneg[i];                    /* r_i = s_i (J d)_i - c */
		int sgn, col;
		for (col = 0; col < 2; col++) {                /* against t, then against a_i */
			for (sgn = 1; sgn >= -1; sgn -= 2) {       /* +r_i <= .. and -r_i <= .. */
				r++;
				glp_set_row_bnds(lp, r, GLP_UP, 0.0, sgn * c);
				for (j = 0; j < n; j++) {
					double v = gsl_matrix_get(Jc, i, j);
					if (v != 0.0) { ne++; ia[ne] = r; ja[ne] = j + 1; ar[ne] = sgn * sc[i] * v; }
				}
				ne++; ia[ne] = r; ja[ne] = (col == 0) ? ct : ca + i; ar[ne] = -1.0;
			}
		}
	}
	r++;
	if (*phase == 1) glp_set_row_bnds(lp, r, GLP_FR, 0.0, 0.0);
	else             glp_set_row_bnds(lp, r, GLP_UP, 0.0, avgcap * navg);
	for (k = 0; k < navg; k++) { ne++; ia[ne] = r; ja[ne] = ca + k; ar[ne] = 1.0; }
	for (j = 1; j < n - 1; j++) {
		int sgn;
		for (sgn = 1; sgn >= -1; sgn -= 2) {
			r++;
			glp_set_row_bnds(lp, r, GLP_UP, 0.0, smooth);
			ne++; ia[ne] = r; ja[ne] = j;     ar[ne] =  1.0 * sgn;
			ne++; ia[ne] = r; ja[ne] = j + 1; ar[ne] = -2.0 * sgn;
			ne++; ia[ne] = r; ja[ne] = j + 2; ar[ne] =  1.0 * sgn;
		}
	}
	glp_load_matrix(lp, (int)ne, ia, ja, ar);
	glp_scale_prob(lp, GLP_SF_AUTO);
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_OFF;
	parm.presolve = GLP_ON;
	ok = (glp_simplex(lp, &parm) == 0 && glp_get_status(lp) == GLP_OPT);
	if (ok) {                                          /* predicted residuals from J d itself (in phase A */
		double sa = 0.0;                               /* t is only a bound, and a_i >= |r_i| is not tight) */
		for (j = 0; j < n; j++) d[j] = glp_get_col_prim(lp, j + 1);
		*tpred = 0.0;
		for (i = 0; i < n; i++) {
			double ri = -sc[i] * Fneg[i];
			for (j = 0; j < n; j++) ri += sc[i] * gsl_matrix_get(Jc, i, j) * d[j];
			*tpred = fmax(*tpred, fabs(ri));
			if (i < navg) sa += fabs(ri);
		}
		*apred = sa / navg;
	}
	glp_delete_prob(lp);
	free(ia); free(ja); free(ar);
	return ok;
}
#endif

void newguess_NR(double *x_grid, double *ne_grid, double *ni_grid, double *phi_grid,
                 int size_phigrid, int size_ngridin, double invgammasq, double v_cutDS,
                 double pfac, double weight,
                 double *ne_corr_delta, double *ne_corr_chiM, double *ni_corr,
                 double *jac_y, double *jac_h, int jac_K, double *dir_h, double *dir_y, double *ls_v, double *ls_tmodel)
{
	int i, j, k, s;
	int size_ngrid = size_ngridin;
	int ninner = size_ngrid - 1;  /* unknowns: phi_grid[1] .. phi_grid[size_ngrid-1] */
	double gamma2   = 1.0 / invgammasq;
	double deltaxsq = x_grid[1] * x_grid[1];
	double pdec     = 2.0 / (1.0 - pfac);
	double phiW_impose = -(0.5 * v_cutDS * v_cutDS);
	double phi0, CC, temp;

	double *F_vec      = malloc(ninner * sizeof(double));
	double *ne_corr_total = malloc(ninner * sizeof(double));
	gsl_vector   *dphi_gsl = gsl_vector_alloc(ninner);
	gsl_matrix   *J        = gsl_matrix_alloc(ninner, ninner);

	printf("weight= %f, phi_grid[0] = %f\n", weight, phi_grid[0]);
	printf("0.5*v_cutDS*v_cutDS = %f\n", 0.5 * v_cutDS * v_cutDS);

	/* Build ne_corr_total at grid point i+1 (same index as residual row i).
	 * ne_corr_chiM[size_ngrid-1] is NaN; substitute [size_ngrid-2] there. */
	for (i = 0; i < ninner; i++) {
		int idx = i + 1;
		double chiM_i;
		// if (idx == size_ngrid - 1)
		// 	chiM_i = ne_corr_chiM[size_ngrid - 2];
		// else
		// 	chiM_i = ne_corr_chiM[idx];
		chiM_i = ne_corr_chiM[idx];
		ne_corr_total[i] = ne_grid[idx] + ne_corr_delta[idx] + chiM_i;
		//ne_corr_total[i] = ne_grid[idx];
	}

	/* Normalize so ne_corr_total[ninner-1] == ni_corr[ninner-1] */
	double scale = ni_corr[ninner-1] / ne_corr_total[ninner-1];
	for (i = 0; i < ninner; i++)
		ne_corr_total[i] *= scale;

	/* Save phi_grid[0] so the Apply step can weight the BC change by alpha,
	 * matching what is done for all interior points.  The residual and
	 * Jacobian are still built with the target phiW_impose so that the NR
	 * step (dphi[0] ~ delta/2) already accounts for the BC shift. */
	double phi0_before = phi_grid[0];

	/* Enforce wall boundary condition for residual/Jacobian computation */
	phi_grid[0] = phiW_impose;

	/* Right end of the unknowns (the n_e grid end, e = size_ngrid-1): the next point belongs to the
	 * asymptotic tail phi0 (x + CC)^pdec. Rather than holding it fixed, which puts a corner at the join as
	 * soon as phi_e moves, tie it to phi_e through the tail: phi_{e+1} = rho phi_e,
	 * rho = ((x_{e+1} + CC)/(x_e + CC))^pdec, with CC from the last valid tail or, before the first one,
	 * from the tail points just beyond e. rho = 0 (no usable CC) keeps the fixed value. */
	static double CC_prev = 0.0;
	static int have_CC = 0;
	int e = size_ngrid - 1;
	double CCb = CC_prev, rho = 0.0;
	if (!have_CC) {
		CCb = -1e30;
		if (e + 2 < size_phigrid && phi_grid[e+1] < 0.0 && phi_grid[e+2] > phi_grid[e+1])
			CCb = pdec * phi_grid[e+1] * (x_grid[e+2] - x_grid[e+1]) / (phi_grid[e+2] - phi_grid[e+1]) - x_grid[e+1];
	}
	if (x_grid[e] + CCb > 0.0) rho = pow((x_grid[e+1] + CCb) / (x_grid[e] + CCb), pdec);

	/* Build Jacobian J (tridiagonal) */
	for (i = 0; i < ninner; i++) {
		double jac_dens = gamma2 * (ne_corr_total[i] - ni_corr[i]);
		for (j = 0; j < ninner; j++) {
			if (i == j)
				gsl_matrix_set(J, i, j, -2.0 / deltaxsq - jac_dens);
			else if ((i == j + 1) || (i == j - 1))
				gsl_matrix_set(J, i, j, 1.0 / deltaxsq);
			else
				gsl_matrix_set(J, i, j, 0.0);
		}
	}
	gsl_matrix_set(J, ninner-1, ninner-1, gsl_matrix_get(J, ninner-1, ninner-1) + rho / deltaxsq);

	/* Nonlocal electron response (NONLOCAL_JAC > 0): the loop above took dne/dphi to be the local
	 * value ne_corr_total on the diagonal, but n_e at one point depends on phi along the whole orbit.
	 * GYRAZE.c supplies, for each localized perturbation h_k of phi, the measured response
	 * y_k = dn_e/deps.  Correct J so that it reproduces those measurements, i.e. so the electron
	 * block A satisfies A h_k = y_k, while keeping the local model on everything orthogonal to the
	 * h_k.  The h_k have disjoint support, so the rank-1 corrections simply add. */
	if (jac_y != NULL && jac_K > 0) printf("NR: nonlocal Jacobian from %d measured directions\n", jac_K);
	for (k = 0; k < jac_K && jac_y != NULL; k++) {
		double hh = 0.0;
		for (j = 0; j < ninner; j++) hh += jac_h[k*size_phigrid + j+1] * jac_h[k*size_phigrid + j+1];
		if (hh < TINY) continue;
		for (i = 0; i < ninner; i++) {
			/* with NONLOCAL_JAC_IONS the measurement is d(ne - ni)/dphi, so the local model it replaces is
			 * ne_corr_total - ni_corr; otherwise it is the electron response alone */
			double dens_loc = ne_corr_total[i] - (NONLOCAL_JAC_IONS ? ni_corr[i] : 0.0);
			double c = gamma2 * (jac_y[k*size_phigrid + i+1] - dens_loc*jac_h[k*size_phigrid + i+1]) / hh;
			if (!isfinite(c)) continue;   /* no valid measurement in this row: keep the local model */
			for (j = 0; j < ninner; j++)
				gsl_matrix_set(J, i, j, gsl_matrix_get(J, i, j) - c * jac_h[k*size_phigrid + j+1]);
		}
	}

	/* Measured response along the previous Newton direction (NONLOCAL_JAC_STEPDIR): dir_y is dn_e along
	 * dir_h.  Near a stall the step follows one smooth mode that the bumps above resolve worst, and there
	 * even a 10% error in the density response can exceed the Jacobian's (small) singular value and turn
	 * the step uphill.  Correct J by a rank-1 update so that J dir_h reproduces the measured change. */
	if (dir_h != NULL && dir_y != NULL) {
		double hh = 0.0, err = 0.0, tt = 0.0;
		for (j = 0; j < ninner; j++) hh += dir_h[j+1] * dir_h[j+1];
		for (i = 0; i < ninner && hh > TINY; i++) {
			double lap = (dir_h[i] - 2.0*dir_h[i+1] + (i < ninner-1 ? dir_h[i+2] : 0.0)) / deltaxsq;
			double target = lap - gamma2 * (dir_y[i+1] - (NONLOCAL_JAC_IONS ? 0.0 : ni_corr[i]*dir_h[i+1]));
			double Jh = 0.0;
			if (!isfinite(target)) continue;   /* no valid measurement in this row: keep J */
			for (j = 0; j < ninner; j++) Jh += gsl_matrix_get(J, i, j) * dir_h[j+1];
			err += (target - Jh) * (target - Jh); tt += target * target;
			for (j = 0; j < ninner; j++)
				gsl_matrix_set(J, i, j, gsl_matrix_get(J, i, j) + (target - Jh) * dir_h[j+1] / hh);
		}
		printf("NR: step-direction correction, |J h - measured| / |measured| = %.3f\n", sqrt(err / fmax(tt, TINY)));
	}

	/* Build -F[i] = -(phi''[i+1]/dx^2 - gamma^2*(ne[i+1] - ni[i+1]))
	 * phi_grid[0] = phiW_impose is already set; phi_grid[size_ngrid] is from
	 * the previous iteration's asymptotic extension (right boundary of last row). */
	for (i = 0; i < ninner; i++) {
		double right = (i == ninner-1 && rho > 0.0) ? rho * phi_grid[e] : phi_grid[i+2];
		double phipp = (right - 2.0*phi_grid[i+1] + phi_grid[i]) / deltaxsq;
		F_vec[i] = -(phipp - gamma2 * (ne_grid[i+1] - ni_grid[i+1]));  /* = -F */
	}

	/* Max relative Poisson error at the current phi (the error_DS[1] criterion in GYRAZE.c).
	 * The densities were evaluated at this phi, so this is the true residual. */
	double E0 = 0.0;
	for (i = 0; i < ninner; i++) {
		double dev = fabs(F_vec[i]) * invgammasq / ni_grid[i+1];
		if (dev > E0) E0 = dev;
	}
	double alpha = weight;

#if NR_SAFESTEP == 1
	/* Safeguards: compare the true residual with the previous call's and, if it rose, go back to
	 * that phi and halve the step. The step is also trust-region limited (see below), since the
	 * local Jacobian is near-singular for smooth modes where phi approaches 0. */
	static double *phi_prev = NULL, *F_prev = NULL, E_prev, alpha_prev, phiW_prev, phi0_prev;
	static int n_prev = -1, np_alloc = -1, nF_alloc = -1;
	static double bscale = 1.0;   /* fraction of the wall step allowed; halved on rejection like alpha */
	int hold_step = 0;            /* set when a rejected step cannot be retried at this grid size */
	/* phi_prev holds the whole profile, not just the inner points, so it stays valid when the n_e grid
	 * resizes (the step also rewrites the asymptotic tail) and the comparison below survives a resize. */
	if (phi_prev == NULL || np_alloc != size_phigrid) {
		free(phi_prev); phi_prev = malloc(size_phigrid * sizeof(double));
		np_alloc = size_phigrid; n_prev = -1;
	}
	/* Error used to accept/reject steps: the RMS of the relative Poisson residual. A Newton step is a
	 * descent direction for the 2-norm, not the max-norm, so judging steps by the max rejects good
	 * steps that lower the residual overall but raise it at one point. (Convergence is still tested
	 * on the max, error_DS[1], in GYRAZE.c.) Row 0 uses the wall at its current value rather than the
	 * imposed one, so a pending wall-potential change doesn't dominate it. */
	double r0 = ((phi_grid[2] - 2.0*phi_grid[1] + phi0_before)/deltaxsq - gamma2*(ne_grid[1] - ni_grid[1]))
	            * invgammasq / ni_grid[1];
	double E_act = r0*r0;
	for (i = 1; i < ninner; i++) { double r = F_vec[i] * invgammasq / ni_grid[i+1]; E_act += r*r; }
	E_act = sqrt(E_act / ninner);
#if DS_MINIMAX || DS_SLP
	/* merit on the max, as the convergence test is (see DS_MINIMAX, DS_SLP in mps.h) */
	E_act = fabs(r0);
	for (i = 1; i < ninner; i++) E_act = fmax(E_act, fabs(F_vec[i] * invgammasq / ni_grid[i+1]));
	const double merit_tol = DS_MINIMAX_MERIT_TOL;
#else
	const double merit_tol = DS_MERIT_TOL;
#endif
	if (n_prev < 0 || fabs(phiW_prev - phiW_impose) > DS_WALL_TOL) {
		/* First call, or the wall target moved enough that the residuals before and after belong to
		 * different problems. The target is compared with DS_WALL_TOL rather than exactly: it drifts by
		 * ~1e-4 every iteration as the MP updates, and an exact test threw the comparison away every
		 * call, leaving the step unguarded. A change in the n_e grid size is NOT a reason to drop the
		 * comparison: E is an rms, normalized by ninner, so it stays comparable, and phi_prev covers the
		 * whole profile. Only the stored residual cannot be reused across a resize, handled below. */
		bscale = 1.0;
	}
	/* A held step (alpha_prev = 0) left phi unchanged, so there is nothing to test: resume stepping.
	 * A step that moved the wall by more than DS_WALL_TOL changed the boundary condition, so the residuals
	 * before and after it belong to different problems and are not compared either: while the wall
	 * approaches its target each step is taken (it is still bounded by the trust region) and becomes the
	 * new baseline. Smaller moves (the target drifts a little as the MP updates) still get the test.
	 * Otherwise allow for the run-to-run noise in the recomputed densities and for the n_e grid
	 * resizing under the rms (DS_MERIT_TOL), either of which is enough to make an unchanged error look
	 * like a rise and hold the step at zero forever. */
	else if (alpha_prev > 0.0 && fabs(phi0_before - phi0_prev) < DS_WALL_TOL && E_act > E_prev * (1.0 + merit_tol)) {
		double E_rose = E_act;
		for (i = 0; i < size_phigrid; i++) phi_grid[i] = phi_prev[i];
		E_act = E_prev;
		phi0_before = phi0_prev;
		bscale *= 0.5;
		if (n_prev == ninner && F_prev != NULL) {
			for (i = 0; i < ninner; i++) F_vec[i] = F_prev[i];
			alpha = 0.5 * alpha_prev;
			/* At the smallest step the Newton direction still does not reduce the error: keep the best
			 * phi and stop stepping, rather than rejecting the same step forever. */
			if (alpha < weight / 64.0) alpha = 0.0;
			printf("NR: merit rose (%f > %f), back to previous phi with alpha = %f\n", E_rose, E_prev, alpha);
		}
		else {
			/* The n_e grid resized, so the stored residual does not fit F_vec and the densities in hand
			 * belong to the phi we just threw away. Take no step at all this call: the next one evaluates
			 * the restored phi on its own grid and everything is consistent again. Retrying with a stale
			 * F_vec is what let a bad step double the error and blow the domain out from 70 to 180. */
			alpha = 0.0;
			hold_step = 1;
			printf("NR: merit rose (%f > %f) and the n_e grid resized %d -> %d; restoring phi, no step\n",
			       E_rose, E_prev, n_prev, ninner);
		}
	}
	else
	{
		alpha = fmin(weight, fmax(2.0 * alpha_prev, weight / 64.0));   /* fmax: resume stepping after a stall */
		bscale = fmin(1.0, 2.0 * bscale);
	}
#endif

	/* DS_SLP: in the endgame the step comes from an LP on this J (slp_solve). Keep a copy, since the SVD
	 * below overwrites J with U. alpha as the accept/reject logic left it scales the LP's trust region. */
	double alpha_sg = alpha;
	int slp_step = 0;
#if DS_SLP
	int use_slp = (fabs(phiW_impose - phi0_before) < DS_WALL_TOL);
	gsl_matrix *Jslp = NULL;
	if (use_slp) { Jslp = gsl_matrix_alloc(ninner, ninner); gsl_matrix_memcpy(Jslp, J); }
#endif
	gsl_vector_view rhs = gsl_vector_view_array(F_vec, ninner);
	gsl_permutation *p  = gsl_permutation_alloc(ninner);
	gsl_matrix *Vsv = NULL;
	gsl_vector *Ssv = NULL;
	/* Levenberg-Marquardt (LM_LAMBDA > 0): d = V diag(s/(s^2 + lambda^2)) U^T F rather than J^{-1} F.
	 * Directions with singular value s >> lambda get the Newton step unchanged; the near-singular smooth
	 * mode (s ~ 0.02-0.16 in the cases looked at, the next one ~ 1) is damped instead of dominating the
	 * step, so the trust region no longer has to shrink the whole step to accommodate it. */
	/* Once the wall sits on its target the continuation is over and the remaining residual is the part that
	 * lives in the near-singular mode, so damp it less (LM_LAMBDA_END): the trust region and the accept/
	 * reject on the rms still guard the step. */
	double lm_lam = LM_LAMBDA;
	if (LM_LAMBDA_END > 0.0 && fabs(phiW_impose - phi0_before) < DS_WALL_TOL) lm_lam = LM_LAMBDA_END;
	/* Row weights for the max-targeted solve (DS_MINIMAX). Fsolve is the right-hand side the damped solve
	 * actually uses, so the diagnostics below project the same vector the step was computed from. */
	double *sw = malloc(ninner * sizeof(double));
	gsl_vector *Fw = gsl_vector_alloc(ninner);
	for (i = 0; i < ninner; i++) { sw[i] = 1.0; gsl_vector_set(Fw, i, F_vec[i]); }
#if DS_MINIMAX
	if (lm_lam > 0.0) {
		double mxr = 0.0;
		for (i = 0; i < ninner; i++) mxr = fmax(mxr, fabs(F_vec[i] * invgammasq / ni_grid[i+1]));
		int nlow = 0;
		for (i = 0; i < ninner && mxr > 0.0; i++) {
			double w = pow(fabs(F_vec[i] * invgammasq / ni_grid[i+1]) / mxr, DS_MINIMAX_P - 2.0);
			if (w < DS_MINIMAX_WFLOOR) { w = DS_MINIMAX_WFLOOR; nlow++; }
			sw[i] = sqrt(w);
			for (j = 0; j < ninner; j++) gsl_matrix_set(J, i, j, sw[i] * gsl_matrix_get(J, i, j));
			gsl_vector_set(Fw, i, sw[i] * F_vec[i]);
		}
		/* Normalise to unit rms. Unnormalised, most rows sit at sqrt(floor) ~ 0.03, the weighted matrix's
		 * singular values drop ~30x, and LM_LAMBDA (calibrated on the unweighted matrix) then damps nearly
		 * everything: steps of 1e-4 that could not outrun the MP drift. Normalised, the worst rows get
		 * weight > 1 and the rest < 1, and lambda keeps its meaning. */
		double ss = 0.0;
		for (i = 0; i < ninner; i++) ss += sw[i]*sw[i];
		ss = sqrt(ss / ninner);
		for (i = 0; i < ninner && ss > 0.0; i++) {
			for (j = 0; j < ninner; j++) gsl_matrix_set(J, i, j, gsl_matrix_get(J, i, j) / ss);
			gsl_vector_set(Fw, i, gsl_vector_get(Fw, i) / ss);
			sw[i] /= ss;
		}
		printf("NR: minimax row weights, %d of %d rows at the floor %g, rms-normalised (max weight %.2f)\n",
		       nlow, ninner, DS_MINIMAX_WFLOOR, 1.0/ss);
	}
#endif
	const double *Fsolve = Fw->data;
	if (lm_lam > 0.0) {
		int nd = 0;
		Vsv = gsl_matrix_alloc(ninner, ninner);
		Ssv = gsl_vector_alloc(ninner);
		gsl_vector *work = gsl_vector_alloc(ninner);
		gsl_linalg_SV_decomp(J, Vsv, Ssv, work);   /* J now holds U */
		gsl_vector_free(work);
		lm_solve(J, Vsv, Ssv, lm_lam, Fw, dphi_gsl);
		for (i = 0; i < ninner; i++) if (gsl_vector_get(Ssv, i) < lm_lam) nd++;
		printf("NR: Levenberg-Marquardt, lambda = %g, smallest singular values %.3g %.3g, %d damped below lambda\n",
		       lm_lam, gsl_vector_get(Ssv, ninner-1), gsl_vector_get(Ssv, ninner-2), nd);
		/* Where the residual lives. A small singular value only stalls the solve if the residual has a
		 * component along the matching LEFT singular vector: removing u_k's component needs a step of
		 * |u_k^T F| / s_k along v_k, so a large projection on a near-null mode is a residual the linear
		 * model cannot remove at all, while a small one means the residual sits in the well-conditioned
		 * subspace and the stall is not a conditioning problem. */
		{
			double Fn = 0.0;
			for (i = 0; i < ninner; i++) Fn += Fsolve[i]*Fsolve[i];
			Fn = sqrt(Fn);
			if (Fn > 0.0) {
				int k, kk;
				double tail = 0.0;
				for (k = 0; k < 3 && k < ninner; k++) {
					int col = ninner-1-k;
					double pk = 0.0;
					for (kk = 0; kk < ninner; kk++) pk += gsl_matrix_get(J, kk, col)*Fsolve[kk];
					printf("NR:   mode %d: s = %.3g, |u^T F|/|F| = %.3g, step needed = %.3g\n",
					       col, gsl_vector_get(Ssv, col), fabs(pk)/Fn, fabs(pk)/fmax(gsl_vector_get(Ssv, col), 1e-300));
				}
				for (k = 0; k < 5 && k < ninner; k++) {
					int col = ninner-1-k;
					double pk = 0.0;
					for (kk = 0; kk < ninner; kk++) pk += gsl_matrix_get(J, kk, col)*Fsolve[kk];
					tail += pk*pk;
				}
				printf("NR:   bottom 5 modes hold %.1f%% of |F|^2\n", 100.0*tail/(Fn*Fn));
			}
			/* and where the max relative residual (the convergence test) actually sits */
			double rmx = 0.0; int imx = 0;
			for (i = 1; i < ninner; i++) {
				double r = fabs(F_vec[i]*invgammasq/ni_grid[i+1]);
				if (r > rmx) { rmx = r; imx = i; }
			}
			printf("NR:   max relative residual %.5f at x = %.3f (i = %d of %d)\n", rmx, x_grid[imx+1], imx, ninner);
		}

		/* Hand back the weakest right singular vector for the line search in GYRAZE.c (DS_LINESEARCH).
		 * SV_decomp orders the singular values descending, so this is the last column of V. It is the
		 * direction the LM solve suppresses most, and the one the residual at a stall lives in. Scaled to
		 * max 1 so the amplitude the caller probes is directly a change in phi, and tapered to zero over
		 * the last points so the asymptotic tail attached below is left alone. */
		if (ls_v != NULL) {
			double vm = 0.0; int imax = 0;
			for (i = 0; i < size_phigrid; i++) ls_v[i] = 0.0;
			for (i = 0; i < ninner; i++) {
				double a = fabs(gsl_matrix_get(Vsv, i, ninner-1));
				if (a > vm) { vm = a; imax = i; }
			}
			if (ls_tmodel != NULL) *ls_tmodel = 0.0;
			if (vm > 0.0) {
				double sg = (gsl_matrix_get(Vsv, imax, ninner-1) < 0.0) ? -1.0 : 1.0;   /* fix the sign */
				/* Amplitude the linear model wants along this mode, in the same units as ls_v (max 1).
				 * The Newton step is sum_i (u_i^T F / s_i) v_i, and ls_v = sg*v/vm, so the component
				 * along it is (u^T F / s)*sg*vm. */
				if (ls_tmodel != NULL) {
					double pk = 0.0;
					for (i = 0; i < ninner; i++) pk += gsl_matrix_get(J, i, ninner-1)*Fsolve[i];
					*ls_tmodel = sg * vm * pk / fmax(gsl_vector_get(Ssv, ninner-1), 1e-300);
				}
				int ntap = DS_LS_TAPER;
				if (ntap > ninner/2) ntap = ninner/2;
				for (i = 0; i < ninner; i++) {
					double w = 1.0;
					if (ntap > 0 && i > ninner-1-ntap)
						w = 0.5*(1.0 - cos(M_PI*(double)(ninner-1-i)/(double)ntap));
					ls_v[i+1] = w * sg * gsl_matrix_get(Vsv, i, ninner-1) / vm;
				}
			}
		}
	}
	else {
		gsl_linalg_LU_decomp(J, p, &s);
		gsl_linalg_LU_solve(J, p, &rhs.vector, dphi_gsl);
		if (ls_v != NULL) for (i = 0; i < size_phigrid; i++) ls_v[i] = 0.0;   /* no SVD, no direction */
		if (ls_tmodel != NULL) *ls_tmodel = 0.0;
	}

	/* The right-hand side is linear in the wall value, so split the step into the response to the wall
	 * change (d_bc, from row 0 alone) and the rest (dphi_gsl becomes d_rest). The wall and its own
	 * response are applied together (beta below), so the wall reaches its target without the phi''
	 * spike that moving the wall ahead of the interior would cause; alpha and the trust region then
	 * act on d_rest only. */
	gsl_vector *d_bc = gsl_vector_calloc(ninner);
	gsl_vector_set(d_bc, 0, -(phiW_impose - phi0_before) / deltaxsq);
	if (lm_lam > 0.0) {
		gsl_vector *b = gsl_vector_alloc(ninner);
		gsl_vector_memcpy(b, d_bc);
		for (i = 0; i < ninner; i++) gsl_vector_set(b, i, sw[i] * gsl_vector_get(b, i));   /* same row weights */
		lm_solve(J, Vsv, Ssv, lm_lam, b, d_bc);
		gsl_vector_free(b);
	}
	else gsl_linalg_LU_svx(J, p, d_bc);
	gsl_vector_sub(dphi_gsl, d_bc);
	free(sw); gsl_vector_free(Fw);
	double bmax = 0.0, dW = phiW_impose - phi0_before;
	for (i = 0; i < ninner; i++) bmax = fmax(bmax, fabs(gsl_vector_get(d_bc, i)));
	/* The same trust region as the rest of the step: the response to the wall change reaches into the
	 * interior (through the near-singular mode), where even a change the size of the wall step is far too
	 * large. A large gap (after a restart) closes over several steps, with the wall and its response still
	 * moving together; the small per-iteration drift of the target closes in one (beta = 1). */
	double beta = (bmax > DS_DPHIMAX) ? DS_DPHIMAX / bmax : 1.0;
#if NR_SAFESTEP == 1
	beta *= bscale;
	if (hold_step) beta = 0.0;   /* restored phi: leave the wall alone too, so nothing moves this call */
#endif

	/* Hand back the Newton direction (scaled to max 1) for GYRAZE.c to measure along next call */
	if (dir_h != NULL) {
		double dm = 0.0;
		for (i = 0; i < ninner; i++) dm = fmax(dm, fabs(gsl_vector_get(dphi_gsl, i)));
		for (i = 0; i < size_phigrid; i++) dir_h[i] = 0.0;
		for (i = 0; i < ninner && dm > 0.0; i++) dir_h[i+1] = gsl_vector_get(dphi_gsl, i) / dm;
	}

#if NR_SAFESTEP == 1
	/* Trust region: limit the largest change in phi per step (smooth modes can be near-singular) */
	double dmax = 0.0;
	for (i = 0; i < ninner; i++) dmax = fmax(dmax, fabs(gsl_vector_get(dphi_gsl, i)));
	if (alpha * dmax > DS_DPHIMAX) alpha = DS_DPHIMAX / dmax;
	/* ...and on the combined step, wall part included */
	double jmax = 0.0;
	for (i = 0; i < ninner; i++) jmax = fmax(jmax, fabs(beta * gsl_vector_get(d_bc, i) + alpha * gsl_vector_get(dphi_gsl, i)));
	if (jmax > DS_DPHIMAX) { alpha *= DS_DPHIMAX / jmax; beta *= DS_DPHIMAX / jmax; }
	printf("NR: alpha = %f, max|dphi| = %f (merit %f); wall %f -> %f (beta = %f)\n", alpha, alpha * dmax, E_act, phi0_before, phi0_before + beta*dW, beta);

#if DS_SLP
	if (use_slp && !hold_step && alpha_sg > 0.0) {
		double Delta = DS_SLP_DELTA * fmin(1.0, alpha_sg / weight), tpred = 0.0, apred = 0.0, dm = 0.0;
		int phase = 0;
		double *sv = malloc(ninner * sizeof(double)), *bd = malloc(ninner * sizeof(double));
		double *dl = malloc(ninner * sizeof(double));
		for (i = 0; i < ninner; i++) {
			sv[i] = invgammasq / ni_grid[i+1];
			bd[i] = Delta * fmin(1.0, (x_grid[e] - x_grid[i+1]) / DS_SLP_TAPER);   /* 0 at the last unknown */
		}
		/* mean over rows 0..ninner-2 = nodes 1..size_ngrid-2, the ones error_Poisson averages */
		if (slp_solve(Jslp, F_vec, sv, ninner, ninner - 1, bd, DS_SLP_SMOOTH, DS_SLP_AVGCAP, dl, &tpred, &apred, &phase)) {
			for (i = 0; i < ninner; i++) { gsl_vector_set(dphi_gsl, i, dl[i]); dm = fmax(dm, fabs(dl[i])); }
			gsl_vector_set_zero(d_bc);   /* F_vec was built with the wall at its target: the LP already
			                              * accounts for moving it there, so the wall goes all the way */
			beta = 1.0;
			slp_step = 1;
			alpha = alpha_sg;            /* the accept/reject bookkeeping tracks the trust-region scale */
			printf("NR: SLP step (%s), Delta = %.2e, max|dphi| = %.2e, predicted max|r| %.6f mean %.6f (now %.6f)\n",
			       phase == 1 ? "A: mean, max held" : "B: max, mean capped", Delta, dm, tpred, apred, E_act);
		}
		else printf("NR: SLP LP not solved, keeping the LM step\n");
		free(sv); free(bd); free(dl);
	}
	if (Jslp) gsl_matrix_free(Jslp);
#endif
	if (F_prev == NULL || nF_alloc != ninner) {
		free(F_prev); F_prev = malloc(ninner * sizeof(double)); nF_alloc = ninner;
	}
	for (i = 0; i < size_phigrid; i++) phi_prev[i] = phi_grid[i];
	for (i = 0; i < ninner; i++) F_prev[i] = F_vec[i];
	E_prev = E_act; alpha_prev = alpha; n_prev = ninner; phiW_prev = phiW_impose; phi0_prev = phi0_before;
#else
	/* Backtracking: start at alpha=weight, halve until max relative Poisson error decreases,
	 * with the densities frozen at the current phi. Uses the whole step, as originally. */
	gsl_vector_add(dphi_gsl, d_bc); gsl_vector_set_zero(d_bc);
	int bt;
	for (bt = 0; bt < 30 && alpha > weight / 10.; bt++) {
		double Enew = 0.0;
		for (i = 0; i < ninner; i++) {
			double phi_l = (i == 0) ? phi0_before + alpha*(phiW_impose - phi0_before)
			                        : phi_grid[i] + alpha*gsl_vector_get(dphi_gsl, i-1);
			double phi_c =                                         phi_grid[i+1] + alpha*gsl_vector_get(dphi_gsl, i);
			double phi_r = (i == ninner-1) ? phi_grid[ninner+1] : phi_grid[i+2] + alpha*gsl_vector_get(dphi_gsl, i+1);
			double phipp_trial = (phi_r - 2.0*phi_c + phi_l) / deltaxsq;
			double dev = fabs((-ne_grid[i+1] + phipp_trial*invgammasq)/ni_grid[i+1] + 1.0);
			if (dev > Enew) Enew = dev;
		}
		if (Enew <= E0) break;
		alpha *= 0.5;
	}
	printf("NR backtracking: alpha = %f after %d halvings (E0=%f)\n", alpha, bt, E0);
	beta = alpha;   /* original behaviour: the wall and the whole step move by alpha */
#endif

	/* Apply the step: the wall and its response by beta, the rest by alpha */
	phi_grid[0] = phi0_before + beta * dW;
	for (i = 0; i < ninner; i++) {
		temp = phi_grid[i+1] + beta * gsl_vector_get(d_bc, i) + (slp_step ? 1.0 : alpha) * gsl_vector_get(dphi_gsl, i);
		if (temp > 0.0)
			printf("WARNING: phi > 0.0 at x = %f, non-monotonic in Debye sheath\n", x_grid[i+1]);
		if (temp < phi_grid[i])
			printf("WARNING: phi non-monotonic at x = %f\n", x_grid[i+1]);
		phi_grid[i+1] = temp;
	}

	/* Asymptotic extension from size_ngrid to size_phigrid.
	 * phi_grid[size_ngrid-1] is now solved by NR; anchor the power law there. */
	printf("In Debye sheath NR, asymptotic result starts at x = %f\n\n", x_grid[size_ngrid]);
	/* phi = phi0 (x + CC)^pdec, matched in value and slope at the edge, so the join stays smooth when the
	 * edge point later becomes interior (the n_e grid end moves by a point from one iteration to the next);
	 * the slope is second-order one-sided. If that is not a valid decaying tail (slope <= 0, or x + CC <= 0),
	 * use the CC of the boundary condition (value-matched), rather than clamping CC to 0, which put a
	 * corner in phi at the join. */
	int extend = 1;
	double slope = (3.0*phi_grid[e] - 4.0*phi_grid[e-1] + phi_grid[e-2]) / (2.0*(x_grid[e] - x_grid[e-1]));
	if (phi_grid[e] < 0.0 && slope > 0.0 && pdec*phi_grid[e]/slope > 0.0) {
		CC = pdec*phi_grid[e]/slope - x_grid[e];
		CC_prev = CC; have_CC = 1;
	}
	else if (rho > 0.0) {   /* the tail used as boundary condition above */
		CC = CCb;
		printf("Warning: no valid tail slope at the edge, keeping CC = %f\n", CC);
	}
	else {
		extend = 0;
		printf("Warning: no valid tail at the edge yet, leaving the tail unchanged\n");
	}
	if (extend) {
		printf("In Debye sheath NR CC = %f\n\n", CC);
		phi0 = phi_grid[e] / pow(x_grid[e] + CC, pdec);
		printf("pdec = %f\nphi0 = %f\n", pdec, phi0);
		for (i = size_ngrid; i < size_phigrid; i++)
			phi_grid[i] = phi0 * pow(x_grid[i] + CC, pdec);
	}

	gsl_permutation_free(p);
	gsl_matrix_free(J);
	gsl_vector_free(dphi_gsl);
	gsl_vector_free(d_bc);
	if (Vsv != NULL) { gsl_matrix_free(Vsv); gsl_vector_free(Ssv); }
	free(F_vec);
	free(ne_corr_total);
}

/* Linearized BVP correction for DS phi when restarting from a different phi_wall.
 * Solves: delta_phi'' = gamma2*(f(x)*delta_phi - g(x)*delta_phi_wall)
 * where f = ne_corr_total, g = ni_corr (sumni_DS_corr).
 * BC: delta_phi[0] = phiW_target - phi_DSgrid[0], delta_phi[N_bvp] = 0.
 * Adds the correction to phi_DSgrid[0..N_bvp]. */
void correct_phi_DS_restart(double *x_DSgrid, double *phi_DSgrid, int size_phiDSgrid,
                             double invgammasq, double v_cutDS,
                             double *ne_DSgrid, double *ne_DSgrid_corr_delta, double *ne_DSgrid_corr_chiM,
                             double *ni_corr, int N_bvp)
{
	double phiW_target = -0.5 * v_cutDS * v_cutDS;
	double delta_phi_wall = phiW_target - phi_DSgrid[0];
	if (fabs(delta_phi_wall) < 1e-2) return;
	printf("Restart DS BVP: delta_phi_wall = %.6f (phi_DSgrid[0]=%.6f -> phiW_target=%.6f)\n",
	       delta_phi_wall, phi_DSgrid[0], phiW_target);
	double gamma2 = 1.0 / invgammasq;
	double dx  = x_DSgrid[1] - x_DSgrid[0];
	double dx2 = dx * dx;
	int n = N_bvp - 1;
	if (n <= 0) { printf("Restart DS BVP: domain too small, skipping.\n"); return; }
	gsl_matrix      *A    = gsl_matrix_alloc(n, n);
	gsl_vector      *b    = gsl_vector_alloc(n);
	gsl_vector      *dphi = gsl_vector_alloc(n);
	gsl_permutation *p    = gsl_permutation_alloc(n);
	/* ne_DSgrid_corr_delta[0] is NaN; set it from index 1 */
	ne_DSgrid_corr_chiM[0] = ne_DSgrid_corr_chiM[1];

	/* Normalize ne_corr_total so it matches ni_corr at the outer boundary (idx = N_bvp-1) */
	double ne_total_outer = ne_DSgrid[N_bvp-1] + ne_DSgrid_corr_delta[N_bvp-1] + ne_DSgrid_corr_chiM[N_bvp-1];
	//double ne_total_outer = ne_DSgrid[N_bvp-1];
	printf("ne_total_outer is %f\n", ne_total_outer);
	double scale = (fabs(ne_total_outer) > 1e-14) ? ni_corr[N_bvp-1] / ne_total_outer : 1.0;
	printf("ni_corr[-1] = %f\n", ni_corr[N_bvp-1]);
	printf("scale is %f\n", scale);
	gsl_matrix_set_zero(A);
	for (int j = 0; j < n; j++) {
		int idx = j + 1;
		double f_j = scale * (ne_DSgrid[idx] + ne_DSgrid_corr_delta[idx] + ne_DSgrid_corr_chiM[idx]) - ni_corr[idx];
		//double f_j = scale * (ne_DSgrid[idx]) - ni_corr[idx];
		double g_j = ne_DSgrid_corr_chiM[idx];
		//double g_j = 0.0;
		gsl_matrix_set(A, j, j, -2.0/dx2 - gamma2*f_j);
		if (j > 0)   gsl_matrix_set(A, j, j-1, 1.0/dx2);
		if (j < n-1) gsl_matrix_set(A, j, j+1, 1.0/dx2);
		double rhs = -gamma2 * g_j * delta_phi_wall;
		if (j == 0) rhs -= delta_phi_wall / dx2;
		gsl_vector_set(b, j, rhs);
		printf("ne_DSgrid[%d] = %f, ni_corr[%d] = %f, rhs = %f\n", j, ne_DSgrid[idx], j, ni_corr[idx], rhs);
	}
	int signum;
	gsl_linalg_LU_decomp(A, p, &signum);
	gsl_linalg_LU_solve(A, p, b, dphi);
	phi_DSgrid[0] = phiW_target;
	for (int j = 0; j < n; j++)
		phi_DSgrid[j + 1] += gsl_vector_get(dphi, j);
	printf("Restart DS BVP correction applied (N_bvp=%d).\n", N_bvp);
	gsl_matrix_free(A);
	gsl_vector_free(b);
	gsl_vector_free(dphi);
	gsl_permutation_free(p);
}
