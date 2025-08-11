/* 
   AUTHORS: Alessandro Geraldini and Robbie Ewart

   DESCRIPTION: Main file of GYRAZE
		solves the magnetised sheath = magnetic presheath + Debye sheath 
		using a grazing-angle gyrokinetic approach
 
   MODIFIED: 3 May 2025 by Alessandro Geraldini
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_bessel.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "mps.h"

#define TEST_EL 0

const char *strqty[8] = {"alpha=","gamma=","nspec=", "ni:ne=", "Ti:Te=","mi:me=","jwall=","pwall="};
const int lenstrqty = 6;

void densionDS(double alpha, double TiovTe, double *Bohm, double *ni_DS, double *phi_DS, double phi0, double **FF, double *mu, double *Uminmu, double *vy, double *mu_op, double *chiM, double *twopidmudvy, int size_phi, int size_mu, int size_U, int size_op_i) {
/* 
This function calculates the ion density in the Debye sheath, which exploits a 1D acceleration of ions in the direction normal to the wall.
It is not a particularly accurate density calculation; a trapezium rule is used.
Two different alternative methods are includeded (method 1 and method 2 are the same, method 2 will be deleted, method 3 is different but less accurate)
INPUTS: alpha, temperature ratio TiovTe, potential in Debye sheath phi_DS, potential at the Debye sheath entrance relative to the magnetic presheath phi0 = phi_mp(0),
        distribution function entering the magnetised sheath, its arguments mu and Uminmu, vy at x=0 in the magnetic presheath which is effectively xbar, effective potential maximum chiM, 
        open orbit integral twopidmudvy, and the size of the phi grid, the size of the mu grid, the size of the Uminmu grid.
OUTPUT: density profile ni_DS
*/
	int i, j, k, count;
	double Bohm1, deltaUperp, halfVx0sq, intgrd, intgrdold, vzk, vzkm, n_inf;
	double intgrdBohmold, intgrdBohm;
	double Fk, Fkm1;
	n_inf = 0.0;
	Bohm1 = 0.0;
	intgrd = 0.0;
	intgrdBohm = 0.0;
	for (j=0; j< size_op_i; j++) { // the open orbit integral in vy = Omega * xbar 
		//printf("vy = %f\tmuop = %f\tchiM = %f\tdmudvy=%f\n", vy[j], mu_op[j], chiM[j], twopidmudvy[j]);
		intgrdold = intgrd;
		intgrdBohmold = intgrdBohm;
		intgrd=0.0;
		intgrdBohm=0.0;
		for (k=1; k < size_U; k++) { // the open orbit integral in vz or U
			halfVx0sq = chiM[j] - 0.5*vy[j]*vy[j] - phi0/TiovTe;
			deltaUperp = mu_op[j] - chiM[j] ;
			vzk = sqrt(2.0*(deltaUperp + Uminmu[k]));
			vzkm = sqrt(2.0*(deltaUperp + Uminmu[k-1]));
			Fk = bilin_interp(mu_op[j], Uminmu[k], FF, mu, Uminmu, size_mu, size_U, -1, -1);
			Fkm1 = bilin_interp(mu_op[j], Uminmu[k-1], FF, mu, Uminmu, size_mu, size_U, -1, -1);
			// intgrd is adding contributions of F * Delta vx
			intgrd += ( (sqrt(2.0*(halfVx0sq + alpha*vzk*twopidmudvy[j])) - sqrt(2.0*halfVx0sq)) * Fk + (sqrt(2.0*(halfVx0sq + alpha*vzkm*twopidmudvy[j])) - sqrt(2.0*halfVx0sq)) * Fkm1 ) * 0.5 * ( vzk - vzkm );
			// intgrdBohm is adding contributions of F * Delta (1/vx) 
			intgrdBohm += ( (-1.0/sqrt(2.0*(halfVx0sq + alpha*vzk*twopidmudvy[j])) + 1.0/sqrt(2.0*halfVx0sq)) * Fk + (-1.0/sqrt(2.0*(halfVx0sq + alpha*vzkm*twopidmudvy[j])) + 1.0/sqrt(2.0*halfVx0sq)) * Fkm1 ) * 0.5 * ( vzk - vzkm );
		}
		if (j > 0) {
			n_inf += (intgrd + intgrdold)*0.5*(vy[j] - vy[j-1]);
			Bohm1 += (intgrdBohm + intgrdBohmold)*0.5*(vy[j] - vy[j-1]);
		}
		else if (j==0) {
			intgrd = 0.0;
			intgrdBohm = 0.0;
		}
	}
	printf("n_inf = %f\n", n_inf);
	*Bohm = Bohm1/n_inf;
	printf("Bohm = %f\n", *Bohm);
	for (i=0; i < size_phi; i++) { 
		ni_DS[i] = 0.0;
		for (j=0; j< size_op_i; j++) {
			intgrdold = intgrd;
			intgrd = 0.0;
			count = 0;
			for (k=1; k < size_U; k++) {
				halfVx0sq = chiM[j] - 0.5*vy[j]*vy[j] - phi_DS[i]/TiovTe - phi0/TiovTe;
				deltaUperp = mu_op[j] - chiM[j];
				if (phi0 == 0.0) {
					halfVx0sq = -phi_DS[i]; 
					deltaUperp = 0.0;
				}
				vzk = sqrt(2.0*(deltaUperp + Uminmu[k]));
				vzkm = sqrt(2.0*(deltaUperp + Uminmu[k-1]));
				Fk = bilin_interp(mu_op[j], Uminmu[k], FF, mu, Uminmu, size_mu, size_U, -1, -1);
				Fkm1 = bilin_interp(mu_op[j], Uminmu[k-1], FF, mu, Uminmu, size_mu, size_U, -1, -1);
				intgrd += ( (sqrt(2.0*(halfVx0sq + alpha*vzk*twopidmudvy[j])) - sqrt(2.0*halfVx0sq)) * Fk + (sqrt(2.0*(halfVx0sq + alpha*vzkm*twopidmudvy[j])) - sqrt(2.0*halfVx0sq)) * Fkm1 ) * 0.5 * ( vzk - vzkm );
				if ((count == 0) && (intgrd != intgrd) ) {
					count = 1;
					//printf("???%f %f %f %f\n", halfVx0sq, deltaUperp, vzk, vzkm);
				}
			}
			if (j > 0) {
				ni_DS[i] += (intgrd + intgrdold)*0.5*(vy[j] - vy[j-1]);
			}
			else if (j==0) {
				intgrd = 0.0;
			}
		}
		ni_DS[i] /= n_inf;
		//printf("ni_DS[i] = %f, phi_DS[i] = %f\n", ni_DS[i], phi_DS[i]);
		//printf("ni_DS = %f\tphi_DS = %f\n", ni_DS[i], phi_DS[i]);
		if (i == size_phi-1) 
			printf("in densionDS: derivative wrt phi is dndphi = %f\n", (ni_DS[i] - ni_DS[i-1])/(phi_DS[i] - phi_DS[i-1]));
	}
	return;
}

void Bohmeval(double *Bohm, double alpha, double TiovTe, double phi0, double **FF, double *mu, double *Uminmu, double *vy, double *mu_op, double *chiM, double *twopidmudvy, int size_mu, int size_U, int sizemuop) {
/*
This function evaluates the Bohm integral for the ions: int f d^3v/v_x^2
*/
	int i, j, k, sizevx=250;
	double deltaUperp, halfVx0sq, intgrdn=0.0, intgrdnold, intgrdB=0.0, intgrdBold, vzk, vzkm;
	double Bohmval=0.0, Bohmval2=0.0, Bohmval3=0.0,  nval = 0.0, fvx[sizevx], vx[sizevx];
	double maxf = 0.0, extraphi = 0.0; //extraphi = any non-zero value is just a TEST;
	double Fk, Fkm1;
	for (j=1; j< sizemuop; j++) {
		intgrdBold = intgrdB;
		intgrdnold = intgrdn;
		intgrdn = 0.0;
		intgrdB = 0.0;
		for (k=1; k < size_U; k++) {
			halfVx0sq = chiM[j] - 0.5*vy[j]*vy[j] - phi0/TiovTe + extraphi;
			deltaUperp = mu_op[j] - chiM[j];
			vzk = sqrt(2.0*(deltaUperp + Uminmu[k]));
			vzkm = sqrt(2.0*(deltaUperp + Uminmu[k-1]));
			Fk = bilin_interp(mu_op[j], Uminmu[k], FF, mu, Uminmu, size_mu, size_U, -1, -1);
			Fkm1 = bilin_interp(mu_op[j], Uminmu[k-1], FF, mu, Uminmu, size_mu, size_U, -1, -1);
			intgrdB += ( (-1.0/sqrt(2.0*(halfVx0sq + alpha*vzk*twopidmudvy[j])) + 1.0/sqrt(2.0*halfVx0sq)) * Fk + (-1.0/sqrt(2.0*(halfVx0sq + alpha*vzkm*twopidmudvy[j])) + 1.0/sqrt(2.0*halfVx0sq)) * Fkm1 ) * 0.5 * ( vzk - vzkm );
			intgrdn += ( ( (sqrt(2.0*(halfVx0sq + alpha*vzk*twopidmudvy[j])) - sqrt(2.0*halfVx0sq)) * Fk + (sqrt(2.0*(halfVx0sq + alpha*vzkm*twopidmudvy[j])) - sqrt(2.0*halfVx0sq)) * Fkm1 ) * 0.5 * ( vzk - vzkm ) );
		}
		if (j > 1) {
			nval += (intgrdn + intgrdnold)*0.5*(vy[j] - vy[j-1]);
			Bohmval += (intgrdB + intgrdBold)*0.5*(vy[j] - vy[j-1]);
		}
		else if (j==1) {
			intgrdn = 0.0;
			intgrdB = 0.0;
		}

	}
	//printf("density = %f\n", nval);
	//printf("true Bohm = %f\n", Bohmval/nval);
	*Bohm= Bohmval/nval;
	Bohmval = 0.0;
	Bohmval2 = 0.0;
	Bohmval3 = 0.0;
	printf("n = %f, ", nval);
	nval = 0.0;
	for (i=0; i<sizevx; i+=1) {
		vx[i] = i*5.0/sizevx;
		fvx[i] = 0.0;
		for (j=0; j< sizemuop; j++) {
			intgrdBold = intgrdB;
			intgrdnold = intgrdn;
			intgrdn = 0.0;
			intgrdB = 0.0;
			for (k=1; k < size_U; k++) {
				halfVx0sq = chiM[j] - 0.5*vy[j]*vy[j] - phi0/TiovTe;
				deltaUperp = mu_op[j] - chiM[j] - 0.5*vx[i]*vx[i] + halfVx0sq;
				deltaUperp = mu_op[j] - chiM[j];
				vzk = sqrt(2.0*(deltaUperp + Uminmu[k]));
				vzkm = sqrt(2.0*(deltaUperp + Uminmu[k-1]));
				Fk = bilin_interp(mu_op[j], Uminmu[k], FF, mu, Uminmu, size_mu, size_U, -1, -1);
				Fkm1 = bilin_interp(mu_op[j], Uminmu[k-1], FF, mu, Uminmu, size_mu, size_U, -1, -1);
				if ( (0.5*vx[i]*vx[i] > halfVx0sq) && (0.5*vx[i]*vx[i] < halfVx0sq +  alpha*vzk*twopidmudvy[j]) && (0.5*vx[i]*vx[i] < halfVx0sq +  alpha*vzkm*twopidmudvy[j]) ) 
					intgrdn += 0.5*( (Fk + Fkm1)*(vzk - vzkm) );
				
				else if ( (0.5*vx[i]*vx[i] > halfVx0sq) && (0.5*vx[i]*vx[i] < halfVx0sq +  alpha*vzk*twopidmudvy[j]) && (0.5*vx[i]*vx[i] > halfVx0sq +  alpha*vzkm*twopidmudvy[j]) ) 
					intgrdn += 0.5*( (Fk + 0.0)*(vzk - sqrt((0.5*vx[i]*vx[i]-halfVx0sq)/(alpha*twopidmudvy[j])) ) );
			}
			if (j > 1) 
				fvx[i] += (intgrdn + intgrdnold)*0.5*(vy[j] - vy[j-1]);
			else if (j==1)
				intgrdn = 0.0;
			if (fvx[i] > maxf)
				maxf = fvx[i];
		}
		//printf("fvx[%d] = %f\n", i, fvx[i]);
		if (i != 0) {
			Bohmval += ( 2.0*(fvx[i] - fvx[i-1])/(vx[i] + vx[i-1])) ;
			//Bohmval2 += ( (vx[i] - vx[i-1])*0.5*(fvx[i]/(vx[i] * vx[i]) + fvx[i-1]/(vx[i-1] * vx[i-1]+TINY) ) );
			Bohmval2 += ( (vx[i] - vx[i-1])*(fvx[i] - fvx[0] - vx[i]*(fvx[1]-fvx[0])/vx[1])/(vx[i] * vx[i]) );
			nval += ( 0.5*(fvx[i] + fvx[i-1])*(vx[i] - vx[i-1]) );
		}
		if (i>1) {
			Bohmval3 += (1.0 - log(vx[i-1]))*(fvx[i] - 2.0*fvx[i-1] + fvx[i-2])/(vx[i]-vx[i-1]);
		}
	}
	printf("maxf = %f\n", maxf);
	printf("nval = %f\n", nval);
	printf("reconstructed Bohm = (%f, %f, %f) same value calculated 3 times using slightly different numerical method for sanity check\n", Bohmval/nval, Bohmval2/nval, Bohmval3/nval);
	return;
}

//void Bohmshouldbe(double *Bohmshouldbe, double phi0, double *n_grid, int p_size, double *Phi_point, double **distfunc, double *vpar, double *mu, int size_vpar, int size_mu, double *vpar_cut_lookup, double gamma, double *x_grid) {
void evalBohmshouldbe(double *Bohmshouldbe, double phi0, double **distfunc, double *vpar, double *mu, int size_vpar, int size_mu, double *vpar_cut_lookup) {
	//define variables
	int vi, vi_1, p, len_F; //vi is a counting variable that will be saved for sums over velocity space, vi_1 is a special value in velocity space devoted to the first point, p is a counting variable that will be saved for counting over phi space, len_F saves the number of entries in the distribution;
	double phi;
	double v_max, v_s, v_min;//v_max is the maximum velocity the distribution function will go up to before effectively just reading 0 from there on, v_s is the separation in velocity space between ajacent points. v_cut is the velocity that would be require in order to just reach the plasma wall boundary, v_min is the minimum velocity (the entrance to the presheath) required to reach a given x
	double sqrt_up, sqrt_lo, theta_up, theta_lo, sqrt_cut, theta_cut; //sqrt_up/lo reprresent the upper and lower limits of one strip integral, theta_up/lo are related to the hyperbolic arcsinh of some specific values (see document that hopefully exists)

	double *F, *Fp, *Fpp, **ddistdvpar, **ddistdvpartwo, **ddistdvparthree; // the zeroth first and second derivatives of the distribution function respectively

	//double* n_pre; //n_pre[x] is an array that will contain the values that we expect the integral to give (for a maxwellian input)

	int mu_ind;
	double *nepart, vpar_cut, n_grid[2], dphi = 0.000001;

	ddistdvpar = malloc(size_mu * sizeof(double));
	ddistdvpartwo = malloc(size_mu * sizeof(double));
	ddistdvparthree = malloc(size_mu * sizeof(double));
	len_F = size_vpar;
	v_max = vpar[len_F - 1];
	v_s = v_max / (len_F - 1);

	F = malloc(len_F * sizeof(double));
	Fp = malloc(len_F * sizeof(double));
	Fpp = malloc(len_F * sizeof(double));
	//Generate the first and second derivative of the distribution function
	for (mu_ind = 0; mu_ind < size_mu; mu_ind++) {
		ddistdvpar[mu_ind] = malloc(len_F * sizeof(double));
		ddistdvpartwo[mu_ind] = malloc(len_F * sizeof(double));
		ddistdvparthree[mu_ind] = malloc(len_F * sizeof(double));
		for (vi = 0; vi < len_F; vi++) {
			if ((vi < len_F - 1) && (vi > 0)) {
				ddistdvpar[mu_ind][vi] = (distfunc[mu_ind][vi + 1] - distfunc[mu_ind][vi - 1]) / (2.0 * v_s);
				ddistdvpartwo[mu_ind][vi] = (distfunc[mu_ind][vi + 1] + distfunc[mu_ind][vi - 1] - (2.0 * distfunc[mu_ind][vi])) / (pow(v_s, 2.0));
			}
			else
			{
				if (vi < len_F - 1)
				{
					ddistdvpar[mu_ind][vi] = (distfunc[mu_ind][vi + 1] - distfunc[mu_ind][vi]) / (v_s);
					ddistdvpartwo[mu_ind][vi] = (distfunc[mu_ind][vi + 2] + distfunc[mu_ind][vi] - (2.0 * distfunc[mu_ind][vi + 1])) / (pow(v_s, 2.0));
				}
				else
				{
					if (vi > 0)
					{
						ddistdvpar[mu_ind][vi] = (distfunc[mu_ind][vi] - distfunc[mu_ind][vi - 1]) / (v_s);
						ddistdvpartwo[mu_ind][vi] = (distfunc[mu_ind][vi - 2] + distfunc[mu_ind][vi] - (2.0 * distfunc[mu_ind][vi - 1])) / (pow(v_s, 2.0));
					}
					else
					{
						printf("Fuck");
					}
				}

			}
		}
		for (vi = 0; vi < len_F; vi ++) {
			if ((vi < len_F - 1) && (vi > 0)) 
				ddistdvparthree[mu_ind][vi] = (ddistdvpartwo[mu_ind][vi + 1] - ddistdvpartwo[mu_ind][vi - 1]) / (2.0 * v_s);
			else if (vi==0) 
				ddistdvparthree[mu_ind][vi] = (ddistdvpartwo[mu_ind][vi + 1] - ddistdvpartwo[mu_ind][vi]) / (v_s);
			else if (vi==len_F-1)
				ddistdvparthree[mu_ind][vi] = (ddistdvpartwo[mu_ind][vi] - ddistdvpartwo[mu_ind][vi-1]) / (v_s);
		}
	}

	nepart = malloc(size_mu*sizeof(double));
	// for finite electron gyroradius effects, we need to integrate in mu as well
	for (p=0;p<2;p++) {
		n_grid[p] = 0.0;
		phi = phi0 + dphi*p;
		v_min = sqrt(-2.0 * phi);
		for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
			nepart[mu_ind] = 0.0;
			//vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
			vpar_cut = vpar_cut_lookup[mu_ind] ; 
			//printf("vpar_cut before = %f\n", vpar_cut);
			vpar_cut = sqrt(vpar_cut*vpar_cut +  (-2.0*phi0) + TINY);
			//printf("vpar_cut = %f\n", vpar_cut);
			for (vi = 0; vi < len_F - 1; vi++) {
				F[vi]     = distfunc[mu_ind][vi];
				Fp[vi]    = ddistdvpar[mu_ind][vi];
				Fpp[vi]   = ddistdvpartwo[mu_ind][vi];
			}
			if (vpar_cut >= v_max) // effectively no cut off
			{
				if ((int)floor(v_min / v_s) >= len_F - 1)
				{
					//ne[p] = 0;
					nepart[mu_ind] = 0.0;
					printf("v_min outside velocity grid");
				}
				else
				{
					if (fabs(v_min) < 2.0*TINY)//special case where the integration becomes simple
					{
						for (vi = 0; vi < len_F - 2; vi++)
						{
							nepart[mu_ind] += (F[vi] + F[vi + 1]) * v_s;
						}
					}
					else //regular case
					{
						vi_1 = (int)floor(v_min / v_s);
						sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi));
						theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi)) / v_min);

						nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
						nepart[mu_ind] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));

						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
						{
							sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi);
							sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi);
							theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi) / (v_min));
							theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi) / (v_min));

							nepart[mu_ind] += 2.0 * (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							nepart[mu_ind] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
						}
					}
				}
			}
			else { // now there is a cut-off
				v_min = sqrt(-2.0 * phi);
				if ((int)floor(v_min / v_s) >= len_F - 1)
				{
					nepart[mu_ind] = 0.0;
				}
				else
				{
					if (fabs(v_min) < 1e-9)
					{
						vi_1 = 0;
						sqrt_up = ((vi_1 + 1) * v_s);
						nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
						nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));

						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
						{
							sqrt_up = (vi + 1) * v_s;
							sqrt_lo = (vi * v_s);

							nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
							nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
						}

						for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
						{
							sqrt_up = (vi + 1) * v_s;
							sqrt_lo = (vi * v_s);

							nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
							nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
						}
						vi_1 = (int)floor(vpar_cut / v_s);
						sqrt_up = vpar_cut;
						sqrt_lo = vi_1 * v_s;

						nepart[mu_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
						nepart[mu_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
					}
					else {
						if ((int)floor(v_min / v_s) == (int)floor(vpar_cut / v_s))
						{
							vi_1 = (int)floor(v_min / v_s);
							sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi));
							theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi)) / (v_min));

							nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
							nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

							sqrt_cut = sqrt(pow(vpar_cut, 2.0) + 2.0 * phi);
							theta_cut = asinh(sqrt(pow(vpar_cut, 2.0) + 2.0 * phi) / (v_min));

							nepart[mu_ind] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * vpar_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
							nepart[mu_ind] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (vpar_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));

							for (vi = vi_1 + 1; vi < len_F - 1; vi++)
							{
								sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi);
								sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi);
								theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi) / (v_min));
								theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi) / (v_min));

								nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
								nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							}
						}
						else //the most likely case
						{
							vi_1 = (int) (floor(v_min / v_s));
							sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phi);
							theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phi) / (v_min));
							nepart[mu_ind] = 0.0;
							

							nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
							nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));


							for (vi = vi_1 + 1; vi < len_F - 1; vi++)
							{
								sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi);
								sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi);
								theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi) / (v_min));
								theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi) / (v_min));

								nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
								nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							}

							for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
							{
								sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi);
								sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi);
								theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi) / (v_min));
								theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi) / (v_min));

								nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
								nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

							}

							vi_1 = (int)floor(vpar_cut / v_s);
							sqrt_up = sqrt(pow(vpar_cut, 2.0) + 2.0 * phi);
							sqrt_lo = sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phi);
							theta_up = asinh(sqrt(pow(vpar_cut, 2.0) + 2.0 * phi) / (v_min));
							theta_lo = asinh(sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phi) / (v_min));

							nepart[mu_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							nepart[mu_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));

						}
					}
				}
			}
			if (mu_ind != 0) {
				n_grid[p] += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1])*(mu[mu_ind] - mu[mu_ind-1]);
				//printf("n_grid %f\n", n_grid[p]);
			}
		}
	}
	//printf("n_grid[0,1] = %f,%f\n", n_grid[0], n_grid[1]);
	*Bohmshouldbe = (n_grid[1]/n_grid[0] - 1.0)/dphi;
	//printf("Bohmshouldbe = %f\n", *Bohmshouldbe);

	free(F); free(Fp); free(Fpp);
	for (mu_ind=0; mu_ind< size_mu; mu_ind++){
		 free(ddistdvpar[mu_ind]);
		 free(ddistdvpartwo[mu_ind]);
		 free(ddistdvparthree[mu_ind]);
	}
	free(ddistdvpar);
	free(ddistdvpartwo);
	free(ddistdvparthree);
	free(nepart);
}

void Figen(double ***ffarr, double **Uminmuarr, double **muarr, int num_spec, double *nioverne, double *mioverme, double *TioverTe, int *sizevpar, int *sizevperp, double dvpar, double dvperp) {
	int n, i, j, coldelectrons;
	double mu, Uminmu, ff;
	double u, condition, chodura, normalization;// flow;

	for (n=0; n<num_spec; n++) {
		if (TioverTe[n] > 10.0) coldelectrons = 1;
		else coldelectrons = 0;
		printf("are electrons cold (1=yes, 0=no): %d\n", coldelectrons);
		u = 0.0;
		condition = TioverTe[n]; 
		if (coldelectrons == 0) {
			if (TioverTe[n]<=1.0) {
				chodura = 9999999.0;
				u = 0.0;
				while ( (chodura > condition) || (chodura < condition - 0.05) ) {
					if (chodura > condition)
						u += 0.0001;
					else 
						u -= 0.0001;
					normalization =  sqrt(2.0)*(1 + erf(u/sqrt(2.0)))*(1.0+u*u) + sqrt(2.0/M_PI)*u*exp(-0.5*u*u);
					//printf("normalization = %f\n", normalization);
					chodura = sqrt(2.0)*(1+erf(u/sqrt(2.0)))/(normalization);
					//flow = ( 0.5*(1+0.5*u*u)*exp(-0.5*u*u) + (1.5 + 0.5*u*u)*(sqrt(M_PI/2.0)*u/2.0)*(1 + erf(u/sqrt(2.0))))*4.0/sqrt(2.0*M_PI)/normalization;
					//printf("%f, %f, %f, %f %f\n", normalization, u, flow, chodura, condition);
				}
			}
			else {
				chodura = 0.0;
				u = 0.1;
				while  ( (chodura > condition) || (chodura < condition - 0.05) ) {
					if (chodura > condition)
						u -= 0.001;
					else
						u += 0.001;
					normalization = (1.0/sqrt(2.0))*2.0*(sqrt(M_PI*u) - M_PI*exp(1/u)*(1-erf(1.0/sqrt(u))))/(2.0*u*sqrt(u));
					//printf("normalization = %f\n", normalization);
					chodura = (1.0/sqrt(2.0))*0.5*2.0*M_PI*exp(1.0/u)*(1.0-erf(1.0/sqrt(u)))/(2.0*sqrt(u)*normalization);
					//printf("%f, %f, %f, %f %f\n", normalization, u, flow, chodura, condition);
				}
			}
		}

		if (coldelectrons == 0) {
			if ( TioverTe[n] <= 1.0 ) 
				//normalization =  (1 + erf(u/sqrt(2.0)))*(1.0+u*u) + sqrt(2.0/M_PI)*u*exp(-0.5*u*u);
				normalization =  sqrt(2.0)*(1 + erf(u/sqrt(2.0)))*(1.0+u*u) + sqrt(2.0/M_PI)*u*exp(-0.5*u*u);
				//normalization =  (1 + erf(u))*(1.0+2.0*u*u) + (2.0/sqrt(M_PI))*u*exp(-u*u);
			else
				//normalization =  2.0*(sqrt(M_PI*u) - M_PI*exp(1/u)*(1-erf(1/sqrt(u))))/(2.0*u*sqrt(u));
				normalization = (1.0/sqrt(2.0))*2.0*(sqrt(M_PI*u) - M_PI*exp(1/u)*(1-erf(1.0/sqrt(u))))/(2.0*u*sqrt(u));
		}

		if (coldelectrons == 0) {
			if (TioverTe[n]<=1.0) {
				//printf("TioverTe[%d] = %f <= 1\n", n, TioverTe[n]);
				printf("normalization = %f\n", normalization);
				for (i=0; i<sizevperp[n]; i++) {
					//mu = 0.5*i*dvperp*i*dvperp;
					mu = i*dvperp;
					muarr[n][i] = mu;
					for (j=0; j<sizevpar[n]; j++) { 
						Uminmu = 0.5*j*dvpar*j*dvpar;
						if (i==0) {
							Uminmuarr[n][j] = Uminmu;
							//printf("Uminmu = %f\n", Uminmu);
						}
						ff = (2.0/(M_PI*sqrt(M_PI)))*(1.0/normalization)*(Uminmu)*exp(- Uminmu - mu + u*sqrt(2.0*Uminmu) - 0.5*u*u);
						ffarr[n][i][j] = ff;
						//printf("%f ", ff);
						//printf("%f %d", ffarr[n][i][j], n);
					}
					//printf("\n", ff);
				}
			}
			else {
				//printf("TioverTe[%d] = %f > 1\n", n, TioverTe[n]);
				for (i=0; i<sizevperp[n]; i++) {
					//mu = 0.5*i*dvperp*i*dvperp;
					mu = i*dvperp;
					muarr[n][i] = mu;
					for (j=0; j<sizevpar[n]; j++) { 
						Uminmu = 0.5*j*dvpar*j*dvpar;
						if (i==0)
							Uminmuarr[n][j] = Uminmu;

						ff = (1.0/M_PI)*(1.0/normalization)*((Uminmu)/(1+u*(Uminmu)))*exp(- Uminmu - mu );
						ffarr[n][i][j] = ff;
					}
				}
			}
		}
		else { // coldelectrons == 1:
			printf("TioverTe[%d] = infinite: ELECTRONS ARE COLD\n", n);
			normalization = 1.0;
			for (i=0; i<sizevperp[n]; i++) {
				//mu = 0.5*i*dvperp*i*dvperp;
				mu = i*dvperp;
				muarr[n][i] = mu;
				for (j=0; j<sizevpar[n]; j++) { 
					Uminmu = 0.5*j*dvpar*j*dvpar;
					if (i==0)
						Uminmuarr[n][j] = Uminmu;
					ff = (1.0/(sqrt(2.0)*(M_PI)*sqrt(M_PI)))*exp(- Uminmu - mu );
					ffarr[n][i][j] = ff;
				}
			}
		}
	}

	//for (i=0; i<sizevperp[0]; i++) {
	//	for (j=0; j<sizevpar[0]; j++) { 
	//		printf("%f ", ffarr[0][i][j]);
	//	}
	//	printf("\n");
	//}
	return;
}

void Figenerate(struct distfuncDKGK *FiGK, int num_spec, double *nioverne, double *mioverme, double *TioverTe, double dvpar, double dvperp) {
	int n, i, j, coldelectrons;
	double mu, Uminmu, ff;
	double u, condition, chodura, normalization;// flow;

	for (n=0; n<num_spec; n++) {
		if (TioverTe[n] > 10.0) coldelectrons = 1;
		else coldelectrons = 0;
		printf("are electrons cold (1=yes, 0=no): %d\n", coldelectrons);
		u = 0.0;
		condition = TioverTe[n]; 
		if (coldelectrons == 0) {
			if (TioverTe[n]<=1.0) {
				chodura = 9999999.0;
				u = 0.0;
				while ( (chodura > condition) || (chodura < condition - 0.05) ) {
					if (chodura > condition)
						u += 0.0001;
					else 
						u -= 0.0001;
					normalization =  sqrt(2.0)*(1 + erf(u/sqrt(2.0)))*(1.0+u*u) + sqrt(2.0/M_PI)*u*exp(-0.5*u*u);
					//printf("normalization = %f\n", normalization);
					chodura = sqrt(2.0)*(1+erf(u/sqrt(2.0)))/(normalization);
					//flow = ( 0.5*(1+0.5*u*u)*exp(-0.5*u*u) + (1.5 + 0.5*u*u)*(sqrt(M_PI/2.0)*u/2.0)*(1 + erf(u/sqrt(2.0))))*4.0/sqrt(2.0*M_PI)/normalization;
					//printf("%f, %f, %f, %f %f\n", normalization, u, flow, chodura, condition);
				}
			}
			else {
				chodura = 0.0;
				u = 0.1;
				while  ( (chodura > condition) || (chodura < condition - 0.05) ) {
					if (chodura > condition)
						u -= 0.001;
					else
						u += 0.001;
					normalization = (1.0/sqrt(2.0))*2.0*(sqrt(M_PI*u) - M_PI*exp(1/u)*(1-erf(1.0/sqrt(u))))/(2.0*u*sqrt(u));
					//printf("normalization = %f\n", normalization);
					chodura = (1.0/sqrt(2.0))*0.5*2.0*M_PI*exp(1.0/u)*(1.0-erf(1.0/sqrt(u)))/(2.0*sqrt(u)*normalization);
					//printf("%f, %f, %f, %f %f\n", normalization, u, flow, chodura, condition);
				}
			}
		}

		if (coldelectrons == 0) {
			if ( TioverTe[n] <= 1.0 ) 
				//normalization =  (1 + erf(u/sqrt(2.0)))*(1.0+u*u) + sqrt(2.0/M_PI)*u*exp(-0.5*u*u);
				normalization =  sqrt(2.0)*(1 + erf(u/sqrt(2.0)))*(1.0+u*u) + sqrt(2.0/M_PI)*u*exp(-0.5*u*u);
				//normalization =  (1 + erf(u))*(1.0+2.0*u*u) + (2.0/sqrt(M_PI))*u*exp(-u*u);
			else
				//normalization =  2.0*(sqrt(M_PI*u) - M_PI*exp(1/u)*(1-erf(1/sqrt(u))))/(2.0*u*sqrt(u));
				normalization = (1.0/sqrt(2.0))*2.0*(sqrt(M_PI*u) - M_PI*exp(1/u)*(1-erf(1.0/sqrt(u))))/(2.0*u*sqrt(u));
		}

		if (coldelectrons == 0) {
			if (TioverTe[n]<=1.0) {
				//printf("TioverTe[%d] = %f <= 1\n", n, TioverTe[n]);
				printf("normalization = %f\n", normalization);
				for (i=0; i<FiGK[n].len_perp; i++) {
					mu = 0.5*i*dvperp*i*dvperp;
					mu = i*dvperp;
					FiGK[n].perp[i] = mu;
					for (j=0; j<FiGK[n].len_par; j++) { 
						Uminmu = 0.5*j*dvpar*j*dvpar;
						if (i==0)
							FiGK[n].par[j] = Uminmu;
						ff = (2.0/(M_PI*sqrt(M_PI)))*(1.0/normalization)*(Uminmu)*exp(- Uminmu - mu + u*sqrt(2.0*Uminmu) - 0.5*u*u);
						FiGK[n].F[i][j] = ff;
					}
				}
			}
			else {
				//printf("TioverTe[%d] = %f > 1\n", n, TioverTe[n]);
				for (i=0; i<FiGK[n].len_perp; i++) {
					mu = 0.5*i*dvperp*i*dvperp;
					mu = i*dvperp;
					FiGK[n].perp[i] = mu;
					for (j=0; j<FiGK[n].len_par; j++) { 
						Uminmu = 0.5*j*dvpar*j*dvpar;
						if (i==0)
							FiGK[n].par[j] = Uminmu;

						ff = (1.0/M_PI)*(1.0/normalization)*((Uminmu)/(1+u*(Uminmu)))*exp(- Uminmu - mu );
						FiGK[n].F[i][j] = ff;
					}
				}
			}
		}
		else { // coldelectrons == 1:
			printf("TioverTe[%d] = infinite: ELECTRONS ARE COLD\n", n);
			normalization = 1.0;
			for (i=0; i<FiGK[n].len_perp; i++) {
				mu = 0.5*i*dvperp*i*dvperp;
				mu = i*dvperp;
				FiGK[n].perp[i] = mu;
				for (j=0; j<FiGK[n].len_par; j++) { 
					Uminmu = 0.5*j*dvpar*j*dvpar;
					if (i==0)
						FiGK[n].par[j] = Uminmu;
					ff = (1.0/(sqrt(2.0)*(M_PI)*sqrt(M_PI)))*exp(- Uminmu - mu );
					FiGK[n].F[i][j] = ff;
				}
			}
		}
	}
	return;
}

void Fegen(double **ffarr, double *vpar_e, double *mue, int sizevpar_e, int sizemu_e, double dvpar, double dvperp) {
	double mu, vpar, ff, Uminmu;
	int i,j;
	double halfdvperpsq = dvperp;
	for (i=0; i<sizemu_e; i++) { 
		//mu = 0.5*i*dvperp*i*dvperp;
		mu = i*halfdvperpsq;
		mue[i] = mu;
		for (j=0; j<sizevpar_e; j++) { 
			vpar = j*dvpar;
			Uminmu = 0.5*vpar*vpar;
			if (i==0) vpar_e[j] = vpar;
			ff = (1.0/pow(2.0*M_PI, 1.5))*exp(- Uminmu - mu );
			ffarr[i][j] = ff;
		}
	}
	return;
}

//double xifp(double phi) {
//	double y;
//	y = 1.0/sqrt( -3.0 -2.0*phi + 4.0*exp(phi) - exp(2.0*phi) );
//	return y;
//}

void phi_fluidMPS(double *phigrid, double *xgrid, int size_x, double alpha) {
	double *x_for_interp, *phi_for_interp, *psi_for_interp, phi0=log(alpha), dxdphi=0.0, dxdphiold;
	int i, size_for_interp = 500;
	double dphi = -phi0/size_for_interp;
	x_for_interp = malloc(size_for_interp*sizeof(double));
	phi_for_interp = malloc(size_for_interp*sizeof(double));
	psi_for_interp = malloc(size_for_interp*sizeof(double));
	for (i=0; i<size_for_interp; i++) {
		dxdphiold = dxdphi;
		phi_for_interp[i] = log(alpha) + i*dphi;
		psi_for_interp[i] = phi_for_interp[i] + 0.5*alpha*alpha*(exp(-2.0*phi_for_interp[i])-1.0);
		dxdphi = 1.0/sqrt( -3.0 -2.0*psi_for_interp[i] + 4.0*exp(psi_for_interp[i]) - exp(2.0*psi_for_interp[i]) ) ; //xifp(psi_for_interp[i]);
		if (i!=0) x_for_interp[i] = x_for_interp[i-1] + (dxdphi + dxdphiold)*0.5*dphi;
		else x_for_interp[0] = 0.0;
		//printf(" x,phi = %f, %f\n", x_for_interp[i], phi_for_interp[i]);
	}
	for (i=0; i<size_x; i++) {
		phigrid[i] = lin_interp(x_for_interp, phi_for_interp, xgrid[i], size_for_interp, 1228);
	}
	return;
}


void make_phigrid(double *x_grid, double *phi_grid, int size_phigrid, double grid_parameter, double deltax, int initial, double phi_jump, double len_scale, double alpha) {
	// improved = 0 is only working case
	int i, improved = 0;//, size_phigrid_improved;
	double xi, *ff, *new_phi, *new_x, power = 0.5, deltaxmin = 0.2, deltaphi=0.1;
	new_phi = malloc(size_phigrid*sizeof(double));
	new_x = malloc(size_phigrid*sizeof(double));
  	ff = malloc(size_phigrid*sizeof(double));
	printf("in make_phigrid: size_phigrid = %d\n", size_phigrid);
	if (initial != 0) {
		if (improved == 0) {
			//printf("grid_parameter = %f\n", grid_parameter);
			for (i=0; i<size_phigrid; i++) {
			    //ff[i] = pow(sqrt(grid_parameter)+sqrt(x_grid[i]), 2.0) - grid_parameter;
			    ff[i] = pow(pow(grid_parameter, power)+pow(x_grid[i], power), 1.0/power) - grid_parameter;
				//printf("ff[%d] = %f\n", i, ff[i]);
			}
			gsl_interp_accel *acc
			= gsl_interp_accel_alloc ();
			gsl_spline *spline
			= gsl_spline_alloc (gsl_interp_cspline, size_phigrid);

			gsl_spline_init (spline, x_grid, phi_grid, size_phigrid);

			for (i=0; i < size_phigrid; i++) {
				xi = i*deltax;
				new_phi[i] = phi_grid[i];
				new_x[i] = pow( pow(grid_parameter+xi, power) - pow(grid_parameter, power), 1.0/power);
				//new_phi[i] = lin_interp(x_grid, phi_grid, new_x[i], size_phigrid, 1);
				//new_phi[i] = gsl_spline_eval( spline, new_x[i], acc) ;
				//printf("%f ", new_x[i]);
				//g = sqrt(new_x[i]);
			}
			gsl_spline_free (spline);
			gsl_interp_accel_free (acc);
			for (i=0; i < size_phigrid; i++) {
				phi_grid[i] = new_phi[i];
				x_grid[i] = new_x[i];
			}
		}
		else { // improved =1 STILL NOT WORKING
			gsl_interp_accel *acc
			= gsl_interp_accel_alloc ();
			gsl_spline *spline
			= gsl_spline_alloc (gsl_interp_cspline, size_phigrid);

			gsl_spline_init (spline, phi_grid, x_grid, size_phigrid);

			gsl_interp_accel *acc2
			= gsl_interp_accel_alloc ();
			gsl_spline *spline2
			= gsl_spline_alloc (gsl_interp_cspline, size_phigrid);

			gsl_spline_init (spline2, x_grid, phi_grid, size_phigrid);

			new_x[0] = 0.0;
			new_phi[0] = phi_grid[0];
			i=0;
			while (i < size_phigrid) {
			//for (i=1; i < size_phigrid; i++) 
				if (phi_grid[0] + i*deltaphi < phi_grid[size_phigrid-1])
					xi = gsl_spline_eval (spline, phi_grid[0] + deltaphi*i, acc);
				else xi = x_grid[i-1] + deltaxmin + TINY;
				printf("xi = %f\n", xi);
				if (xi < x_grid[i-1] + deltaxmin) { 
					new_x[i] = xi;
					new_phi[i] = phi_grid[0] + deltaphi*i;
					printf("new_phi = %f\n", new_phi[i]);
				}
				else if (new_x[i-1] + deltaxmin > x_grid[size_phigrid-1]) {
					//i = size_phigrid-1;
					new_x[i] = new_x[i-1] + deltaxmin; 
					new_phi[i] = new_phi[i-1]  + TINY;
				}
				else {
					new_x[i] = new_x[i-1] + deltaxmin;
					new_phi[i] = gsl_spline_eval (spline2, new_x[i], acc);
				}
				printf("i=%d\tnew_x = %f\tnew_phi = %f\n", i, new_x[i], new_phi[i]);
				i++;
			}
			gsl_spline_free (spline);
			gsl_interp_accel_free (acc);
			gsl_spline_free (spline2);
			gsl_interp_accel_free (acc2);
			for (i=0; i < size_phigrid; i++) {
				phi_grid[i] = new_phi[i];
				x_grid[i] = new_x[i];
			}
		}
	}
	else {
		for (i=0; i < size_phigrid; i++) {
			xi = i*deltax;
			ff[i] = xi;
			//new_x[i] = pow( pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), 2.0);
			new_x[i] = pow( pow(grid_parameter+xi, power) - pow(grid_parameter, power), 2.0);
			//g = sqrt(new_x[i]);
			new_phi[i] = phi_jump*pow(len_scale, 2.0)/pow(new_x[i] + len_scale, 2.0);
			new_phi[i] = phi_jump*exp(-new_x[i]/len_scale);
			//if (new_x[i] < 10.0/sqrt(2)) 
			//	new_phi[i] = 3.0*pow((new_x[i]*sqrt(2.0)/10.0 - 1.0), 5.0);
			//else 
			//	new_phi[i] = 0.0;
		}
		//phi_fluidMPS(new_phi, new_x, size_phigrid, alpha);
		for (i=0; i < size_phigrid; i++) {
			phi_grid[i] = new_phi[i];
			x_grid[i] = new_x[i];
		}
	}
	free(new_phi);
	free(new_x);
	free(ff);
	return;
}

// this function is no longer used, but keep it just in case
void rescale_array(double *array, int size_array, double jump) {
	int i;
	double old_jump = array[0];
	for (i=0; i < size_array; i++) {
		array[i] *= (jump/old_jump);
	}
	return;
}

/* the set of functions below are necessary to model electron reflection from the infinitely thin Debye sheath obtain parallel velocity cutoff as a function of magnetic moment with assumption rho_e >> lambda_D
*/
double FF(double beta) {
	double integrand, integral=0.0, beta_s = 0.005, betaval;
	int ind, num_beta = (int) (beta/beta_s);
	for (ind = 0; ind < num_beta+1; ind++) {
		betaval = beta*ind/num_beta;
		integrand = 0.5*(1.0-cos(2.0*betaval)) / ( M_PI - betaval + 0.5*sin(2.0*betaval) );
		integral += integrand*beta/num_beta;
	}
        return integral;
}
double FFp(double beta) {
	double integrand = 0.5*(1.0-cos(2.0*beta)) / ( M_PI - beta + 0.5*sin(2.0*beta) );
        return integrand;
}
double vparcut(double beta, double vcut) {
	double vparcutval = sqrt((1.0-exp(-2.0*FF(beta))))*vcut/sin(beta);
	return vparcutval;
}
double mucut(double beta, double vcut) {
	double mucutval = 0.5*exp(-2.0*FF(beta))*pow(vcut/sin(beta), 2.0);
	return mucutval;
}

double vparcut_mu(double mu, double vcut) {
	double beta, musearch=-10000.0, vparcutn;
	double betalow = 0.0, betaup = M_PI;
	do {
		beta = (betaup + betalow)/2.0;
		musearch = 0.5*exp(-2.0*FF(beta))*pow(vcut/sin(beta), 2.0);
		if (musearch < mu) {
			betaup = beta;
		}
		else {
			betalow = beta;
		}
		//printf("%f %f\n", musearch, mu);
	} while (fabs(musearch-mu) > TINY) ;
	vparcutn = sqrt((1.0-exp(-2.0*FF(beta))))*vcut/sin(beta);
	return vparcutn;
}

// The main function of MAGSHEATH
int main(void) {
// computation time
	int MAX_IT, ZOOM_DS, ZOOM_MP;
	double INITIAL_GRID_PARAMETER, SYS_SIZ, MAXMU, DMU, MAXVPAR, MAXVPAR_I, DVPAR, DVPAR_I, SMALLGAMMA, tol_MP[2], tol_DS[2], tol_current, WEIGHT_MP, WEIGHT_DS, WEIGHT_j, MARGIN_MP, MARGIN_DS, GRIDSIZE_MP, GRIDSIZE_DS; // DXMIN

	clock_t begin_it = clock(); // Finds the start time of the computation
	double tot_time;
// strings for reading files
	char line_million[1000000], line_hundred[100], dirname[200], dirname_it[200], dirnameit[150];
// random
	int type_distfunc_entrance;
	int i, j, n=0, nrows_distfile = 0, ncols=0, ndirname, ind; //s
	double phi0_init_MP = -0.0, ioncharge= 1.0, EW;
// pointer to double to which strings above is converted
	double *storevals;
// velocity grid parameters universal to all species
// spatial grid parameters in magnetic presheath / Chodura sheath
	double grid_parameter, deltax, system_size, DS_size;
// spatial grid parameters in Debye sheath
	double deltaxDS; 
// input parameters set in inputfile.txt
	int num_spec, fix_current=0;
	double alpha, *nioverne, *TioverTe, *mioverme, gamma_ref, target_current;
// other parameters derived from input ones
	double alpha_deg, factor_small_grid_parameter=1.0;
// quantities related to overall current, potential drop or any bump in potential in Debye sheath (not implemented yet)
	double v_cut, current, v_cutDS = 0.5;
// quantities related to electrostatic potential in the magnetic preasheath / Chodura sheath
	double *x_grid, *phi_grid;
// quantities related to electrostatic potential the Debye sheath
	double *phi_DSgrid, *x_DSgrid;
	int size_phigrid, size_phiDSgrid;
// electron quantities
	double *ne_grid, *ne_DSgrid, ne_inf = 0.0;
	double *vy_e_wall, *mu_e_op, *chiM_e, *twopidmudvy_e, *vpar_e_cut_lookup;
	double **dist_e_DK, **dist_e_GK, *vpar_e, *mu_e, *U_e_DS, *vpar_e_DS, *vpar_e_cut;
	double flux_e=0.0, flux_eDS = 0.0, Q_e=0.0, Q_eDS=0.0;
	double Uminmu_MPE, garbage = 0.0;
	int size_mu_e, size_vpar_e, size_neDSgrid = 0, size_op_e=0;
// quantities related to individual ion species: dim1 is for species
	double **ni_grid, **ni_DSgrid;
	double **vy_i_wall, **mu_i_op, **chiM_i, **twopidmudvy_i;
	double *flux_i, *lenfactor, Bohm=0.0, lenMP, *Q_i;
	int *size_ngrid, *size_mu_i, *size_U_i, *sizevxopen, *size_op_i;
	struct distfuncDKGK *FiGK;
	double ***dist_i_GK, **mu_i, **U_i;
// quantities related to overall ion properties
	double *sumni_grid, *sumni_DSgrid, sumflux_i=0.0, sumni_norm=0.0, sumQ_i=0.0;
	int size_sumnigrid;
// quantities related to the iteration
	int zoomfactor;
	int N=0, N_DS=0, convergence_MP = 0, convergence_DS = 0, convergence_j = 0;
	double error_MP[2], error_DS[2];
	double weight_j, weight_MP, weight_DS;
	double Bohmshouldbe=0.0;

	//printf("%f\n", acos(-1.0));
	mkdir("OUTPUT", S_IRWXU);
	FILE *fout = fopen("OUTPUT/output_MAGSHEATH.txt", "w");
	FILE *fp; 
	char fpstr[150];
	char fBohmstr[150];

	if (fout == NULL) {
		printf("Problem opening output file for writing\n");
		exit(1);
	}

	/* READ INPUT FILE 
	Contains the values of: 
	     alpha
	     number of ion species
	     ni/ne for all species (line to array)
	     Ti/Te for all species (line to array)
	     mi/me for all ion species (line to array)
	     fix_current = 1 or 0 (fixes potential),
	     fixed current/potential 
	*/

	// Note: multiple species not yet implemented
	FILE *input, *fBohm, *numinput;

	if ((input = fopen("inputfile.txt", "r")) == NULL) { 
		printf("Cannot open inputfile.txt. Try input_physparams.txt");
		fprintf(fout, "Cannot open inputfile.txt. Try input_physparams.txt");
		if ((input = fopen("input_physparams.txt", "r")) == NULL) { 
			printf("Cannot open %s\n", "input_physparams.txt");
			fprintf(fout, "Cannot open %s\n", "input_physparams.txt");
			exit(EXIT_FAILURE); 
		}
	}
	i=0; // for counting rows (lines) of file
	j=0;
	ncols = 0; // for counting numbers in each row (line) of file
	ndirname=0;
	char *outputstr = "OUTPUT/";
	for (ndirname=0; ndirname<strlen(outputstr);ndirname++) 
		dirname[ndirname] = outputstr[ndirname];
	while (fgets(line_hundred, 100, input) != NULL) {	
		if ( ( (line_hundred[0] != '#') && (line_hundred[0] != ' ') ) && (line_hundred[0] != '\n') ) {
		storevals = linetodata(line_hundred, strlen(line_hundred), &ncols);
		//printf("ndirname = %d\n", ndirname);
		if ( (i!= 7) && (i!=0) ) {
			for (j=0; j < lenstrqty-1; j++) {
				dirname[ndirname+j] = strqty[i-fix_current-1][j];
			}
			ndirname += (lenstrqty -1);
			//dirname[0] = 'a'; dirname[1] = 'l'; dirname[2] 
			for (j=0;j<strlen(line_hundred)-1;j++) {
				//printf("%c\n", line_hundred[j]);
				dirname[j+ndirname] = line_hundred[j];
				printf("%s\n", dirname);
			}
		}
		else if (i==0) {
			for (j=0;j<strlen(line_hundred)-1;j++) {
				//printf("%c\n", line_hundred[j]);
				dirname[j+ndirname] = line_hundred[j];
				printf("%s\n", dirname);
			}
		}
		if (i==0){
			type_distfunc_entrance = strncmp(line_hundred, "ADHOC", 5); // type_distfunc_entrance = 0 if the functions are ADHOC
		}
		if (i==1)
			// input in degrees only for convenience
			// should be below 5 degrees (0.1 rad) for asymptotic theory to be valid
			alpha_deg = *storevals; 
		else if (i==2)
			gamma_ref = *storevals; 
		else if (i==3)
			num_spec = (int) (*storevals); 
		else if (i==4)
			nioverne = storevals; 
		else if (i==5)
			TioverTe = storevals; //linetodata(line_hundred, strlen(line_hundred), &ncols); 
		else if (i==6)
			mioverme = storevals; 
		else if (i==7)
			fix_current = (int) (*storevals); 
		else if (i==8) {
			if (fix_current != 0) {
				target_current = current = *storevals; // in units of thermal electron velocity 
				v_cut = 3.0; // set to a reasonable value
				// we will set v_cut to some value later
			}
			else {
				v_cut = sqrt(2.0*(*storevals)); 
				target_current = current = 0.0; // gets calculated afterwords
			}
		}
		if (i!=7) {
			ndirname += strlen(line_hundred)-1;
			dirname[ndirname] = '/';
			ndirname += 1;
		}
		i += 1; // count the rows in the file
		//printf("dirname = %s\n", dirname);
		printf("%s\n", dirname);
		mkdir(dirname, S_IRWXU);
		}
	}

	snprintf(fBohmstr, 150, "%s/Bohmvalues.txt", dirname);
	fBohm = fopen(fBohmstr, "w");
	if (fBohm == NULL) { 
		printf("Cannot open %s\n", fBohmstr);
		fprintf(fout, "Cannot open %s\n", fBohmstr);
		exit(EXIT_FAILURE); 
	}
	snprintf(fpstr, 150, "%s/inputfile.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file\n");
	}
	if ((input = fopen("inputfile.txt", "r")) == NULL) { 
		printf("Cannot open inputfile.txt. Try input_physparams.txt");
		fprintf(fout, "Cannot open inputfile.txt. Try input_physparams.txt");
		if ((input = fopen("input_physparams.txt", "r")) == NULL) { 
			printf("Cannot open %s\n", "input_physparams.txt");
			fprintf(fout, "Cannot open %s\n", "input_physparams.txt");
			exit(EXIT_FAILURE); 
		}
	}
	while (fgets(line_hundred, 100, input) != NULL) {
		if ( (line_hundred[0] != '#') && (isalpha(line_hundred[0]) == 0) )
			fprintf(fp, "%s", line_hundred);
	}
	
	fclose(input);
	fclose(fp);


	for (i=0; i< ndirname; i++) {
		dirname_it[i] = dirname[i];
	}
	char *strit = "iteration";
	int nstrit;
	nstrit = strlen(strit);
	for (i=0; i< nstrit; i++) {
		dirname_it[ndirname+i] = strit[i];
	}
	printf("dirname for iteration = %s\n", dirname_it);
	mkdir(dirname, S_IRWXU);
	for (n=0; n<num_spec; n++) sumni_norm += nioverne[n];
	for (n=0; n<num_spec; n++) {
		nioverne[n] /= sumni_norm;
		printf("nioverne[%d] = %f\n", n, nioverne[n]);
	}
	//dirname[ndirname+1] = '\0';
	//ndirname+=1;
//<<<<<<< HEAD
//=======
//	/*fclose(input);*/
//>>>>>>> refs/remotes/origin/main
	if ( (fabs(gamma_ref) < TINY) && (alpha_deg < 1.0) ) deltax = 0.6;
	//if (TioverTe[0] < 0.4) weight_MP = WEIGHT_MP/6.0;
	printf("directory where output will be stored is %s\n", dirname);
	alpha = alpha_deg*M_PI/180; // alpha used in radians from now on
	// initial iteration assumes flat potential profile in magnetic presheath
	// therefore, the parameter v_cutDS is equal to v_cut

/* 
	FINISHED READING PHYSICAL INPUT FILE 
	AND GENERATING THE NEW (IF NOT ALREADY EXISTING) 
	DIRECTORIES WHERE FILES WILL BE STORED
	INCLUDING ITERATION SUBDIRECTORIES
*/

/* 
	READ NUMERICAL INPUT FILE
	AND GENERATE ARRAYS OF THE CORRECT SIZES
*/

	//Now open numerical input file
	if ((numinput = fopen("input_numparams.txt", "r")) == NULL) { 
		printf("Cannot open %s\n", "input_numparams.txt");
		fprintf(fout, "Cannot open %s\n", "input_numparams.txt");
		exit(EXIT_FAILURE); 
	}

	ncols = 0; // for counting numbers in each row (line) of file
	i=0;
	while (fgets(line_hundred, 100, numinput) != NULL) {
		if ( ( (line_hundred[0] != '#') && (line_hundred[0] != ' ') ) && (line_hundred[0] != '\n') ) {
			storevals = linetodata(line_hundred, strlen(line_hundred), &ncols);
			if (i==0) MAX_IT = *storevals;
			else if (i==1) {
				INITIAL_GRID_PARAMETER = storevals[0];
				printf("INITIAL GP = %f\n", INITIAL_GRID_PARAMETER);
				SYS_SIZ = storevals[1];
				GRIDSIZE_MP =  storevals[2]; 
				GRIDSIZE_DS =  storevals[3]; 
				grid_parameter=INITIAL_GRID_PARAMETER;
				deltax = GRIDSIZE_MP;
			}
			else if (i==2) {
				MAXMU = storevals[0]; 
				MAXVPAR = storevals[1]; 
				MAXVPAR_I = storevals[2]; 
				DMU =  storevals[3]; 
				DVPAR =  storevals[4]; 
				DVPAR_I =  storevals[5]; 
			}
			if (i==3) SMALLGAMMA = *storevals; //linetodata(line_hundred, strlen(line_hundred), &ncols); 
			if (i==4) {
				tol_MP[0] = storevals[0]; tol_MP[1] = storevals[1]; tol_DS[0] = storevals[2]; tol_DS[1] = storevals[3]; tol_current = storevals[4];
			}
			if (i==5) {
				WEIGHT_MP =  storevals[0]; WEIGHT_DS =  storevals[1]; WEIGHT_j =  storevals[2]; 
				weight_j=WEIGHT_j; weight_MP=WEIGHT_MP; weight_DS=WEIGHT_DS;
			}
			if (i==6) {
				MARGIN_MP =  storevals[0]; MARGIN_DS =  storevals[1]; 
			}
			if (i==7) {
				ZOOM_MP =  storevals[0]; ZOOM_DS =  storevals[1]; 
			}
			i+=1;
		}
	}
	fclose(numinput);
	printf("SMALLGAMMA = %f\n", SMALLGAMMA);

	lenfactor = malloc(num_spec*sizeof(double));
	lenMP = 0.0;
	for (n=0; n<num_spec; n++) {
		lenMP += (nioverne[n]*mioverme[n]) ;
	}
	lenMP = sqrt(lenMP);
	for (n=0; n<num_spec; n++)
		lenfactor[n] = lenMP/sqrt(mioverme[n]*TioverTe[n]); 
		
	n=0; 
	lenMP = 0.0;
	printf("lenMP = %f\tlenfactor[0] = %f\n", lenMP, lenfactor[0]);
	while ( (n < num_spec) && (lenMP < 1.0/lenfactor[n]) ) {
		lenMP = 1.0/lenfactor[n];
		n++;
	}
	printf("lenMP = %f\n", lenMP);
	if (TioverTe[0] > 10.0) lenMP = 1.0;

	v_cutDS = v_cut;
	if (gamma_ref > 1.0)
		deltaxDS = GRIDSIZE_DS/gamma_ref;
	else
		deltaxDS = GRIDSIZE_DS;
	printf("deltaxDS = %f\n", deltaxDS);
	if (lenMP > 1.0) {
		system_size = SYS_SIZ*lenMP;
		//deltax *= lenMP;
		//zoomfactor = (int) (deltax / DXMIN);
		zoomfactor = ZOOM_MP; // *( (int) lenMP );
	}
	else {
		system_size = SYS_SIZ;
		//zoomfactor = (int) ( deltax / (DXMIN*lenMP) ); 
		zoomfactor = (int) ( ZOOM_MP*(1.0/lenMP) );
	}
	printf("lenMP = %f\n", lenMP);
	printf("zoomfactor = %d\n", zoomfactor);
	
/* 
	FINISHED READING NUMERICAL INPUT FILE 
*/


/* 
	PRINT SOME INPUTS, ALSO TO A FILE
*/

i=0;
	printf("INPUT PARAMETERS:\n\tmagnetic field angle = α = %f deg (%f rad)\n\telectron gyroradius over Debye length (reference value at MP entrance) = ρ_e/λ_D = %f\n\tnumber of species = %d\n", alpha_deg, alpha, gamma_ref, num_spec);
	fprintf(fout, "INPUT PARAMETERS:\n\tmagnetic field angle = α = %f deg (%f rad)\n\telectron gyroradius over Debye length  (reference value at MP entrance) = ρ_e/λ_D = %f\n\tnumber of species = %d\n", alpha_deg, alpha, gamma_ref, num_spec);
	for (i=0; i<num_spec; i++) {
		printf("\tion temperature = τ = ZT_i/T_e for species %d = %f\n\tmass ratio = m_i/m_e for species ? = %f\n", i+1, TioverTe[i], mioverme[i]);
		fprintf(fout, "\tion temperature = τ = ZT_i/T_e for species %d = %f\n\tmass ratio = m_i/m_e for species ? = %f\n", i+1, TioverTe[i], mioverme[i]);
	}
	printf("\tfix_current = %d\n", fix_current);
	fprintf(fout, "\tfix_current = %d\n", fix_current);
	if (fix_current != 0) {
		printf("\tcurrent = j/(e*n_MPE*v_t,i) = %f\n", target_current);
		fprintf(fout, "\tcurrent = j/(e*n_MPE*v_t,i) = %f\n", target_current);
	}
	else {
		printf("\twall potential = eφ_W/T_e = %f\n", -0.5*v_cut*v_cut);
		fprintf(fout, "\twall potential = eφ_W/T_e = %f\n", -0.5*v_cut*v_cut);
	}

	printf("magnetic presheath size in simulation (in ion gyroradii) = %f\n", system_size);
	fprintf(fout, "magnetic presheath size in simulation (in ion gyroradii) = %f\n", system_size);
	size_phigrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size), 2.0) - grid_parameter ) / deltax );
	 
	printf("size of coarse potential grid in magnetic presheath = %d\n", size_phigrid);
	fprintf(fout, "size of coarse potential grid in magnetic presheath = %d\n", size_phigrid);
	printf("grid parameter = %f\n", grid_parameter);

/* 
	FINISHED PRINTING INPUTS
*/

/* 
	BEGIN ALLOCATING MEMORY
*/

	// form x, phi, ne and ni grids for magnetic presheath
	x_grid = malloc(size_phigrid*sizeof(double));
	phi_grid = malloc(size_phigrid*sizeof(double));
	ne_grid = malloc(size_phigrid*sizeof(double));
	sumni_grid = malloc(size_phigrid*sizeof(double));
	ni_grid = malloc(num_spec*sizeof(double));
	// WALL DISTRIBUSION FUNCTION BELOW
// gives size of the above arrays for each species
	size_op_i = malloc(num_spec*sizeof(int));
// where values of vy will be stored
	vy_i_wall = malloc(num_spec*sizeof(double));
// where values of mu corresponding to a given value of vy at the wall will be stored
	mu_i_op = malloc(num_spec*sizeof(double)); 
// gives minimum value of U_perp = (1/2)vx^2 + (1/2)vy^2 + phi(x) corresponding to a given value of vy at the wall
	chiM_i = malloc(num_spec*sizeof(double)); 
// gives range of values of (1/2)vx^2, upon multiplying by alpha*vz
	twopidmudvy_i = malloc(num_spec*sizeof(double)); 
	// WALL DISTRIBUTION FUNCTION ABOVE
	flux_i = malloc(num_spec*sizeof(double));
	Q_i = malloc(num_spec*sizeof(double));
	mu_i = malloc(num_spec*sizeof(double));
	U_i = malloc(num_spec*sizeof(double));
	dist_i_GK = malloc(num_spec*sizeof(double));
	FiGK = malloc(num_spec*sizeof(struct distfuncDKGK));
	size_ngrid = calloc(num_spec,sizeof(int));
	size_mu_i = malloc(num_spec*sizeof(int));
	size_U_i = malloc(num_spec*sizeof(int));
	sizevxopen = malloc(num_spec*sizeof(int));
	for (n=0; n<num_spec; n++)
		ni_grid[n] = malloc(size_phigrid*sizeof(double));

	if (type_distfunc_entrance == 0) {
		size_mu_e = (int) (MAXMU/DMU);
		size_vpar_e = (int) (MAXVPAR/DVPAR);
		mu_e = malloc(size_mu_e*sizeof(double));
		vpar_e = malloc(size_vpar_e*sizeof(double));
		dist_e_DK = malloc(size_mu_e*sizeof(double));
		dist_e_GK = malloc(size_mu_e*sizeof(double));
		for (i=0; i < size_mu_e; i++) {
			dist_e_DK[i] = malloc(size_vpar_e*sizeof(double));
			dist_e_GK[i] = malloc(size_vpar_e*sizeof(double));
		}
		Fegen(dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, DVPAR, DMU);
		for (n=0; n< num_spec; n++) {
			size_mu_i[n] = (int) (MAXMU/DMU);
			size_U_i[n] = (int)  ((MAXVPAR_I + 1.0/TioverTe[n])/DVPAR_I);
			mu_i[n] = malloc(size_mu_i[n]*sizeof(double));
			U_i[n] = malloc(size_U_i[n]*sizeof(double));
			dist_i_GK[n] = malloc(size_mu_i[n]*sizeof(double));
			FiGK[n].F = malloc(size_mu_i[n]*sizeof(double));
			FiGK[n].len_perp = size_mu_i[n];
			FiGK[n].len_par = size_U_i[n];
			FiGK[n].par = malloc(size_U_i[n]*sizeof(double));
			FiGK[n].perp = malloc(size_mu_i[n]*sizeof(double));
			for (i=0; i < size_mu_i[n]; i++) {
				dist_i_GK[n][i] = malloc(size_U_i[n]*sizeof(double));
				FiGK[n].F[i] = malloc(size_U_i[n]*sizeof(double));
			}
		}
		Figen(dist_i_GK, U_i, mu_i, num_spec, nioverne, mioverme, TioverTe, size_U_i, size_mu_i, DVPAR_I, DMU);
		//Figenerate(FiGK, num_spec, nioverne, mioverme, TioverTe, DVPAR_I, DMU);
		//for (i=0; i < size_mu_i[0]; i++) {
		//	printf("mu = %f\n", mu_i[0][i]);
		//	for (j=0; j < size_U_i[0]; j++) {
		//		if (i==0)
		//			printf("Uminmu = %f\n", U_i[0][j]);
		//		printf("%f ", dist_i_GK[0][i][j]);
		//	}
		//	printf("\n");
		//}
	}
	else { 
		// import the distribution function from a file into a 2 dimensional array, F(mu,U)
		FILE *file;
		if ((file = fopen("Fi_mpe.txt", "r")) == NULL)
		{	
			printf("cannot open file %s\n", "Fi_mpe.txt");
			fprintf(fout, "cannot open file %s\n", "Fi_mpe.txt");
			exit(-1); 
		}
		/* Count the number of rows in the distfuncin.txt file */
		while(fgets(line_million, 1000000, file) != NULL) 
			nrows_distfile += 1; 
		// Allocate the right amount of memory to the distribution function pointer
		dist_i_GK[0] = (double**) calloc(nrows_distfile,sizeof(double*));
		FiGK[0].F = (double**) calloc(nrows_distfile,sizeof(double*));
		// The number of rows is also the the size of the array in U_i 
		nrows_distfile = 0; // Set number of rows counter to zero again
		rewind(file); // rewind file
		// Read each line of file, extract the data (ie the numbers in the line) and count the columns
		// Once columns are counted, allocate memory to dist_i_GK and assign each number to a different column of dist_i_GK (for a fixed row)
		while (fgets(line_million, 1000000, file) != NULL) {
			storevals = linetodata(line_million, strlen(line_million), &ncols);
			dist_i_GK[0][nrows_distfile] = storevals;
			FiGK[0].F[nrows_distfile] = storevals;
			nrows_distfile +=1; 
		}
		fclose(file);
		//printf("~~~~~The second element of dist_i_GK is %f~~~~~\n",dist_i_GK[0][1]);
		// Extract two 1D arrays representing the grid points on which dist_i_GK is defined
		// first open file and check for error
		if ((file = fopen("Fi_mpe_args.txt", "r")) == NULL) {	
			printf("cannot open file %s\n", "Fi_mpe_args.txt");
			fprintf(fout, "cannot open file %s\n", "Fi_mpe_args.txt");
			exit(-1); 
		}
		// Read each line of file, extract the data and assign the values of the first line to mu_i, second to U_i
		i=0;
		while (fgets(line_million, 1000000, file) != NULL) {	
			storevals = linetodata(line_million, strlen(line_million), &ncols);
			if (i == 0) {
				size_mu_i[0] = ncols;
				mu_i[0] = storevals;
				printf("size_mu_i[0] = %d\n", size_mu_i[0]);
			}
			else if (i==1) {
				size_U_i[0] = ncols;
				//printf("size_U_i = %d\n", size_U_i);
				U_i[0] = storevals;
				//for (i=0; i< size_U_i; i++) printf("U_i[%d] = %f\n", i, U_i[i]);
			}
			i += 1; 
		}
		fclose(file);

		/*
		EXTRACT ELECTRON DISTRIBUTION FUNCTION
		Import the distribution function into a 2 dimensional array, F(mu,vpar). 
		*/
		file = fopen("Fe_mpe.txt", "r");
		if (file == NULL)
		{	
			printf("Cannot open file Fe_mpe.txt");
			fprintf(fout, "Cannot open file Fe_mpe.txt");
			exit(-1); 
		}
		/* Count the number of rows in the distfuncin.txt file */
		while(fgets(line_million, 1000000, file) != NULL) {	
			nrows_distfile += 1; 
		}
		/* Allocate the right amount of memory to the distribution function pointer */
		dist_e_DK = (double**) calloc(nrows_distfile,sizeof(double*));
		dist_e_GK = (double**) calloc(nrows_distfile,sizeof(double*));
		nrows_distfile = 0; // Set number of rows counter to zero again
		rewind(file); // rewind to first line of file
		while (fgets(line_million, 1000000, file) != NULL) {
			// assign each number to a different column of dist_e_DK (for a fixed row)
			dist_e_DK[nrows_distfile] = linetodata(line_million, strlen(line_million), &ncols);
			// allocate memory to dist_e_GK which might be used to solve Debye sheath
			// dist_e_GK depends on magnetic presheath and Debye sheath solutions
			// will assign different values to dist_e_GK at each iteration
			dist_e_GK[nrows_distfile] = (double*) calloc(ncols,sizeof(double*));
			nrows_distfile +=1; 
			// jump to a new line and repeat the process
		}
		fclose(file);
		//printf("~~~~~The second element of dist_e_DK is %f~~~~~\n",dist_e_DK[0][1]);
		// Extract the two 1D arrays representing the grid points on which dist_i_GK is defined
		// first open file and check for error
		if ((file = fopen("Fe_mpe_args.txt", "r")) == NULL) {	
			printf("cannot open file Fe_mpe_args.txt");
			fprintf(fout, "Cannot open file Fe_mpe_args.txt");
			exit(-1); 
		}
		// Read each line of the file, extract the data (ie numbers in a line) 
		// Assign values of first line to mu_e, second to vpar_e
		i=0;
		while (fgets(line_million, 1000000, file) != NULL) {	
			storevals = linetodata(line_million, strlen(line_million), &ncols);
			if (i == 0) {
				size_mu_e = ncols;
				mu_e = storevals;
			}
			else {
				size_vpar_e = ncols;
				vpar_e = storevals;
			}
			i += 1; 
		}
		printf("size_mu_e = %d\nsize_vpar_e=%d\n", size_mu_e, size_vpar_e);
		fclose(file);
	}
	// WALL ION DISTRIBUTION FUNCTION BELOW
	for (n=0; n<num_spec; n++) {
// gives size of the above arrays for each species
		sizevxopen[n] = (int) 50*sqrt((1.0+1.0/TioverTe[n])*(1.0+TioverTe[n]));
// where values of vy will be stored
		vy_i_wall[n]  = malloc((zoomfactor*size_phigrid+1)*sizeof(double));
// where values of mu corresponding to a given value of vy at the wall will be stored
		mu_i_op[n]  = malloc((zoomfactor*size_phigrid+1)*sizeof(double));
// gives minimum value of U_perp = (1/2)vx^2 + (1/2)vy^2 + phi(x) corresponding to a given value of vy at the wall
		chiM_i[n]  = malloc((zoomfactor*size_phigrid+1)*sizeof(double));
// gives range of values of (1/2)vx^2, upon multiplying by alpha*vz
		twopidmudvy_i[n]  = malloc((zoomfactor*size_phigrid+1)*sizeof(double));
	}
	// WALL DISTRIBUTION FUNCTION ABOVE
	U_e_DS  = malloc(size_vpar_e*sizeof(double));
	vpar_e_DS  = malloc(size_vpar_e*sizeof(double));
	vpar_e_cut  = malloc(size_mu_e*sizeof(double));

	/*
	CARRY OUT A SINGLE ION DENSITY CALCULATION
	*/
	if (TioverTe[0] > 10.0) {	
		printf("WARNING: ion temperature too high. Proceed by using hot ion limit i.e. flat potential profile (otherwise code will be too slow for little gain)\n");
		make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, 0, phi0_init_MP, 1.0, alpha);
		//size_ngrid = (int) ( ( pow(sqrt(grid_parameter) + sqrt(system_size - 5.0), 2.0) - grid_parameter ) / deltax );
		for (n=0; n<num_spec; n++) 
			densfinorb(TioverTe[n], lenfactor[n], alpha, size_phigrid, size_ngrid+n, ne_grid, x_grid, phi_grid, ioncharge, dist_i_GK[n], mu_i[n], U_i[n], size_mu_i[n], size_U_i[n], grid_parameter, flux_i+n, Q_i+n, zoomfactor, MARGIN_MP, -999.9, vy_i_wall[n], mu_i_op[n], chiM_i[n], twopidmudvy_i[n], size_op_i+n);
		FILE *fout; 
		if ((fout = fopen("TESTS/densfinorbflat.txt", "w")) == NULL) {	
			printf("Cannot open TESTS/densfinorbflat_out.txt");
			exit(EXIT_FAILURE);
		}
		for (i=0; i < size_ngrid[0]; i++)
			fprintf(fout, "%f %f\n", x_grid[i], ne_grid[i]);
		exit(0);
	}

	/* 
		THIS PART IS JUST A TEST FOR THE FINITE ELCTRON DENSITY CALCULATION
		IT SHOULD ALWAYS BE IGNORED
		TEST_EL = 0 SHOULD BE USED
	*/
	if (TEST_EL==1) { //////// IGNORE
		printf("WARNING: This is just a test for the finite electron orbits. Set flag TEST_EL = 0 to avoid being here...\n"); 
		if (gamma_ref < TINY) DS_size = SYS_SIZ;
		else if (gamma_ref < 1.0) DS_size = SYS_SIZ/gamma_ref;
		else DS_size = SYS_SIZ;
		size_phiDSgrid = (int) ( DS_size/deltaxDS );
		printf("size of coarse potential grid in Debye sheath = %d\n", size_phiDSgrid);
		fprintf(fout, "size of coarse potential grid in Debye sheath = %d\n", size_phiDSgrid);
		x_DSgrid = malloc(size_phiDSgrid*sizeof(double)); phi_DSgrid = malloc(size_phiDSgrid*sizeof(double));
		ne_DSgrid = malloc(size_phiDSgrid*sizeof(double)); ni_DSgrid = malloc(num_spec*sizeof(double));
		sumni_DSgrid = malloc(size_phiDSgrid*sizeof(double)); vy_e_wall  = malloc((ZOOM_DS*size_phiDSgrid+1)*size_phiDSgrid*sizeof(double));
		mu_e_op  = malloc((ZOOM_DS*size_phiDSgrid+1)*sizeof(double)); chiM_e  = malloc((ZOOM_DS*size_phiDSgrid+1)*sizeof(double));
		twopidmudvy_e  = malloc((ZOOM_DS*size_phiDSgrid+1)*sizeof(double));
		v_cut = 1.0;
		make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, 0, -0.0, 1.0, alpha);
		for (i=0; i < size_phigrid; i++)
			printf("phi_grid(%f)=%f\n", x_grid[i], phi_grid[i]);
		v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]) ;
		make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, 0.0, deltaxDS, 0, -0.5*v_cutDS*v_cutDS, 1.0/gamma_ref, alpha);
		printf("v_cutDS = %f\nphi_DSgrid[0] = %f\n", v_cutDS, phi_DSgrid[0]);
		printf("now evaluate electron density in MPS\n");
		denszeroorb(-1.0, 1.0, phi_grid, ne_grid, size_phigrid, &flux_e, &Q_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, vpar_e_cut, 0.0, x_grid, &ne_inf);
		printf("now evaluate ion density in MPS\n");
		for (n=0; n<num_spec; n++) 
			densfinorb(TioverTe[n], lenfactor[n], alpha, size_phigrid, size_ngrid+n, ne_grid, x_grid, phi_grid, ioncharge, dist_i_GK[n], mu_i[n], U_i[n], size_mu_i[n], size_U_i[n], grid_parameter, flux_i+n, Q_i+n, zoomfactor, -MARGIN_MP, -999.9, vy_i_wall[n], mu_i_op[n], chiM_i[n], twopidmudvy_i[n], size_op_i+n);
		for (ncols=0; ncols<size_mu_e; ncols+=1) {
			for (ind=0; ind<size_vpar_e; ind+=1) {
				U_e_DS[ind] = 0.5*ind*DVPAR*ind*DVPAR;
				Uminmu_MPE = sqrt(2.0*U_e_DS[ind] - 2.0*phi_grid[0]);
				if (ncols == 0) vpar_e_DS[ind] = sqrt(2.0*U_e_DS[ind] - 2.0*phi_grid[0]);
				dist_e_GK[ncols][ind] = bilin_interp(mu_e[ncols], Uminmu_MPE, dist_e_DK, mu_e, vpar_e, size_mu_e, size_vpar_e, -1, -1)/ne_grid[0];
			}
		}
		printf("size grid = %d\n", size_phiDSgrid);
		densfinorb(1.0, 1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e, U_e_DS, size_mu_e, size_vpar_e, 0.0, &flux_eDS, &garbage, ZOOM_DS, -0.5, -999.9, vy_e_wall, mu_e_op, chiM_e, twopidmudvy_e, &size_op_e); 
		printf("size_neDSgrid = %d\n", size_neDSgrid);
		for (i=0; i< size_phiDSgrid; i++) printf("ne = %f at x = %f\n", ne_DSgrid[i], x_DSgrid[i]);
		exit(1);
	}
	/* 
		END OF TEST FOR THE FINITE ELECTRON DENSITY CALCULATION
	*/


/*
		ITERATE MAGNETIC PRESHEATH POTENTIAL TO FIND SOLUTION
		FIRST WITH SIMPLIFIED ELECTRON REFLECTION (CUTOFF) MODEL
*/

	N=0; // set iteration number to zero
	make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, phi0_init_MP, 1.0, alpha); // construct the grid of x and phi (initial guess) values in the magnetic presheath
	while ( ( (convergence_MP <= 1) || ( convergence_j <= 1 && fix_current == 1 ) ) && (N<MAX_IT) ) {
		printf("weight_MP = %f\nITERATION # = %d\n", weight_MP, N);
		fprintf(fout, "weight_MP = %f\nITERATION # = %d\n", weight_MP, N);
		current = target_current;
		if (phi_grid[0] + 0.5*v_cut*v_cut < 0.0) 
			grid_parameter = 0.0;
		else if ( (2.0*factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut) < INITIAL_GRID_PARAMETER) )  {
			grid_parameter = 2.0*factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut);
			
		}
		make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, phi0_init_MP, 1.0, alpha);
		printf("grid_parameter = %f\n", grid_parameter);
		fprintf(fout, "grid parameter = %f\n", grid_parameter);
		printf("\t(phi_mp0, phi_ds0, phi_wall) = (%f, %f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut - phi_grid[0], -0.5*v_cut*v_cut);
		fprintf(fout, "\t(phi_mp0, phi_ds0, phi_wall) = (%f, %f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut - phi_grid[0], -0.5*v_cut*v_cut);

/* 	
		NOW EVALUATE THE ION DENSITY IN THE MAGNETIC PRESHEATH
*/

		printf("evaluate ion density in MPS\n");
		for (i=0; i<size_phigrid; i++) sumni_grid[i] = 0.0;
		sumflux_i = 0.0;
		sumQ_i = 0.0;
		for (n=0; n<num_spec; n++) { // n is ion species index
/*	
		THIS FUNCTION (densfinorb) IS WHERE ALL THE COMPLICATION IS
		THE DENSITY IS EVALUATED ONLY ASSUMING SMALL MAGNETIC FIELD ANGLE AT THE WALL
		ARBITRARY ORBIT DISTORTION (NON-CIRCULAR) IS INCLUDED
*/
			densfinorb(TioverTe[n], lenfactor[n], alpha, size_phigrid, size_ngrid+n, ni_grid[n], x_grid, phi_grid, ioncharge, dist_i_GK[n], mu_i[n], U_i[n], size_mu_i[n], size_U_i[n], grid_parameter, flux_i+n, Q_i+n, zoomfactor, MARGIN_MP, -999.9, vy_i_wall[n], mu_i_op[n], chiM_i[n], twopidmudvy_i[n], size_op_i+n);
			do
				size_sumnigrid = size_ngrid[n] ;
			while (size_sumnigrid > size_ngrid[n]) ;
			sumflux_i += (nioverne[n]*sqrt(TioverTe[n]/mioverme[n])*flux_i[n]);
			sumQ_i += (nioverne[n]*sqrt(TioverTe[n]/mioverme[n])*Q_i[n]);
			for (i=0; i<size_sumnigrid; i++) { // i is position index
				sumni_grid[i] += (nioverne[n]*ni_grid[n][i]);
				//printf("sumni_grid[%d] = %f\n", i, sumni_grid[i]);
			}
		}
		
		printf("densfinorb module ran for all ion species\n");
		printf("size of total ion density grid = %d\n", size_sumnigrid);
		fprintf(fout, "size of total ion density grid = %d\n", size_sumnigrid);

/*
		ION DENSITY CALCULATION FINISHED
*/

		// calculate new electron velocity cutoff (square root of 2 times potential drop across Debye sheath)
		v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
		// if the magnetic presheath potential drop phi_grid[0] exceeds the total potential drop (1/2)vcut^2
		// set the Debye sheath potential drop to zero
		if (v_cut*v_cut < - 2.0*phi_grid[0]) v_cutDS = 0.0;
		printf("v_cutDS = %f\n", v_cutDS);
		//if ( 0.5*v_cutDS*v_cutDS < 0.02) {
		//	weight_DS = WEIGHT_DS/30.0;
		//	weight_MP = WEIGHT_MP/30.0;
		//	weight_j = WEIGHT_j/30.0;
		//}
		//else if ( 0.5*v_cutDS*v_cutDS < 0.1) {
		//	weight_DS = WEIGHT_DS/10.0;
		//	weight_MP = WEIGHT_MP/10.0;
		//	weight_j = WEIGHT_j/10.0;
		//}
		//else if (0.5*v_cutDS*v_cutDS < 0.25) {
		//	weight_DS = WEIGHT_DS/3.0;
		//	weight_MP = WEIGHT_MP/3.0;
		//	weight_j = WEIGHT_j/3.0;
		//}
		if (0.5*v_cutDS*v_cutDS < 0.5) { // adjust weights in the iteration
			//weight_DS = WEIGHT_DS/3.0;
			weight_MP = pow(1.0*v_cutDS*v_cutDS, 0.5)*WEIGHT_MP;
			weight_j = pow(1.0*v_cutDS*v_cutDS, 0.5)*WEIGHT_j;
			weight_j = v_cutDS*WEIGHT_j;
			if (weight_MP < 0.00001) weight_MP = 0.00001;
			if (weight_j < 0.00001) weight_j = 0.00001;
		}
		if ( (N%100 == 0) && (N!=0) ) {	 // adjust weights
			weight_MP /= 2.0;
			weight_j /= 2.0;
		}
		if (gamma_ref >= TINY) { // ELECTRON CUTOFF MODEL
			// gamma = infinity model
			//vpar_cut_lookup[0] = 1e10; 
			//mue_cut_lookup[0] = 0.0;
			//for (i=1; i< size_cut; i++) {
			//	vpar_cut_lookup[i] = vparcut(M_PI - i*M_PI/size_cut, v_cutDS);
			//	mue_cut_lookup[i] = mucut(M_PI - i*M_PI/size_cut, v_cutDS); }
			//vpar_cut_lookup[size_cut] = 0.0;
			//mue_cut_lookup[size_cut] = 1e10;
			// finite gamma model
			for (i=0; i< size_mu_e; i++) {
				vpar_e_cut[i] = vparcut_mu(mu_e[i], v_cutDS);
				if (0.5*v_cutDS*v_cutDS - 0.54*gamma_ref*sqrt(2.0*mu_e[i])*pow(v_cutDS*v_cutDS*0.5,0.7) > 0.001)
					vpar_e_cut[i] = sqrt( v_cutDS*v_cutDS - 2.0*0.54*gamma_ref*sqrt(2.0*mu_e[i])*pow(v_cutDS*v_cutDS*0.5,0.7) );
				else vpar_e_cut[i] = 0.0;
				vpar_e_cut[i] = ( v_cutDS - (1/v_cutDS)*0.54*gamma_ref*sqrt(2.0*mu_e[i])*pow(v_cutDS*v_cutDS*0.5,0.7) );
				if (DEBUG == 1)
					printf("vpar = %f for mu = %f\n", vpar_e_cut[i], mu_e[i]);
			}
		}
		else if (gamma_ref < TINY) { // first solve MPS w/ simplified e- reflection
		      for (i=0; i< size_mu_e; i++) 
			vpar_e_cut[i] = v_cutDS;
		}
/* 
		ROBBIE's ELECTRON DENSITY CALCULATION (ADAPTED TO 2D)
*/
		denszeroorb(-1.0, 1.0, phi_grid, ne_grid, size_phigrid, &flux_e, &Q_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, vpar_e_cut, 0.0, x_grid, &ne_inf);

		// TARGET BOHM CONDITION 
		Bohmshouldbe = TioverTe[0]*(ne_grid[1] - ne_grid[0])/( 0.5*(ne_grid[0]+ne_grid[1])*v_cutDS*(sqrt(2.0*(phi_grid[1] - phi_grid[0]) + v_cutDS*v_cutDS) - v_cutDS) ); // phi_DSgrid[0])*ne_grid[0]);

/* 		
		SAVE DATA TO FILES IN ITERATION DIRECTORIES
		THIS IS JUST TO CREATE GIFS THAT SHOW THE 
		EVOLUTION AS THE ITERATION GOES
*/
		snprintf(dirnameit, 150, "%s%d", dirname_it, N);
		mkdir(dirnameit, S_IRWXU);
		printf("%s%d/phi_n_MP.txt\n", dirname_it, N);
		snprintf(fpstr, 150, "%s%d/phi_n_MP.txt", dirname_it, N);
		fp = fopen(fpstr, "w");
		if (fp == NULL) {
			printf("error when opening file %s\n", fpstr);
		}
		for (i=0; i<size_phigrid; i++) {
			fprintf(fp, "%f %f %f %f\n", x_grid[i], phi_grid[i], sumni_grid[i], ne_grid[i]);
		}
		fclose(fp);
		snprintf(fpstr, 150, "%s%d/Fi_W.txt", dirname_it, N);
		fp = fopen(fpstr, "w");
		if (fp == NULL) {
			printf("error when opening file\n");
		}
		// SAVE TARGET DISTRIBUTION FUNCTION INFORMATION
		for (i=0; i<size_op_i[0]; i++) { // at the moment assumes 1 ion species, hence index 0
			fprintf(fp, "%f %f %f %f\n", vy_i_wall[0][i], mu_i_op[0][i], chiM_i[0][i], twopidmudvy_i[0][i]);
		}
		fclose(fp);
		// SAVE IMPORTANT MISCELLANEOUS OUTPUT
		snprintf(fpstr, 150, "%s%d/misc_output.txt", dirname_it, N);
		fp = fopen(fpstr, "w");
		if (fp == NULL) {
			printf("error when opening file\n");
		}
		current = sumflux_i - flux_e;
		fprintf(fp, "%f\n", current); // current at target
		fprintf(fp, "%f\n", 0.5*v_cut*v_cut); // potential drop across magnetic presheath + Debye sheath
		fprintf(fp, "%f\n", Q_e); // electron heat flux
		fprintf(fp, "%f\n", sumQ_i); // ion heat flux
		fprintf(fp, "%f\n", flux_e); // electron particle flux/current
		fprintf(fp, "%f\n", sumflux_i); // total ion particle flux/current (although must include Z in sumflux for multiple species)

		/* 
			HAVE ION AND ELECTRON DENSITY
			CALCULATE THE ERROR IN QUASINEUTRALITY
			BASED ON ERROR, DECIDE WHAT TO DO
		*/
		error_Poisson(error_MP, x_grid, ne_grid, sumni_grid, nioverne, phi_grid, size_phigrid, size_sumnigrid, 0.0);
		printf("error_av = %f\terror_max = %f\n", error_MP[0], error_MP[1]);
		// if error thresholds are satisfied, increase convergence_MP flag by 1
		if ( (error_MP[0] < tol_MP[0]) && (error_MP[1] < tol_MP[1]) ) convergence_MP += 1 ; 
		else convergence_MP = 0;
		if (convergence_MP == 0) { // calculate new guess
			printf("MP not converged --> calculate new MP potential guess\n");
			fprintf(fout, "MP not converged --> calculate new MP potential guess\n");
			newguess(x_grid, ne_grid, sumni_grid, phi_grid, size_phigrid, size_sumnigrid, 0.0, 0.0, 1.5, weight_MP);// p, m);
		}
		else { // don't calculate new guess, potential is converged
			printf("MP converged --> no iteration needed\n");
			fprintf(fout, "MP converged --> no iteration needed\n");
		}
		if (fix_current == 1) { // calculate new guess for total potential drop
			if (fabs(target_current - current) > tol_current*sumflux_i) convergence_j = 0;
			else convergence_j += 1;
			printf("target current = %f +/- %f\n", target_current, tol_current*sumflux_i);
			fprintf(fout, "target current = %f +/- %f\n", target_current, tol_current*sumflux_i);
			printf("current = %f = ion current (%f) - electron current (%f) = %f x electron current\n", current, sumflux_i, flux_e, current/flux_e);
			fprintf(fout, "current = %f = ion current (%f) - electron current (%f) = %f x electron current\n", current, sumflux_i, flux_e, current/flux_e);
			if (convergence_j == 0) {
				fprintf(fout, "new WALL potential guess\n");
				printf("new WALL potential guess\n");
				newvcut(&v_cut, v_cutDS, sumflux_i, flux_e, target_current, tol_current, weight_j);
			}
			else {
				printf("current converged, NO wall potential iteration\n");
				fprintf(fout, "current converged, NO wall potential iteration\n");
			}
		}
		else {
			printf("current = %f = ion current (%f) - electron current (%f) = %f x electron current\n", current, sumflux_i, flux_e, current/flux_e);
			fprintf(fout, "current = %f = ion current (%f) - electron current (%f) = %f x electron current\n", current, sumflux_i, flux_e, current/flux_e);
		}
		printf("\tcurrent = %f (target = %f)\n", current, target_current);
		fprintf(fout, "\tcurrent = %f (target = %f)\n", current, target_current);
		printf("convergence (MP, j) = (%d, %d)\n", convergence_MP, convergence_j);

		//if ( (convergence_MP > 1) && ( convergence_j > 1 || fix_current == 0) ) {
		//	for (i=0; i<size_sumnigrid; i++) sumni_grid[i] = 0.0;
		//	sumflux_i = 0.0;
		//	sumQ_i = 0.0;
		//	for (n=0; n<num_spec; n++) {
		//		densfinorb(TioverTe[n], lenfactor[n], alpha, size_phigrid, size_ngrid+n, ni_grid[n], x_grid, phi_grid, ioncharge, dist_i_GK[n], mu_i[n], U_i[n], size_mu_i[n], size_U_i[n], grid_parameter, flux_i+n, Q_i+n, zoomfactor, MARGIN_MP, -999.9, vy_i_wall[n], mu_i_op[n], chiM_i[n], twopidmudvy_i[n], size_op_i+n);
		//		do
		//			size_sumnigrid = size_ngrid[n] ;
		//		while (size_sumnigrid > size_ngrid[n]) ;
		//		sumflux_i += (nioverne[n]*sqrt(TioverTe[n]/mioverme[n])*flux_i[n]);
		//		sumQ_i += (nioverne[n]*sqrt(TioverTe[n]/mioverme[n])*Q_i[n]);
		//		for (i=0; i<size_sumnigrid; i++)
		//			sumni_grid[i] += (nioverne[n]*ni_grid[n][i]);
		//	}
		//	for (n=0; n<num_spec; n++) 
		//		densfinorb(TioverTe[n], lenfactor[n], alpha, size_phigrid, size_ngrid+n, ni_grid[n], x_grid, phi_grid, ioncharge, dist_i_GK[n], mu_i[n], U_i[n], size_mu_i[n], size_U_i[n], grid_parameter, flux_i+n, Q_i+n, zoomfactor, MARGIN_MP, -999.9, vy_i_wall[n], mu_i_op[n], chiM_i[n], twopidmudvy_i[n], size_op_i+n);
		//	error_Poisson(error_MP, x_grid, ne_grid, sumni_grid, nioverne, phi_grid, size_phigrid, size_sumnigrid, 0.0);
		//	printf("error_av = %f\terror_max = %f\n", error_MP[0], error_MP[1]);
		//	Bohmeval(&Bohm, alpha, TioverTe[0], phi_grid[0], dist_i_GK[0], mu_i[0], U_i[0], vy_i_wall[0], mu_i_op[0], chiM_i[0], twopidmudvy_i[0], size_mu_i[0], size_U_i[0],  size_op_i[0]); 
		//	printf("Bohm integral should be %f\n", Bohmshouldbe);
		//	evalBohmshouldbe(&Bohmshouldbe, phi_grid[0], dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, vpar_e_cut);
		//	printf("Bohm integral should be %f\n", Bohmshouldbe);
		//	fprintf(fBohm, "%f %f\n", Bohm, Bohmshouldbe);
		//	fprintf(fout, "Bohm integral of converged MP solution should be %f\n", Bohmshouldbe);
		//} 
		//else N++;
		N++;
		//Bohmeval(&Bohm, alpha, TioverTe[0], phi_grid[0], dist_i_GK[0], mu_i[0], U_i[0], vy_i_wall[0], mu_i_op[0], chiM_i[0], twopidmudvy_i[0], size_mu_i[0], size_U_i[0], size_op_i[0]); 
		printf("Bohm integral is %f\n", Bohm);
		evalBohmshouldbe(&Bohmshouldbe, phi_grid[0], dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, vpar_e_cut);
		printf("Bohm integral should be %f\n", Bohmshouldbe);
		printf("electron heat flux Q_e = %f\n", Q_e);
		printf("electron particle flux Phi = %f\n", flux_e);
		printf("electron sheath heat transmission coefficient??? = %f\n", Q_e/flux_e);
	}
	clock_t end_it = clock(); // finds end time of last iteration
	tot_time = (double) (end_it - begin_it) / CLOCKS_PER_SEC;
	printf("At %dth iteration MP iteration with simplified DS model converged successfully in %f seconds\n", N, tot_time);
	fprintf(fout, "At %dth iteration MP iteration with simplified DS model converged successfully in %f seconds\n", N, tot_time);
	printf("FINAL CHECK on accuracy of electrostatic potential solution in MP with simplified electron reflection model in DS\n");
	fprintf(fout, "FINAL CHECK on accuracy of electrostatic potential solution in MP with simplified electron reflection model in DS\n");
	for (i=0; i<size_sumnigrid; i++) sumni_grid[i] = 0.0;
	sumflux_i = 0.0;
	sumQ_i = 0.0;
	for (n=0; n<num_spec; n++) {
		// CALCULATE DENSITY OVER ENTIRE DOMAIN 
		// UNTIL 6 or 7 rho_i FROM THE DOMAIN CEILING
		densfinorb(TioverTe[n], lenfactor[n], alpha, size_phigrid, size_ngrid+n, ni_grid[n], x_grid, phi_grid, ioncharge, dist_i_GK[n], mu_i[n], U_i[n], size_mu_i[n], size_U_i[n], grid_parameter, flux_i+n, Q_i+n, zoomfactor, -1.0, -999.9, vy_i_wall[n], mu_i_op[n], chiM_i[n], twopidmudvy_i[n], size_op_i+n);
		do
			size_sumnigrid = size_ngrid[n] ;
		while (size_sumnigrid > size_ngrid[n]) ;
		sumflux_i += (nioverne[n]*sqrt(TioverTe[n]/mioverme[n])*flux_i[n]);
		sumQ_i += (nioverne[n]*sqrt(TioverTe[n]/mioverme[n])*Q_i[n]);
		for (i=0; i<size_sumnigrid; i++)
			sumni_grid[i] += (nioverne[n]*ni_grid[n][i]);
	}
	error_Poisson(error_MP, x_grid, ne_grid, sumni_grid, nioverne,  phi_grid, size_phigrid, size_sumnigrid, 0.0);
	printf("error_av = %f\terror_max = %f\n", error_MP[0], error_MP[1]);
	//if ( error_MP[1] > tol_MP[1] ) {
	//	printf("ERROR: MP solution rejected because it does not satisfy Poisson's equation accurately enough on the extended domain\n");
	//	fprintf(fout, "ERROR: MP solution rejected because it does not satisfy Poisson's equation accurately enough on the extended domain\n");
	//	exit(-1);
	//}
	//printf("FINAL CHECK passed. HURRAY!\n");
	//}

	if ( (gamma_ref <= 5.0) && (gamma_ref >= 0.09) ) {
	// FULL DEBYE SHEATH SOLUTION CALCULATED WITH FINITE (DISTORTED) ELECTRON GYROORBITS
		// form x, phi, ne and ni grids for Debye sheath
		if (gamma_ref < TINY) DS_size = SYS_SIZ;
		else if (gamma_ref < 1.0) DS_size = SYS_SIZ/gamma_ref;
		else DS_size = SYS_SIZ;
		size_phiDSgrid = (int) ( DS_size/deltaxDS );
		printf("size of coarse potential grid in Debye sheath = %d\n", size_phiDSgrid);
		fprintf(fout, "size of coarse potential grid in Debye sheath = %d\n", size_phiDSgrid);
		x_DSgrid = malloc(size_phiDSgrid*sizeof(double));
		phi_DSgrid = malloc(size_phiDSgrid*sizeof(double));
		ne_DSgrid = malloc(size_phiDSgrid*sizeof(double));
		ni_DSgrid = malloc(num_spec*sizeof(double));
		sumni_DSgrid = malloc(size_phiDSgrid*sizeof(double));
		vy_e_wall  = malloc(size_phiDSgrid*sizeof(double));
		mu_e_op  = malloc(size_phiDSgrid*sizeof(double));
		chiM_e  = malloc(size_phiDSgrid*sizeof(double));
		twopidmudvy_e  = malloc(size_phiDSgrid*sizeof(double));
		vpar_e_cut_lookup  = malloc(size_phiDSgrid*sizeof(double));
		for (n=0; n<num_spec; n++) ni_DSgrid[n] = malloc(size_phiDSgrid*sizeof(double));
		//make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, phi0_init_MP, 1.0, alpha);
		if (gamma_ref < 1.0) 
			make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, 0.0, deltaxDS, 0, -phi_grid[0] - 0.5*v_cut*v_cut, 1.0/gamma_ref, alpha);
		else 
			make_phigrid(x_DSgrid, phi_DSgrid, size_phiDSgrid, 0.0, deltaxDS, 0, -phi_grid[0] - 0.5*v_cut*v_cut, 1.0, alpha);
		v_cutDS = sqrt(v_cut*v_cut + 2.0*phi_grid[0]);
		if (gamma_ref >= TINY) {
			for (i=0; i< size_mu_e; i++) {
				vpar_e_cut[i] = vparcut_mu(mu_e[i], v_cutDS);
				if (0.5*v_cutDS*v_cutDS - 0.54*gamma_ref*sqrt(2.0*mu_e[i])*pow(v_cutDS*v_cutDS*0.5,0.7) > 0.001)
					vpar_e_cut[i] = sqrt( v_cutDS*v_cutDS - 2.0*0.54*gamma_ref*sqrt(2.0*mu_e[i])*pow(v_cutDS*v_cutDS*0.5,0.7) );
				else vpar_e_cut[i] = 0.0;
				vpar_e_cut[i] = ( v_cutDS - (1/v_cutDS)*0.54*gamma_ref*sqrt(2.0*mu_e[i])*pow(v_cutDS*v_cutDS*0.5,0.7) );
			}
		}
		else if (gamma_ref < TINY) { // first solve MPS w/ simplified e- reflection
		      for (i=0; i< size_mu_e; i++) 
			vpar_e_cut[i] = v_cutDS;
		}
		N_DS = 0;
		convergence_MP = convergence_j = 0;
		weight_MP = WEIGHT_MP;
		weight_j = WEIGHT_j;//*2.0/3.0;
		error_MP[0] = error_DS[0] = 100000.0;
		if (fix_current == 0) convergence_j=2;
		while ( ( ( convergence_MP <= 1) || (convergence_DS <= 1) || (convergence_j <= 1) ) && (N_DS < MAX_IT) ) {
			printf("weight_(MP,DS,j) = (%f, %f, %f)\n", weight_MP, weight_DS, weight_j);
			make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, phi0_init_MP, 1.0, alpha);
			printf("in while loop for combined DS+MP iteration\n");
			printf("v_cutDS = %f\n", v_cutDS);
			printf("evaluate electron density in MPS\n");
			denszeroorb(-1.0, 1.0, phi_grid, ne_grid, size_phigrid, &flux_e, &Q_e, dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, vpar_e_cut, 0.0, x_grid, &ne_inf);
			//Bohmshouldbe = TioverTe[0]*(ne_grid[1] - ne_grid[0])/((phi_grid[1] - phi_grid[0])*ne_grid[0]);
			Bohmshouldbe = TioverTe[0]*(ne_grid[1] - ne_grid[0])/( 0.5*(ne_grid[0]+ne_grid[1])*v_cutDS*(sqrt(2.0*(phi_grid[1] - phi_grid[0]) + v_cutDS*v_cutDS) - v_cutDS) ); // phi_DSgrid[0])*ne_grid[0]);
			printf("evaluate ion density in MP\n");
			for (i=0; i<size_sumnigrid; i++) sumni_grid[i] = 0.0;
			sumflux_i = 0.0;
			sumQ_i = 0.0;
			for (n=0; n<num_spec; n++) {
				densfinorb(TioverTe[n], lenfactor[n], alpha, size_phigrid, size_ngrid+n, ni_grid[n], x_grid, phi_grid, ioncharge, dist_i_GK[n], mu_i[n], U_i[n], size_mu_i[n], size_U_i[n], grid_parameter, flux_i+n, Q_i+n, zoomfactor, MARGIN_MP, -999.9, vy_i_wall[n], mu_i_op[n], chiM_i[n], twopidmudvy_i[n], size_op_i+n);
				do
					size_sumnigrid = size_ngrid[n] ;
				while (size_sumnigrid > size_ngrid[n]) ;
				sumflux_i += (nioverne[n]*sqrt(TioverTe[n]/mioverme[n])*flux_i[n]);
				sumQ_i += (nioverne[n]*sqrt(TioverTe[n]/mioverme[n])*Q_i[n]);
				for (i=0; i<size_sumnigrid; i++)
					sumni_grid[i] += (nioverne[n]*ni_grid[n][i]);
			}
			printf("flux_eDS = (%f, %f)\tflux_i = %f\n", flux_eDS*ne_grid[0], flux_e, sumflux_i);
			if (phi_grid[0] + 0.5*v_cut*v_cut < 0.0) 
				grid_parameter = 0.0;
			else if ( (2.0*factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut) < INITIAL_GRID_PARAMETER) )  {// && (phi_grid[0] + 0.5*v_cut*v_cut >1.0e-13 ) )   
				grid_parameter = 2.0*factor_small_grid_parameter*(phi_grid[0] + 0.5*v_cut*v_cut);
			}
			make_phigrid(x_grid, phi_grid, size_phigrid, grid_parameter, deltax, N, phi0_init_MP, 1.0, alpha);
			printf("grid_parameter = %f\n", grid_parameter);
			printf("evaluate ion density in DS\n");
			for (i=0; i<size_phiDSgrid; i++) sumni_DSgrid[i] = 0.0;
			for (n=0; n<num_spec; n++) {
				//printf("n (species index) = %d\n", n);
				densionDS(alpha, TioverTe[n], &Bohm, ni_DSgrid[n], phi_DSgrid, phi_grid[0], dist_i_GK[n], mu_i[n], U_i[n], vy_i_wall[n], mu_i_op[n], chiM_i[n], twopidmudvy_i[n], size_phiDSgrid, size_mu_i[n], size_U_i[n], size_op_i[n]);
				//printf("nioverne[0]=%f\n", nioverne[n]);
				for (i=0; i<size_phiDSgrid; i++){
					sumni_DSgrid[i] += (nioverne[n]*ni_DSgrid[n][i]);
					//printf("niDS[%d] = %f, sumni = %f\n", i, ni_DSgrid[n][i], sumni_DSgrid[i]);
				}
			}
			printf("evaluate electron density in DS\n");
			for (ncols=0; ncols<size_mu_e; ncols+=1) {
				for (ind=0; ind<size_vpar_e; ind+=1) {
					U_e_DS[ind] = 0.5*ind*DVPAR*ind*DVPAR;
					Uminmu_MPE = sqrt(2.0*U_e_DS[ind] - 2.0*phi_grid[0]);
					if (ncols == 0) vpar_e_DS[ind] = Uminmu_MPE;
					dist_e_GK[ncols][ind] = bilin_interp(mu_e[ncols], Uminmu_MPE, dist_e_DK, mu_e, vpar_e, size_mu_e, size_vpar_e, -1, -1)/(ne_inf*ne_grid[0]);
				}
			}
			printf("size_phiDSgrid = %d\n", size_phiDSgrid);
			EW = phi_DSgrid[1] - phi_DSgrid[0]; EW/=(x_DSgrid[1] - x_DSgrid[0]);
			if (gamma_ref > SMALLGAMMA) {
				densfinorb(1.0, 1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e, U_e_DS, size_mu_e, size_vpar_e, 0.0, &flux_eDS, &garbage, ZOOM_DS, MARGIN_DS, -999.9, vy_e_wall, mu_e_op, chiM_e, twopidmudvy_e, &size_op_e); 
				//for (i=0; i<size_neDSgrid; i++)
				//	ne_DSgrid[i] *= (sumni_DSgrid[size_neDSgrid-1]/ne_DSgrid[size_neDSgrid-1]);
				//printf("chiM mu\n");
				for (i=0; i<size_op_e; i++) {
					vpar_e_cut_lookup[i] = sqrt(2.0*(chiM_e[i] - mu_e_op[i]));
					//printf("%f %f \n", chiM_e[i], mu_e[i]);
				} 
				for (i=0; i<size_mu_e; i++) {
					vpar_e_cut[i] = lin_interp(mu_e_op, vpar_e_cut_lookup, mu_e[i], size_op_e, -1);
					//vpar_e_cut[i] = sqrt(2.0*(chiM_e[i] - mu_e[i]));
					printf("%f %f %f\n", chiM_e[i], mu_e[i], vpar_e_cut[i]);
				} 
				printf("flux_eDS = %f\n", flux_eDS);
			}
			else {
				for (i=0; i< size_mu_e; i++) 
					vpar_e_cut[i] = 0.0; //sqrt(-2.0*lin_interp(x_DSgrid, phi_DSgrid, sqrt(mue_cut_lookup[i]), size_phiDSgrid, 1234)); 
				denszeroorb(-1.0, 1.0, phi_DSgrid, ne_DSgrid, size_phiDSgrid, &flux_eDS, &Q_eDS, dist_e_GK, vpar_e_DS, mu_e, size_vpar_e, size_mu_e, vpar_e_cut, gamma_ref, x_DSgrid, &ne_inf); // needs dist_e_GK which in this case acts like dist_e_DK_DS
				//for (i=0; i< size_mu_e; i++) 
				//	vpar_e_cut[i] = sqrt(-2.0*lin_interp(x_DSgrid, phi_DSgrid, sqrt(2.0*mu_e[i]), size_phiDSgrid, 1234)); 
				for (i=0; i<size_mu_e; i++) {
					//vpar_e_cut[i] = sqrt(2.0*(-phi_DSgrid[0]*exp(((phi_DSgrid[1]-phi_DSgrid[0])/(phi_DSgrid[0]*deltaxDS))*sqrt(2.0*mu_e[i]))));
					//if (-phi_DSgrid[0] - 0.54*gamma_ref*sqrt(mu_e[i])*pow(-phi_DSgrid[0],0.7) < 0.0)
					//	vpar_e_cut[i] = sqrt(2.0*(-phi_DSgrid[0] - 0.54*gamma_ref*sqrt(mu_e[i])*pow(-phi_DSgrid[0],0.7)));
					//else vpar_e_cut[i] = 0.0;
					vpar_e_cut[i] = sqrt(2.0*(-phi_DSgrid[0] + 0.5*EW*EW)*exp(-EW*sqrt(2.0*mu_e[i])/(-phi_DSgrid[0]))*gsl_sf_bessel_I0(EW*sqrt(2.0*mu_e[i])/(-phi_DSgrid[0])));
					//if (-phi_DSgrid[0] < ((phi_DSgrid[1]-phi_DSgrid[0])/deltaxDS)*sqrt(2.0*mu_e[i])) 
					//	vpar_e_cut[i] = 0.0;
					//else 
					//	vpar_e_cut[i] = sqrt(2.0*(-phi_DSgrid[0] - ((phi_DSgrid[1]-phi_DSgrid[0])/deltaxDS)*sqrt(2.0*mu_e[i])));
				}
				i=0;
				while (ne_DSgrid[i] < 1.0-MARGIN_DS) i++;
				size_neDSgrid = i;
				printf("size_neDSgrid = %d\n", size_neDSgrid);
			}
			error_Poisson(error_MP, x_grid, ne_grid, sumni_grid, nioverne, phi_grid, size_phigrid, size_sumnigrid, 0.0);
			printf("error_av = %f\terror_max = %f\n", error_MP[0], error_MP[1]);
			if ( (error_MP[0] < tol_MP[0]) && (error_MP[1] < tol_MP[1]) ) convergence_MP += 1 ;
			else 
				convergence_MP = 0;
			if (convergence_MP == 0) { // || (convergence_j == 0) ) 
				newguess(x_grid, ne_grid, sumni_grid, phi_grid, size_phigrid, size_sumnigrid, 0.0, 0.0, 1.5, weight_MP);
				printf("MP not converged\n");
				fprintf(fout, "MP not converged\n");
			}
			else if (convergence_MP == 1) {
				//newguess(x_grid, ne_grid, sumni_grid, phi_grid, size_phigrid, size_sumnigrid, 0.0, 0.0, 1.5, weight_MP);
				//weight_MP/=2.0;
				printf("MP converged\n");
				fprintf(fout, "MP converged\n");
			}
			else {
				//newguess(x_grid, ne_grid, sumni_grid, phi_grid, size_phigrid, size_sumnigrid, 0.0, 0.0, 1.5, weight_MP);
				printf("MP converged\n");
				fprintf(fout, "MP converged\n");
			}
			
			current = sumflux_i - flux_e; 
			if (fix_current == 1) {
				if (fabs(target_current - current) > tol_current*sumflux_i) 
					convergence_j = 0;
				else convergence_j += 1;
				printf("target current = %f +/- %f\n", target_current, tol_current*sumflux_i);
				fprintf(fout, "target current = %f +/- %f\n", target_current, tol_current*sumflux_i);
				printf("current = %f = ion current (%f) - electron current (%f) = %f x electron current\n", current, sumflux_i, flux_e, current/flux_e);
				fprintf(fout, "current = %f = ion current (%f) - electron current (%f) = %f x electron current\n", current, sumflux_i, flux_e, current/flux_e);
				//if ( (convergence_MP >= 1) && (convergence_DS >= 1) ) {
				if (convergence_j == 0) {
					printf("new WALL potential guess\n");
					fprintf(fout, "new WALL potential guess\n");
					newvcut(&v_cut, v_cutDS, sumflux_i, flux_e, target_current, tol_current, weight_j);
				}
				else if (convergence_j == 1) {
					//if ( (N_DS != 0) && (weight_j > WEIGHT_j/10.0) ) weight_j *= 0.95;
					//weight_j/=2.0;
					printf("current converged, NO wall potential iteration\n");
					fprintf(fout, "current converged, NO wall potential iteration\n");
				}
				else {
					printf("current converged, NO wall potential iteration\n");
					fprintf(fout, "current converged, NO wall potential iteration\n");
				}
				//}
			}
			else {
				printf("current = %f = %f x ion current", current, current/flux_e);
				fprintf(fout, "current = %f = %f x ion current", current, current/flux_e);
			}

			printf("new DS potential guess\n");
// Calculating v_cutDS after new phiwall and new phi_MP[0] have been calculated, to ensure that the new imposed phi_DS[0] is consistent
			if (v_cut*v_cut < - 2.0*phi_grid[0]) {
				printf("WARNING: The total sheath and presheath potential drop is smaller than the presheath potential drop\n");
				printf("\tTo avoid a non-monotonic potential, setting the total potential drop to be just above the presheath potential drop\n");
				v_cut = sqrt(-2.0*phi_grid[0] + 0.01);
				v_cutDS = sqrt(0.01);
			}
			else v_cutDS = sqrt(2.0*phi_grid[0] + v_cut*v_cut);

			//if ( 0.5*v_cutDS*v_cutDS < 0.1) {
			//	weight_DS = WEIGHT_DS/10.0;
			//	weight_MP = WEIGHT_MP/10.0;
			//	weight_j = WEIGHT_j/10.0;
			//}
			//else if (0.5*v_cutDS*v_cutDS < 0.25) {
			//	weight_DS = WEIGHT_DS/3.0;
			//	weight_MP = WEIGHT_MP/3.0;
			//	weight_j = WEIGHT_j/3.0;
			//}
			printf("0.5*v_cutDS*v_cutDS = %f\n", 0.5*v_cutDS*v_cutDS);
			if (0.5*v_cutDS*v_cutDS < 0.5) {
				weight_MP = v_cutDS*WEIGHT_MP;
				weight_j = v_cutDS*WEIGHT_j;
				if (weight_MP < 0.00001) weight_MP = 0.00001;
				if (weight_DS < 0.00001) weight_DS = 0.00001;
				if (weight_j < 0.00001) weight_j = 0.00001;
			}
			if ( (N_DS%100 == 0) && (N_DS!=0) ) {	
				//weight_MP /= 2.0;
				//weight_j /= 2.0;
				//weight_DS /= 2.0;
			}
			error_Poisson(error_DS, x_DSgrid, ne_DSgrid, sumni_DSgrid, nioverne, phi_DSgrid, size_phiDSgrid, size_neDSgrid, 1.0/(gamma_ref*gamma_ref));
			printf("error_av = %f\terror_max = %f\n", error_DS[0], error_DS[1]);


			snprintf(dirnameit, 150, "%s%d", dirname_it, N);
			mkdir(dirnameit, S_IRWXU);
			printf("%s%d/phi_n_DS.txt\n", dirname_it, N);
			snprintf(fpstr, 150, "%s%d/phi_n_DS.txt", dirname_it, N);
			fp = fopen(fpstr, "w");
			if (fp == NULL)  
				printf("error when opening file %s\n", fpstr);
			for (i=0; i<size_phiDSgrid; i++) {
				fprintf(fp, "%f %f %f %f\n", x_DSgrid[i], phi_DSgrid[i], sumni_DSgrid[i], ne_DSgrid[i]);
			}
			fclose(fp);
			snprintf(fpstr, 150, "%s%d/vparcut.txt", dirname_it, N);
			fp = fopen(fpstr, "w");
			if (fp == NULL) {
				printf("error when opening file %s\n", fpstr);
			}
			for (i=0; i<=size_mu_e; i++) {
				fprintf(fp, "%f %f\n", mu_e[i], vpar_e_cut[i]);
			}
			fclose(fp);
			snprintf(fpstr, 150, "%s%d/phi_n_MP.txt", dirname_it, N);
			fp = fopen(fpstr, "w");
			if (fp == NULL) {
				printf("error when opening file %s\n", fpstr);
			}
			for (i=0; i<size_phigrid; i++) {
				fprintf(fp, "%f %f %f %f\n", x_grid[i], phi_grid[i], sumni_grid[i], ne_grid[i]);
			}
			fclose(fp);
			snprintf(fpstr, 150, "%s%d/Fi_W.txt", dirname_it, N);
			fp = fopen(fpstr, "w");
			if (fp == NULL) {
				printf("error when opening file %s\n", fpstr);
			}
			for (i=0; i<size_op_i[0]; i++) {
				fprintf(fp, "%f %f %f %f\n", vy_i_wall[0][i], mu_i_op[0][i], chiM_i[0][i], twopidmudvy_i[0][i]);
			}
			fclose(fp);
			snprintf(fpstr, 150, "%s%d/misc_output.txt", dirname_it, N);
			fp = fopen(fpstr, "w");
			if (fp == NULL) {
				printf("error when opening file %s\n", fpstr);
			}
			fprintf(fp, "%f\n", sumflux_i-flux_e);
			fprintf(fp, "%f\n", 0.5*v_cut*v_cut);
			fprintf(fp, "%f\n", Q_e);
			fprintf(fp, "%f\n", sumQ_i);
			fprintf(fp, "%f\n", flux_e);
			fprintf(fp, "%f\n", sumflux_i);

			if ( (error_DS[0] < tol_DS[0]) && (error_DS[1] < tol_DS[1]) ) convergence_DS += 1 ;
			else convergence_DS = 0;
			if (convergence_DS == 0 || convergence_MP == 0 || convergence_j == 0) { 
				printf("phi_DSgrid[0] = %f\n", phi_DSgrid[0]);
				newguess(x_DSgrid, ne_DSgrid, sumni_DSgrid, phi_DSgrid, size_phiDSgrid, size_neDSgrid, 1.0/(gamma_ref*gamma_ref), v_cutDS, 2.0, weight_DS);// p, m);
				printf("DS not converged\n");
				fprintf(fout, "DS not converged\n");
			}
			else if (convergence_DS == 1) {
				//newguess(x_DSgrid, ne_DSgrid, sumni_DSgrid, phi_DSgrid, size_phiDSgrid, size_neDSgrid, 1.0/(gamma_ref*gamma_ref), v_cutDS, 2.0, weight_DS);// p, m);
				//weight_DS /= 2.0;
				printf("DS converged\n");
				fprintf(fout, "DS converged\n");
			}
			else {
				//newguess(x_DSgrid, ne_DSgrid, sumni_DSgrid, phi_DSgrid, size_phiDSgrid, size_neDSgrid, 1.0/(gamma_ref*gamma_ref), v_cutDS, 2.0, weight_DS);// p, m);
				printf("DS converged\n");
				fprintf(fout, "DS converged\n");
			}
			printf("\t(phi_mp0, phi_ds0, phi_wall) = (%f, %f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut - phi_grid[0], -0.5*v_cut*v_cut);
			fprintf(fout, "\t(phi_mp0, phi_ds0, phi_wall) = (%f, %f, %f)\n", phi_grid[0], -0.5*v_cut*v_cut - phi_grid[0], -0.5*v_cut*v_cut);
			printf("ITERATION N = %d of which N_DS = %d\n\n", N, N_DS);
			fprintf(fout, "ITERATION N = %d of which N_DS = %d\n\n", N, N_DS);
			printf("convergence (MP,DS, j) = (%d, %d, %d)\n\tfull convergence achieved when all of these numbers are 2 or above\n", convergence_MP, convergence_DS, convergence_j);
			fprintf(fout, "convergence (MP,DS, j) = (%d, %d, %d)\n\tfull convergence achieved when all of these numbers are 2 or above\n", convergence_MP, convergence_DS, convergence_j);

			N_DS++; N++;
			//gsl_permutation_free (p);
			//gsl_matrix_free (m);
		}
		printf("FINAL CHECK on accuracy of electrostatic potential solutions in MP and DS\n");
		fprintf(fout, "FINAL CHECK on accuracy of electrostatic potential solutions in MP and DS\n");
		for (i=0; i<size_sumnigrid; i++) sumni_grid[i] = 0.0;
		sumflux_i = 0.0;
		sumQ_i = 0.0;
		for (n=0; n<num_spec; n++) {
			densfinorb(TioverTe[n], lenfactor[n], alpha, size_phigrid, size_ngrid+n, ni_grid[n], x_grid, phi_grid, ioncharge, dist_i_GK[n], mu_i[n], U_i[n], size_mu_i[n], size_U_i[n], grid_parameter, flux_i+n, Q_i+n, zoomfactor, -1.0, -999.9, vy_i_wall[n], mu_i_op[n], chiM_i[n], twopidmudvy_i[n], size_op_i+n);
			do
				size_sumnigrid = size_ngrid[n] ;
			while (size_sumnigrid > size_ngrid[n]) ;
			sumflux_i += (nioverne[n]*sqrt(TioverTe[n]/mioverme[n])*flux_i[n]);
			sumQ_i += (nioverne[n]*sqrt(TioverTe[n]/mioverme[n])*Q_i[n]);
			for (i=0; i<size_sumnigrid; i++)
				sumni_grid[i] += (nioverne[n]*ni_grid[n][i]);
		}
		if (gamma_ref > SMALLGAMMA) {
			densfinorb(1.0, 1.0, alpha, size_phiDSgrid, &size_neDSgrid, ne_DSgrid, x_DSgrid, phi_DSgrid, -1.0, dist_e_GK, mu_e, U_e_DS, size_mu_e, size_vpar_e, 0.0, &flux_eDS, &garbage, ZOOM_DS, 0.0, -999.9, vy_e_wall, mu_e_op, chiM_e, twopidmudvy_e, &size_op_e); 
			//printf("chiM mu\n");
			for (i=0; i<size_op_e; i++) {
				vpar_e_cut_lookup[i] = sqrt(2.0*(chiM_e[i] - mu_e_op[i]));
				//printf("%f %f  %f\n", chiM_e[i], mu_e[i], vpar_e_cut_lookup[i]);
			} 
			for (i=0; i<size_mu_e; i++) {
				vpar_e_cut[i] = lin_interp(mu_e_op, vpar_e_cut_lookup, mu_e[i], size_op_e, -1);
				//vpar_e_cut[i] = sqrt(2.0*(chiM_e[i] - mu_e[i]));
				//printf("%f %f %f\n", chiM_e[i], mu_e[i], vpar_e_cut[i]);
			} 
			printf("flux_eDS = %f\n", flux_eDS);
		}
		else {
			for (i=0; i< size_mu_e; i++) 
				vpar_e_cut[i] = 0.0; //sqrt(-2.0*lin_interp(x_DSgrid, phi_DSgrid, sqrt(mue_cut_lookup[i]), size_phiDSgrid, 1234)); 
			denszeroorb(-1.0, 1.0, phi_DSgrid, ne_DSgrid, size_phiDSgrid, &flux_eDS, &Q_eDS, dist_e_GK, vpar_e_DS, mu_e, size_vpar_e, size_mu_e, vpar_e_cut, gamma_ref, x_DSgrid, &ne_inf);
			//for (i=0; i< size_mu_e; i++) 
			//	vpar_e_cut[i] = sqrt(-2.0*lin_interp(x_DSgrid, phi_DSgrid, sqrt(2.0*mu_e[i]), size_phiDSgrid, 1234)); 
			for (i=0; i<size_mu_e; i++) {
				//vpar_e_cut[i] = sqrt(2.0*(-phi_DSgrid[0]*exp(((phi_DSgrid[1]-phi_DSgrid[0])/(phi_DSgrid[0]*deltaxDS))*sqrt(2.0*mu_e[i]))));
				//if (-phi_DSgrid[0] - 0.54*gamma_ref*sqrt(mu_e[i])*pow(-phi_DSgrid[0],0.7) < 0.0)
					//vpar_e_cut[i] = sqrt(2.0*(-phi_DSgrid[0] - 0.54*gamma_ref*sqrt(mu_e[i])*pow(-phi_DSgrid[0],0.7)));
				//else vpar_e_cut[i] = 0.0;
				vpar_e_cut[i] = sqrt(2.0*(-phi_DSgrid[0] + 0.5*EW*EW)*exp(-EW*sqrt(2.0*mu_e[i])/(-phi_DSgrid[0]))*gsl_sf_bessel_I0(EW*sqrt(2.0*mu_e[i])/(-phi_DSgrid[0])));
			}
			
			i=0;
			while (ne_DSgrid[i] < 1.0 - MARGIN_DS) i++;
			size_neDSgrid = i;
			printf("size_neDSgrid = %d\n", size_neDSgrid);
		}
		error_Poisson(error_MP, x_grid, ne_grid, sumni_grid, nioverne,  phi_grid, size_phigrid, size_sumnigrid, 0.0);
		printf("error_av = %f\terror_max = %f\n", error_MP[0], error_MP[1]);
		//if ( error_MP[1] > tol_MP[1] ) {
		//	printf("ERROR: MP solution rejected because it does not satisfy Poisson's equation accurately enough on the extended domain\n");
		//	fprintf(fout, "ERROR: MP solution rejected because it does not satisfy Poisson's equation accurately enough on the extended domain\n");
		//	exit(-1);
		//}
		error_Poisson(error_DS, x_DSgrid, ne_DSgrid, sumni_DSgrid, nioverne, phi_DSgrid, size_phiDSgrid, size_neDSgrid, 1.0/(gamma_ref*gamma_ref));
		printf("error_av = %f\terror_max = %f\n", error_DS[0], error_DS[1]);
		i=1;
		if ( error_DS[1] > tol_DS[1]) {
			printf("WARNING: DS solution should be rejected because it does not satisfy Poisson's equation accurately enough on the extended domain\n");
			fprintf(fout, "WARNING: DS solution should be rejected because it does not satisfy Poisson's equation accurately enough on the extended domain\n");
			exit(-1);
		}
		//printf("FINAL CHECK passed. HURRAY!\n");
		//Bohmeval(&Bohm, alpha, TioverTe[0], phi_grid[0], dist_i_GK[0], mu_i[0], U_i[0], vy_i_wall[0], mu_i_op[0], chiM_i[0], twopidmudvy_i[0], size_mu_i[0], size_U_i[0], size_op_i[0]); 
		evalBohmshouldbe(&Bohmshouldbe, phi_grid[0], dist_e_DK, vpar_e, mu_e, size_vpar_e, size_mu_e, vpar_e_cut);
		printf("?Bohm integral should be %f\n", Bohmshouldbe);
		fprintf(fout, "Bohm integral should be %f\n", Bohmshouldbe);
		fprintf(fBohm, "%f %f\n", Bohm, Bohmshouldbe);
		printf("electron heat flux Q_e = %f\n", Q_e);
		printf("electron particle flux Phi = %f\n", flux_e);
		printf("electron sheath heat transmission coefficient??? = %f\n", Q_e/flux_e);
	// Print potential and density profiles to files
	}

	if (N== MAX_IT) {
		printf("No convergence after %d iterations :(\n", MAX_IT);
		fprintf(fout, "No convergence after %d iterations :(\n", MAX_IT);
	}
	else {
		clock_t end_itDS = clock(); // finds end time of last iteration
		tot_time = (double) (end_itDS - begin_it) / CLOCKS_PER_SEC;
		printf("At %dth iteration MP+DS combined iteration converged successfully in %f seconds\n", N, tot_time);
		fprintf(fout, "At %dth iteration MP+DS combined iteration converged successfully in %f seconds\n", N, tot_time);
		printf("\t(phi_DSE, phi_wall) = (%f, %f)\n\tcurrent = %f\n", phi_grid[0], -0.5*v_cut*v_cut, current);
		fprintf(fout, "\t(phi_DSE, phi_wall) = (%f, %f)\n\t current = %f\n", phi_grid[0], -0.5*v_cut*v_cut, current);
	}
	fclose(fout);

	// Print potential and density profiles to files
	snprintf(fpstr, 150, "%s/phi_n_DS.txt", dirname);
	//fp = fopen("OUTPUT/phi_n_DS.txt", "w");
	fp = fopen(fpstr, "w");
	if (fp == NULL)  
		printf("error when opening file %s\n", fpstr);

	if (N_DS != 0) {
		for (i=0; i<size_phiDSgrid; i++) {
			fprintf(fp, "%f %f %f %f\n", x_DSgrid[i], phi_DSgrid[i], sumni_DSgrid[i], ne_DSgrid[i]);
		}
	}
	else 
		fprintf(fp, "%f\n", -0.5*v_cutDS*v_cutDS);
	fclose(fp);
	snprintf(fpstr, 150, "%s/phi_n_MP.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file\n");
	}
	for (i=0; i<size_phigrid; i++) {
		fprintf(fp, "%f %f %f %f\n", x_grid[i], phi_grid[i], sumni_grid[i], ne_grid[i]);
	}
	fclose(fp);
	snprintf(fpstr, 150, "%s/vparcut.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file\n");
	}
	for (i=0; i<size_mu_e; i++) {
		fprintf(fp, "%f %f\n", mu_e[i], vpar_e_cut[i]);
	}
	fclose(fp);
	snprintf(fpstr, 150, "%s/Fi_W.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file\n");
	}
	for (i=0; i<size_op_i[0]; i++) {
		fprintf(fp, "%f %f %f %f\n", vy_i_wall[0][i], mu_i_op[0][i], chiM_i[0][i], twopidmudvy_i[0][i]);
	}
	fclose(fp);
	snprintf(fpstr, 150, "%s/misc_output.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file\n");
	}
	fprintf(fp, "%f\n", sumflux_i-flux_e);
	fprintf(fp, "%f\n", 0.5*v_cut*v_cut);
	fprintf(fp, "%f\n", Q_e);
	fprintf(fp, "%f\n", sumQ_i);
	fprintf(fp, "%f\n", flux_e);
	fprintf(fp, "%f\n", sumflux_i);
	fclose(fp);
	snprintf(fpstr, 150, "%s/number_iterations.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file %s\n", dirname);
	}
	fprintf(fp, "%d %d\n", N, N_DS);
	fclose(fp);
	// save also input files in the output folder
	snprintf(fpstr, 150, "%s/Fi_mpe_args.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file\n");
	}
	for (i=0; i<size_mu_i[0]; i++) {
		fprintf(fp, "%f ", mu_i[0][i]);
	}
	fprintf(fp, "\n");
	for (j=0; j<size_U_i[0]; j++) {
		fprintf(fp, "%f ", U_i[0][j]);
	}
	fclose(fp);
	snprintf(fpstr, 150, "%s/Fi_mpe.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file\n");
	}
	for (i=0; i<size_mu_i[0]; i++) {
		for (j=0; j<size_U_i[0]; j++) {
			fprintf(fp, "%f ", dist_i_GK[0][i][j]);
			if (j == size_U_i[0] - 1)	
				fprintf(fp, "\n");
		}
	}
	fclose(fp);
	snprintf(fpstr, 150, "%s/Fe_mpe_args.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file\n");
	}
	for (i=0; i<size_mu_e; i++) {
		fprintf(fp, "%f ", mu_e[i]);
	}
	fprintf(fp, "\n");
	for (j=0; j<size_vpar_e; j++) {
		fprintf(fp, "%f ", vpar_e[j]);
	}
	fclose(fp);
	snprintf(fpstr, 150, "%s/Fe_mpe.txt", dirname);
	fp = fopen(fpstr, "w");
	if (fp == NULL) {
		printf("error when opening file\n");
	}
	for (i=0; i<size_mu_e; i++) {
		for (j=0; j<size_vpar_e; j++) {
			fprintf(fp, "%f ", dist_e_DK[i][j]);
			if (j == size_vpar_e - 1)	
				fprintf(fp, "\n");
		}
	}
	fclose(fp);
	fclose(fBohm); //??????
	free(Q_i);
	free(vpar_e_cut);
	free(x_grid);
	free(phi_grid);
	free(ne_grid);
	free(ni_grid);
	free(x_DSgrid);
	free(phi_DSgrid);
	free(ne_DSgrid);
	free(ni_DSgrid);
	free(U_e_DS);
	free(vpar_e_DS);
	exit(0);
}
