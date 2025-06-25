// Authors: Aessandro Geraldini and Robbie Ewart
/* This script calculates the density profiles of charged particles in the magnetised sheath for a given potential profile and entrance distribution function: denszeroorb neglects gyro-orbit size; densfinorb includes distorted gyro-orbits for magnetic field angle α<<1. 
 * */
/* Authors: 
   Robbie Ewart 
wrote the function denszeroorb (originally makelookup) which calculated the electron density from 1D distribution of parallel velocities
   Alessandro Geraldini 
wrote the function densfinorb and upgraded denszeroorb by including integration over mu and tweaking to allow calculation of density of particles accelerated to wall

MODIFIED on 15 JUL 2022 by Alessandro Geraldini
*/

#define TESTELL 0
#define APPROXMUFORSMALLORBIT 0

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "mps.h"
#define SMALLGAMMA 0.21
#define numb 0.00000001

double tophat(double x1, double x2, double x) {
	double y;
	if ((x1 <= x ) && (x <= x2 )) 
		y = 1.0;
	else
		y = 0.0;
	return y;
}



double smallrhoefrac(double x, double mu) {
	double frac;
	// calculates the fraction of ions present for a given magnetic moment mu at a position x
	// gyrophase angles corresponding to gyro-orbits touching the wall are excluded
	if (x== 0.0) 
		frac = 0.0;
	else if (x<2.0*sqrt(2.0*mu))
		frac = 1.0 - (1.0/M_PI)*acos(x/(sqrt(2.0*mu)) - 1.0) ;
	else
		frac = 1.0;
	return frac;
}

// denszeroorb (previously makelookup) calculates the density integral of a species when this is only a function of the potential at a given point 

/* AG
   to use for ions in the Debye sheath, phi must be a list of POSITIVE numbers 
   equal to minus the potential relative to the Debye sheath entrance
   set number pointed at by vpar_cut_lookup to be equal to 0.0
   set size_mu = 0 while size_cut, size_vpar and mue_cut_lookup are redundant (can set to zero/NULL)
*/

void denszeroorb(double charge, double TeovTs, double *phi_real, double *n_grid, int p_size, double *Phi_point, double *Qe_point, double **distfunc, double *vpar, double *mu, int size_vpar, int size_mu, double *vpar_cut_lookup, double gamma, double *x_grid, double *n_infp) {
	//define variables
	int count = 0;
	int vi, vi_1, p, len_F; //vi is a counting variable that will be saved for sums over velocity space, vi_1 is a special value in velocity space devoted to the first point, p is a counting variable that will be saved for counting over phi space, len_F saves the number of entries in the distribution;
	double *phi, *vparacc, phip = -0.000000001 ;
	double Phi=0.0, Qe=0.0, n_inf = 0.0;
	double v_max, v_s, v_min;//v_max is the maximum velocity the distribution function will go up to before effectively just reading 0 from there on, v_s is the separation in velocity space between ajacent points. v_cut is the velocity that would be require in order to just reach the plasma wall boundary, v_min is the minimum velocity (the entrance to the presheath) required to reach a given x
	double sqrt_up, sqrt_lo, theta_up, theta_lo, sqrt_cut, theta_cut; //sqrt_up/lo reprresent the upper and lower limits of one strip integral, theta_up/lo are related to the hyperbolic arcsinh of some specific values (see document that hopefully exists)

	double *F, *Fp, *Fpp, **ddistdvpar, **ddistdvpartwo; // the zeroth first and second derivatives of the distribution function respectively

	double* in_err, *af_err, * n_res, w, w_up, w_lo; // for calculating errors in different parts of the code

	double* n_pre; //n_pre[x] is an array that will contain the values that we expect the integral to give (for a maxwellian input)

	int show_err;// change show_err to 1 or 0 if you do/don't want the error output to be printed to a file (note that the error output will only be meaningful for maxwellian inputs). Change p_option to control how the phi values are distributed with 1 being in an attempt to make n_grid roughly linear and 0 making phi linear;
	int mu_ind;
	double *nepart, vpar_cut;
	FILE* fptr1;
	show_err = 0;

	vparacc = malloc(size_vpar*sizeof(double));
	phi = malloc(p_size*sizeof(double));
	for (p=0; p<p_size; p++) {
		phi[p] = -phi_real[p]*charge*TeovTs;
	}


	ddistdvpar = malloc(size_mu * sizeof(double));
	ddistdvpartwo = malloc(size_mu * sizeof(double));
	in_err = malloc(p_size * sizeof(double));
	af_err = malloc(p_size * sizeof(double));
	n_res = malloc(p_size * sizeof(double));
	n_pre = malloc(p_size * sizeof(double));
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
		for (vi = 0; vi < len_F; vi++)
		{
			if ((vi < len_F - 1) && (vi > 0))
			{
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
	}

	/* comment by AG
	the integration has two options 
	if size_mu is zero then the distribution function is 1D 
	e.g. rho_e << lambda_D model of electrons in magnetic presheath x ~ rho_i (velocity coordinate = v_parallel) or ions in the Debye sheath x ~ lambda_D << rho_i (velocity coordinate = v_x)
	if size_mu is non-zero then the distribution function is 2D but, differently from the integration in densfinorb.c, is still a local function of x
	e.g. rho_e >~ lambda_D for electrons in the magnetic presheath x ~ rho_i
	*/
	if (size_mu > 1) {
		printf("In 2D integration\n");
		nepart = malloc(size_mu*sizeof(double));
		// for finite electron gyroradius effects, we need to integrate in mu as well
		// cutoff function for large rhoe

		// Normalization at infinity
		n_inf = 0.0;
		v_min = 0.0;
		for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
			count = 0;
			nepart[mu_ind] = 0.0;
			vpar_cut = vpar_cut_lookup[mu_ind]; 
			//lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
			//printf("in makelookup.c: mu = %f and vpar_cut_DSE = %f\n", mu[mu_ind], vpar_cut);
			vpar_cut = sqrt(vpar_cut*vpar_cut + 2.0*(-phi[0]) + TINY);
			for (vi = 0; vi < len_F - 1; vi++) {
				F[vi]     = distfunc[mu_ind][vi];
				Fp[vi]    = ddistdvpar[mu_ind][vi];
				Fpp[vi]   = ddistdvpartwo[mu_ind][vi];
			}
			if (vpar_cut >= v_max) // effectively no cut off
			{
				for (vi = 0; vi < len_F - 2; vi++)
					nepart[mu_ind] += (F[vi] + F[vi + 1]) * v_s;
				count++;
			}
			else { // now there is a cut-off
				vi_1 = 0;
				sqrt_up = ((vi_1 + 1) * v_s);
				nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
				nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));

				w = 0.5 * pow(v_s, 2.0);

				for (vi = vi_1 + 1; vi < len_F - 1; vi++) {
					sqrt_up = (vi + 1) * v_s;
					sqrt_lo = (vi * v_s);

					nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
					nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

					w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
					w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
				}

				for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++) {
					sqrt_up = (vi + 1) * v_s;
					sqrt_lo = (vi * v_s);

					nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
					nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

					w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
					w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
				}
				vi_1 = (int)floor(vpar_cut / v_s);
				sqrt_up = vpar_cut;
				sqrt_lo = vi_1 * v_s;

				nepart[mu_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
				nepart[mu_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));

				w_up = 0.5 * pow(vpar_cut, 2.0);
				w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
			}
			//printf("in makelookup.c: mu = %f and vpar_cut_DSE = %f, nepart[mu_ind] = %f\n", mu[mu_ind], vpar_cut, nepart[mu_ind]);
			if (mu_ind != 0) {
				n_inf += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1])*(mu[mu_ind] - mu[mu_ind-1]);
			}
			//printf("count=%d\n", count);
		}
		printf("n_inf = %f\n", n_inf);
		
		n_inf = 0.0;
		v_min = sqrt(-2.0 * phip);
		for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
			nepart[mu_ind] = 0.0;
			//vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
			vpar_cut = vpar_cut_lookup[mu_ind] ; 
			//vpar_cut = sqrt(vpar_cut*vpar_cut + TINY);
			vpar_cut = sqrt(vpar_cut*vpar_cut + 2.0*(-phi[0]) + TINY);
			//printf("vpar_cut = %f\n", vpar_cut);
			for (vi = 0; vi < len_F - 1; vi++) {
				//F[vi]   = exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
				//Fp[vi]   = -vi*v_s*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
				//Fpp[vi]   = (vi*v_s*vi*v_s - 1.0)*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
				F[vi]     = distfunc[mu_ind][vi];
				Fp[vi]    = ddistdvpar[mu_ind][vi];
				Fpp[vi]   = ddistdvpartwo[mu_ind][vi];
			}
			if (vpar_cut >= v_max) // effectively no cut off
			{
				if ((int)floor(v_min / v_s) >= len_F - 1) {
					//ne[p] = 0;
					nepart[mu_ind] = 0.0;
					printf("v_min outside velocity grid");
				}
				else {
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
						sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phip));
						theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phip)) / v_min);

						nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
						nepart[mu_ind] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));

						w = ((0.5 * pow((((vi_1 + 1) * v_s) - v_min), 2.0)) + (v_min * (((vi_1 + 1) * v_s) - v_min)));

						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
						{
							sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip);
							sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip);
							theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip) / (v_min));
							theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip ) / (v_min));

							nepart[mu_ind] += 2.0 * (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							nepart[mu_ind] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
							w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
						}
					}
				}
			//n_res[p] = ne[p] - n_pre[p];
			}
			else { // now there is a cut-off
				v_min = sqrt(-2.0 * phip);
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

						w = 0.5 * pow(v_s, 2.0);

						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
						{
							sqrt_up = (vi + 1) * v_s;
							sqrt_lo = (vi * v_s);

							nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
							nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
							w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
						}

						for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
						{
							sqrt_up = (vi + 1) * v_s;
							sqrt_lo = (vi * v_s);

							nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
							nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
							w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
						}
						vi_1 = (int)floor(vpar_cut / v_s);
						sqrt_up = vpar_cut;
						sqrt_lo = vi_1 * v_s;

						nepart[mu_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
						nepart[mu_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));

						w_up = 0.5 * pow(vpar_cut, 2.0);
						w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
					}
					else
					{
						if ((int)floor(v_min / v_s) == (int)floor(vpar_cut / v_s))
						{
							vi_1 = (int)floor(v_min / v_s);
							sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phip));
							theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phip)) / (v_min));

							nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
							nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

							w_up = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));

							sqrt_cut = sqrt(pow(vpar_cut, 2.0) + 2.0 * phip);
							theta_cut = asinh(sqrt(pow(vpar_cut, 2.0) + 2.0 * phip) / (v_min));

							nepart[mu_ind] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * vpar_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
							nepart[mu_ind] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (vpar_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));

							w_lo = 0.5 * (pow(vpar_cut, 2.0) - pow(v_min, 2.0));

							for (vi = vi_1 + 1; vi < len_F - 1; vi++)
							{
								sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip);
								sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip);
								theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip) / (v_min));
								theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip) / (v_min));

								nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
								nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

								w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
								w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
							}
						}
						else //the most likely case
						{
							vi_1 = (int) (floor(v_min / v_s));
							sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phip);
							theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phip) / (v_min));
							nepart[mu_ind] = 0.0;
							

							nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
							nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

							w = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));

							for (vi = vi_1 + 1; vi < len_F - 1; vi++)
							{
								sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip);
								sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip);
								theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip) / (v_min));
								theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip) / (v_min));

								nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
								nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

								w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
								w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
							}

							for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
							{
								sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip);
								sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip);
								theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip) / (v_min));
								theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip) / (v_min));

								nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
								nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

								w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
								w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
							}

							vi_1 = (int)floor(vpar_cut / v_s);
							sqrt_up = sqrt(pow(vpar_cut, 2.0) + 2.0 * phip);
							sqrt_lo = sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phip);
							theta_up = asinh(sqrt(pow(vpar_cut, 2.0) + 2.0 * phip) / (v_min));
							theta_lo = asinh(sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phip) / (v_min));

							nepart[mu_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							nepart[mu_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));

							w_up = 0.5 * (pow(vpar_cut, 2.0) - pow(v_min, 2.0));
							w_lo = 0.5 * (pow((vi_1 * v_s), 2.0) - pow(v_min, 2.0));
						}
					}
				//n_res[p] = ne[p] - n_pre[p];
				}
			}
			if (mu_ind != 0) {
				n_inf += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1])*(mu[mu_ind] - mu[mu_ind-1]);
				//printf("n_grid %f\n", n_grid[p]);
			}
		}
	//n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] -phi[0] )));
	//printf("n_grid[%d] = %f, n_pre = %f\n", p, n_grid[p], n_pre[p]);
	//}

		// ELECTRON FLUX electron flux
		for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
			nepart[mu_ind] = 0.0;
			//vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
			vpar_cut = vpar_cut_lookup[mu_ind] ; 
			//vpar_cut = sqrt(vpar_cut*vpar_cut + TINY);
			vpar_cut = sqrt(vpar_cut*vpar_cut + 2.0*(-phi[0]) + TINY);
			for (vi = 0; vi < len_F - 1; vi++) {
				//F[vi]   = vi*v_s*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
				//Fp[vi]  = (1.0 - vi*v_s*vi*v_s)*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
				//Fpp[vi] = (-3.0*vi*v_s + vi*v_s*vi*v_s*vi*v_s )*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
				F[vi]   = vi*v_s*distfunc[mu_ind][vi];
				Fp[vi]  = vi*v_s*ddistdvpar[mu_ind][vi] + distfunc[mu_ind][vi];
				Fpp[vi] = vi*v_s*ddistdvpartwo[mu_ind][vi] + 2.0*ddistdvpar[mu_ind][vi];
				//printf("F[%d]=%f (should be %f)\n", vi, F[vi], vi*v_s*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s));
			}
			//printf("v_max = %f\tvpar_cut = %f\n", v_max, vpar_cut);
			if (vpar_cut >= v_max) // effectively no cut off
				Phi = 0.0;
			else { // now there is a cut-off
				v_min = 0.0;
				if ((int)floor(v_min / v_s) >= len_F - 1)
				{
					nepart[mu_ind] = 0.0;
				}
				else
				{
					vi_1 = -1;
					for (vi = vi_1 + 1; vi < len_F - 1; vi++)
					{
						sqrt_up = (vi + 1) * v_s;
						sqrt_lo = (vi * v_s);

						nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
						nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
						w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
					}

					for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
					{
						sqrt_up = (vi + 1) * v_s;
						sqrt_lo = (vi * v_s);

						nepart[mu_ind] -= (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
						nepart[mu_ind] -= ( ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))) );

						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
						w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
					}
					vi_1 = (int)floor(vpar_cut / v_s);
					sqrt_up = vpar_cut;
					sqrt_lo = vi_1 * v_s;

					nepart[mu_ind] -= (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
					nepart[mu_ind] -= ( ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))) );

					w_up = 0.5 * pow(vpar_cut, 2.0);
					w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
				}
			}
			if (mu_ind != 0) 
				Phi += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1 ])*(mu[mu_ind]-mu[mu_ind-1]);
		}

		// ELECTRON HEAT FLUX
		for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
			nepart[mu_ind] = 0.0;
			//vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
			vpar_cut = vpar_cut_lookup[mu_ind] ; 
			//vpar_cut = sqrt(vpar_cut*vpar_cut + TINY);
			vpar_cut = sqrt(vpar_cut*vpar_cut + 2.0*(-phi[0]) + TINY);
			for (vi = 0; vi < len_F - 1; vi++) {
				F[vi]   = (0.5*vi*v_s*vi*v_s + mu[mu_ind])*vi*v_s*distfunc[mu_ind][vi];
				Fp[vi]  = (0.5*vi*v_s*vi*v_s + mu[mu_ind])*(vi*v_s*ddistdvpar[mu_ind][vi] + distfunc[mu_ind][vi]) + vi*v_s*vi*v_s*distfunc[mu_ind][vi];
				Fpp[vi] = (0.5*vi*v_s*vi*v_s + mu[mu_ind])*(vi*v_s*ddistdvpartwo[mu_ind][vi] + 2.0*ddistdvpar[mu_ind][vi]) + (vi*v_s)*(vi*v_s*ddistdvpar[mu_ind][vi] + distfunc[mu_ind][vi]) + 2.0*vi*v_s*distfunc[mu_ind][vi] + vi*v_s*vi*v_s*ddistdvpar[mu_ind][vi];
				//F[vi]   = (0.5*vi*v_s*vi*v_s)*vi*v_s*distfunc[mu_ind][vi];
				////Fp[vi]  = 1.5*vi*v_s*vi*v_s*distfunc[mu_ind][vi] + 0.5*vi*v_s*vi*v_s*vi*v_s*ddistdvpar[mu_ind][vi];
				//Fp[vi]  = (0.5*vi*v_s*vi*v_s)*(vi*v_s*ddistdvpar[mu_ind][vi] + distfunc[mu_ind][vi]) + vi*v_s*vi*v_s*distfunc[mu_ind][vi];
				////Fpp[vi] = 3.0*vi*v_s*distfunc[mu_ind][vi] + 3.0*vi*v_s*vi*v_s*ddistdvpar[mu_ind][vi] + 0.5*vi*v_s*vi*v_s*vi*v_s*ddistdvpartwo[mu_ind][vi]; 
				//Fpp[vi] = (0.5*vi*v_s*vi*v_s)*(vi*v_s*ddistdvpartwo[mu_ind][vi] + 2.0*ddistdvpar[mu_ind][vi]) + (vi*v_s)*(vi*v_s*ddistdvpar[mu_ind][vi] + distfunc[mu_ind][vi]) + 2.0*vi*v_s*distfunc[mu_ind][vi] + vi*v_s*vi*v_s*ddistdvpar[mu_ind][vi];
			}
			//printf("v_max = %f\tvpar_cut = %f\n", v_max, vpar_cut);
			if (vpar_cut >= v_max) // effectively no cut off
				Qe = 0.0;
			else { // now there is a cut-off
				v_min = 0.0;
				if ((int)floor(v_min / v_s) >= len_F - 1)
				{
					nepart[mu_ind] = 0.0;
				}
				else
				{
					vi_1 = -1;
					for (vi = vi_1 + 1; vi < len_F - 1; vi++)
					{
						sqrt_up = (vi + 1) * v_s;
						sqrt_lo = (vi * v_s);

						nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
						nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
						w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
					}

					for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
					{
						sqrt_up = (vi + 1) * v_s;
						sqrt_lo = (vi * v_s);

						nepart[mu_ind] -= (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
						nepart[mu_ind] -= ( ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))) );

						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
						w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
					}
					vi_1 = (int)floor(vpar_cut / v_s);
					sqrt_up = vpar_cut;
					sqrt_lo = vi_1 * v_s;

					nepart[mu_ind] -= (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
					nepart[mu_ind] -= ( ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))) );

					w_up = 0.5 * pow(vpar_cut, 2.0);
					w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
				}
			}
			if (mu_ind != 0) 
				Qe += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1 ])*(mu[mu_ind]-mu[mu_ind-1]);
		}

		// NOW CALCULATE THE DENSITY PROFILE
		for (p = 0; p < p_size; p++)
		{
			//printf("phi[%d] = %f\n", p, phi[p]);
			if (fabs(phi[p]) < TINY) {
				n_grid[p] = n_inf;
			}
			else {
				n_grid[p] = 0.0;
				v_min = sqrt(-2.0 * phi[p]);
				for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
					nepart[mu_ind] = 0.0;
					//vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
					vpar_cut = vpar_cut_lookup[mu_ind] ; 
					vpar_cut = sqrt(vpar_cut*vpar_cut +  (-2.0*phi[0]) + TINY);
					//printf("vpar_cut = %f\n", vpar_cut);
					for (vi = 0; vi < len_F - 1; vi++) {
						//F[vi]   = exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
						//Fp[vi]   = -vi*v_s*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
						//Fpp[vi]   = (vi*v_s*vi*v_s - 1.0)*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
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
								sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p]));
								theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p])) / v_min);

								nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
								nepart[mu_ind] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));

								w = ((0.5 * pow((((vi_1 + 1) * v_s) - v_min), 2.0)) + (v_min * (((vi_1 + 1) * v_s) - v_min)));
								in_err[p] += 2.0 * ((((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w))));
								in_err[p] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));

								for (vi = vi_1 + 1; vi < len_F - 1; vi++)
								{
									sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
									sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
									theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
									theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

									nepart[mu_ind] += 2.0 * (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
									nepart[mu_ind] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

									w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
									w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
									af_err[p] += 2.0 * ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
									af_err[p] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
								}
							}
						}
						//n_res[p] = ne[p] - n_pre[p];
					}
					else { // now there is a cut-off
						v_min = sqrt(-2.0 * phi[p]);
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

								w = 0.5 * pow(v_s, 2.0);
								af_err[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s)))) - erf(sqrt(w));
								af_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));

								for (vi = vi_1 + 1; vi < len_F - 1; vi++)
								{
									sqrt_up = (vi + 1) * v_s;
									sqrt_lo = (vi * v_s);

									nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
									nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

									w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
									w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
									af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
									af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
								}

								for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
								{
									sqrt_up = (vi + 1) * v_s;
									sqrt_lo = (vi * v_s);

									nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
									nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

									w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
									w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
									af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
									af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
								}
								vi_1 = (int)floor(vpar_cut / v_s);
								sqrt_up = vpar_cut;
								sqrt_lo = vi_1 * v_s;

								nepart[mu_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
								nepart[mu_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));

								w_up = 0.5 * pow(vpar_cut, 2.0);
								w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
								af_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
								af_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
							}
							else
							{
								if ((int)floor(v_min / v_s) == (int)floor(vpar_cut / v_s))
								{
									vi_1 = (int)floor(v_min / v_s);
									sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p]));
									theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p])) / (v_min));

									nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
									nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

									w_up = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
									in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w_up))));
									in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

									sqrt_cut = sqrt(pow(vpar_cut, 2.0) + 2.0 * phi[p]);
									theta_cut = asinh(sqrt(pow(vpar_cut, 2.0) + 2.0 * phi[p]) / (v_min));

									nepart[mu_ind] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * vpar_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
									nepart[mu_ind] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (vpar_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));

									w_lo = 0.5 * (pow(vpar_cut, 2.0) - pow(v_min, 2.0));
									in_err[p] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * vpar_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
									in_err[p] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (vpar_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));

									for (vi = vi_1 + 1; vi < len_F - 1; vi++)
									{
										sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
										sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
										theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
										theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

										nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
										nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

										w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
										w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
										af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
										af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
									}
								}
								else //the most likely case
								{
									vi_1 = (int) (floor(v_min / v_s));
									sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phi[p]);
									theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
									nepart[mu_ind] = 0.0;
									

									nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
									nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

									w = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
									in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w))));
									in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

									for (vi = vi_1 + 1; vi < len_F - 1; vi++)
									{
										sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
										sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
										theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
										theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

										nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
										nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

										w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
										w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
										af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
										af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
									}

									for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
									{
										sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
										sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
										theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
										theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

										nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
										nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

										w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
										w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
										af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
										af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
									}

									vi_1 = (int)floor(vpar_cut / v_s);
									sqrt_up = sqrt(pow(vpar_cut, 2.0) + 2.0 * phi[p]);
									sqrt_lo = sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phi[p]);
									theta_up = asinh(sqrt(pow(vpar_cut, 2.0) + 2.0 * phi[p]) / (v_min));
									theta_lo = asinh(sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

									nepart[mu_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
									nepart[mu_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));

									w_up = 0.5 * (pow(vpar_cut, 2.0) - pow(v_min, 2.0));
									w_lo = 0.5 * (pow((vi_1 * v_s), 2.0) - pow(v_min, 2.0));
									in_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
									in_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));
								}
							}
						//n_res[p] = ne[p] - n_pre[p];
						}
					}
					//if (nepart[mu_ind] != nepart[mu_ind]) printf("mu_ind = %d, p = %d\n", mu_ind, p);

//piece below is attempt at including first order in gamme effect of density variation --> does not work so well, so leave commented out
					if ( (gamma < SMALLGAMMA) && (gamma > TINY) ) nepart[mu_ind] *= smallrhoefrac(x_grid[p], mu[mu_ind]);
//instead, use wall electric field to calculate cutoff correction in vparallel --> works better
					if (mu_ind != 0) {
						n_grid[p] += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1])*(mu[mu_ind] - mu[mu_ind-1]);
						//printf("n_grid %f\n", n_grid[p]);
					}
				}
			}
			//n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] -phi[0] )));
			//printf("n_grid[%d] = %f, n_pre = %f\n", p, n_grid[p], n_pre[p]);
			if (p== p_size-1)
				printf("in denszeroorb: charge = %f, derivative wrt phi is dndphi = %f\n", charge, (n_grid[p] - n_grid[p-1])/(phi[p] - phi[p-1]));
		}
		free(nepart);
	}
	else {
		printf("In 1D ion integration\n");
		n_inf = 0.0;
		for (vi = 0; vi < len_F - 1; vi++) {
			F[vi] = distfunc[0][vi];
			if (vi != 0)
				n_inf += 0.5*(F[vi] + F[vi-1])*(vpar[vi] - vpar[vi-1]);
		}
		for (p=0; p<p_size;p++) {
			n_grid[p] = 0.0;
			//if (phi[p] < 0.0) n_grid[p] = 1.0;
			//else {
				for (vi = 0; vi < len_F - 1; vi++) {
					vparacc[vi] = sqrt(vpar[vi]*vpar[vi] + 2.0*phi[p]);
					if (vi != 0)
						n_grid[p] += 0.5*(F[vi] + F[vi-1])*(vparacc[vi] - vparacc[vi-1]);
				}
			//}
		}
	}

	Phi =  Phi / n_inf;
	Qe /= n_inf;
	printf("Phi = %f\n", Phi);
	printf("n_inf = %f\n", n_inf);
	*n_infp = n_inf;
	for (p = 0;p < p_size; p++)
	{
		//n_grid[p] = n_grid[p] / n_grid[p_size - 1];
		n_grid[p] = n_grid[p] / n_inf;
		//if (n_grid[p] != n_grid[p]) n_grid[p] = 1.0;
		n_res[p] = n_grid[p] - n_pre[p];
		if ( (DEBUG == 1) && (x_grid[p] < 10.0) )
			printf("n_denszeroorb[%d/%d] = %f, x_grid = %f, phi = %f\n", p, p_size, n_grid[p], x_grid[p], phi[p]);
	}

	//only works if size_mu = 1 (or for the final value of mu)
	if (show_err == 1)
	{
		fptr1 = fopen("OUTPUT/total_error.txt", "w");
		for (p = 0; p < p_size; p++)
		{
			fprintf(fptr1, "%.15f\n", n_res[p]);
		}
		fclose(fptr1);
		fptr1 = fopen("OUTPUT/regular_point_error.txt", "w");
		for (p = 0; p < p_size; p++)
		{
			fprintf(fptr1, "%.15f\n", af_err[p]);
		}
		fclose(fptr1);
		fptr1 = fopen("OUTPUT/initial_point_error.txt", "w");
		for (p = 0; p < p_size; p++)
		{
			fprintf(fptr1, "%.15f\n", in_err[p]);
		}
		fclose(fptr1);
		//fptr1 = fopen("phi.txt", "w");
		//for (p = 0; p < p_size; p++)
		//{
		//	fprintf(fptr1, "%.15f\n", phi[p]);
		//}
		//fclose(fptr1);
		fptr1 = fopen("OUTPUT/Edens.txt", "w");
		for (p = 0; p < p_size; p++)
		{
			fprintf(fptr1, "%.15f\n", in_err[p]);
		}
		fclose(fptr1);
	}

	for (mu_ind = 0; mu_ind < size_mu; mu_ind++) {
		free(ddistdvpar[mu_ind]);
		free(ddistdvpartwo[mu_ind]);
	}
	free(ddistdvpar);
	free(ddistdvpartwo);
	free(vparacc);
	free(F);
	free(Fp);
	free(Fpp);
	free(in_err);
	free(af_err);
	free(n_res);
	free(n_pre);
	free(phi);

	*Phi_point = Phi;
	*Qe_point = Qe;

	printf("end of denszeroorb\n");
	return;
}


void densfinorb(double Ti, double lenfactor, double alpha, int size_phigrid, int *size_ngrid, double* n_grid, double *x_grid, double* phi_grid, double charge, double **FF, double *mumu, double *UU, int sizemumu, int sizeUU, double grid_parameter, double *flux, double *Qflux, int zoomfactor, double margin, double phi_DSbump, double *vy_op, double *mu_op, double *chiMax_op, double *dmudvy_op, int *size_op) {
	// declare variables
	clock_t begin = clock(); // Finds the start time of the computation
	double limit_rho = 8.0;
	double n_inf=0.0;
	double Ucrit = 0.0, frac_reflected=0.0, frac = 0.0, vzcrit;
	double deltax, deltax_inf, deltaE = 0.1, phibar;
	double *chiinf, *muinf, **vxinf;
	double *xx, *phi, *phip, *phipp, **chi;
	/* xx is the array containing the position (distance from the wall) on the fine grid; phi contains the values of phi(x), extracted from a file. phip is phi prime, first derivative of phi. n; newphi is the new electrostatic potential guess; chi is the effective potential;*/ 
	int s, w, ind, muapproxtype = 1;
	int stop = 0, sizexbar, maxj, j_inf;
	int reflected = 0;
	int icrit;
	/* n is the domain of the position x (the largest value of the index i, so e.g. x_grid[n] = L_2 in the paper); sizexbar is the domain of the position xbar (the largest value of the index j); */ 
	double *openorbit, openorbitantycal, **mu, *muopen, *chiMopen, *xbaropen, *openorbitopen, *Ucritf, *xbar, xbarcrit, chiMcrit, *xifunction;
	/* openorbit is the Delta_M = 2*pi*dmu/dxbar; openorbitantycal is the analytical value of openorbit for a flat potential; mu is the array containing values of mu(xbar, Uperp), index j for values of xbar, k for values of Uperp; xbar is the grid of values used in the closed orbit integral; FF contains the distribution function, read from the file distfile.txt. UU and mumu contain the values of U and mu corresponding to the function FF (which is F(mu, U)); FFprime is the numerical first derivative of F with respect to U; xifunction is the function of x which defines the grid of values of xbar by finding a chi whose minimum lies exactly at each grid point x. 
	Note: first index of FF and FFprime is mu, second one is U; */
	int i=0, ic=0, j=0, k=0, l=0;
	int *jmclosed, *jmopen,  sizeU, size_finegrid, size_xlim;
	/* jmclosed represent minimum values of xbar above which we integrate open and closed orbit density integrals respectively (xbar_m,o and xbar_m in the paper); i is an index usually representing the positin x  j is an index usually representing the orbit position xbar; k is an index usually representing the energy Uperp (or velocity vx). It's always used in conjunction with j (and sometimes i); l is an index (used in for loops) usually representing the total energy (or velocity vz). It's used only in the DENSITY INTEGRALS part of the code; sizeU is the size of the integration range over U (or velocity vz). It is set later on in the code */
	int *crossed_max, *crossed_min, *kdrop;
	/* crossed_max and crossed_min is non-zero when a minimum (or maximum) of chi is found for some xbar[j]; kdrop is an integer which is non-zero only if there is an additional potential barrier for the particle at x << rho e.g. for an ion if the sheath reverses Uperp is allowed to be above chiMax and therefore the k = 0 of Uperp[j][k] starts for Uperp[j][0] = chiMax[j] + phi_barrier instead of Uperp[j][0] = chiMax[j] which is the conventional way */
	int *lowerlimit, *upperlimit, **upper, *imax, *imin;
	/* lowerlimit represents the lower limit of k in the integrals over Uperp (or vx). It's needed because some of the earlies energies; (which are the largest because thy are values of chi stored after the maximum is found); may be so large that they are associated with very small values of the distribution function. This avoids integrating in an empty portion of phase space; upperlimit[j] represents the largest value of k (the smallest stored energy Uperp = chi_minimum) associated with some value of j; upper[j][i] represents the value of k associated with the smallest value of vx when integrating over Uperp. Going above upperlimit[j][i] makes Uperp < chi so velocities imaginary; imax/imin[j] stores the position of the maximum/minimum of the effective potential chi (It's x_M/x_m in the paper, which depends on xbar). */
	double **Uperp, ***vx, *chiMax, *chimpp, *chimin, oorbintgrd, oorbintgrdantycal;
	/* Uperp stores the possible values of Uperp associated with closed orbits, and so does vx; chiMax and chimin store the local maxima and minima of the effective potential maximum, oorbintgrd is the value of the integrand in the first open orbit integral (oorbintgrdantycal is the analytical result for flat potential) */
	double vz, U, dvz = 0.1, dvzopen = 0.1, dvx, dxbar, intdU=0.0, intdUopen=0.0;
		/* vz used in the density integral; U is the total energy, used in the density integral; dvz is the thickness of the vz grid used to take the integral over U (which is taken over vz in practice), dvzopen is the same for the open orbit piece; dvx is the thickness of the vx grid used to take the integral over Uperp ( which is taken over vx in practice). It must be evaluated because it depends on stored values of vx[j][i][k]; dxbar is the thickness of the xbar grid; intdU is the value of the integral over U in the closed orbit density integration process; intdUopen is the same as above, for the open orbit integral  */
	double intdUold=0.0, intdvx=0.0, intdvxold = 0.0, intdxbar=0.0, intdxbaropen=0.0, F, Fold=0.0, Fold_ref=0.0, Ucap;
		/* intdUold is a variable which stores the old intdU, so that the trapezium rule of integration can be applied (intdUold + intdU)*dvz; intdvx stores the integral over Uperp (hence over vx) in the closed orbit integral; intdxbar stores the value of the integral over xbar (which is the final result!), intdxbaropen does the same in the open orbit density integral; intdxbaropenBohm does the same for the Bohm integral; idealBohm is what the Bohm integral shoult be if Bohm condition is marginally satisfied; F is the value of the distribution function evaluated in the density integrals by interpolating FF, and Fold is the `old' needed to apply the trapezium rule; Fprime is the bilinearly interpolated value of FFprime, and Fprimeold is the same at the previous grid point (needed for trapezium rule); used in INTEGRALS OF DISTRIBUTION FUNCTION AT INFINITY; Ucap is the topmost total energy integrated to */
	double intdUopenflow = 0.0, intdUopenflowold = 0.0, intdxbaropenflow = 0.0, oorbintgrdflow = 0.0, oorbintgrdflowold = 0.0;
	// values of various integrals
	double oorbintgrdold=0.0, Fopen=0.0, intdUopenold=0.0;
	double vx0open; 
	double intdUantycal=0.0, intdvxantycal=0.0, vxnew=0.0, vxold = 0.0, Uperpnew = 0.0, *xtop, intdUopenantycal=0.0;
	double openorbitnew, chinew, munew = 0.0;
	/* intdUantycal is the integral over U (or v_z) for a flat potential profile (phi =0) for some value of xbar and Uperp; intdvxantycal  is the integral over Uperp (or vx) for a flat potential profile for some value of xbar; vxnew is the value of vx at the 'new' grid point, used in the vx integral (taken using the trapezium rule); vxold is the value of vx at the 'old' grid point, used in the vx integral; Uperpnew is the value of Uperp (used in the closed orbit density integral); munew is the valye of mu (used in the closed orbit density integral); xtop is the top bounce point x_t of the last closed orbit; intdUopenantycal is the analytical value of the integral over U  in the open orbit density integral */
	double xi, *gg, *ff;
	double flux0, du, fluxinf1old, fluxinfintgrdold, fluxinfintgrd, Qfluxinf1old, Qfluxinfintgrdold, Qfluxinfintgrd, u, Chodura2, Chodura2old, Chodura1old, Chodura1, Chodura;
	double fluxinf, Qfluxinf, fluxinf1, Qfluxinf1, densinf1, densinf, densinf1old; 
	double musmall, muell, Omegaell, minphiformucalc = 0.3;

	//printf("charge = %f\n", charge);
	//for (i=0; i<size_phigrid;i++) printf("index %d\tx = %f\tphi = %f\n", i, x_grid[i], phi_grid[i]);

	// Make grids in g = sqrt(x), x, and phi
	gg = malloc(size_phigrid*sizeof(double));
	ff = malloc(size_phigrid*sizeof(double));
	for (i=0; i<size_phigrid; i++) {
		gg[i] = sqrt(x_grid[i]);
		ff[i] = pow(sqrt(grid_parameter)+gg[i], 2.0) - grid_parameter;
		//printf("ff[%d] = %f\n", i, ff[i]);
	}
	deltax = ff[1];
	//deltax = gg[1];
	//size_finegrid = (int) (zoomfactor*ff[size_phigrid-1]/deltax) ;
	size_finegrid = zoomfactor*size_phigrid-zoomfactor;
	printf("size_finegrid = %d\t it is %d x %d - %d\n", size_finegrid, zoomfactor, size_phigrid, zoomfactor);
	//printf("size_finegrid = %d\n", size_finegrid);
	xx = (double*)calloc(size_finegrid,sizeof(double)); // xx = x has correct size
	phi = (double*)calloc(size_finegrid,sizeof(double)); // phi now has correct size
	FILE *fp, *filellip;
	filellip = fopen("checkellip.txt", "w");
	if (filellip == NULL) {	
		printf("Cannot open filellip.txt");
		exit(EXIT_FAILURE);
	}
	if (zoomfactor != 1) {
		i=0;
		gsl_interp_accel *acc = gsl_interp_accel_alloc ();
		gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, size_phigrid);
		gsl_spline_init (spline, ff, phi_grid, size_phigrid);
		//gsl_spline_init (spline, gg, phi_grid, size_phigrid);
		fp = fopen("OUTPUT/phispline.txt", "w");
		if (fp == NULL) printf("Error: phispline not created\n");
		for (i = 0; i < size_finegrid; i += 1)
		{
			xi = i*deltax/zoomfactor;
			if (i == 0) xi += 0.00001; //TINY;
			if (i == size_finegrid-1) xi -= 0.00001; //TINY;
			xx[i] = pow( pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), 2.0);
			//xx[i] = xi*xi;
			phi[i] = gsl_spline_eval (spline, xi, acc);
			//phi[i] = lin_interp(gg, phi_grid, sqrt(xx[i]), size_phigrid, 923);
			fprintf(fp, "%f %f\n", xx[i], phi[i]);
			xx[i] *= lenfactor;
			phi[i] *= (charge/Ti);
		}
		fclose(fp);
		gsl_spline_free (spline);
		gsl_interp_accel_free (acc);

	}
	else {
		i=0;
		fp = fopen("OUTPUT/phispline.txt", "w");
		if (fp == NULL) printf("Error: phispline not created\n");
		for (i = 0; i < size_finegrid; i += 1) {
			xi = i*deltax/zoomfactor;
			xx[i] = x_grid[i]*lenfactor;
			//if (i== size_finegrid-1) xx[i] -= TINY;
			phi[i] = phi_grid[i];
			fprintf(fp, "%f %f\n", xx[i], phi[i]);
			phi[i] *= (charge/Ti);
		}
		fclose(fp);
	}
	printf("lenfactor = %f\n", lenfactor);
	printf("Ti= %f\tcharge = %f\n", Ti, charge);
	printf("phi[0] = %f\tphi_grid[0] = %f\n", phi[0], phi_grid[0]);

	// Introduce a cap in energy (U, Uperp) high enough that we can safely assume F = 0
	Ucap = 10.0 + 10.0/Ti;
	if (bilin_interp(0.0, Ucap, FF, mumu, UU, sizemumu, sizeUU, -1, -1) > 1e-6) printf("ERROR in densfinorb_renorm.c: increase Ucap please\n");
	else if (bilin_interp(Ucap, 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1) > 1e-6) printf("ERROR in densfinorb_renorm.c: increase Ucap please\n");
	else printf("Ucap has been checked to be large enough\n");

	// Evaluate integrals of distribution function at infinity
	Qfluxinfintgrd = Qfluxinf1 = Qfluxinf = fluxinfintgrd = fluxinf1 = fluxinf = densinf = densinf1 = Chodura2 = Chodura1 = Chodura = 0.0;
	du = 0.01;
	printf("sizemumu = %d\tsizeUU=%d\n", sizemumu, sizeUU);
	for (i=0; i<sizemumu; i++) {
		fluxinf1old = fluxinf1;
		Qfluxinf1old = Qfluxinf1;
		Chodura1old = Chodura1;
		Chodura1 = 0.0;
		fluxinf1 = 0.0;
		Qfluxinf1 = 0.0;
		densinf1old = densinf1;
		densinf1 = 0.0;
		//Fprimeold = Fprime = 0.0;
		Fold = F = 0.0;
		fluxinfintgrdold = fluxinfintgrd = 0.0;
		for (j=0; j< sizeUU; j++) {
			u = sqrt(2.0*UU[j]);
			fluxinfintgrdold = fluxinfintgrd;
			Qfluxinfintgrdold = Qfluxinfintgrd;
			Chodura2old = Chodura2;
			Fold = F;
			F = FF[i][j];
			fluxinfintgrd = F*u;
			Qfluxinfintgrd = F*(mumu[i] + 0.5*u*u)*u; // 
			if (j>0) Chodura2 = F/(u*u);
			else Chodura2 = 0.0;
			if (j != 0) {
				du = u - sqrt(2.0*UU[j-1]); 
				fluxinf1 += 0.5*(fluxinfintgrd + fluxinfintgrdold)*du;
				Qfluxinf1 += 0.5*(Qfluxinfintgrd + Qfluxinfintgrdold)*du;
				Chodura1 += 0.5*(Chodura2 + Chodura2old)*du;
				densinf1 += 0.5*(F + Fold)*du;
			}
		}
		if (i!=0) {
			Qfluxinf += 0.5*(mumu[i]-mumu[i-1])*(Qfluxinf1 + Qfluxinf1old);
			fluxinf += 0.5*(mumu[i]-mumu[i-1])*(fluxinf1 + fluxinf1old);
			Chodura += 0.5*(mumu[i]-mumu[i-1])*(Chodura1 + Chodura1old);
			densinf += 0.5*(mumu[i]-mumu[i-1])*(densinf1 + densinf1old);
		}
	}
	densinf *= (4.0*M_PI);
	fluxinf *= (4.0*M_PI);
	Qfluxinf *= (4.0*M_PI);
	Chodura *= (4.0*M_PI);
	*flux = fluxinf/densinf;
	*Qflux = Qfluxinf/densinf;
	Chodura /= densinf;
	printf("densinf = %f\tfluxinf = %f\tQfluxinf = %f\tChodura = %f\n", densinf, *flux, *Qflux, Chodura);
	printf("sizemumu = %d\tsizeUU = %d\n", sizemumu, sizeUU);

	/* We initialize all arrays that contain functions of position x with the correct size n */
	phip       = malloc(size_finegrid*sizeof(double)); 
	phipp      = malloc(size_finegrid*sizeof(double)); 
	jmclosed   = malloc(size_finegrid*sizeof(int));
	jmopen     = malloc(size_finegrid*sizeof(int));
	xifunction = malloc(size_finegrid*sizeof(double));
	/* FORM XBAR GRIDS 
	Take derivatives of phi and use them to obtain two grids for xbar, one to be used for closed orbits and one to be used for open orbits. */
	xbar = (double*)calloc(size_finegrid,sizeof(double));	
	for (i=0; i<size_finegrid; i++) {	
		//if (i!=0)	xiprime[i-1] = (xiprime[i]-xiprime[i-1])/(xx[i]-xx[i-1]); 
		// Evaluate derivative of phi
		jmopen[i] = jmclosed[i] = 0; //*
		if (i == 0)
		{	
			phip[0] = (phi[1] - phi[0])/(xx[1]-xx[0]);
			//phip[0] = (phi[1] - phi[0])/(xx[1]-xx[0]) + (xx[0] - xx[1])*((phi[2] - phi[1])/(xx[2]-xx[1]) - (phi[1] - phi[0])/(xx[1]-xx[0]))/(xx[2]-xx[1]) ; 
		}
		else if (i == size_finegrid-1)
		{	
			//phip[i] = 0.5*(phi[i] - phi[i-1])/(xx[i] - xx[i-1]); 
			phip[i] = phip[i-1];
			//phip[0] = phip[1] + (xx[0] - xx[1])*(phip[2] - phip[1])/(xx[2]-xx[1]); 
		}
		else {	
			//phip[i] = (phi[i+1] - phi[i])/(xx[i+1] - xx[i]);
			phip[i] = ((xx[i] - xx[i-1])/(xx[i+1]- xx[i-1]))*(phi[i+1] - phi[i])/(xx[i+1] - xx[i]) + ((xx[i+1] - xx[i])/(xx[i+1]- xx[i-1]))*(phi[i] - phi[i-1])/(xx[i] - xx[i-1]); 
		}
		// xifunction is xbar corresponding to given position x being a stationary point
		xifunction[i] = xx[i] + phip[i];
		if (DEBUG == 1)
			printf("xifunction[%d] = %f\n", i, xifunction[i]); 
		if (i == 1) {	
			if  (xifunction[i] > xifunction[i-1]) { // immediately found that xi is increasing at x=0 telling us x_c = 0
				icrit = i-1;//=0;
			}
		}
		else if (i > 1) {	
			if ( (xifunction[i] > xifunction[i-1]) && (xifunction[i-1] < xifunction[i-2] ) ) {
			//if ( (xifunction[i] + TINY > xifunction[i-1]) && (xifunction[i-1] < xifunction[i-2] + TINY) ) 
				icrit = i-1;
			}
			//else if ( (xifunction[i] - TINY < xifunction[i-1]) && (xifunction[i-1] + TINY > xifunction[i-2]) ) 
			else if ( (xifunction[i] < xifunction[i-1]) && (xifunction[i-1]  > xifunction[i-2]) ) {
			// found maximum of xi: this only happens when phi has noise in second derivative
				//phip[i] = 0.5*(phi[i] - phi[i-1])/(xx[i] - xx[i-1]); 

				printf("ERROR in densfinorb: too much noise in second derivative\n");
				printf("xifunction = %f\n", xifunction[i-1]);
				printf("i = %d\n", i);
			}
		} 
	}
	//printf("icrit = %d\n", icrit);

	for (i=0; i<size_finegrid; i++) {	
		//printf("xx[%d] = %f\n", i, xx[i]);
		//if (i < size_phigrid)
		//	printf("x_grid[%d] = %f\n", i, x_grid[i]);
		if (i == 0)
			phipp[0] = (phip[1] - phip[0])/(xx[1]-xx[0]);
		else if (i == size_finegrid-1)
			phipp[i] = phipp[i-1];
		else 	
			phipp[i] = ((xx[i] - xx[i-1])/(xx[i+1]- xx[i-1]))*(phip[i+1] - phip[i])/(xx[i+1] - xx[i]) + ((xx[i+1] - xx[i])/(xx[i+1]- xx[i-1]))*(phip[i] - phip[i-1])/(xx[i] - xx[i-1]); 
	}

	// Whole commented section below perhaps was overkill in xbar resolution
	//////// but keep in case I ever change my mind
	//j=0;
	//k = icrit+1; 
	//// the smallest value of xbar is not calculated at the critical point, but one step ahead
	//for (i=icrit-1;i>=0;i--) {
	//	while ( (xifunction[k] < xifunction[i]) && (k<size_finegrid) ) {	
	//		xbar[j] = xifunction[k]; 
	//		j++; k++; 
	//	}	
	//	xbar[j] = xifunction[i]; 
	//	j++; 
	//}
	//while (k<size_finegrid) { 	
	//	xbar[j] = xifunction[k]; 
	//	k++; j++; 
	//}	
	////////

	j=0; // set counting index to zero
	for (k=icrit+1; k < size_finegrid; k++) { //*
		//if (j!=0) xbar[2*j-1] = xbar[2*j-2] + 0.5*(xifunction[k]-xbar[2*j-2]); 
		//xbar[2*j] = xifunction[k]; 
		xbar[j] = xifunction[k]; 
		//if (xbar[j] > xx[size_finegrid-1] - 0.5*limit_rho) maxj = j;
		//j+=2;
		j++;
	}
	xbarcrit = xifunction[icrit];
	chiMcrit = 0.5*phip[icrit]*phip[icrit] + phi[icrit];
	sizexbar = j;
	// ensure last xbar is within bounds //xbar[sizexbar-1] = xx[size_finegrid-1];
	//sizexbar = j-1;
	//for (j = 1; j < sizexbar; j+=2) xbar[j] = (xbar[j-1] + xbar[j+1])*0.5;
	if (DEBUG == 1) for (j = 0; j < sizexbar; j++) printf("xbar[%d] = %f\n", j, xbar[j]); 
	//printf("sizexbar = %d, size_finegrid (x) =%d\n", sizexbar, size_finegrid);
	//
	// Lots of array allocations now that size of xbar (vy + x)
	chi = (double **)calloc(sizexbar,sizeof(double*)); // chi(xbar, x) indices j and i 
	Uperp = (double**)calloc(sizexbar,sizeof(double*)); 
	mu = (double**)calloc(sizexbar,sizeof(double*)); 
	vx = (double***)calloc(sizexbar,sizeof(double**)); 
	upper = (int**)calloc(sizexbar,sizeof(int*)); 
	chimin = (double*)calloc(sizexbar,sizeof(double)); 
	chiMax = (double*)calloc(sizexbar,sizeof(double));
	chimpp = (double*)calloc(sizexbar,sizeof(double));
	crossed_min = (int*)calloc(sizexbar,sizeof(int)); 
	crossed_max = (int*)calloc(sizexbar,sizeof(int));
	kdrop = (int*)calloc(sizexbar,sizeof(int));
	openorbit = (double*)calloc(sizexbar,sizeof(double));
	upperlimit = (int*)calloc(sizexbar,sizeof(int));
	lowerlimit = (int*)calloc(sizexbar,sizeof(int));
	xtop = (double*)calloc(sizexbar,sizeof(double));
	imax = (int*)calloc(sizexbar,sizeof(int));
	imin = (int*)calloc(sizexbar,sizeof(int));

	/////////////////////////////////////////////////////
	/* CLOSED ORBIT ARRAY FILLING
	Set up the grid in xbar and also initialize all arrays that contain a different number at different values of xbar, indexed j */
	for (j=0;j<sizexbar;j++)  {
		imax[j] = imin[j] = -1;
		openorbit[j] = 0.0;
		crossed_min[j] = 0;
		crossed_max[j] = 0;
		xtop[j] = 0.0;
		chimin[j] = 0.0;
		chiMax[j] = 0.0;
		upperlimit[j] = -1;
		lowerlimit[j] = 0; 
	}
	/* Loop below initializes all 2d arrays which are functions of xbar and x. It allocates the right amount of memory to arrays of pointers of size xbar. The result is a 2D array indexed j (size sizexbar) and i or k (size n, see above) */
	for (j=0; j < sizexbar; j++) {
		chi[j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing chi(xbar, x) indexes j and i 
		Uperp[j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing Uperp for closed orbits indexes j and k
		mu[j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing adiab invariant mu(xbar, Uperp), indexes j and k
		upper[j] = (int*)calloc(size_finegrid,sizeof(int)); // array containing the index of energy Uperp corresponding to chi(x) 
		vx[j] = (double**)calloc(size_finegrid,sizeof(double*)); // see below
		/* The loop below initializes the 3d array containing the velocity of a particle at a given orbit position xbar, with particle position x and with energy Uperp, indexed i, j and k. */
		for (i = 0; i < size_finegrid; i++) {
			vx[j][i] = (double*)calloc(size_finegrid,sizeof(double)); 
		} 
	} // 3D array containing value of vx at different xbar, x and Uperp
	for (j=0; j<sizexbar; j++) {	
		chi[j][0] =  0.5*pow((xx[0] - xbar[j]), 2.0) + phi[0];
	}
	//chiMax[0] = chiMcrit;
	//mu[0][0] = 0.0;
	/* This for loop fills in the arrays, calculating the integrals where necessary */
	for (i=1; i<size_finegrid; i++) {
		for (j=0; j<sizexbar; j++)
		{	
			chi[j][i] =  0.5*pow((xx[i] - xbar[j]), 2.0) + phi[i];
	/* 
	FINDING MAXIMA/MINIMA
	*/
	/* Below, we use the array elements to define arrays for the effective potential maxima and minima that exist for every xbar (index j). We search through the chi curve from x=0 to the largest value of x but we search through every chi curve first (so first scan in xbar at fixed x, then move to the next x). We create arrays of mu, Uperp and vx.
	Because I am comparing neighbouring values of x to find a maximum, I need to consider two distinct cases.The first one is treated in the if loop below, which the program should enter only if chi is decreasing at x=0. Implying that chi(x=0) is an effective potential maximum. The second one is the else if loop after that, which finds maxima of chi that are stationary points by comparing the value of the function before and after each point. Note that because we compare the function at a point with the function at points before and after, the point we consider at every iteration step in the index i is indexed (i-1), compared with (i-2) and i.  */
			if ( (i==1) && (chi[j][i] < chi[j][i-1]) ) {	
				crossed_max[j] += 1;		
				imax[j] = 0;
				//printf("phiDSbump = 
				if (phi_DSbump > 1e-6) { // phi_DSbump is an attempt, still not implemented, to include a potential bump that reflects ions from the Debye sheath
					chiMax[j] = chi[j][0] + phi_DSbump;
					kdrop[j] = ceil(phi_DSbump/deltaE); 
					deltaE = phi_DSbump/kdrop[j];
				}
				else {
					chiMax[j] = chi[j][0];
					kdrop[j] = 0;
					deltaE=0.0;
				}
				for (k=0; k<=kdrop[j]; k++) {
					Uperp[j][k] = chiMax[j] - k*deltaE;
					mu[j][k] = 0.0;
					vx[j][0][k] = sqrt(2.0*Uperp[j][k] - chi[j][0]); 
					//printf("vx[%d][%d][%d] = %f, deltaE = %f, Uperp = %f\n", j, i-1, k, vx[j][0][k], deltaE, Uperp[j][k]);
						
				}
				upper[j][0] = kdrop[j];
				//printf("upper[%d][0] = %d\n", j, kdrop[j]);
			} 
			else if ( (i > 1) && ((chi[j][i] < chi[j][i-1]) && (chi[j][i-1] > chi[j][i-2])) ) {
				//printf("i = %d/%d\tj = %d/%d\n", i, size_finegrid, j, sizexbar);
				crossed_max[j] += 1;
	// crossed_min and crossed_max let the code know whether we go up or down the effective potential well. We are starting to go down!
				if (crossed_max[j] > 1)
				{	
					printf("***WARNING*** There is more than one maximum!\n");
					printf("j is %d and maxima at %d and %d for first and second\n",  j, imax[j], i-1);
					for (k=imax[j]-1;k<i+1;k++) {
					printf("chi[%d][%d] = %f\n", j, k, chi[j][k]);}
					crossed_max[j] = 1; 
					printf("exit code now\n");
					exit(-1);
				}
				imax[j] = i-1; // Store the index of the maximum 
				chiMax[j] = chi[j][i-1]; //Store chi Maximum itself, a function of xbar (index j)
			}
	/* We store the index corresponding to the position of a minimum for a given value of xbar.*/
			//watch out HERE
			else if ( (i > 1) && ((chi[j][i] > chi[j][i-1]) && (chi[j][i-1] < chi[j][i-2])) ) {	
				crossed_min[j] += 1;
				//printf("i-1 = %d, imax[%d] = %d, crossed_min[%d] = %d, crossed_max[%d] = %d\n", i-1, j,imax[j], j, crossed_min[j], j, crossed_max[j]);
				if (crossed_min[j] > 1) {	
					printf("***WARNING*** There is more than one minimum!\n");
					printf("j is %d and minima at %d and %d for first and second\n",  j, imin[j], i-1);
					for (k=imin[j]-1;k<i+1;k++) {
					printf("chi[%d][%d] = %f\n", j, k, chi[j][k]);}
					//if (i-1 - imin[j] > 2) 	
					//{	exit(1); } 
					crossed_min[j] = 1; 
				}
				imin[j] = i-1; // Store the index of the maximum 
				//printf("imin[%d] = %d\n", j, imin[j]);
				upperlimit[j] = imin[j] - imax[j] + kdrop[j];
				chimpp[j] = ( ( chi[j][i] - chi[j][i-1] ) / (xx[i] - xx[i-1]) - (chi[j][i-1] - chi[j][i-2]) / (xx[i-1] - xx[i-2]) ) *2.0/ (xx[i] - xx[i-2] ) ;  
				//printf("chimpp[%d] = %f\n", chimpp[j]);
				chimin[j] = chi[j][i-1];  //Store chi minimum itself, a function of xbar (index j)
				crossed_max[j] = 1;
	/* We will now start going up the effective potential well! The temporary flag below is used because we still have to store the effective potential minimum as a possible value of Uperp. If we don't have this flag we miss the region near the minimum of chi. */
			}
	/* FILLING IN ARRAYS */
	/* Once we cross a maximum, we start filling in arrays for mu, Uperp and vx related to the orbit we are in. As we go down the maximum we store the values of Uperp = chi(x) we encounter which will form a grid of (unevenly spaced) allowed values of Uperp. We also store the value of the small deltamu associated with every value of Uperp above the current value of chi(x), and add it to the previous values */                                                                               	
			if ( ( (crossed_max[j] == 1  && crossed_min[j] == 0) || (i-1 == imin[j]) )  ) {
				Uperp[j][i-1-imax[j]+kdrop[j]] = chi[j][i-1];
				if (Uperp[j][i-1-imax[j]] > Ucap && lowerlimit[j] != 0)
					lowerlimit[j] = i-1-imax[j]+kdrop[j]; 
				mu[j][i-1-imax[j]+kdrop[j]] = 0.0;
				//printf("kdrop[j] = %d\n", kdrop[j]);
				upper[j][i-1] = i-1-imax[j]+kdrop[j];
	/* Note that the size of the dimension of the array with values of Uperp is set to n, which is larger than the size it will turn out to be. This is because in C there is no way to append elements to arrays as I go along, enough memory has to be given to the array from the start. n is the largest possible size the array could have. */
				for (k=0;k<=upper[j][i-1]; k++) 
				{	
					//printf("k=%d/%d\n", k, upper[j][i-1]);
					// replaced k with upper below
					if ( (upper[j][i-1] == 0) || (i-1 == 0) ) {
						vx[j][i-1][k] = sqrt((2.0*Uperp[j][k] - chi[j][i-1])); // equiv 0
						//vx[j][i-1][k] = 0.0;
						//printf("j = %d, i-1 = %d, %f\n", j, i-1, chi[j][i-2]-chi[j][i-1]);
						mu[j][k] += 0.0; 

					}
					else if ( (k == imin[j] - imax[j] + kdrop[j]) && (crossed_min[j] == 1) ) {
						vx[j][i-1][k] = 0.0;
						mu[j][k] = 0.0; 
					}
					else if ( (k == imin[j] - imax[j] - 1 + kdrop[j]) && (crossed_min[j] == 1) ) {
						vx[j][i-1][k] = sqrt(2.0*(Uperp[j][k] - chi[j][i-1])); // equiv 0
						//mu[j][k] += sqrt(0.5)*(1.0/M_PI)*sqrt(0.5*chimpp[j])*pow(xx[imin[j]-1] - xx[imin[j]], 2.0);
						mu[j][k] = 0.5*pow(xx[imin[j]-1] - xx[imin[j]], 2.0)*pow(chimpp[j], 0.5);
						//printf("1..%f\n", sqrt(0.5*chimpp[j])*pow(xx[imin[j]-1] - xx[imin[j]], 2.0));
						// should uncomment above and comment below (spurious)
						//mu[j][k] += (2.0/M_PI)*sqrt(chi[j][i-2]-chi[j][i-1])*(2.0/3.0)*(xx[i-1] - xx[i-2]);
					}
					else if ( k == upper[j][i-1] - 1 )  {	
						vx[j][i-1][k] = sqrt(2.0*(Uperp[j][k] - chi[j][i-1]));// equiv 0
						mu[j][k] += (sqrt(2.0)/M_PI)*sqrt(chi[j][i-2]-chi[j][i-1])*(2.0/3.0)*(xx[i-1] - xx[i-2]);
						//mu[j][k] += (2.0/M_PI)*sqrt(chi[j][i-2]-chi[j][i-1])*(2.0/3.0)*(xx[i-1] - xx[i-2]); // with previous normalization
						//printf("2..i=%d,j=%dterm=%f\n", i,j, sqrt(chi[j][i-2]-chi[j][i-1])*(2.0/3.0)*(xx[i-1] - xx[i-2]));
					}
					else if (k == upper[j][i-1]) {
						vx[j][i-1][k] = 0.0;
					}
					else {
						vx[j][i-1][k] = sqrt(2.0*(Uperp[j][k] - chi[j][i-1]));
						mu[j][k] += (1.0/M_PI)*0.5*(vx[j][i-1][k] + vx[j][i-2][k])*(xx[i-1] - xx[i-2]); 
						//printf("3..%f\n", (vx[j][i-1][k] + vx[j][i-2][k])*(xx[i-1] - xx[i-2]));
					}
					if (mu[j][k] != mu[j][k]) {
						printf("BEFORE: mu[%d][%d] is NAN, kdrop[%d] = %d\n", j, k, j, kdrop[j]); 
						exit(-1);
					}  
				}
			}
	/* Once we cross the minimum, we stop creating array elements with values of Uperp. However, we keep storing the value of vx associated with any given point x on an effective potential curve with xbar, with energy Uperp and using this value to finish performing the mu integral. This should happen as long the effective potential at the point under consideration is smaller than the effective potential maximum. */
			else if ( ( crossed_min[j] == 1 && crossed_max[j] == 1 && chi[j][i-1] < chiMax[j] && ( i-1 != imin[j] ) ) ) {
				for (k=0;k <= upperlimit[j] ;k++)
				{	
					//if (i-1 == imin[j] +1) {
					//	vx[j][i-1][k] = sqrt(2.0*(Uperp[j][k] - chi[j][i-1]));
					//	mu[j][k] += (1.0/M_PI)*0.5*(vx[j][i-1][k] + vx[j][i-2][k])*(xx[i-1] - xx[i-2]); 
					//}
					//else if (chi[j][i-1] < Uperp[j][k]) {	
					if ( (chi[j][i-1] < Uperp[j][k]) && (chi[j][i-2] < Uperp[j][k]) ) {	
						vx[j][i-1][k] = sqrt(2.0*(Uperp[j][k] - chi[j][i-1]));
						mu[j][k] += (1.0/M_PI)*0.5*(vx[j][i-1][k] + vx[j][i-2][k])*(xx[i-1] - xx[i-2]); 
						//printf("vx[j][i-2][k] = %f\n", vx[j][i-2][k]);
					}
					else if (Uperp[j][k] <= chi[j][i-1] && Uperp[j][k-1] > chi[j][i-1]) {
						upper[j][i-1] = k;
				//mu[j][k] += (2.0/M_PI)*(vx[j][i-2][k])*(xx[i-1] - xx[i-2])*(Uperp[j][k] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]);
						ind = 0;
						while (Uperp[j][k] < chi[j][i-2-ind]) 
							ind++;
						
						mu[j][k] += (sqrt(2.0)/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[j][k] - chi[j][i-2-ind], 1.5)/(chi[j][i-1-ind] - chi[j][i-2-ind]); // double check normalization
						//}
						//printf("DEBUG: %f and mu = %f: i = %d, j = %d, k = %d \nchi = %f, %f\n", mu[j][k], (2.0/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[j][k] - chi[j][i-2-ind], 1.5)/(chi[j][i-1-ind] - chi[j][i-2-ind]), i, j, k, chi[j][i-1], chi[j][i-2]);
						//printf("MORE DEBUG: %f\n", (2.0/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[j][k] - chi[j][i-2-ind], 1.5)/(chi[j][i-1-ind] - chi[j][i-2-ind]));
			
					}
					else if (Uperp[j][k] <= chi[j][i-1] && Uperp[j][k] > chi[j][i-2]) {
						mu[j][k] += (sqrt(2.0)/M_PI)*(2.0/3.0)*(xx[i-1] - xx[i-2])*pow(Uperp[j][k] - chi[j][i-2], 1.5)/(chi[j][i-1] - chi[j][i-2]); 
					}
					if (mu[j][k] != mu[j][k]) {
						printf("mu[%d][%d] is NAN, kdrop[%d] = %d\n", j, k, j, kdrop[j]); 
						exit(-1);
					}  
				}
			}
	/* When the effective potential at the iteration point (i-1) under consideration becomes larger than the effective potential maximum, we finish performing the mu integral. We also store the position of the top of the orbit which has chi = chiMax, in order to perform the open orbit integral. If the loop below is accessed, a switch it turned off to signify the no more closed orbits can be present */
			//else if ( ( crossed_min[j] == 1 ) && ( crossed_max[j] == 1) && ( chi[j][i-1] > chiMax[j] - TINY) ) 	
			else if ( ( crossed_min[j] == 1 ) && ( crossed_max[j] == 1) && ( chi[j][i-1] > chiMax[j] ) ) {	
				//itop[j] = i-2;
				xtop[j] = xx[i-2] + ((chiMax[j] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]))*(xx[i-1] - xx[i-2]);
				for (k=0; k<upper[j][i-2]; k++) {
					ind = 0;
					//while (Uperp[j][0] < chi[j][i-2-ind]) 
					//	ind++; // find top bounce point index i-2-ind for a given xbar[j] and for Uperp[j][0] = chiMax
					//mu[j][k] += (sqrt(2.0)/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[j][0] - chi[j][i-2-ind], 1.5)/(chi[j][i-1-ind] - chi[j][i-2-ind]);  // CHECK NORMALIZATION
					while (Uperp[j][k] < chi[j][i-2-ind]) 
						ind++; // find top bounce point index i-2-ind, with i the smallest index such that chi[j][i-1] > chiMax[j] for a given xbar[j], for orbits with Uperp[j][k] > chi[j][i-2] 
					mu[j][k] += (sqrt(2.0)/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[j][k] - chi[j][i-2-ind], 1.5)/(chi[j][i-1-ind] - chi[j][i-2-ind]);  // CHECK NORMALIZATION
					//printf("mu[%d][%d] = %f\n", j, k, mu[j][k]);
					//mu[j][k] += (1.0/M_PI)*(vx[j][i-2][k])*(xx[i-1] - xx[i-2])*(Uperp[j][k] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]);
				}
				crossed_max[j] = 0; 
			}
			if (j!=0) {	
				if ( ((chiMax[j-1] < chi[j-1][i-1] + TINY ) && (chiMax[j] > chi[j][i-1] - TINY)) ) { // || (imax[j-1] == -1 && imax[j] != -1) )
					jmclosed[i-1] = j-1;
					if (i - 1 > icrit) { 	
						jmopen[i-1] = j-1; 
						if (DEBUG == 1) {
							printf("i = %d, j = %d, icrit = %d\n", i-1, j-1, icrit); 
							printf("jmopen[%d] = %d\n", i-1, jmopen[i-1]); 
						}
					} 
				}
			} 
			//printf("upper = %d, upperlimit = %d\n", upper[j][i-1], upperlimit[j]);
		} 
	}	


	upperlimit[sizexbar-1] = upperlimit[sizexbar-2] +1;
	//printf("xbar[maxj=%d] = %f\n", maxj, xbar[maxj]);	
	// OPEN ORBIT INTEGRAL
	/* Now we perform the open orbit integral. We use a change of variables which makes the integrand smooth at the top bounce point. The change of variables is to some var = sqrt(x_t - x) */
	maxj = sizexbar +1 ; //*
	muopen = calloc(maxj,sizeof(double)); //*
	chiMopen = calloc(maxj,sizeof(double));//*	
	xbaropen = calloc(maxj,sizeof(double));//*
	openorbitopen = calloc(maxj,sizeof(double));
	Ucritf = calloc(maxj,sizeof(double));	

	xbaropen[0] = xbarcrit; //*
	chiMopen[0] = chiMcrit; //*
	muopen[0] = 0.0; //*
	openorbitopen[0] = 0.0; //*
	//Ucritf[0] = phi[0] + 0.5*xbarcrit*xbarcrit;
	Ucritf[0] = chiMcrit ;
	printf("sizexbar = %d\n", sizexbar);


	for (j=0;j<maxj-1;j++) { //*
		//printf("lowerlimit[%d] = %d and upper[%d][0] = %d, kdrop = %d, upperlimit = %d\n", j, lowerlimit[j],j, upper[j][0], kdrop[j], upperlimit[j]);
		mu[j][upperlimit[j]] = 0.0;
		//mu[j][upperlimit[j]-1] = 0.5*pow(xx[imin[j]] - xx[imin[j]-1], 2.0)*pow(chimpp[j], 0.5);
		chiMopen[j+1] = chiMax[j];
		k=0;
		if (APPROXMUFORSMALLORBIT == 1) {
			if (fabs(phi[0]) < minphiformucalc) {
				if (xbar[j] > 0.0) { //otherwise lin_interp fails
					phibar = lin_interp(xx, phi, xbar[j], size_finegrid, 1);
					//printf("APPROXIMATE MU if potential drop across region is small\n");
					//while (fabs(phibar) < minphiformucalc && k < upperlimit[j]) {
					//while (fabs(lin_interp(xx, phi, xbar[j]+0.1, size_finegrid, 1) - phibar)/0.1 < minphiformucalc && k <= upperlimit[j]) {
					Omegaell = sqrt(1.0 + lin_interp(xx, phipp, xbar[j], size_finegrid, 1));
					//printf("Omegaell = %f\n", Omegaell);
					muell = Uperp[j][upperlimit[j]-k] - phibar + 0.5*pow(lin_interp(xx, phip, xbar[j], size_finegrid, 1)/Omegaell, 2.0);
					//printf("Uperp[j][k] = %f\n", muell);
					muell /= Omegaell;
					//printf("muell = %f\n", muell);
					//mu[j][upperlimit[j]-k] = 0.5*pow(xx[imin[j]] - xx[imin[j]-k], 2.0)*pow(chimpp[j], 0.5);
					musmall = 0.5*pow(xx[imin[j]] - xx[imin[j]-k], 2.0)*pow(chimpp[j], 0.5);
					//if ( (xbar[j] < 8.0) ) // (muell < 1.0) ) //&& 
					//	printf("at xbar[j] = %f\tmu[j][k] = %f\tmuell = %f\tfractional error = %f\tmuexp = %f\tfractional error = %f\n", xbar[j], mu[j][upperlimit[j]-k], muell, muell/mu[j][upperlimit[j]-k] - 1.0, muexp, muexp/mu[j][upperlimit[j]-k] - 1.0);
					if (muapproxtype == 1) 
						mu[j][upperlimit[j]-k] = muell;
					else
						mu[j][upperlimit[j]-k] = musmall;
					k++;
				}
			}
		}
		muopen[j+1] = mu[j][0]; //*
		xbaropen[j+1] = xbar[j];
		Ucritf[j+1] = chiMax[j] - mu[j][0] ; //*
		//printf("mucrit , Ucrit = %f, %f\n", muopen[j+1], Ucritf[j+1]);
		if (DEBUG == 1)
			printf("jmopen[%d] = %d\n", j, jmopen[j]); 
		//if (j==0) twopimuprime[0] = 0.0; //*
		if (j==0) //*
			openorbit[j] = (2.0*M_PI) * ( ((xbar[j] - xbarcrit)/(xbar[j+1] - xbarcrit)) *(mu[j+1][0] - mu[j][0])/(xbar[j+1] - xbar[j]) + ((xbar[j+1] - xbar[j])/(xbar[j+1] - xbarcrit)) * (mu[j][0] - 0.0)/(xbar[j] - xbarcrit) );	
		else if (j==maxj-2) 
			openorbit[j] = openorbit[j-1];	
		else //* 
			openorbit[j] = (2.0*M_PI) * (((xbar[j] - xbar[j-1])/(xbar[j+1] - xbar[j-1])) *(mu[j+1][0] - mu[j][0])/(xbar[j+1] - xbar[j]) + ((xbar[j+1] - xbar[j])/(xbar[j+1] - xbar[j-1])) * (mu[j][0] - mu[j-1][0])/(xbar[j] - xbar[j-1]) );	
			//openorbit[j-1] = 2.0*M_PI*(mu[j][0] - mu[j-1][0])/(xbar[j] - xbar[j-1]);	
		openorbitopen[j+1] = openorbit[j];
		openorbitantycal = 2.0*M_PI*xbar[j];
		if (DEBUG == 1)
			printf("%f %f %f %f %f\n", xbar[j], mu[j][0], Uperp[j][0], openorbit[j], openorbitantycal); 
		if (TESTELL == 1) { // test elliptical orbit approximation // ignore
			printf("TESTELL entered\n");
			if (xbar[j] < xx[size_finegrid-1]) {
				printf("xbar condition fulfilled\n");
				for (k=0; k < upperlimit[j]; k++) {
					Omegaell = sqrt(1.0 + lin_interp(xx, phipp, xbar[j], size_finegrid, 1));
					muell = Uperp[j][k] - lin_interp(xx, phi, xbar[j], size_finegrid, 1) + 0.5*pow(lin_interp(xx, phip, xbar[j], size_finegrid, 1)/Omegaell, 2.0);
					muell /= Omegaell;
					if ( (mu[j][k] < 2.5) && (mu[j][k] > 1.5) ) {
						printf("xbar = %f\tOmegaell = %f\tmu (actual, ellipmodel) = (%f, %f)\n", xbar[j], Omegaell, mu[j][k], muell);
						fprintf(filellip, "%f %f %f %f\n", xbar[j], Omegaell, mu[j][k], muell);
					}
				}
			}
		}
	}
	// temporary
	muopen[maxj-1] = 999999.0;
	chiMopen[maxj-1] = 999999.0; 
	xbaropen[maxj-1] = 100.0;
	openorbitopen[maxj-1] = 100.0;
	Ucritf[maxj-1] = 0.0;
	//
	//Ucritf[0] = Ucritf[1] - muopen[1]*(Ucritf[2] - Ucritf[1])/(muopen[2] - muopen[1]);
	//Ucritf[0] = Ucritf[1];
	//Ucritf[0] = sqrt(phi[0] + xbarcrit*xbarcrit);
	//if (Ucritf[0] < 0.0) Ucritf[0] = phi[0];
	if (DEBUG == 1) {
		printf("~~~~~The second element of FF is %f~~~~~\n",FF[0][1]);
		printf("~~~~~The second element of UU is %f~~~~~\n",UU[1]);
		printf("~~~~~The fourth element of mu is %f~~~~~\n",mumu[3]);
	}

	i=0;
	clock_t int1 = clock(); // Finds the time of the computation so far
	double inttime  = (double)(int1 - begin) / CLOCKS_PER_SEC;
	if (DEBUG == 1) 
		printf("in densfinorb: Array filling DONE: time is %f\n", inttime);
	/* DENSITY INTEGRALS 
	This part calculates the density integrals and outputs the result of the integration to a file fout and also the yz distribution function to three files one containing the distribution function the other two containing the velocity grid */
	FILE *fout; 
	if (Ti < 10.1) 
		fout = fopen("OUTPUT/densfinorb_out.txt", "w");
	else 
		fout = fopen("TESTS/densfinorb_out.txt", "w");
	if (fout == NULL) {	
		printf("Cannot open densfinorb_out.txt");
		exit(EXIT_FAILURE);
	}

	////////////////////////////////////////////
	// calculate normalization at infinity
	ic = 0;
	while(x_grid[ic+1] < x_grid[size_phigrid-1] - limit_rho) ic++;
	size_xlim = ic;
	printf("size_xlim = %d/%d, limit_rho = %f\n", size_xlim, size_phigrid, limit_rho);
	//deltax_inf = (xx[size_finegrid-1] - xx[size_finegrid-2])*0.8; // overkill resolution to get n_inf correct
	deltax_inf = (xx[size_xlim*zoomfactor] - xx[size_xlim*zoomfactor-1]); // overkill resolution to get n_inf correct
	printf("deltax_inf = %f\n", deltax_inf);
	//deltax_inf = (xx[1] - xx[0]); // higher resolution than at infinity
	j_inf = (int) sqrt(2.0*Ucap)/deltax_inf ;
	printf("j_inf = %d, deltax_inf = %f\n", j_inf, deltax_inf);
	muinf = (double*)calloc(j_inf,sizeof(double*)); 
	vxinf = (double**)calloc(j_inf,sizeof(double**)); 
	chiinf = (double*)calloc(j_inf,sizeof(double*));
	for (j=j_inf-1; j>=0; j--) {	
		vxinf[j] = (double*)calloc(j_inf,sizeof(double*)); 
		chiinf[j] =  0.5*deltax_inf*j*deltax_inf*j;
		for (k=j; k<j_inf; k++) {
			vxinf[j][k] = sqrt(2.0*(chiinf[k] - chiinf[j]));
			//if (j!=j_inf-1) {
			//	if (k>j+1)
			//		muinf[k] += (2.0/M_PI)*(vxinf[j+1][k] + vxinf[j][k])*deltax_inf;
			//	else if (k==j+1)
			//		muinf[k] += (4.0/M_PI)*sqrt(chiinf[j+1]-chiinf[j])*(2.0/3.0)*deltax_inf;
			//}
		}
		muinf[j] = chiinf[j];
	}

	intdxbar = 0.0;
	for (j=0; j<j_inf; j++) {
		vxnew = 0.0;
		intdvxold = intdvx;
		intdvx = 0.0;
		intdU = 0.0;
		for (k=j; k<j_inf; k++) {
			vxold = vxnew;
			intdUold = intdU;
			intdU = 0.0;
			munew = muinf[k]; 
			vxnew = vxinf[j][k];
			Ucrit = lin_interp(muopen, Ucritf, munew, maxj, 790);
			//printf("vxnew = %f\tmunew = %f\tUcrit = %f\n", vxnew, munew, Ucrit);
			Uperpnew = chiinf[k];
			sizeU = (int) sqrt(2.0*Ucap - 2.0*Uperpnew)/dvz;
			if (phi[0] > 0.0) 
				reflected = 1;
			else
				reflected = 0;
			for (l=0; l < sizeU; l++) {	
				if (l!=0) {	
					Fold = F;
					Fold_ref = F;
					vz = dvz*l;
					U = Uperpnew + 0.5*pow(vz, 2.0);
					//if ( (U > munew) && (U - 0.5*vz*vz + 0.5*(vz-dvz)*(vz-dvz) < munew) ) {
					//	frac = (vz - sqrt(2.0*(munew - Uperpnew)))/dvz;
					//	Fold = bilin_interp(munew, 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					//	F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					//}
					//else if (U > munew){
					//	frac = 1.0;
					//	F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					//}
					//else {
					//	frac = 1.0;
					//	F = 0.0;
					//}
					frac = 1.0;
					F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 

					if ( (U-munew < Ucrit - numb) && (reflected == 1) ) 
						frac_reflected = 1.0;
					else if ( (U-munew > Ucrit - numb) && (reflected == 1) ) {
						reflected = 0;
						if (Ucrit < Uperpnew - munew + numb) frac_reflected = 0.0;
						else {
							vzcrit = sqrt(2.0*(Ucrit + munew - Uperpnew));
							Fold_ref = bilin_interp(munew, Ucrit, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
							frac_reflected = (vzcrit - (vz - dvz))/dvz;
						}
					}
					else frac_reflected = 0.0;

					//if ( (charge < 0) && (frac != 1.0) ) printf("WARNING\n\n\n\n\n\n\n\n");
					//if ( (charge < 0) && (l == 1) ) printf("frac_reflected = %f\n", frac_reflected);
					
					//if (fabs(frac_reflected) > TINY) printf("frac_reflected = %f, charge = %f\n", frac_reflected, charge);
					if (charge > 0.0) frac_reflected = 0.0;
					//frac_reflected = 0.0;
					intdU += 0.5*frac*dvz*((F+Fold) + frac_reflected*(F + Fold_ref));
					//if (reflected == 1) printf("particles are being reflected\n");
				}
				else {	
					U = Uperpnew ;
					F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					intdU += 0.0;
				} 
			}

			intdUantycal = exp(-chiinf[k])*(1.0/(2.0*M_PI));// result with phi =0
			if (DEBUG == 1) printf("Analytical intdU is %f, numerical one is %f\n", intdUantycal, intdU);
			if (k!=j) {
				dvx = vxnew - vxold; 
				//printf("intdU = %f\tintdUold = %f\nvx = %f\tvxold = %f\n", intdU, intdUold, vxnew, vxold);
				intdvx += 2.0*0.5*dvx*(intdU+intdUold);
			}
			intdvxantycal = (2.0/(2.0*sqrt(M_PI)))*exp(-0.5*deltax_inf*j*deltax_inf*j);//*erf(sqrt(xbar[j]*xbar[j]-(pos-xbar[j])*(pos-xbar[j])));
			intdvxantycal = sqrt(1.0/(2.0*M_PI))*exp(-0.5*deltax_inf*j*deltax_inf*j);//*erf(sqrt(xbar[j]*xbar[j]-(pos-xbar[j])*(pos-xbar[j])));
			//intdvx = intdvxantycal;
		}
		if (DEBUG == 1)	printf("intdvx is %f, analytical one is %f\n", intdvx, intdvxantycal); 
		//printf("intdvx is %f, while the analytical one is %f\n", intdvx, intdvxantycal);
		//intdvx = intdvxantycal;
		
		dxbar = deltax_inf;
		if (j!=0) {
			intdxbar += 0.5*(intdvx+intdvxold)*dxbar; 
			//printf("intdxbar = %f\n", intdxbar);
		}
	}
	n_inf = 2.0*intdxbar; // only half of the xbar domain is considered
	printf("n_inf = %f\n", n_inf);
	if (DEBUG == 1) {	
		printf("for charge = %f\n", charge);
		printf("n_inf = %f\n", n_inf);
		if ( (n_inf  != n_inf) || (n_inf < TINY) )
		{ 	
			printf("n_inf = %f\n", n_inf);
			exit(-1);
		} 
	} 

	// density profile
	stop = 0; // set stop index to zero; it turns to 1 if density exceeds threshold in the input (expressed as fraction of density at infinity)
	ic = 0;
	while (stop == 0) {
		i = ic*zoomfactor;
		intdxbar = 0.0;
		intdxbaropen = 0.0;
		intdxbaropenflow = 0.0;
		for (j=0; j<sizexbar; j++) { //*
			vxnew = 0.0;
			intdUopenold = intdUopen;
			intdUopenflowold = intdUopenflow; //intdUopenBohmold = intdUopenBohm;
			intdUopen = 0.0;
			intdUopenflow = 0.0; //intdUopenBohm = 0.0;
			//intdUopenxbar = 0.0; //intdUopensquare = 0.0;
			if (j == jmopen[i]) {
				oorbintgrd = oorbintgrdold = 0.0;
				oorbintgrdflow = oorbintgrdflowold = 0.0; 
				sizeU = (int) sqrt(2.0*(Ucap - chiMax[j]))/dvzopen;
				if (j==0) {
					munew = 0.0;
					openorbitnew = 0.0;
				}
				else {
					munew = mu[j+1][0] + ( (chiMax[j+1] - chi[j+1][i]) / (chiMax[j+1] - chi[j+1][i] + chi[j][i] - chiMax[j]) ) * (mu[j][0] - mu[j+1][0]);
					openorbitnew = openorbit[j+1] + ( (chiMax[j+1] - chi[j+1][i]) / (chiMax[j+1] - chi[j+1][i] + chi[j][i] - chiMax[j]) ) * (openorbit[j] - openorbit[j+1]);
				}
				for (l=0; l < sizeU; l++) {
					oorbintgrdold = oorbintgrd;
					oorbintgrdflowold = oorbintgrdflow; 
					vz = dvzopen*l;
					if (j!=0) { 
						chinew = chi[j+1][i] + ( (chiMax[j+1] - chi[j+1][i]) / (chiMax[j+1] - chi[j+1][i] + chi[j][i] - chiMax[j]) ) * (chi[j][i] - chi[j+1][i]) ;
						vx0open = 0.0;
					}
					else {
						chinew = chi[j][i];
						if (chiMax[j] > chinew) 
							vx0open = sqrt(2.0*(chiMax[j] - chinew));
						else vx0open = 0.0;
					}
					U = chinew + 0.5*pow(vz, 2.0); 
					if ( (U > munew) && (U - 0.5*vz*vz + 0.5*(vz-dvzopen)*(vz-dvzopen) < munew) ) {
						frac = (vz - sqrt(2.0*(munew - chinew)))/dvzopen;
						Fopen = bilin_interp(munew, 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
						//printf("frac = %f\n", frac);
						oorbintgrdold = ( sqrt(vx0open*vx0open + 2.0*alpha*sqrt(2.0*(munew - chinew))*openorbitnew) - vx0open )*Fopen;
						oorbintgrdflowold = 0.5*alpha*sqrt(2.0*(munew - chinew))*openorbitnew*Fopen;
						Fopen = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
						//****
						oorbintgrd = ( sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbitnew) - vx0open )*Fopen;
						oorbintgrdflow = 0.5*alpha*vz*openorbitnew*Fopen;
					}
					else {
						frac = 1.0;
						Fopen = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
						//printf("U = %f, munew = %f, Fopen = %f\n", U, munew, Fopen);
						oorbintgrd = ( sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbitnew + TINY) - vx0open )*Fopen;
						oorbintgrdflow = 0.5*alpha*vz*openorbitnew*Fopen;

					}
					if (oorbintgrd != oorbintgrd) {	
						printf("vx0open = %f, openorbit[%d] = %f, imaginary oorbintgrd in initial piece of integral due to negative value of openorbit?\n", vx0open, j, openorbit[j]); 
						exit(-1);
					}
					if (i==0) 
						//****
						oorbintgrdantycal = sqrt(2.0*alpha*vz*2.0*M_PI*xbar[j])*(2.0*(U-chiMax[j]))*exp(-U)/pow(M_PI, 1.5);  // CHECK NORMALIZATION
						// for a flat potential, oorbintgrd can be calculated analytically
						//oorbintgrd = oorbintgrdantycal;
						//printf("oorbintgrd is %f, analytical is %F\n", oorbintgrd, oorbintgrdantycal);
					if (l!=0) {	
						intdUopen += 0.5*frac*dvzopen*(oorbintgrd+oorbintgrdold);
						intdUopenflow += 0.5*frac*dvzopen*(oorbintgrdflow+oorbintgrdflowold);
					}
					else {	
						intdUopen += 0.0;
						intdUopenflow += 0.0;
					} 
				}
				if (j == 0) { //*
					intdUopenold = 0.0;
					//printf("intdUopenold = %f\n", intdUopenold);
					dxbar = xbar[0] - xbarcrit;
					dxbar = 0.0;
					intdxbaropen += 0.5*(intdUopen+intdUopenold)*dxbar;
					intdxbaropenflow += 0.5*(intdUopenflow+intdUopenflowold)*dxbar;
				}
			}
			else if (j > jmopen[i]) {
				oorbintgrd = oorbintgrdold = 0.0;
				oorbintgrdflow = oorbintgrdflowold = 0.0;
				sizeU = (int) sqrt(2.0*(Ucap - chiMax[j]))/dvzopen;
				for (l=0; l < sizeU; l++) {
					oorbintgrdold = oorbintgrd;
					oorbintgrdflowold = oorbintgrdflow;
					vz = dvzopen*l;
					U = chiMax[j] + 0.5*pow(vz, 2.0);
					//vx0open = sqrt(TINY + chiMax[j] - chi[j][i]);
					vx0open = sqrt(2.0*(chiMax[j] - chi[j][i])); // SEARCH HERE TO FIND CURRENT POSITION OF NEW CHANGES TO NORMALIZATION
					if (vx0open != vx0open) {
						printf("HERE imaginary vx0open, j = %d, i is %d, chi[j][i] = %f, chiMax[j] = %f\n", j, i, chi[j][i], chiMax[j]); 
						exit(-1);
					}
					if ( (U >= mu[j][0]) && (U - 0.5*vz*vz + 0.5*(vz-dvzopen)*(vz-dvzopen) < mu[j][0]) ) {
						frac = (vz - sqrt(2.0*(mu[j][0] - chiMax[j])))/dvzopen;
						Fopen = bilin_interp(mu[j][0], 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
						oorbintgrdold = sqrt(2.0*alpha*sqrt(2.0*(mu[j][0] - chiMax[j]))*openorbit[j])*Fopen;
						oorbintgrdflowold = alpha*sqrt(2.0*(mu[j][0] - chiMax[j]))*openorbit[j]*Fopen;
						Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
						oorbintgrd = sqrt(2.0*alpha*vz*openorbit[j])*Fopen;
						oorbintgrdflow = alpha*vz*openorbit[j]*Fopen;
					}
					else {
						frac = 1.0;
						Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
						oorbintgrd = (sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbit[j]) - vx0open)*Fopen;
						oorbintgrdflow = alpha*vz*openorbit[j]*Fopen;
						//if (oorbintgrd != oorbintgrd) 
						if (frac != frac) 
							printf("AHA\n\n\n\n\n");
					}
					//Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1);
					if (oorbintgrd != oorbintgrd) {
						printf("vx0open = %f, openorbit[%d] = %f, HERE imaginary oorbintgrd in initial piece of integral due to negative value of openorbit?\n", vx0open, j, openorbit[j]); 
						exit(-1);
					}
					if (i==0) {
						oorbintgrdantycal = sqrt(2.0*alpha*vz*M_PI*xbar[j])*(U-chiMax[j])*exp(-U)/pow(M_PI, 1.5); 
						if (DEBUG == 1) 
							printf("openorbit = %f (should be %f with flat potential profile)\n", oorbintgrd, oorbintgrdantycal);
						// check not working
					}
	// for a flat potential, oorbintgrd can be calculated analytically
					//oorbintgrd = oorbintgrdantycal;
	//printf("oorbintgrd is %f, analytical is %F\n", oorbintgrd, oorbintgrdantycal);
						//if (oorbintgrdold < TINY ) 
					//if ( (l == 1) && (vx0open < TINY) ) 
					//	intdUopen += 4.0* (2.0/3.0) *pow(dvzopen, 1.5) * sqrt(alpha*openorbit[j])*bilin_interp(mu[j][0], 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					if (l!=0) {
						intdUopen += 0.5*frac*dvzopen*(oorbintgrd+oorbintgrdold);
						intdUopenflow += 0.5*frac*dvzopen*(oorbintgrdflow+oorbintgrdflowold);
					}
					else {
						intdUopen += 0.0;
						intdUopenflow += 0.0;
					} 
				}
				if (intdUopen != intdUopen)  {
					printf("HERE, j is %d\n", j); 
					exit(-1);
				}
				if (i==0) {	
					intdUopenantycal = 2.0*0.919*sqrt(2.0*alpha*xbar[j])*exp(-xbar[j]*xbar[j])/M_PI; 
					if (DEBUG == 1) { 
						printf("xbar[%d] = %f, analytical is %f, numerical is %f\n", j, xbar[j], intdUopenantycal, intdUopen); printf("vx0open is %f (should be zero)\n", vx0open); 
					}
				} 
				dxbar = xbar[j] - xbar[j-1];
				if ( (j == jmopen[i]+1) &&  (j != 1) ) { // in this case dxbar is different //*
					dxbar = (xbar[j] - xbar[j-1])*(chiMax[j] - chi[j][i])/ (chiMax[j] - chi[j][i] + chi[j-1][i] - chiMax[j-1]);// open orbit density does not need to be so accurate at this point
					//intdUopenold = intdUopen + ( (chiMax[j] - chi[j][i])/ (chiMax[j] - chi[j][i] + chi[j-1][i] - chiMax[j-1]) ) * (intdUopenold - intdUopen);// open orbit density does not need to be so accurate at this point
					if (DEBUG == 1)
						printf("dxbar = %f\n", dxbar);
				}
	//(chi[j][k]-chi[j][i])/(2.0*(xx[i]-xx[k])); 
				intdxbaropen += 0.5*(intdUopen+intdUopenold)*dxbar;
				intdxbaropenflow += 0.5*(intdUopenflow+intdUopenflowold)*dxbar;
				//intdxbaropenBohm += 0.5*(intdUopenBohm+intdUopenBohmold)*dxbar;
				if (intdxbaropen != intdxbaropen) {
					printf("PROBLEM HERE, j is %d\n", j); 
					exit(-1);
				} 
			}
			if ( j>=jmclosed[i] ) {	
				/* We have entered the closed orbit integral */
				intdvxold = intdvx;
				//intdvxflowold = intdvxflow;
				//intdvxflow = 0.0;
				intdvx = 0.0;
				if (j == jmclosed[i]) {
					intdvx = 0.0;
					//intdxbar += 0.0; 
				}
				else if (j > jmclosed[i]) {	
					if (DEBUG ==1) {
						printf("lowerlimit[%d] = %d\tupper[%d][%d] = %d\n", j, lowerlimit[j], j, i, upper[j][i]);
						printf("mu = %f\n", mu[j][lowerlimit[j]]);
						printf("imax[%d] = %d\n", j, imax[j]);
					}
					for (k=lowerlimit[j]; k<upper[j][i]+1; k++) {	
						vxold = vxnew;
						intdUold = intdU;
						intdU = 0.0;
						if (k == upper[j][i]) {
							Uperpnew = chi[j][i];
							if (lowerlimit[j] == upper[j][i] ) {
								//munew = 0.0;
								munew = mu[j][k]; 
								//printf("why would I ever be here? i = %d mu = %f\n\n\n", i, munew);
							}
							else if (k == upperlimit[j]) // correct but unnecessary as taken ino account below
								munew = 0.0;
							else {
								munew = ((chi[j][i] - Uperp[j][k])*mu[j][k-1] + (Uperp[j][k-1] - chi[j][i])*mu[j][k])/(Uperp[j][k-1] - Uperp[j][k]);
							}
							//munew = mu[j][k]; 
							vxnew = 0.0;
							dvx = vxold - vxnew; 
						}
						else {	
							Uperpnew = Uperp[j][k];
							vxnew = vx[j][i][k];
							dvx = vxold - vxnew;
							munew = mu[j][k]; 
						}
						//printf("muopen[j] = %f, munew = %f, mu[j][k] = %f\t j = %d\tk=%d\n", muopen[maxj-1], munew, mu[j][k], j, k);
						Ucrit = lin_interp(muopen, Ucritf, munew, maxj, 1879);
						//if (i == 20) printf("muopen = %f\tUcrit = %f\n", munew, Ucrit);
						//printf("munew, Ucrit = %f, %f\n", munew, Ucrit);
						//printf("j, k, upperlimit = %d, %d, %d\n", j, k, upperlimit[j]);
						sizeU = (int) sqrt(2.0*(Ucap - Uperpnew))/dvz;
						reflected = 1;
						for (l=0; l < sizeU; l++)
						{	
							if (l!=0)
							{	
								Fold = F;
								Fold_ref = F;
								vz = dvz*l;
								U = Uperpnew + 0.5*vz*vz;
								if ( (U > munew) && (U - 0.5*vz*vz + 0.5*(vz-dvz)*(vz-dvz) < munew) ) {
									frac = (vz - sqrt(2.0*(munew - Uperpnew+TINY)))/dvz;
									Fold = bilin_interp(munew, 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
									F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
								}
								else if (U > munew){
									frac = 1.0;
									F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
								}
								else {
									frac = 1.0;
									F = 0.0;
								}

								if ( (U-munew < Ucrit - numb) && (reflected == 1) ) 
									frac_reflected = 1.0;
								else if ( (U-munew > Ucrit - numb) && (reflected == 1) ) {
									reflected = 0;
									if (Ucrit < Uperpnew - munew + numb) frac_reflected = 0.0;
									else {
										vzcrit = sqrt(2.0*(Ucrit + munew - Uperpnew));
										Fold_ref = bilin_interp(munew, Ucrit, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
										frac_reflected = (vzcrit - (vz - dvz))/dvz;
									}
									//if ( i == 40 && munew < 0.7 && munew > 0.5 ) {
									//printf("frac_reflected = %f, l=%d, vzcrit = %f, vz - dvz = %f\n", frac_reflected, l, vzcrit, vz - dvz);
									//printf("U = %f, Ucrit = %f, munew = %f, Uperpnew = %f\n", U, Ucrit, munew, Uperpnew);
									//}
								}
								else { 
									frac_reflected = 0.0;
								}
								//frac_reflected = 0.0;
								intdU += 0.5*frac*dvz*((F+Fold) + frac_reflected*(F + Fold_ref));
							}
							else {	
								//vz = dvz*l;

								U = Uperpnew ;//+ 0.5*vz*vz;

								//if (phi[0] < 0.0) reflected = 0.0;
								//if (U-munew > Ucrit) 
								//	reflected = 0;

								//else printf("not here\n");
								//reflected = 1;
								F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
								intdU += 0.0;
							} 
						}
						intdUantycal = exp(-Uperpnew)*(1.0/(2.0*M_PI));// result with phi =0
						//intdU = intdUantycal;
						if (DEBUG == 1) 	
							printf("Analytical intdU is %f, numerical one is %f\n", intdUantycal, intdU);
						if (k==lowerlimit[j]) {	
							intdvx += 0.0; 
						}
						else {	
							intdvx += 2.0*0.5*dvx*(intdU+intdUold);
						}
						if (intdvx != intdvx) {	
							printf("intdvx is NAN, j=%d, i=%d\n", j, i); 
							exit(-1);
						} 
					}
					intdvxantycal = (2.0/(2.0*sqrt(M_PI)))*exp(-(xx[i]-xbar[j])*(xx[i]-xbar[j]))*erf(sqrt(xbar[j]*xbar[j]-(xx[i]-xbar[j])*(xx[i]-xbar[j])));
					if (DEBUG == 1) {
						printf("pos=%f, i=%d, j=%d, Uperp is %f, chi is %f\nvxnew and vxold are %f and %f and and vx is %f, dvx is %f\nintdvx is %f, analytical one is %f, upper is %d, upperlimit is %d\n", xx[i], i, j, Uperpnew, chi[j][i], vxnew, vxold, vx[j][i][k], dvx, intdvx, intdvxantycal, upper[j][i], upperlimit[j]); 
					}
					//printf("intdvx is %f, while the analytical one is %f\n", intdvx, intdvxantycal);
					//intdvx = intdvxantycal;
					if (j==jmclosed[i]+1 && i!=0) {
						//dxbar = (xx[imax[j]] - xx[imax[j]])*(xbar[j] - xbar[j-1])/(xx[imax[j]] - xx[i]);
						//dxbar = (chiMax[j]-chi[j][i])/(2.0*(xx[i]-xx[imax[j]]));
						if (i==imax[j]) {	
							dxbar = xbar[j] - xbar[j-1]; 
						}
						else {	
							//dxbar = xbar[j] - (xbar[j-1]*(chiMax[j] - chi[j][i]) - xbar[j]*(chiMax[j-1] - chi[j-1][i]))/(-chiMax[j-1] + chi[j-1][i] + chiMax[j] - chi[j][i]); 
							dxbar = (xbar[j] - xbar[j-1])*(chiMax[j] - chi[j][i])/ (chiMax[j] - chi[j][i] + chi[j-1][i] - chiMax[j-1]);// open orbit density does not need to be so accurate at this point
						}
						intdxbar += 0.5*(intdvx+intdvxold)*dxbar;
						if (intdxbar != intdxbar) {	
							printf("intdxbar is NAN, j=%d, i=%d\n", j, i);  
							exit(-1);
						} 
					}
					else
					{
						dxbar = xbar[j] - xbar[j-1];
						intdxbar += 0.5*(intdvx+intdvxold)*dxbar;
						if (intdxbar != intdxbar) {
							printf("intdxbar is NAN, j=%d, i=%d\n", j, i); 
							exit(-1);
						} 
					}
				} 
			} 
		}
		n_grid[ic] = intdxbar + intdxbaropen;
		if (ic == 0) {
			flux0 = intdxbaropenflow/(n_inf*alpha);
			if (charge < 0.0) 
				*flux = flux0;  
			printf("flow velocity at x=0 = %f\n", flux0/n_grid[0]);
			printf("flux evaluated at x=0 is %f\n", flux0);
		}
		if ( (ic != size_xlim) && ( ( 1.0 - n_grid[ic]/n_inf <= margin ) || ( ( (charge < 0.0) && ( - phi_grid[ic] < margin ) ) && ( - phi_grid[0] < 0.2 ) ) ) ) {	
			printf("ic = %d/%d\n", ic, size_phigrid);
			stop = 1;
			*size_ngrid = ic; 
			printf("stopping density evaluation at x = %f, density = %f*n_inf, because either the potential or the charge density perturbation are too small\n", x_grid[ic], n_grid[ic]);
		}
		//else if (xx[i] > xx[size_finegrid-1] - limit_rho) {
		//	printf("ic = %d/%d\n", ic, size_phigrid);
		//	stop = 1;
		//	*size_ngrid = ic; 
		//	printf("WARNING in densfinorb.c: stopping for positive density\n");

		//}
		if ( (DEBUG == 1) ) {	
			//printf("for charge = %f\n", charge);
			printf("%f, %f, %f is TOTAL, CLOSED and OPEN orbit density at position index %d, position %f, potential = %f, n_inf %f, ic = %d, size_xlim = %d\n", n_grid[ic], intdxbar, intdxbaropen, ic, xx[i], phi_grid[ic], n_inf, ic, size_xlim);
			if ( (n_grid[ic] != n_grid[ic]) || (n_grid[ic] < TINY) ) { 	
				printf("n_finorb[ic] = %f\n", n_grid[ic]/n_inf);
				//exit(-1);
			} 
		} 
		if (ic == size_xlim - 1) {
			printf("ic = %d/%d\n", ic, size_phigrid);
			stop = 1;
			*size_ngrid = ic; 
			printf("In densfinorb.c: reached maximum distance from wall at which can calculate density\n");
		}
		if (ic == size_xlim) {
			n_inf = n_grid[ic];
			ic = 0; //should never reach here if set ic=0 before beginning of while (stop == 0) loop.
		}
		else {
			fprintf(fout, "%f %f %f %f\n", xx[i], n_grid[ic]/n_inf, intdxbar/n_inf, intdxbaropen/n_inf);
			//fprintf(fout, "%f %f %f %f\n", xx[i], n_grid[ic], intdxbar, intdxbaropen);
			ic += 1;
		}
	} 
	printf("NINF = %f\n", n_inf);
	printf("in densfinorb: charge = %f, dndphi = %f\n", charge, (n_grid[*size_ngrid-1] - n_grid[*size_ngrid-2])/(phi_grid[*size_ngrid-1] - phi_grid[*size_ngrid-2]));
	fclose(fout);
	if (stop == 0) { 
		printf("ERROR: the density never reached stopdens*n_inf\n"); 
		exit(-1); 
	}
	for (ic = 0; ic < size_phigrid; ic++) {
		if (DEBUG == 1) 
			printf("before renormalizing n_finorb[%d] = %f, phi_grid[%d] = %f\n", ic, n_grid[ic], ic, phi_grid[ic]);
		if (ic < *size_ngrid) {
			n_grid[ic] /= n_inf;
			if (DEBUG ==1) {
				printf("x_grid[%d] = %f\tn_finorb[%d] = %f\tphi_grid[%d] = %f\n", ic, x_grid[ic], ic, n_grid[ic], ic, phi_grid[ic]);
			}
		}
		else n_grid[ic] = 0.0;
	}

	//for (j=0; j<sizemumu; j++) {
	//	vy_op[j] = lin_interp(muopen, xbaropen, mumu[j], maxj, 2700);
	//	//vy_op[j] = lin_interp(muopen, xbar, mumu[j], sizexbar, 2700);
	//	chiMax_op[j] = lin_interp(muopen, chiMopen, mumu[j], maxj, 2700);
	//	dmudvy_op[j] = lin_interp(muopen, openorbitopen, mumu[j], maxj, 2700);
	//	if (DEBUG == 1) 
	//		printf("%d/%d %f %f %f %f\n", j, sizemumu, mumu[j], vy_op[j], chiMax_op[j], dmudvy_op[j]);
	//}
	for (j=0; j<maxj/2; j++) {
		vy_op[j] = xbaropen[j];
		mu_op[j] = muopen[j];
		chiMax_op[j] = chiMopen[j];
		dmudvy_op[j] = openorbitopen[j];
		if (DEBUG == 0) 
			printf("%d/%d %f %f %f %f\n", j, sizemumu, mu_op[j], vy_op[j], chiMax_op[j], dmudvy_op[j]);
	}
	*size_op = maxj/2;

	fclose(filellip);

	// If you love your variables (and your memory) set them free // wise words Robbie
	free(chiMopen);//
	printf("1\n");
	free(openorbitopen);//
	printf("2\n");
	free(xbaropen);//
	printf("3\n");
	free(Ucritf);//
	printf("4\n");
	free(muopen);//
	printf("5\n");
	free(xx);//
	printf("6\n");
	free(jmclosed);//
	printf("7\n");
	free(jmopen);//
	printf("8\n");
	free(xifunction);//
	printf("9\n");
	free(xbar);//
	printf("10\n");
	free(chimin);//
	printf("11\n");
	free(chiMax);//
	printf("12\n");
	free(chimpp);//
	printf("13\n");
	free(crossed_min);//
	printf("13\n");
	free(crossed_max);//
	printf("13\n");
	free(openorbit);//
	printf("13\n");
	free(upperlimit);//
	printf("13\n");
	free(lowerlimit);//
	printf("14\n");
	free(xtop);//
	printf("14\n");
	free(imax);//
	printf("14\n");
	free(imin);//
	printf("14\n");

	for (w = 0; w < sizexbar; w++) {
		free(chi[w]);//
		free(Uperp[w]);//
		free(mu[w]);//
		free(upper[w]);//
		for (s = 0; s < size_finegrid; s++)
			free(vx[w][s]);//
		free(vx[w]);//
	}
	printf("14\n");
	free(chi);//
	printf("15\n");
	free(Uperp);//
	printf("15\n");
	free(mu);//
	printf("15\n");
	free(upper);//
	printf("15\n");
	free(vx);//
	printf("15\n");

	free(phi);//
	free(gg);//
	free(ff);//
	free(phip);//
	free(kdrop);//
	printf("16\n");

	for (j=0; j<j_inf; j++) {
		free(vxinf[j]);//
	}
	free(vxinf); //
	printf("17\n");
	free(chiinf); //
	printf("18\n");
	free(muinf);//

	clock_t end = clock(); // finds the end time of the computation
	double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("in densfinorb: module ran in %f seconds\n", jobtime);
	return;
}
// closed densfinorb function 

