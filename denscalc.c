// Authors: Aessandro Geraldini and Robbie Ewart
// MODIFIED 16 JAN 2022
/* This script calculates the density profiles of charged particles in the magnetised sheath for a given potential profile and entrance distribution function: denszeroorb neglects gyro-orbit size; densfinorb includes distorted gyro-orbits for magnetic field angle α<<1. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "mps.h"

double tophat(double x1, double x2, double x) {
	double y;
	if ((x1 <= x ) && (x <= x2 )) 
		y = 1.0;
	else
		y = 0.0;
	return y;
}


/* Authors: 
   Robbie Ewart 
(wrote the script makelookup.c which calculated the electron density from 1D distribution of parallel velocities, from which this function has been adapted)
   Alessandro Geraldini 
(upgraded the function makelookup into denszeroorb by including integration over mu and tweaking to alllow calculation of density of particles accelerated to wall)

MODIFIED on 15/01/2022 by Alessandro Geraldini

denszeroorb (previously makelookup) calculates the density integral of a species when this is only a function of the potential at a given point */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mps.h"

/* AG
   to use for ions in the Debye sheath, phi must be a list of POSITIVE numbers 
   equal to minus the potential relative to the Debye sheath entrance
   set number pointed at by vpar_cut_lookup to be equal to 0.0
   set size_mu = 0 while size_cut, size_vpar and mue_cut_lookup are redundant (can set to zero/NULL)
*/

void denszeroorb(double charge, double *phi_real, double *n_grid, int p_size, double *Phi_point, double **distfunc, double *vpar, double *mu, int size_vpar, int size_mu, double *mue_cut_lookup, double *vpar_cut_lookup, int size_cut) {// double power) {
//defining des variables
int count = 0;
int vi, vi_1, p, len_F; //vi is a counting variable that will be saved for sums over velocity space, vi_1 is a special value in velocity space devoted to the first point, p is a counting variable that will be saved for counting over phi space, len_F saves the number of entries in the distribution;
double *phi, *vparacc, phip = -0.000000001 ;
double Phi=0.0, n_inf = 0.0;
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
	phi[p] = -phi_real[p]*charge;
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
		vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
		//printf("in makelookup.c: mu = %f and vpar_cut_DSE = %f\n", mu[mu_ind], vpar_cut);
		vpar_cut = sqrt(vpar_cut*vpar_cut + 2.0*(-phi[0]) + TINY);
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
			for (vi = 0; vi < len_F - 2; vi++)
			{
				nepart[mu_ind] += (F[vi] + F[vi + 1]) * v_s;
			}
			//printf("nepart[%d] = %f\n", mu_ind, nepart[mu_ind]);
			count++;
		}
		else { // now there is a cut-off
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
			//printf("nepart[%d] = %f\n", mu_ind, nepart[mu_ind]);
		}
		if (mu_ind != 0) {
			n_inf += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1])*(mu[mu_ind] - mu[mu_ind-1]);
		}
		//printf("count=%d\n", count);
	}
	printf("n_inf = %f\n", n_inf);
	
//else {
	n_inf = 0.0;
	v_min = sqrt(-2.0 * phip);
	for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
		nepart[mu_ind] = 0.0;
		vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
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
				else //regualar case
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


	// ELECTRON CURRENT
	for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
		nepart[mu_ind] = 0.0;
		vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
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
				vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
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
						else //regualar case
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
				if (mu_ind != 0) {
					n_grid[p] += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1])*(mu[mu_ind] - mu[mu_ind-1]);
					//printf("n_grid %f\n", n_grid[p]);
				}
			}
		}
		//n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] -phi[0] )));
		//printf("n_grid[%d] = %f, n_pre = %f\n", p, n_grid[p], n_pre[p]);
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
		for (vi = 0; vi < len_F - 1; vi++) {
			vparacc[vi] = sqrt(vpar[vi]*vpar[vi] + 2.0*phi[p]);
			if (vi != 0)
				n_grid[p] += 0.5*(F[vi] + F[vi-1])*(vparacc[vi] - vparacc[vi-1]);
		}
	}
}
	
//else {
//	for (vi = 0; vi < len_F - 1; vi++) {
//		F[vi]     = distfunc[0][vi];
//		Fp[vi]    = ddistdvpar[0][vi];
//		Fpp[vi]   = ddistdvpartwo[0][vi];
//		//printf("F = %f, Fp = %f, Fpp = %f\n", F[vi], Fp[vi], Fpp[vi]);
//	}
//	vpar_cut = *vpar_cut_lookup;
//	if (phi[0] < 0.0) { // the case for electrons in MPS, not the case for ions in DS
//		v_cut = sqrt(vpar_cut*vpar_cut + 2.0*(-phi[0]));
//	}
//	else v_cut = 0.0;
//
//
//
//	// Normalization at infinity
//	n_inf = 0.0;
//	v_min = 0.0;
//
//	if (phi[0] > 0.0) {
//		n_inf = 0.0;	
//		for (vi = 0; vi < len_F - 2; vi++)
//		{
//			n_inf += 0.5*(F[vi] + F[vi + 1]) * v_s;
//		}
//	}
//	else {
//	if (v_cut >= v_max) { //effectively no cut off
//		v_min = 0.0;
//		if ((int)floor(v_min / v_s) >= len_F - 1)
//		{
//			n_inf = 0;
//			printf("v_min outside velocity grid");
//		}
//		else
//		{
//			if (fabs(v_min) < TINY)//special case where the integration becomes simple
//			{
//				for (vi = 0; vi < len_F - 2; vi++)
//				{
//					n_inf += (F[vi] + F[vi + 1]) * v_s;
//				}
//			}
//			else //regualar case
//			{
//				vi_1 = (int)floor(v_min / v_s);
//				sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0));
//				theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0)) / v_min);
//
//				n_inf += 2.0 * (((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//				n_inf += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//
//				w = ((0.5 * pow((((vi_1 + 1) * v_s) - v_min), 2.0)) + (v_min * (((vi_1 + 1) * v_s) - v_min)));
//
//				for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//				{
//					sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0));
//					sqrt_lo = sqrt(pow(((vi)* v_s), 2.0));
//					theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0)) / (v_min));
//					theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0)) / (v_min));
//
//					n_inf += 2.0 * (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//					n_inf += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//					w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//					w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//				}
//			}
//		}
//	}
//	else {
//		// for ions in DS always end up here
//		//now instead there is a cutoff
//			v_min = 0.0;
//			vi_1 = 0;
//			sqrt_up = ((vi_1 + 1) * v_s);
//
//			n_inf += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
//			n_inf += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));
//			w = 0.5 * pow(v_s, 2.0);
//
//			for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//			{
//				sqrt_up = (vi + 1) * v_s;
//				sqrt_lo = (vi * v_s);
//
//				n_inf += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
//				n_inf += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//
//				w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//				w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
//			}
//
//			for (vi = vi_1 + 1; vi < (int)floor(v_cut / v_s); vi++)
//			{
//				sqrt_up = (vi + 1) * v_s;
//				sqrt_lo = (vi * v_s);
//
//				n_inf += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
//				n_inf += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//
//				w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//				w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
//			}
//			vi_1 = (int)floor(v_cut / v_s);
//			sqrt_up = v_cut;
//			sqrt_lo = vi_1 * v_s;
//
//			n_inf += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
//			n_inf += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
//
//			w_up = 0.5 * pow(v_cut, 2.0);
//			w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
//	}
//	}
//
//	//printf("n_inf = %f\n", n_inf);
//
//	// ELECTRON CURRENT
//	// fill in
//
//
//	if (phi[0] > 0.0) {
//		for (p = 0;p < p_size; p++) {
//			n_grid[p] = 0.0;
//			for (vi = 0; vi < len_F - 2; vi++) {
//				n_grid[p] += 0.5*(F[vi] + F[vi + 1]) * ( sqrt((vi+1)*v_s*(vi+1)*v_s + 2.0*phi[p]) - sqrt((vi)*v_s*(vi)*v_s + 2.0*phi[p]) );
//			}
//		}
//	}
//	else {
//	if (v_cut >= v_max) { //effectively no cut off
//		for (p = 0; p < p_size; p++) {
//			n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] + (0.5 * pow(v_cut, 2.0)))));
//			v_min = sqrt(-2.0 * phi[p]);
//			//printf("v_min = %f\n", v_min);
//			if ((int)floor(v_min / v_s) >= len_F - 1)
//			{
//				n_grid[p] = 0.0;
//				printf("v_min outside velocity grid");
//			}
//			else
//			{
//				n_grid[p] = 0.0;
//				if (fabs(phi[p]) < TINY) {
//					n_grid[p] = n_inf;
//				}
//				else if (fabs(v_min) < TINY)//special case where the integration becomes simple
//				{
//					for (vi = 0; vi < len_F - 2; vi++)
//					{
//						n_grid[p] += (F[vi] + F[vi + 1]) * v_s;
//					}
//				}
//				else //regualar case
//				{
//					vi_1 = (int)floor(v_min / v_s);
//					sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p]));
//					theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p])) / v_min);
//
//					n_grid[p] += 2.0 * (((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//					n_grid[p] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//
//					w = ((0.5 * pow((((vi_1 + 1) * v_s) - v_min), 2.0)) + (v_min * (((vi_1 + 1) * v_s) - v_min)));
//					in_err[p] += 2.0 * ((((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w))));
//					in_err[p] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//
//					for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//					{
//						sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
//						sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
//						theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//						theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//
//						n_grid[p] += 2.0 * (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//						n_grid[p] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//						w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//						af_err[p] += 2.0 * ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
//						af_err[p] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//					}
//				}
//
//			}
//			n_res[p] = n_grid[p] - n_pre[p];
//			//printf("n_grid[%d/%d] = %f\n", p, p_size, n_grid[p]);
//		}
//	}
//	else {
//		// for ions in DS always end up here
//		//now instead there is a cutoff
//		for (p = 0;p < p_size; p++)
//		{
//			n_grid[p] = 0.0;
//			if ((phi[p] + (0.5 * pow(v_cut, 2.0))) > 0)
//			{
//				n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] + (0.5 * pow(v_cut, 2.0)))));
//				if (DEBUG == 1)
//					printf("n_pre[%d] = %f\n", p, n_pre[p]);
//			}
//			else
//			{
//				n_pre[p] = 0.5 * exp(phi[p]);
//			}
//			//if (phi[0] < 0.0) {
//			v_min = sqrt(-2.0 * phi[p]);
//			//}
//			//else { 
//			//	v_min = 2.0*TINY;
//			//}
//			if ((int)floor(v_min / v_s) >= len_F - 1)
//			{
//				n_grid[p] = 0.0;
//			}
//			else
//			{
//				if (fabs(v_min) < TINY)
//				{
//					vi_1 = 0;
//					sqrt_up = ((vi_1 + 1) * v_s);
//
//					n_grid[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
//					n_grid[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));
//					w = 0.5 * pow(v_s, 2.0);
//					af_err[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s)))) - erf(sqrt(w));
//					af_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));
//
//					for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//					{
//						sqrt_up = (vi + 1) * v_s;
//						sqrt_lo = (vi * v_s);
//
//						n_grid[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
//						n_grid[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//
//						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//						w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
//						af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
//						af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//					}
//
//					for (vi = vi_1 + 1; vi < (int)floor(v_cut / v_s); vi++)
//					{
//						sqrt_up = (vi + 1) * v_s;
//						sqrt_lo = (vi * v_s);
//
//						n_grid[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
//						n_grid[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//
//						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//						w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
//						af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
//						af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//					}
//					vi_1 = (int)floor(v_cut / v_s);
//					sqrt_up = v_cut;
//					sqrt_lo = vi_1 * v_s;
//
//					n_grid[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
//					n_grid[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
//
//					w_up = 0.5 * pow(v_cut, 2.0);
//					w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
//					af_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
//					af_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
//				}
//				else
//				{
//					//printf("am I always here for %d/%d\n", p, p_size);
//					if ((int)floor(v_min / v_s) == (int)floor(v_cut / v_s))
//					{
//						vi_1 = (int)floor(v_min / v_s);
//						sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p]));
//						theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p])) / (v_min));
//
//						n_grid[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//						n_grid[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));
//
//						w_up = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
//						in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w_up))));
//						in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));
//
//						sqrt_cut = sqrt(pow(v_cut, 2.0) + 2.0 * phi[p]);
//						theta_cut = asinh(sqrt(pow(v_cut, 2.0) + 2.0 * phi[p]) / (v_min));
//
//						n_grid[p] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * v_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
//						n_grid[p] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (v_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
//
//						w_lo = 0.5 * (pow(v_cut, 2.0) - pow(v_min, 2.0));
//						in_err[p] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * v_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
//						in_err[p] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (v_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
//
//						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//						{
//							sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
//							sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
//							theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//							theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//
//							n_grid[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//							n_grid[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//							w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//							af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
//							af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//						}
//					}
//					else //the most likely case
//					{
//						vi_1 = (int)floor(v_min / v_s);
//						sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phi[p]);
//						theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//
//						n_grid[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//						n_grid[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));
//
//						w = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
//						in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w))));
//						in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));
//
//						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//						{
//							sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
//							sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
//							theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//							theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//
//							n_grid[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//							n_grid[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//							w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//							af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
//							af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//						}
//
//						for (vi = vi_1 + 1; vi < (int)floor(v_cut / v_s); vi++)
//						{
//							sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
//							sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
//							theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//							theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//
//							n_grid[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//							n_grid[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//							w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//							af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
//							af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//						}
//
//						vi_1 = (int)floor(v_cut / v_s);
//						sqrt_up = sqrt(pow(v_cut, 2.0) + 2.0 * phi[p]);
//						sqrt_lo = sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phi[p]);
//						theta_up = asinh(sqrt(pow(v_cut, 2.0) + 2.0 * phi[p]) / (v_min));
//						theta_lo = asinh(sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//
//						n_grid[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//						n_grid[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));
//
//						w_up = 0.5 * (pow(v_cut, 2.0) - pow(v_min, 2.0));
//						w_lo = 0.5 * (pow((vi_1 * v_s), 2.0) - pow(v_min, 2.0));
//						in_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
//						in_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));
//
//
//
//					}
//				}
//			}
//			n_res[p] = n_grid[p] - n_pre[p];
//			printf("n_grid[%d] = %f, n_pre = %f\n", p, n_grid[p], n_pre[p]);
//		}
//	}
//	}
//}


Phi =  Phi / n_inf;
//printf("Phi = %f\n", Phi);
printf("n_inf = %f\n", n_inf);
for (p = 0;p < p_size; p++)
{
	//n_grid[p] = n_grid[p] / n_grid[p_size - 1];
	n_grid[p] = n_grid[p] / n_inf;
	//if (n_grid[p] != n_grid[p]) n_grid[p] = 1.0;
	n_res[p] = n_grid[p] - n_pre[p];
	//printf("n_makelookup[%d/%d] = %f, phi = %f\n", p, p_size, n_grid[p], phi[p]);
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

printf("end of makelookup\n");
return;
}


void densfinorb(int output_type, int dim_output, double Te, double alpha, int size_phigrid, int *size_ngrid, double* n_grid, double *x_grid, double* phi_grid, double charge, double **FF, double *mumu, double *UU, int sizemumu, int sizeUU, double grid_parameter, double *vx_i_DS, double *fx_i_DS, double *flux, int zoomfactor, double stopdens) {

// DECLARATIONS
//inpu: double rho_over_unitgrid = for ions in the magnetic presheath this is 1.0, for electrons in the Debye sheath this can be chosen as min(1.0, rho_e/lambda_D)
//double temp_over_unitgrid = 1.0;
clock_t begin = clock(); // Finds the start time of the computation
double limit_rho = 6.0;
double n_inf=0.0, corr;
double phi_cut = -99999.0; // keep phi_cut as large and negative as possible
int reflected = 0;
double Ucrit = 0.0, frac_reflected=0.0, frac = 0.0, vzcrit;
double deltax, deltax_inf;
double *chiinf, *muinf, **vxinf;
double *xx, pos, *phi, *phip, *niclosed, *niopen, *nitotal, **chi, **chiop, Telarge;
/* pos is x, the position (distance from the wall); phi contains the values of phi(x), extracted from a file. phip is phi prime, first derivative of phi. ne is the electron density; niclosed is the closed orbit ion density, niopen is the open orbit ion density, nitotal is the sum of the two; newphi is the new electrostatic potential guess; chi is the effective potential;*/ 
int stop = 0, sizexbar, maxj, j_inf;
int icrit, ind;
/* n is the domain of the position x (the largest value of the index i, so e.g. x_grid[n] = L_2 in the paper); sizexbar is the domain of the position xbar (the largest value of the index j); */ 
double *openorbit, openorbitantycal, **mu, *muopen, *Ucritf, *xbar, xbarcrit, *twopimuprime, *openintegral, *xifunction;
/* openorbit is the integral Delta_M in the paper (an integral over the last closed orbit, representing how much U_perp - chi_M change in one full last orbit); openorbitantycal is the analytical value of openorbit for a flat potential (the first potential guess usually); mu is the array containing values of mu(xbar, Uperp). The index j is for values of xbar, k is for values of Uperp; xbar is the grid of values used in the closed orbit integral; FF contains the distribution function, read from the file distfile.txt. UU and mumu contain the values of U and mu corresponding to the function FF (which is F(mu, U)); FFprime is the numerical first derivative of F with respect to U; Many of these arrays are pointers because we don't initially know the array size, until we read the files (e.g. the file with the potential and the one with the distribution function); Note: first index of FF and FFprime is mu, second one is U; Interesting output variables; xifunction is the function of x which defines the grid of values of xbar by finding a chi whose minimum lies exactly at each grid point x */
int *jmclosed, *jmopen,  i=0, ic=0, j=0, k=0, l=0, sizeU, size_finegrid;
/* jmclosed represent minimum values of xbar above which we integrate open and closed orbit density integrals respectively (xbar_m,o and xbar_m in the paper); i is an index usually representing the positin x  j is an index usually representing the orbit position xbar; k is an index usually representing the energy Uperp (or velocity vx). It's always used in conjunction with j (and sometimes i); l is an index (used in for loops) usually representing the total energy (or velocity vz). It's used only in the DENSITY INTEGRALS part of the code; sizeU is the size of the integration range over U (or velocity vz). It is set later on in the code */
int *crossed_max, *crossed_min, *kdrop;
/* crossed_max and crossed_min is non-zero when a minimum (or maximum) of chi is found for some xbar[j]; kdrop is an integer which is non-zero only if there is an additional potential barrier for the particle at x << rho e.g. for an ion if the sheath reverses Uperp is allowed to be above chiMax and therefore the k = 0 of Uperp[j][k] starts for Uperp[j][0] = chiMax[j] + phi_barrier instead of Uperp[j][0] = chiMax[j] which is the conventional way */
int *lowerlimit, *upperlimit, **upper, *itop, *imax, *imin;
/* lowerlimit represents the lower limit of k in the integrals over Uperp (or vx). It's needed because some of the earlies energies; (which are the largest because thy are values of chi stored after the maximum is found); may be so large that they are associated with very small values of the distribution function. This avoids integrating in an empty portion of phase space; upperlimit[j] represents the largest value of k (the smallest stored energy Uperp = chi_minimum) associated with some value of j; upper[j][i] represents the value of k associated with the smallest value of vx when integrating over Uperp. Going above upperlimit[j][i] makes Uperp < chi so velocities imaginary; imax/imin[j] stores the position of the maximum/minimum of the effective potential chi (It's x_M/x_m in the paper, which depends on xbar). */
double **Uperp, ***vx, *chiMax, *chimpp, *chimin, oorbintgrd, oorbintgrdantycal, oorbintgrdBohm, intdvxopen, intdvxopenBohm=0.0, intdvxopenfluidBohm=0.0, intdxbaropenold;
/* Uperp stores the possible values of Uperp associated with closed orbits, and so does vx; chiMax and chimin store the local maxima and minima of the effective potential maximum, oorbintgrd is the value of the integrand in the first open orbit integral (oorbintgrdantycal is the analytical result for flat potential); oorbintgrdBohm is the value of the integrand in the `Bohm' integral; intdvxopen and similars are where the integral of f_{0x} (v_x) over v_x, and its two important moments <v_x> (fluidBohm) and <1/v_x^2> (Bohm) are stored; They are a check that the extracted distribution function f_{0x} has the same moments it had when we carried out the density integral earlier; */
double vz, U, dvz = 0.2, dvzopen = 0.2, dvx, dvxopen, Deltavx, vxopen, dxbar, intdU=0.0, intdUopen=0.0, intdUopenBohm = 0.0;
	/* vz used in the density integral; U is the total energy, used in the density integral; dvz is the thickness of the vz grid used to take the integral over U (which is taken over vz in practice), dvzopen is the same for the open orbit piece; dvx is the thickness of the vx grid used to take the integral over Uperp ( which is taken over vx in practice). It must be evaluated because it depends on stored values of vx[j][i][k]; dvxopen is the thickness of the vx grid on which f_{0x} (v_x) is defined (and integrated to check consistency of its moments); dxbar is the thickness of the xbar grid; intdU is the value of the integral over U in the closed orbit density integration process; intdUopen is the same as above, for the open orbit integral; intdUopenBohm same, for Bohm integral */
double intdUold=0.0, intdvx=0.0, intdvxold = 0.0, intdxbar=0.0, intdxbaropen=0.0, intdxbaropenBohm = 0.0, Bohm = 0.0, F, Fold=0.0, Fold_ref=0.0, Ucap;
	/* intdUold is a variable which stores the old intdU, so that the trapezium rule of integration can be applied (intdUold + intdU)*dvz; intdvx stores the integral over Uperp (hence over vx) in the closed orbit integral; intdxbar stores the value of the integral over xbar (which is the final result!), intdxbaropen does the same in the open orbit density integral; intdxbaropenBohm does the same for the Bohm integral; idealBohm is what the Bohm integral shoult be if Bohm condition is marginally satisfied; Bohm is the Bohm integral at the Debye sheath entrance x=0; F is the value of the distribution function evaluated in the density integrals by interpolating FF, and Fold is the `old' needed to apply the trapezium rule; Fprime is the bilinearly interpolated value of FFprime, and Fprimeold is the same at the previous grid point (needed for trapezium rule); used in INTEGRALS OF DISTRIBUTION FUNCTION AT INFINITY; Ucap is the topmost total energy integrated to */
double intdUopenflow = 0.0, intdUopenflowold = 0.0, intdxbaropenflow = 0.0, oorbintgrdflow = 0.0, oorbintgrdflowold = 0.0;
// values of various integrals
double oorbintgrdold=0.0, oorbintgrdBohmold=0.0, Fopen=0.0, intdUopenold=0.0, intdUopenBohmold=0.0;
double *chiprimetop, vx0open, aa; 
double intdUantycal=0.0, intdvxantycal=0.0, vxnew=0.0, vxold = 0.0, Uperpnew = 0.0, *xtop, intdUopenantycal=0.0;
double openorbitnew, chinew, munew = 0.0;
/* intdUantycal is the integral over U (or v_z) for a flat potential profile (phi =0) for some value of xbar and Uperp; intdvxantycal  is the integral over Uperp (or vx) for a flat potential profile for some value of xbar; vxnew is the value of vx at the 'new' grid point, used in the vx integral (taken using the trapezium rule); vxold is the value of vx at the 'old' grid point, used in the vx integral; Uperpnew is the value of Uperp (used in the closed orbit density integral); munew is the valye of mu (used in the closed orbit density integral); xtop is the top bounce point x_t of the last closed orbit; intdUopenantycal is the analytical value of the integral over U  in the open orbit density integral */
double oorbintgrdxbar, oorbintgrdsquare, intdUopenxbar = 0.0, intdUopensquare = 0.0, intdxbarxbar=0.0, intdxbarsquare=0.0;// kBohmsquared;
/* oorbintgrdxbar is an integral over the open orbit distribution function at x=0 which is needed to evaluate a coefficient that appears when; expanding quasineutrality near x=0. It was just for playing around and at the moment plays no role in the code; similarly with all other integrals here */
double oorbintgrdxbarold, oorbintgrdsquareold, intdUopenxbarold, intdUopensquareold;
// all used to take integrals above
double xi, *gg, *ff;
double du, fluxinf1old, fluxinfintgrdold, fluxinfintgrd, u, Fprimeold, Fprime;
double **distfunc_iprime, fluxinf, fluxinf1, densinf1, densinf, densinf1old;

//printf("charge = %f\n", charge);
//if (charge < 0.0) dvz = 0.2;
//else dvz = 0.2;
//printf("zoomfactor = %d\n", zoomfactor);

if (DEBUG == 1) {
	for (i=0; i<size_phigrid;i++) {
		printf("index %d\tx = %f\tphi = %f\n", i, x_grid[i], phi_grid[i]);
	}
}

gg = malloc(size_phigrid*sizeof(double));
ff = malloc(size_phigrid*sizeof(double));
for (i=0; i<size_phigrid; i++) {
	gg[i] = sqrt(x_grid[i]);
	ff[i] = pow(sqrt(grid_parameter)+gg[i], 2.0) - grid_parameter;
	//printf("ff[%d] = %f\n", i, ff[i]);
}
deltax = ff[1];
//printf("ff[%d] = %f\n", size_phigrid-1, ff[size_phigrid-1]);
size_finegrid = (int) (zoomfactor*ff[size_phigrid-1]/deltax) ;
//printf("size_finegrid = %d\n", size_finegrid);
xx = (double*)calloc(size_finegrid,sizeof(double)); // xx = x has correct size
phi = (double*)calloc(size_finegrid,sizeof(double)); // phi now has correct size
FILE *fp;

if (zoomfactor != 1) {
	i=0;
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, size_phigrid);
	//gsl_spline *spline = gsl_spline_alloc (gsl_interp_polynomial, size_phigrid);

	//gsl_spline_init (spline, gg, phi_grid, size_phigrid);
	gsl_spline_init (spline, ff, phi_grid, size_phigrid);
	//gsl_spline_init (spline, phi_grid, x_grid, size_phigrid);
	//gsl_spline_init (spline, x_grid, phi_grid, size_phigrid);

	fp = fopen("OUTPUT/phispline.txt", "w");
	if (fp == NULL) printf("Error: phispline not created\n");
	//printf ("%f\n", gg[size_phigrid-1]);
	for (i = 0; i < size_finegrid; i += 1)
	{
		xi = i*deltax/zoomfactor;
		//printf("i=%d/%d\n", i, size_finegrid);
		//deltax = phi_grid[(i+1)/zoomfactor] - phi_grid[i%zoomfactor];
		phi[i] = phi_grid[i/zoomfactor] + (i%zoomfactor)*deltax/zoomfactor;
		if (i == 0) xi += TINY;
		//if (i == size_finegrid-1) xi -= TINY;
		xx[i] = pow( pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), 2.0);
		phi[i] = gsl_spline_eval (spline, xi, acc);
		//xx[i] = gsl_spline_eval (spline, phi[i], acc);
		//fprintf(fp, "%f %f\n", pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), phi[i]);
		fprintf(fp, "%f %f\n", xx[i], phi[i]);
		phi[i] *= charge;
		//xx[i] = xi*xi/rho_over_unitgrid;
	}
	fclose(fp);
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);

}
else {
  //xx = (double*)calloc(size_finegrid,sizeof(double)); // xx = x has correct size
  //phi = (double*)calloc(size_finegrid,sizeof(double)); // phi now has correct size
	i=0;
	fp = fopen("OUTPUT/phispline.txt", "w");
	if (fp == NULL) printf("Error: phispline not created\n");
	//printf ("%f\n", gg[size_phigrid-1]);
	for (i = 0; i < size_finegrid; i += 1) {
		xi = i*deltax/zoomfactor;
		//if (i == 0) xi += TINY;
		//if (i == size_finegrid-1) xi -= TINY;
		//xx[i] = pow( pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), 2.0);
		xx[i] = x_grid[i];
		//printf("xx = %f\n", xx[i]);
		//printf("gg = %f %f\n", xi, gg[i/zoomfactor]);
		//printf("phi = %f\n", phi[i/zoomfactor]);
		phi[i] = phi_grid[i];
		//fprintf(fp, "%f %f\n", pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), phi[i]);
		fprintf(fp, "%f %f\n", xx[i], phi[i]);
		phi[i] *= charge;
		//xx[i] = xi*xi/rho_over_unitgrid;
	}
	fclose(fp);
}

printf("phi[0] = %f\n", phi[0]);


// Introduce a number that equals Te when Te > 1 and 1 when Te<1. When Te is large, this number increases the number of grid points in tot. energy U to account for the fact that vz ~ vB > v_t,i when the electron temperature is large. Moreover, |v| ~ vB > v_t,i at the Debye sheath entrance, and so the number of grid points in xbar must also be increased.
if (Te>1.0) {	
	Telarge = Te;
}
else {
	Telarge = 1.0;
}
//printf("Telarge = %f\n", Telarge);
Ucap = 18.0 + 6.0*Te;



if (phi[0] <= 0.0) {
	distfunc_iprime = calloc(sizemumu,sizeof(double*));

	// evaluate derivative of F with respect to total energy U at U = mu
	for (i=0; i<sizemumu; i++) {
		//printf("i=%d/%d\n", i, sizemumu);
		distfunc_iprime[i] = calloc(sizeUU,sizeof(double));
		for (j=0; j < sizeUU; j++) {
			if ((j==sizeUU-1) || (i==sizemumu-1))
				distfunc_iprime[i][j] = 0.0;	
			else
				distfunc_iprime[i][j] = (FF[i][j+1] - FF[i][j])/(UU[j+1] - UU[j]); 
		}
	}

	fluxinf1old = fluxinf1 = fluxinf = 0.0;
	densinf1old = densinf1 = densinf = 0.0;
	du = 0.01;
	for (i=0; i<sizemumu; i++) {
		//if (i!=0) {
		fluxinf1old = fluxinf1;
		fluxinf1 = 0.0;
		densinf1old = densinf1;
		densinf1 = 0.0;
		//}
		Fprimeold = Fprime = 0.0;
		Fold = F = 0.0;
		fluxinfintgrdold = fluxinfintgrd = 0.0;
		for (j=0; j< sizeUU; j++) {
			u = sqrt(UU[j]);
			fluxinfintgrdold = fluxinfintgrd;
			Fprimeold = Fprime;
			Fold = F;
			Fprime = distfunc_iprime[i][j];
			F = FF[i][j];
			//Fprime = bilin_interp(mumu[i], U, distfunc_iprime, mumu, UU, sizemumu, sizeUU, -1, -1);
			//F = bilin_interp(mumu[i], U, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
			fluxinfintgrd = F*u;
			if (j != 0) {
				du = u - sqrt(UU[j-1]); 
				fluxinf1 += 0.5*(fluxinfintgrd + fluxinfintgrdold)*du;
				densinf1 += 0.5*(F + Fold)*du;
			}
		}
		if (i!=0) {
			fluxinf += 0.5*(mumu[i]-mumu[i-1])*(fluxinf1 + fluxinf1old);
			densinf += 0.5*(mumu[i]-mumu[i-1])*(densinf1 + densinf1old);
		}
	}
	densinf *= (4.0*M_PI);
	fluxinf *= (4.0*M_PI);
	*flux = fluxinf;
	printf("densinf = %f\tfluxinf = %f\n", densinf, *flux);
}


i=0;
/* We initialize all arrays that contain functions of position x with the correct size n */
niclosed = (double*)calloc(size_phigrid,sizeof(double)); // array containing the density of closed orbits
niopen = (double*)calloc(size_phigrid,sizeof(double)); // array containing the density of open orbits
nitotal = (double*)calloc(size_phigrid,sizeof(double)); // array containing the total density of ions 
phip = (double*)calloc(size_finegrid,sizeof(double)); // phi now has correct size
jmclosed = (int*)calloc(size_finegrid,sizeof(int));
jmopen = (int*)calloc(size_finegrid,sizeof(int));
xifunction = (double*)calloc(size_finegrid,sizeof(double));

/* FORM XBAR GRIDS 
Take derivatives of phi and use them to obtain two grids for xbar, one to be used for closed orbits and one to be used for open orbits. */
xbar = (double*)calloc(size_finegrid,sizeof(double));	
//xbarop = (double*)calloc(n,sizeof(double)); // Open orbit grid defined such that every position x has a corresponding chi for which xM = x. This allows better resolution of the places where the open orbit integrand, oorbintgrd, is large (~alpha^(1/2))	
for (i=0; i<size_finegrid; i++)
{	
	//if (i!=0)	xiprime[i-1] = (xiprime[i]-xiprime[i-1])/(xx[i]-xx[i-1]); 
	// Evaluate derivative of phi
	jmopen[i] = jmclosed[i] = 0;
	if (i == 0)
	{	
		phip[0] = (phi[1] - phi[0])/(xx[1]-xx[0]);
		//phip[0] = (phi[1] - phi[0])/(xx[1]-xx[0]) + (xx[0] - xx[1])*((phi[2] - phi[1])/(xx[2]-xx[1]) - (phi[1] - phi[0])/(xx[1]-xx[0]))/(xx[2]-xx[1]) ; 
		//printf("phip[0] = %f\n", phip[0]);
	}
	else if (i == size_finegrid-1)
	{	
		//phip[i] = 0.5*(phi[i] - phi[i-1])/(xx[i] - xx[i-1]); 
		phip[i] = phip[i-1];
		//printf("phip[%d] = %f\n", i, phip[i]);
	}
		//phip[0] = phip[1] + (xx[0] - xx[1])*(phip[2] - phip[1])/(xx[2]-xx[1]); 
	else {	
		//phip[i] = (phi[i+1] - phi[i])/(xx[i+1] - xx[i]);
		phip[i] = ((xx[i] - xx[i-1])/(xx[i+1]- xx[i-1]))*(phi[i+1] - phi[i])/(xx[i+1] - xx[i]) + ((xx[i+1] - xx[i])/(xx[i+1]- xx[i-1]))*(phi[i] - phi[i-1])/(xx[i] - xx[i-1]); 
		//printf("phip[%d] = %f\n", i, phip[i]);
		//printf("phip[%d] = %f\n", i, phip[i]);
	}
	// xifunction is xbar corresponding to given position x being a stationary point
	xifunction[i] = xx[i] + phip[i]/2.0;
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
			if (i<size_finegrid/2) {
				//exit(-1);
			}
		}
	} 
}
//printf("icrit = %d\n", icrit);

// Whole section below perhaps was overkill in xbar resolution
////////
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
j=0;
for (k=icrit+1; k < size_finegrid; k++) {
	//if (j!=0) xbar[2*j-1] = xbar[2*j-2] + 0.5*(xifunction[k]-xbar[2*j-2]); 
	//xbar[2*j] = xifunction[k]; 
	xbar[j] = xifunction[k]; 
	//if (xbar[j] > xx[size_finegrid-1] - 0.5*limit_rho) maxj = j;
	j++;
}
xbarcrit = xifunction[icrit];
sizexbar = j;
//sizexbar = 2*j-1;
if (DEBUG == 1) {
	for (j = 0; j < sizexbar; j++)
		printf("xbar[%d] = %f\n", j, xbar[j]); 
}	
//printf("sizexbar = %d, size_finegrid (x) =%d\n", sizexbar, size_finegrid);
chi = (double**)calloc(sizexbar,sizeof(double*)); // chi(xbar, x) indices j and i 
chiop = (double**)calloc(sizexbar,sizeof(double*)); // chi(xbar, x) indices j and i 
Uperp = (double**)calloc(sizexbar,sizeof(double*)); 
mu = (double**)calloc(sizexbar,sizeof(double*)); 
muopen = (double*)calloc(sizexbar+1,sizeof(double));	
Ucritf = (double*)calloc(sizexbar+1,sizeof(double));	
vx = (double***)calloc(sizexbar,sizeof(double**)); 
upper = (int**)calloc(sizexbar,sizeof(int*)); 
chimin = (double*)calloc(sizexbar,sizeof(double)); 
chiMax = (double*)calloc(sizexbar,sizeof(double));
chimpp = (double*)calloc(sizexbar,sizeof(double));
crossed_min = (int*)calloc(sizexbar,sizeof(int)); 
crossed_max = (int*)calloc(sizexbar,sizeof(int));
kdrop = (int*)calloc(sizexbar,sizeof(int));
twopimuprime = (double*)calloc(sizexbar,sizeof(double));
openintegral = (double*)calloc(sizexbar,sizeof(double));
openorbit = (double*)calloc(sizexbar,sizeof(double));
upperlimit = (int*)calloc(sizexbar,sizeof(int));
lowerlimit = (int*)calloc(sizexbar,sizeof(int));
chiprimetop = (double*)calloc(sizexbar,sizeof(double));
xtop = (double*)calloc(sizexbar,sizeof(double));
itop = (int*)calloc(sizexbar,sizeof(int));
imax = (int*)calloc(sizexbar,sizeof(int));
imin = (int*)calloc(sizexbar,sizeof(int));

/////////////////////////////////////////////////////
/* CLOSED ORBIT ARRAY FILLING
Set up the grid in xbar and also initialize all arrays that contain a different number at different values of xbar, indexed j */
for (j=0;j<sizexbar;j++)  {
	itop[j] = 0;
	imax[j] = imin[j] = -1;
	openorbit[j] = 0.0;
	crossed_min[j] = 0;
	crossed_max[j] = 0;
	xtop[j] = 0.0;
	chiprimetop[j] = 0.0;
	chimin[j] = 0.0;
	chiMax[j] = 0.0;
	upperlimit[j] = -1;
	lowerlimit[j] = 0; 
}
/* Loop below initializes all 2d arrays which are functions of xbar and x. It allocates the right amount of memory to arrays of pointers of size xbar. The result is a 2D array indexed j (size sizexbar) and i or k (size n, see above) */
for (j=0; j < sizexbar; j++) {
	chi[j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing chi(xbar, x) indexes j and i 
	chiop[j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing chi(xbar, x) indexes j and i 
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
	chi[j][0] =  pow((xx[0] - xbar[j]), 2.0) + phi[0];
}
/* This for loop fills in the arrays, calculating the integrals where necessary */
for (i=1; i<size_finegrid; i++) {
	for (j=0; j<sizexbar; j++)
	{	
		chi[j][i] =  pow((xx[i] - xbar[j]), 2.0) + phi[i];
//chiop[j][i] =  pow((xx[i] - xbarop[j]), 2.0) + phi[i];
/* 
FINDING MAXIMA/MINIMA
*/
/* Below, we use the array elements to define arrays for the effective potential maxima and minima that exist for every xbar (index j). We search through the chi curve from x=0 to the largest value of x but we search through every chi curve first (so first scan in xbar at fixed x, then move to the next x). We create arrays of mu, Uperp and vx.
Because I am comparing neighbouring values of x to find a maximum, I need to consider two distinct cases.The first one is treated in the if loop below, which the program should enter only if chi is decreasing at x=0. Implying that chi(x=0) is an effective potential maximum. The second one is the else if loop after that, which finds maxima of chi that are stationary points by comparing the value of the function before and after each point. Note that because we compare the function at a point with the function at points before and after, the point we consider at every iteration step in the index i is indexed (i-1), compared with (i-2) and i.  */
		if ( (i==1) && (chi[j][i] < chi[j][i-1]) )
		{	
			//if (phi_cut > chiMax[j]) {
		//		kdrop[j] = (int)((phi_cut - chiMax[j])/0.2); 
		//		printf("ever here?\n");
		//	}
		//	else kdrop[j] = 0;
			crossed_max[j] += 1;		
			imax[j] = 0;
			chiMax[j] = chi[j][i-1];
			Uperp[j][kdrop[j]] = chi[j][0];
			mu[j][0] = 0.0;
			vx[j][0][0] = 0.0; 
			//printf("charge = %f\n", charge);
			//printf("imax[%d]  = %d\n", j, imax[j]);
		} 
		else if ( (i > 1) && ((chi[j][i] < chi[j][i-1]) && (chi[j][i-1] > chi[j][i-2])) )
		{
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
				exit(EXIT_FAILURE);
			}
			imax[j] = i-1; // Store the index of the maximum 
			chiMax[j] = chi[j][i-1]; //Store chi Maximum itself, a function of xbar (index j)
			//if (phi_cut > chiMax[j]) {
			//	kdrop[j] = (int)((phi_cut - chiMax[j])/0.2); 
			//	printf("ever here?\n");
			//}
			//else kdrop[j] = 0;
		}
/* We store the index corresponding to the position of a minimum for a given value of xbar.*/
		//watch out HERE
		else if ( (i > 1) && ((chi[j][i] > chi[j][i-1]) && (chi[j][i-1] < chi[j][i-2])) )
		{	
			crossed_min[j] += 1;
			//printf("i-1 = %d, imax[%d] = %d, crossed_min[%d] = %d, crossed_max[%d] = %d\n", i-1, j,imax[j], j, crossed_min[j], j, crossed_max[j]);
			if (crossed_min[j] > 1)
			{	
				printf("***WARNING*** There is more than one minimum!\n");
				printf("j is %d and minima at %d and %d for first and second\n",  j, imin[j], i-1);
				for (k=imin[j]-1;k<i+1;k++) {
				printf("chi[%d][%d] = %f\n", j, k, chi[j][k]);}
				//if (i-1 - imin[j] > 2) 	
				//{	exit(1); } 
				crossed_min[j] = 1; 
			}
			imin[j] = i-1; // Store the index of the maximum 
			upperlimit[j] = imin[j] - imax[j];
			chimpp[j] = ( ( chi[j][i] - chi[j][i-1] ) / (xx[i] - xx[i-1]) - (chi[j][i-1] - chi[j][i-2]) / (xx[i-1] - xx[i-2]) ) *2.0/ (xx[i] - xx[i-2] ) ;  
			//printf("chimpp[%d] = %f\n", chimpp[j]);
			chimin[j] = chi[j][i-1];  //Store chi minimum itself, a function of xbar (index j)
			crossed_max[j] = 1;
/* We will now start going up the effective potential well! The temporary flag below is used because we still have to store the effective potential minimum as a possible value of Uperp. If we don't have this flag we miss the region near the minimum of chi. */
		}
/* FILLING IN ARRAYS */
/* Once we cross a maximum, we start filling in arrays for mu, Uperp and vx related to the orbit we are in. As we go down the maximum we store the values of Uperp = chi(x) we encounter which will form a grid of (unevenly spaced) allowed values of Uperp. We also store the value of the small deltamu associated with every value of Uperp above the current value of chi(x), and add it to the previous values */                                                                               	
		if ( (crossed_max[j] == 1  && crossed_min[j] == 0) || (i-1 == imin[j]) )
		{
			Uperp[j][i-1-imax[j]+kdrop[j]] = chi[j][i-1];
			if (Uperp[j][i-1-imax[j]+kdrop[j]] < 24.0*Telarge && Uperp[j][i-2-imax[j]+kdrop[j]] > 24.0*Telarge)
			lowerlimit[j] = i-1-imax[j]+kdrop[j]; 
			mu[j][i-1-imax[j]+kdrop[j]] = 0.0;
			//printf("kdrop[j] = %d\n", kdrop[j]);
			upper[j][i-1] = i-1-imax[j]+kdrop[j];
/* Note that the size of the dimension of the array with values of Uperp is set to n, which is larger than the size it will turn out to be. This is because in C there is no way to append elements to arrays as I go along, enough memory has to be given to the array from the start. n is the largest possible size the array could have. */
			for (k=0;k<=upper[j][i-1]; k++) 
			{	
				// replaced k with upper below
				if (upper[j][i-1] == 0) {
				//if (k==0) 
					vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));// equiv 0
					//printf("j = %d, i-1 = %d, %f\n", j, i-1, chi[j][i-2]-chi[j][i-1]);
					mu[j][k] += 0.0; 

				}
				else if (k == imin[j] - imax[j]) {
					vx[j][i-1][k] = 0.0;
					mu[j][k] = 0.0; 
				}
				else if (k == imin[j] - imax[j] - 1) {
					vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));// equiv 0
					mu[j][k] += (1.0/M_PI)*sqrt(0.5*chimpp[j])*pow(xx[imin[j]-1] - xx[imin[j]], 2.0);
					// should uncomment above and comment below (spurious)
					//mu[j][k] += (2.0/M_PI)*sqrt(chi[j][i-2]-chi[j][i-1])*(2.0/3.0)*(xx[i-1] - xx[i-2]);
					
				}
				else if ( k == upper[j][i-1] - 1 ) {	
					vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));// equiv 0
					mu[j][k] += (2.0/M_PI)*sqrt(chi[j][i-2]-chi[j][i-1])*(2.0/3.0)*(xx[i-1] - xx[i-2]);
				}
				else {
					vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));
					mu[j][k] += (1.0/M_PI)*(vx[j][i-1][k] + vx[j][i-2][k])*(xx[i-1] - xx[i-2]); 
				}
			}
		}
/* Once we cross the minimum, we stop creating array elements with values of Uperp. However, we keep storing the value of vx associated with any given point x on an effective potential curve with xbar, with energy Uperp and using this value to finish performing the mu integral. This should happen as long the effective potential at the point under consideration is smaller than the effective potential maximum. */
		else if ( crossed_min[j] == 1 && crossed_max[j] == 1 && chi[j][i-1] < chiMax[j] && ( i-1 != imin[j] ) ) {
			for (k=0;k <= upperlimit[j] ;k++)
			{	
				//if (i-1 = imin[j] +1) {
				//	vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));
				//	mu[j][k] += (1.0/M_PI)*(vx[j][i-1][k] + vx[j][i-2][k])*(xx[i-1] - xx[i-2]); 
				//}
				//if (Uperp[j][k] <= chi[j][i-1] && (i-2 == imin[j]) && k != upperlimit[j] ) {
				//	mu[j][k] += (1.0/M_PI)*sqrt(0.5*chimpp[j])*pow(xx[imin[j]-1] - xx[imin[j]], 2.0);
				//}
				//else if (chi[j][i-1] < Uperp[j][k]) {	
				if (chi[j][i-1] < Uperp[j][k]) {	
					vx[j][i-1][k] = sqrt((Uperp[j][k] - chi[j][i-1]));
					mu[j][k] += (1.0/M_PI)*(vx[j][i-1][k] + vx[j][i-2][k])*(xx[i-1] - xx[i-2]); 
				}
				else if (Uperp[j][k] <= chi[j][i-1] && Uperp[j][k-1] > chi[j][i-1]) {
					upper[j][i-1] = k;
			//mu[j][k] += (2.0/M_PI)*(vx[j][i-2][k])*(xx[i-1] - xx[i-2])*(Uperp[j][k] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]);
					ind = 0;
					while (Uperp[j][k] < chi[j][i-2-ind]) 
						ind++;
					
					mu[j][k] += (2.0/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[j][k] - chi[j][i-2-ind], 1.5)/(chi[j][i-1-ind] - chi[j][i-2-ind]); 
					//}
					//printf("YO %f and mu = %f: i = %d, j = %d, k = %d \nchi = %f, %f\n", mu[j][k], (2.0/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[j][k] - chi[j][i-2-ind], 1.5)/(chi[j][i-1-ind] - chi[j][i-2-ind]), i, j, k, chi[j][i-1], chi[j][i-2]);
		
				}
				else if (Uperp[j][k] <= chi[j][i-1] && Uperp[j][k] > chi[j][i-2]) {
					//if (i-1 == imin[j] +1) {
					//	//mu[j][k] += (2.0/M_PI)*(2.0/3.0)*(xx[i-1] - xx[i-2])*pow(Uperp[j][k] - chi[j][i-2], 1.5)/(chi[j][i-1] - chi[j][i-2]); 
					//	//mu[j][k] += (2.0/M_PI)*sqrt(0.5*chimpp[j])*pow(xx[imin[j]-1] - xx[imin[j]], 2.0);
					//	mu[j][k] += (1.0/M_PI)*sqrt(0.5*chimpp[j])*pow(xx[imin[j]-1] - xx[imin[j]], 2.0);
					//}
					//else {
				      	mu[j][k] += (2.0/M_PI)*(2.0/3.0)*(xx[i-1] - xx[i-2])*pow(Uperp[j][k] - chi[j][i-2], 1.5)/(chi[j][i-1] - chi[j][i-2]); 
					//}
				}
				if (mu[j][k] != mu[j][k]) {
					printf("mu is NAN\n"); 
					exit(-1);
				}  
			}
		}
/* When the effective potential at the iteration point (i-1) under consideration becomes larger than the effective potential maximum, we finish performing the mu integral. We also store the position of the top of the orbit which has chi = chiMax, in order to perform the open orbit integral. If the loop below is accessed, a switch it turned off to signify the no more closed orbits can be present */
		//else if ( ( crossed_min[j] == 1 ) && ( crossed_max[j] == 1) && ( chi[j][i-1] > chiMax[j] - TINY) ) 	
		else if ( ( crossed_min[j] == 1 ) && ( crossed_max[j] == 1) && ( chi[j][i-1] > chiMax[j] ) ) {	
			itop[j] = i-2;
			xtop[j] = xx[i-2] + ((chiMax[j] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]))*(xx[i-1] - xx[i-2]);
			chiprimetop[j] = (chi[j][i-1] - chi[j][i-2])/(xx[i-1] - xx[i-2]);
			for (k=0; k<upper[j][i-2]; k++) {
				ind = 0;
				while (Uperp[j][0] < chi[j][i-2-ind]) {
					ind++;
				}
				mu[j][k] += (2.0/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[j][0] - chi[j][i-2-ind], 1.5)/(chi[j][i-1-ind] - chi[j][i-2-ind]); 
				//printf("mu[%d][%d] = %f\n", j, k, mu[j][k]);
				//mu[j][k] += (1.0/M_PI)*(vx[j][i-2][k])*(xx[i-1] - xx[i-2])*(Uperp[j][k] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]);
				//if (mu[j][k] != mu[j][k])  printf("1mu is NAN\n");
			}
			crossed_max[j] = 0; 
		}
		if (j!=0)
		{	
			if ( ((chiMax[j-1] < chi[j-1][i-1] ) && (chiMax[j] > chi[j][i-1])) ) { // || (imax[j-1] == -1 && imax[j] != -1) )
				jmclosed[i-1] = j-1;
				if (i - 1 > icrit) { 	
					jmopen[i-1] = j-1; 
					if (DEBUG == 1) {
						printf("i = %d, j = %d, icrit = %d\n", i-1, j-1, icrit); 
						printf("jmopen[%d] = %d\n", i-1, jmopen[i-1]); 
					}
				} 
			}
			//if (i==size_finegrid-1 && crossed_max[j] == 1 && crossed_min[j] == 1) {	
			//	itop[j] = size_finegrid-1;
			//	if (noback == 0) {
			//		maxj = j + 50;	
			//		printf("maxj = %d\n", maxj);
			//		noback = 1; 
			//	} 
			//} 
		} 
		//printf("upper = %d, upperlimit = %d\n", upper[j][i-1], upperlimit[j]);
	} 
}	


//printf("xbar[maxj=%d] = %f\n", maxj, xbar[maxj]);	
// OPEN ORBIT INTEGRAL
/* Now we perform the open orbit integral. We use a change of variables which makes the integrand smooth at the top bounce point. The change of variables is to some var = sqrt(x_t - x) */
maxj = sizexbar - 1;
muopen[0] = 0.0;
for (j=0;j<maxj;j++) {
	//printf("j=%d/%d\n", j, maxj);
	mu[j][upperlimit[j]] = 0.0;
	mu[j][upperlimit[j]-1] = pow(xx[imin[j]] - xx[imin[j]-1], 2.0)*pow(0.5*chimpp[j], 0.5);
	k=0;
	if (charge < 0.0) {
	while (fabs(phi[0]) < 0.3 && k < upperlimit[j]) {
		mu[j][upperlimit[j]-k] = pow(xx[imin[j]] - xx[imin[j]-k], 2.0)*pow(0.5*chimpp[j], 0.5);
		k++;
	}
	}
	//if (j == 20) {
	//	for (k=0; k<upperlimit[j]+1; k++)
	//		printf("mu[%d][%d] = %f\t Uperp[%d][%d] = %f\n", j, k, mu[j][k], j, k, Uperp[j][k]); 
	//		//printf("mu[%d][%d] = %f\t Uperp[%d][%d] = %f\n", j, k, mu[j][k], j, k, Uperp[j][k] - phi[j]); 
	//}
	//printf("mu[%d][%d] = %f (should be %f near infinity)\n", j, upperlimit[j]-1, mu[j][upperlimit[j]-1], Uperp[j][upperlimit[j]-1]);
	//printf("imin = %d\n", imin[j]);
	muopen[j+1] = mu[j][0];
	//Ucritf[j+1] = sqrt(chiMax[j] - mu[j][0]) ;
	Ucritf[j+1] = chiMax[j] - mu[j][0] ;
	//printf("mucrit , Ucrit = %f, %f\n", muopen[j+1], Ucritf[j+1]);
	if (DEBUG == 1)
	{	
		printf("jmopen[%d] = %d\n", j, jmopen[j]); 
	}
	//METHOD 1
	if (j != 0)
	{	
		//twopimuprime[j-1] = 2.0*M_PI*(mu[j][0] - mu[j-1][0])/(xbar[j] - xbar[j-1]);	
		twopimuprime[j] = ((xbar[j] - xbar[j-1])/(xbar[j+1] - xbar[j-1])) *(mu[j+1][0] - mu[j][0])/(xbar[j+1] - xbar[j]) + ((xbar[j+1] - xbar[j])/(xbar[j+1] - xbar[j-1])) * (mu[j][0] - mu[j-1][0])/(xbar[j] - xbar[j-1]);	
		twopimuprime[j] *= (2.0*M_PI);
	
	}
	else
	{
		//twopimuprime[j] = 2.0*M_PI*(mu[j][0] - 0.0)/(xbar[j+1] - xbar[j]);	
		twopimuprime[j] = ((xbar[j] - xbarcrit)/(xbar[j+1] - xbarcrit)) *(mu[j+1][0] - mu[j][0])/(xbar[j+1] - xbar[j]) + ((xbar[j+1] - xbar[j])/(xbar[j+1] - xbarcrit)) * (mu[j][0] - 0.0)/(xbar[j] - xbarcrit);	
		twopimuprime[j] *= (2.0*M_PI);
		//twopimuprime[j] = M_PI*(mu[j+1][0] - mu[j][0])/(xbar[j+1] - xbar[j]);	
	}
// METHOD 2
	//if (imax[j] != 0) {
	//	openintegral[j] += 4.0*(xx[imax[j]+1] - xx[imax[j]])*(xx[imax[j]+1] - xx[imax[j]])/vx[j][i][0];
	//}
	//else if (imax[j] == 0) {
	//	openintegral[j] += 2.0*(xx[imax[j]+1] - xx[imax[j]])*(xx[imax[j]+1] - xx[imax[j]])/vx[j][i][0];
	//}
	for (i=imax[j]+2;i<itop[j];i++)
	{
		openintegral[j] += 2.0*( (xx[i] - xx[imax[j]])/vx[j][i][0] + (xx[i-1] - xx[imax[j]])/vx[j][i-1][0] )*(xx[i] - xx[i-1]);
	}
	//openorbit[j] = openintegral[j];
	openorbit[j] = twopimuprime[j];  
	openorbitantycal = 4.0*M_PI*xbar[j];
	if (DEBUG == 1)
		printf("%f %f %f %f %f\n", xbar[j], mu[j][0], Uperp[j][0], openorbit[j], openorbitantycal); 
}
// temporary
muopen[maxj] = 1000.0;
//
//Ucritf[0] = Ucritf[1] - muopen[1]*(Ucritf[2] - Ucritf[1])/(muopen[2] - muopen[1]);
//Ucritf[0] = Ucritf[1];
//Ucritf[0] = sqrt(phi[0] + xbarcrit*xbarcrit);
Ucritf[0] = phi[0] + xbarcrit*xbarcrit;
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
if ((fout = fopen("OUTPUT/densfinorb_out.txt", "w")) == NULL)
{	
	printf("Cannot open densfinorb_out.txt");
	exit(EXIT_FAILURE);
}

////////////////////////////////////////////
// calculate normalization at infinity
deltax_inf = (xx[size_finegrid-1] - xx[size_finegrid-2]);
j_inf = (int) sqrt(20.0*Telarge)/deltax_inf ;
printf("j_inf = %d, deltax_inf = %f\n", j_inf, deltax_inf);
muinf = (double*)calloc(j_inf,sizeof(double*)); 
vxinf = (double**)calloc(j_inf,sizeof(double**)); 
chiinf = (double*)calloc(j_inf,sizeof(double*));
for (j=j_inf-1; j>=0; j--) {	
	vxinf[j] = (double*)calloc(j_inf,sizeof(double*)); 
	chiinf[j] =  deltax_inf*j*deltax_inf*j;
	for (k=j; k<j_inf; k++) {
		vxinf[j][k] = sqrt(chiinf[k] - chiinf[j]);
		//if (j!=j_inf-1) {
		//	if (k>j+1)
		//		muinf[k] += (2.0/M_PI)*(vxinf[j+1][k] + vxinf[j][k])*deltax_inf;
		//	else if (k==j+1)
		//		muinf[k] += (4.0/M_PI)*sqrt(chiinf[j+1]-chiinf[j])*(2.0/3.0)*deltax_inf;
		//}
	}
	muinf[j] = chiinf[j];
	//chiinf[j] += 0.001;
	//muinf[j] += TINY;
}
//muinf[0] = 0.0;
//muinf[1]


/*
flow_inf = 


*/

intdxbar = 0.0;
intdvx = 0.0;
intdU = 0.0;
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
		//Ucrit = pow(lin_interp(muopen, Ucritf, munew, sizexbar, 791), 2.0);
		Ucrit = lin_interp(muopen, Ucritf, munew, sizexbar, 791);
		//printf("vxnew = %f\tmunew = %f\tUcrit = %f\n", vxnew, munew, Ucrit);
		Uperpnew = chiinf[k];
		sizeU = (int) sqrt(Ucap - Uperpnew)/dvz;
		reflected = 1;
			for (l=0; l < sizeU; l++)
			{	
				if (l!=0)
				{	
					Fold = F;
					Fold_ref = F;
					vz = dvz*l;
					U = Uperpnew + pow(vz, 2.0);
					//if ( (U < munew) && (charge < 0.0) ) printf("???\n\n\n\n");
					if ( (U > munew) && (U - vz*vz + (vz-dvz)*(vz-dvz) < munew) ) {
						frac = (vz - sqrt(munew - Uperpnew))/dvz;
						//if (munew > Uperpnew) {
						//	printf("Uperpnew = %f, munew = %f\n", Uperpnew, munew);
						//}
						//else {
						//	frac = vz/dvz;
						//	//printf("frac = %f\n", frac);
						//	printf("???Uperpnew = %f, munew = %f\n", Uperpnew, munew);
						//}
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


					if ( (U-munew < Ucrit) && (reflected == 1) ) {
						frac_reflected = 1.0;
						//printf("frac_reflected = %f\n", frac_reflected);
						//if (U < munew) frac_reflected = 0.0;
					}
					else if ( (U-munew >= Ucrit) && (reflected == 1) ) {
						reflected = 0;
						if (Ucrit < Uperpnew - munew) frac_reflected = 0;
						else {
							vzcrit = sqrt(Ucrit + munew - Uperpnew);
							Fold_ref = bilin_interp(munew, Ucrit, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
							frac_reflected = (vzcrit - (vz - dvz))/dvz;
						}
						//if (munew < 0.7 && munew > 0.5 ) {
						//printf("frac_reflected = %f, l=%d, vzcrit = %f, vz - dvz = %f\n", frac_reflected, l, vzcrit, vz - dvz);
						//printf("U = %f, Ucrit = %f, munew = %f, Uperpnew = %f\n", U, Ucrit, munew, Uperpnew);
						//}
					}
					else frac_reflected = 0.0;
							
					//if (reflected == 1) printf("particles are being reflected\n");
					//F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					//printf("frac_reflected = %f\n", frac_reflected);
					//printf("F = %f\n", F);
					//intdU += (1.0 + frac_reflected)*frac*dvz*(F+Fold);
					//if ( ( (i==70) && (charge < 0.0) ) && ( l <= 3) )
					//	printf("Fold = %f\tF = %f\n", Fold, F);
					//frac_reflected = 0.0;
					intdU += frac*dvz*(F+Fold) + frac*frac_reflected*dvz*(F + Fold_ref);
					//printf("frac = %f, frac_reflected = %f, F = %f, Fold = %f, Fold_ref = %f\n", frac, frac_reflected, F, Fold, Fold_ref);
				}
				else {	
					//vz = dvz*l;
					U = Uperpnew ;//+ pow(vz, 2.0);
					//if (U-munew  < Ucrit) 
					//	reflected = 1;
					//else 
					//	printf("not here?\n");
					//reflected = 1;
					F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					intdU += 0.0;
				} 
			}




		//for (l=0; l < sizeU; l++) {	
		//	if (l!=0)
		//	{	
		//		Fold = F;
		//		Fold_ref = F;
		//		vz = dvz*l;
		//		U = Uperpnew + pow(vz, 2.0);
		//		if ( (U > munew) && (U - vz*vz + (vz-dvz)*(vz-dvz) < munew) ) {
		//			frac = (vz - sqrt(munew - Uperpnew + TINY))/dvz;
		//			//printf("frac = %f\n", frac);
		//			Fold = bilin_interp(munew, 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
		//			F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
		//		}
		//		else if (U > munew) {
		//			frac = 1.0;
		//			F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
		//		}
		//		else {
		//			frac = 1.0;
		//			F = 0.0;
		//		}
		//			
		//		//U = Uperpnew + pow(vz, 2.0);

		//		if ( (U-munew < Ucrit) && (reflected == 1) ) {
		//			frac_reflected = 1.0;
		//		}
		//		else if ( (U-munew >= Ucrit) && (reflected == 1) ) {
		//			reflected = 0;
		//			vzcrit = sqrt(Ucrit + munew - Uperpnew);
		//			Fold_ref = bilin_interp(munew, Ucrit, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
		//			frac_reflected = (vzcrit - (vz - dvz))/dvz;
		//		}
		//		else frac_reflected = 0.0;
		//				
		//		//frac_reflected = 0.0;
		//		intdU += frac*dvz*(F+Fold) + frac*frac_reflected*dvz*(F + Fold_ref);
		//	}
		//	else {	
		//		//vz = dvz*l;
		//		U = Uperpnew ;//+ pow(vz, 2.0);
		//		if (U-munew < Ucrit) {
		//			reflected = 1;
		//		}
		//		F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
		//		intdU += 0.0;
		//	} 
		//}
		intdUantycal = exp(-chiinf[k])*(1.0/(2.0*M_PI));// result with phi =0
		if (DEBUG == 1) printf("Analytical intdU is %f, numerical one is %f\n", intdUantycal, intdU);
		if (k!=j) {	
			dvx = vxnew - vxold; 
			intdvx += 4.0*0.5*dvx*(intdU+intdUold);
		}
		intdvxantycal = (2.0/(2.0*sqrt(M_PI)))*exp(-deltax_inf*j*deltax_inf*j);//*erf(sqrt(xbar[j]*xbar[j]-(pos-xbar[j])*(pos-xbar[j])));
	}
	if (DEBUG == 1)	printf("intdvx is %f, analytical one is %f\n", intdvx, intdvxantycal); 
	//printf("intdvx is %f, while the analytical one is %f\n", intdvx, intdvxantycal);
	//intdvx = intdvxantycal;
	dxbar = deltax_inf;
	if (j!=0) intdxbar += (intdvx+intdvxold)*dxbar; // multiply by two because only half the xbars are contemplated here
}
n_inf = intdxbar;
if (DEBUG == 0) {	
	printf("for charge = %f\n", charge);
	printf("n_inf = %f\n", n_inf);
	if ( (n_inf  != n_inf) || (n_inf < TINY) )
	{ 	
		printf("n_inf = %f\n", n_inf);
		exit(-1);
	} 
} 


// density profile
i=0;
stop = 0;
ic = 0;
while (stop == 0) {
	i = ic*zoomfactor;
	pos = xx[i];
	intdxbar = 0.0;
	intdxbaropen = 0.0;
	intdxbaropenflow = 0.0;
	for (j=0; j<sizexbar; j++) {
		vxnew = 0.0;
		intdUopenold = intdUopen;
		intdUopenflowold = intdUopenflow; intdUopenBohmold = intdUopenBohm;
		intdUopenxbarold = intdUopenxbar; intdUopensquareold = intdUopensquare;
		intdUopen = 0.0;
		intdUopenflow = 0.0; intdUopenBohm = 0.0;
		intdUopenxbar = 0.0; intdUopensquare = 0.0;
		if (j == jmopen[i]) {
			oorbintgrd = oorbintgrdold = 0.0;
			oorbintgrdflow = oorbintgrdflowold = 0.0; oorbintgrdBohm = oorbintgrdBohmold = 0.0;
			oorbintgrdxbar = oorbintgrdxbarold = 0.0; oorbintgrdsquare = oorbintgrdsquareold = 0.0;
			sizeU = (int) sqrt(Ucap - chiMax[j])/dvzopen;
			if (j==0) {
				munew = 0.0;
				openorbitnew = 0.0;
			}
			else {
				munew = mu[j+1][0] + ( (chiMax[j+1] - chi[j+1][i])/ (chiMax[j+1] - chi[j+1][i] + chi[j][i] - chiMax[j]) ) * (mu[j][0] - mu[j+1][0]);
				openorbitnew = openorbit[j+1] + ( (chiMax[j+1] - chi[j+1][i])/ (chiMax[j+1] - chi[j+1][i] + chi[j][i] - chiMax[j]) ) * (openorbit[j] - openorbit[j+1]);
			}
			for (l=0; l < sizeU; l++) {
				oorbintgrdold = oorbintgrd;
				oorbintgrdflowold = oorbintgrdflow; oorbintgrdBohmold = oorbintgrdBohm;
				oorbintgrdxbarold = oorbintgrdxbar; oorbintgrdsquareold = oorbintgrdsquare;
				vz = dvzopen*l;
				if (j!=0) {
					chinew = chi[j+1][i] + ( (chiMax[j+1] - chi[j+1][i])/ (chiMax[j+1] - chi[j+1][i] + chi[j][i] - chiMax[j]) ) * (chi[j][i] - chi[j+1][i]) ;
					vx0open = 0.0;
				}
				else {
					chinew = chi[j][i];
					if (chiMax[j] > chinew) 
						vx0open = sqrt(chiMax[j] - chinew);
					else vx0open = 0.0;
				}
				U = chinew + pow(vz, 2.0); 
				if ( (U > munew) && (U - vz*vz + (vz-dvzopen)*(vz-dvzopen) < munew) ) {
					frac = (vz - sqrt(munew - chinew))/dvzopen;
					Fopen = bilin_interp(munew, 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					//printf("frac = %f\n", frac);
					oorbintgrdold = ( sqrt(vx0open*vx0open + alpha*sqrt(munew - chinew)*openorbitnew) - vx0open )*Fopen;
					oorbintgrdflowold = 0.5*alpha*sqrt(munew - chinew)*openorbitnew*Fopen;
					Fopen = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					oorbintgrd = ( sqrt(vx0open*vx0open + alpha*vz*openorbitnew) - vx0open )*Fopen;
					oorbintgrdflow = 0.5*alpha*vz*openorbitnew*Fopen;
				}
				else {
					frac = 1.0;
					Fopen = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					oorbintgrd = ( sqrt(vx0open*vx0open + alpha*vz*openorbitnew) - vx0open )*Fopen;
					oorbintgrdflow = 0.5*alpha*vz*openorbitnew*Fopen;
				}
				if (oorbintgrd != oorbintgrd) {	
					if (DEBUG == 1)
					{
						printf("vx0open = %f, openorbit[%d] = %f, imaginary oorbintgrd in initial piece of integral due to negative value of openorbit?\n", vx0open, j, openorbit[j]); 
					}
					exit(-1);
				}
				if (i==0) 
					oorbintgrdantycal = sqrt(2.0*alpha*vz*M_PI*xbar[j])*(U-chiMax[j])*exp(-U)/pow(M_PI, 1.5); 
// for a flat potential, oorbintgrd can be calculated analytically
					//oorbintgrd = oorbintgrdantycal;
					//printf("oorbintgrd is %f, analytical is %F\n", oorbintgrd, oorbintgrdantycal);
				if (l!=0) {	
					intdUopen += 2.0*frac*dvzopen*(oorbintgrd+oorbintgrdold);
					intdUopenflow += 2.0*frac*dvzopen*(oorbintgrdflow+oorbintgrdflowold);
					intdUopenBohm += 2.0*frac*dvzopen*(oorbintgrdBohm+oorbintgrdBohmold);
					intdUopenxbar += 2.0*frac*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
					intdUopensquare += 2.0*frac*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
				}
				else {	
					intdUopen += 0.0;
					intdUopenflow += 0.0;
					intdUopenBohm += 0.0;
					intdUopenxbar += 2.0*frac*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
					intdUopensquare += 2.0*frac*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
				} 
			}
			if (j == 0) {
				intdUopenold = 0.0;
				//printf("intdUopenold = %f\n", intdUopenold);
				dxbar = xbar[0] - xbarcrit;
				dxbar = 0.0;
				intdxbaropen += 0.5*(intdUopen+intdUopenold)*dxbar;
				intdxbaropenflow += 0.5*(intdUopenflow+intdUopenflowold)*dxbar;
				intdxbaropenBohm += 0.5*(intdUopenBohm+intdUopenBohmold)*dxbar;
				intdxbarxbar += 0.5*(intdUopenxbar+intdUopenxbarold)*dxbar;
				intdxbarsquare += 0.5*(intdUopensquare+intdUopensquareold)*dxbar; 
			}
			//else 
			//	//dxbar = (xbar[j] - xbar[j-1])*(chiMax[j] - chi[j][i])/ (chiMax[j] - chi[j][i] + chi[j-1][i] - chiMax[j-1]);// open orbit density does not need to be so accurate at this point
			//	dxbar = 0.0;
			//	//printf("dxbar = %f\n", dxbar);
			// x is already large enough for closed orbit density to dominate if jmopen is not zero
			//printf("dxbar = %f\n\n\n\n\n", dxbar);
		}
		else if (j > jmopen[i]) {
			oorbintgrd = oorbintgrdold = 0.0;
			oorbintgrdflow = oorbintgrdflowold = 0.0;
			oorbintgrdBohm = oorbintgrdBohmold = 0.0;
			oorbintgrdxbar = oorbintgrdxbarold = 0.0;
			oorbintgrdsquare = oorbintgrdsquareold = 0.0;
			sizeU = (int) sqrt(Ucap - chiMax[j])/dvzopen;
			for (l=0; l < sizeU; l++) {
				oorbintgrdold = oorbintgrd;
				oorbintgrdflowold = oorbintgrdflow;
				oorbintgrdBohmold = oorbintgrdBohm;
				oorbintgrdxbarold = oorbintgrdxbar;
				oorbintgrdsquareold = oorbintgrdsquare;
				vz = dvzopen*l;
				U = chiMax[j] + pow(vz, 2.0);
				//vx0open = sqrt(TINY + chiMax[j] - chi[j][i]);
				vx0open = sqrt(chiMax[j] - chi[j][i]);
				if (vx0open != vx0open) {		
					if (DEBUG == 1)
						printf("HERE imaginary vx0open, j = %d, i is %d, chi[j][i] = %f, chiMax[j] = %f\n", j, i, chi[j][i], chiMax[j]); 
					exit(-1);
				}
				if ( (U >= mu[j][0]) && (U - vz*vz + (vz-dvzopen)*(vz-dvzopen) < mu[j][0]) ) {
					frac = (vz - sqrt(mu[j][0] - chiMax[j]))/dvzopen;
					Fopen = bilin_interp(mu[j][0], 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					oorbintgrdold = sqrt(alpha*sqrt(mu[j][0] - chiMax[j])*openorbit[j])*Fopen;
					oorbintgrdflowold = 0.5*alpha*sqrt(mu[j][0] - chiMax[j])*openorbit[j]*Fopen;
					oorbintgrdBohmold = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*sqrt(mu[j][0] - chiMax[j])*openorbit[j]))*Fopen;
					oorbintgrdxbarold = xbar[j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*sqrt(mu[j][0] - chiMax[j])*openorbit[j]))*Fopen;
					oorbintgrdsquareold = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + alpha*sqrt(mu[j][0] - chiMax[j])*openorbit[j]), 1.5))*Fopen;
					Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					oorbintgrd = sqrt(alpha*vz*openorbit[j])*Fopen;
					oorbintgrdflow = 0.5*alpha*vz*openorbit[j]*Fopen;
					oorbintgrdBohm = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*vz*openorbit[j]))*Fopen;
					oorbintgrdxbar = xbar[j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*vz*openorbit[j]))*Fopen;
					oorbintgrdsquare = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + alpha*vz*openorbit[j]), 1.5))*Fopen;

				}
				else {
					frac = 1.0;
					Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
					//if (i ==0) printf("l = %d, Fopen = %f\n", l, Fopen);
					oorbintgrd = (sqrt(vx0open*vx0open + alpha*vz*openorbit[j]) - vx0open)*Fopen;
					oorbintgrdflow = 0.5*alpha*vz*openorbit[j]*Fopen;
					oorbintgrdBohm = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*vz*openorbit[j]))*Fopen;
					oorbintgrdxbar = xbar[j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + alpha*vz*openorbit[j]))*Fopen;
					oorbintgrdsquare = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + alpha*vz*openorbit[j]), 1.5))*Fopen;
					//if (oorbintgrd != oorbintgrd) 
					if (frac != frac) 
						printf("AHA\n\n\n\n\n");
				}
				//Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1);
				if (oorbintgrd != oorbintgrd) {
					if (DEBUG == 1) {
						printf("vx0open = %f, openorbit[%d] = %f, HERE imaginary oorbintgrd in initial piece of integral due to negative value of openorbit?\n", vx0open, j, openorbit[j]); 
					} 
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
					intdUopen += 2.0*frac*dvzopen*(oorbintgrd+oorbintgrdold);
					intdUopenflow += 2.0*frac*dvzopen*(oorbintgrdflow+oorbintgrdflowold);
					intdUopenBohm += 2.0*frac*dvzopen*(oorbintgrdBohm+oorbintgrdBohmold);
					intdUopenxbar += 2.0*frac*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
					intdUopensquare += 2.0*frac*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
				}
				else {
					intdUopen += 0.0;
					intdUopenflow += 0.0;
					intdUopenBohm += 0.0;
					intdUopenxbar += 2.0*frac*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
					intdUopensquare += 2.0*frac*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
				} 
			}
			if (intdUopen != intdUopen)  {
				if (DEBUG == 1)
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
			if ( (j == jmopen[i]+1) &&  (j != 1) ) { //
				//if ( fabs(pos - xx[imax[j]]) < TINY )
				//	intdxbaropen += 0.0; 
				//else { 
				dxbar = (xbar[j] - xbar[j-1])*(chiMax[j] - chi[j][i])/ (chiMax[j] - chi[j][i] + chi[j-1][i] - chiMax[j-1]);// open orbit density does not need to be so accurate at this point
				//printf("intdUopen = %f\t old = %f\n", intdUopen, intdUopenold);
				//printf("intdUopenold = %f\n", intdUopenold);
				//intdUopenold = intdUopen + ( (chiMax[j] - chi[j][i])/ (chiMax[j] - chi[j][i] + chi[j-1][i] - chiMax[j-1]) ) * (intdUopenold - intdUopen);// open orbit density does not need to be so accurate at this point
				
				if (DEBUG == 1)
					printf("dxbar = %f\n", dxbar);
				//}
			}
//(chi[j][k]-chi[j][i])/(2.0*(pos-xx[k])); 
			intdxbaropen += 0.5*(intdUopen+intdUopenold)*dxbar;
			intdxbaropenflow += 0.5*(intdUopenflow+intdUopenflowold)*dxbar;
			intdxbaropenBohm += 0.5*(intdUopenBohm+intdUopenBohmold)*dxbar;
			intdxbarxbar += 0.5*(intdUopenxbar+intdUopenxbarold)*dxbar;
			intdxbarsquare += 0.5*(intdUopensquare+intdUopensquareold)*dxbar; 
			if (intdxbaropen != intdxbaropen) {
				if (DEBUG == 1)
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
						//added spurious extra below
						//if (k == upperlimit[j]) munew = 0.0;
						//if (k == upperlimit[j]-1) munew = Uperp[j][k];
						//if (k == upperlimit[j]-2) munew = Uperp[j][k];
						vxnew = 0.0;
						dvx = vxold - vxnew; 
					}
					else
					{	
						Uperpnew = Uperp[j][k];
						vxnew = vx[j][i][k];
						dvx = vxold - vxnew;
						munew = mu[j][k]; 
						//added spurious extra below
						//if (k == upperlimit[j]) munew = 0.0;
						//if (k == upperlimit[j]-1) munew = Uperp[j][k];
						//if (k == upperlimit[j]-2) munew = Uperp[j][k];
					}
					//Ucrit = pow(lin_interp(muopen, Ucritf, munew, sizexbar, 791), 2.0);
					Ucrit = lin_interp(muopen, Ucritf, munew, sizexbar, 791);
					//if (i == 20) printf("muopen = %f\tUcrit = %f\n", munew, Ucrit);
					//printf("munew, Ucrit = %f, %f\n", munew, Ucrit);
					//printf("j, k, upperlimit = %d, %d, %d\n", j, k, upperlimit[j]);
					sizeU = (int) sqrt(Ucap - Uperpnew)/dvz;
					reflected = 1;
					for (l=0; l < sizeU; l++)
					{	
						if (l!=0)
						{	
							Fold = F;
							Fold_ref = F;
							vz = dvz*l;
							U = Uperpnew + pow(vz, 2.0);
							//if ( (U < munew) && (charge < 0.0) ) printf("???\n\n\n\n");
							if ( (U > munew) && (U - vz*vz + (vz-dvz)*(vz-dvz) < munew) ) {
								frac = (vz - sqrt(munew - Uperpnew+TINY))/dvz;
								//if (munew > Uperpnew) {
								//	printf("Uperpnew = %f, munew = %f\n", Uperpnew, munew);
								//}
								//else {
								//	frac = vz/dvz;
								//	//printf("frac = %f\n", frac);
								//	printf("???Uperpnew = %f, munew = %f\n", Uperpnew, munew);
								//}
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


							if ( (U-munew < Ucrit) && (reflected == 1) ) {
								frac_reflected = 1.0;
								//printf("frac_reflected = %f, l=%d, vzcrit = %f\n", frac_reflected, l, vzcrit);
								//if (U < munew) frac_reflected = 0.0;
							}
							else if ( (U-munew >= Ucrit) && (reflected == 1) ) {
								reflected = 0;
								if (Ucrit <= Uperpnew - munew) frac_reflected = 0.0;
								else {
									vzcrit = sqrt(Ucrit + munew - Uperpnew);
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
								//if ( i == 20 && munew < 0.65 && munew > 0.6 ) {
								//printf("??frac_reflected = %f, l=%d, vzcrit = %f, vz - dvz = %f\n", frac_reflected, l, vzcrit, vz - dvz);
								//printf("U = %f, Ucrit = %f, munew = %f\n", U, Ucrit, munew);
								//}
							}
									
							//if (reflected == 1) printf("particles are being reflected\n");
							//F = bilin_interp(munew, U-munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1); 
							//printf("frac_reflected = %f\n", frac_reflected);
							//printf("F = %f\n", F);
							//intdU += (1.0 + frac_reflected)*frac*dvz*(F+Fold);
							//if ( ( (i==70) && (charge < 0.0) ) && ( l <= 3) )
							//	printf("Fold = %f\tF = %f\n", Fold, F);
							//frac_reflected = 0.0;
							intdU += frac*dvz*(F+Fold) + frac*frac_reflected*dvz*(F + Fold_ref);
							//printf("frac = %f, frac_reflected = %f, F = %f, Fold = %f, Fold_ref = %f\n", frac, frac_reflected, F, Fold, Fold_ref);
						}
						else {	
							//vz = dvz*l;

							U = Uperpnew ;//+ pow(vz, 2.0);

							if (phi[0] < 0.0) reflected = 0.0;
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
						intdvx += 4.0*0.5*dvx*(intdU+intdUold);
					}
					if (intdvx != intdvx) {	
						if (DEBUG ==0)
							printf("intdvx is NAN, j=%d, i=%d\n", j, i); 
						exit(-1);
					} 
				}
				intdvxantycal = (2.0/(2.0*sqrt(M_PI)))*exp(-(pos-xbar[j])*(pos-xbar[j]))*erf(sqrt(xbar[j]*xbar[j]-(pos-xbar[j])*(pos-xbar[j])));
				if (DEBUG == 1) {
					printf("pos=%f, i=%d, j=%d, Uperp is %f, chi is %f\nvxnew and vxold are %f and %f and and vx is %f, dvx is %f\nintdvx is %f, analytical one is %f, upper is %d, upperlimit is %d\n", pos, i, j, Uperpnew, chi[j][i], vxnew, vxold, vx[j][i][k], dvx, intdvx, intdvxantycal, upper[j][i], upperlimit[j]); 
				}
				//printf("intdvx is %f, while the analytical one is %f\n", intdvx, intdvxantycal);
				//intdvx = intdvxantycal;
				if (j==jmclosed[i]+1 && i!=0) {
					//dxbar = (xx[imax[j]] - xx[imax[j]])*(xbar[j] - xbar[j-1])/(xx[imax[j]] - pos);
					//dxbar = (chiMax[j]-chi[j][i])/(2.0*(pos-xx[imax[j]]));
					if (i==imax[j]) {	
						dxbar = xbar[j] - xbar[j-1]; 
					}
					else {	
						//dxbar = xbar[j] - (xbar[j-1]*(chiMax[j] - chi[j][i]) - xbar[j]*(chiMax[j-1] - chi[j-1][i]))/(-chiMax[j-1] + chi[j-1][i] + chiMax[j] - chi[j][i]); 
						dxbar = (xbar[j] - xbar[j-1])*(chiMax[j] - chi[j][i])/ (chiMax[j] - chi[j][i] + chi[j-1][i] - chiMax[j-1]);// open orbit density does not need to be so accurate at this point
					}
					intdxbar += 0.5*(intdvx+intdvxold)*dxbar;
					if (intdxbar != intdxbar) {	
						if (DEBUG == 1) {	
							printf("intdxbar is NAN, j=%d, i=%d\n", j, i);  
						} 
						exit(-1);
					} 
				}
				else
				{
					dxbar = xbar[j] - xbar[j-1];
					intdxbar += 0.5*(intdvx+intdvxold)*dxbar;
					if (intdxbar != intdxbar) {
						if (DEBUG == 1)
							printf("intdxbar is NAN, j=%d, i=%d\n", j, i); 
						exit(-1);
					} 
				}
			} 
		} 
	}
	//ne[ic] = exp(phi[i]/Te); 
	niclosed[ic] = intdxbar;
	niopen[ic] = intdxbaropen;
	nitotal[ic] = niopen[ic] + niclosed[ic];
	//printf("nitotal[%d] = %f = %f + %f\n", ic, niclosed[ic], niopen[ic], nitotal[ic]);
	if (ic == 0) {
		Bohm = intdxbaropenBohm/nitotal[ic];
		if (charge < 0.0) 
			*flux = intdxbaropenflow/alpha; 
		if (DEBUG == 0) 
		//printf("flow velocity at x=0 = %f\n", (*flux)/nitotal[ic]);
		printf("flux evaluated at x=0 is %f\n", intdxbaropenflow/alpha);
		printf("flux evaluated at x=0 is %f\n", *flux);
		//printf("intdxbarsquare = %f\tintdxbarxbar = %f\n", intdxbarsquare, intdxbarxbar);
		//kBohmsquared = intdxbarxbar/(intdxbarsquare - 0.5*nitotal[0]); 
	}
	if ( nitotal[ic] >= stopdens*n_inf) {	
		printf("ic = %d/%d\n", ic, size_phigrid);
		stop = 1;
		*size_ngrid = ic; 
		printf("stopping because density larger than %f times the density at infinity\n", stopdens);
	}
	//else if ( phi_grid[ic] >= - (1.0-stopdens) ) {	
	//	printf("ic = %d/%d\n", ic, size_phigrid);
	//	stop = 1;
	//	*size_ngrid = ic; 
	//	printf("stopping because potential larger than -%fTe\n", 1.0-stopdens);
	//}
	else if (x_grid[ic] > x_grid[size_phigrid-1] - limit_rho) {
		printf("ic = %d/%d\n", ic, size_phigrid);
		stop = 1;
		*size_ngrid = ic; 
		printf("WARNING in densfinorb.c: stopping for positive density\n");

	}
	if (DEBUG == 0) {	
		printf("for charge = %f\n", charge);
		printf("%f, %f, %f is TOTAL, CLOSED and OPEN orbit density at position index %d, position %f, potential %f\n", nitotal[ic], niclosed[ic], niopen[ic], ic, pos, phi_grid[ic]);
		if ( (nitotal[ic] != nitotal[ic]) || (nitotal[ic] < TINY) ) { 	
			printf("nitotal[ic] = %f\n", nitotal[ic]);
			//exit(-1);
		} 
	} 
	fprintf(fout, "%f %f %f %f\n", pos, nitotal[ic]/n_inf, niclosed[ic]/n_inf, niopen[ic]/n_inf);
	//fprintf(fout, "%f %f %f %f\n", pos, nitotal[ic], niclosed[ic], niopen[ic]);
	//i += 1; 
	ic += 1;
} 
if (stop == 0) { 
	printf("ERROR: the density never reached stopdens*n_inf\n"); 
	exit(-1); 
}
for (ic = 0; ic < size_phigrid; ic++) {
	if (DEBUG == 1) 
		printf("before renormalizing nitotal[%d] = %f, phi_grid[%d] = %f\n", ic, nitotal[ic], ic, phi_grid[ic]);
	if (ic < *size_ngrid)
		nitotal[ic] /= n_inf;
	else nitotal[ic] = 0.0;
	if (DEBUG == 1) 
		printf("nitotal[%d] = %f, phi_grid[%d] = %f\n", ic, nitotal[ic], ic, phi_grid[ic]);
}
//if (charge < 0.0) {
//	for (ic = *size_ngrid; ic < size_phigrid; ic++) {
//		//nitotal[ic] = 0.5 * exp(phi_grid[ic]) * (1 + erf(sqrt(phi_grid[ic] -phi_grid[0] )));
//		nitotal[ic] = nitotal[*size_ngrid-1] ;//nitotal[*size_ngrid-1] + (-nitotal[*size_ngrid-1] + n_inf)* (phi_grid[ic] - phi_grid[*size_ngrid-1]) *exp(phi_grid[ic]);
//	}
//	*size_ngrid = size_phigrid;
//}



//ic=0;
//while (nitotal[ic] < stopdens*nitotal[*size_ngrid-1]) {
//	//nitotal[ic] /= nitotal[*size_ngrid-1];
//	ic += 1;
//}
//*size_ngrid = ic;
//if (charge < 0.0) {
//	//for (ic = 0; ic < *size_ngrid; ic++) {
//	//	nitotal[ic] /= n_inf;
//	//}
//	ic=0;
//	while (nitotal[ic] < stopdens*nitotal[*size_ngrid-1]) {
//		//nitotal[ic] /= nitotal[*size_ngrid-1];
//		ic += 1;
//	}
//	*size_ngrid = ic;
//}
fclose(fout);

if (output_type == 1)
{
	dvxopen = (5.0/dim_output)*sqrt(Telarge);
	Ucap = 18.0 + 2.0*Te;
	/* The part below performs the integration over the open orbit distribution function at x=0 in a way that allows to extract the distribution function in v_x (because v_x is the only velocity component that matters in the Debye sheath).*/
	FILE *output, *outputyz, *outputyzvy, *outputyzvz;
	if ((output = fopen("OUTPUT/f0x.txt", "w")) == NULL)
	{        // Check for presence of file
		printf("Cannot open %s\n", "outputfile.txt");
		exit(EXIT_FAILURE);
	}
	if ((outputyz = fopen("OUTPUT/f0yz.txt", "w")) == NULL)
	{
		// Check for presence of file
		printf("Cannot open %s\n", "outputyz.txt");
		exit(EXIT_FAILURE);
	}
	if ((outputyzvy = fopen("OUTPUT/vy0.txt", "w")) == NULL)
	{
		// Check for presence of file
		printf("Cannot open %s\n", "outputyzvy.txt");
		exit(EXIT_FAILURE);
	}
	if ((outputyzvz = fopen("OUTPUT/vz0.txt", "w")) == NULL)
	{
		// Check for presence of file
		printf("Cannot open %s\n", "outputyzvz.txt");
		exit(EXIT_FAILURE);
	}
	intdxbaropen = 0.0;
	intdvxopen = 0.0;
	intdUopen = 0.0;
	F = 0.0;
	dvzopen = 0.1;
	sizeU = sqrt(Ucap)/dvzopen;
	if (DEBUG == 1) printf("maxj = %d, sizeU = %d\n", maxj, sizeU);
	for (i=0; i<dim_output; i++) {
		vxopen = i*dvxopen;
		intdxbaropenold = intdxbaropen;
		intdxbaropen = 0.0;
		for (j=0; j < maxj; j++)
		{	
			if (i==0)
			{
				fprintf(outputyzvy, "%f\n", xbar[j]);
			}
			intdUopenold = intdUopen;	
			intdUopen = 0.0;
			F = 0.0;
			//sizeU = (int) ( sqrt(Ucap - pow(vxopen, 2.0) - chi[j][0])/dvzopen );
			for (l=0; l<sizeU; l++)
			{	
				vz = dvzopen*l;
				if (j == 0 && i == 0)
				{	fprintf(outputyzvz, "%f\n", vz); }
				Fold = F;
				corr = 0.000;
				vx0open = sqrt(corr + chiMax[j] - chi[j][0]);
				//printf("vx0open = %f\n", vx0open);
				U = pow(vz, 2.0) + pow(vx0open, 2.0) + chi[j][0];
				Deltavx = sqrt(pow(vx0open, 2.0) + alpha*vz*openorbit[j]) -  vx0open;
				aa = tophat(vx0open, vx0open + Deltavx, vxopen);
				//Fopen = bilin_interp(mu[j][0], U, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
				Fopen = bilin_interp(mu[j][0], U-mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1);
				F = Fopen*aa;
				if (i==0) // just because these have to be printed only once, and are independent of i
				{
					if (l==sizeU - 1 && i == 0)
					{	
						fprintf(outputyz, "%f\n", Fopen*Deltavx); 
					}
					else 
					{	
						fprintf(outputyz, "%f ", Fopen*Deltavx); 
					}
				}
				if (DEBUG == 1)
					printf("tophat = %f when vx0open = %f, Deltavx = %f and vxopen = %f\n, F = %f\n", aa, vx0open, Deltavx, vxopen, F);
				if (l!=0) {	
					intdUopen += 2.0*(F+Fold)*dvzopen; 
				}
			}
			if (j != 0) {	
				intdxbaropen += 0.5*(intdUopen + intdUopenold)*(xbar[j] - xbar[j-1]); 
			} 
		}
		if (DEBUG == 1)
			printf("FBohm = %f\n", intdxbaropen);
		if (i!=0)
		{
			intdvxopen += 0.5*(intdxbaropen + intdxbaropenold)*dvxopen;
			intdvxopenfluidBohm += 0.5*(intdxbaropen*vxopen + intdxbaropenold*(vxopen-dvxopen))*dvxopen;
			if (i>1)
			{
				intdvxopenBohm += 0.5*(intdxbaropen/pow(vxopen, 2.0) + intdxbaropenold/pow(vxopen-dvxopen, 2.0))*dvxopen;
			}
			else if (i==1)
			{
				intdvxopenBohm += 0.5*intdxbaropen/pow(vxopen, 2.0)*dvxopen;
			}
		}
		vx_i_DS[i] = vxopen;
		fx_i_DS[i] = intdxbaropen;
	}
	for (i=0; i<dim_output; i++) {	
		fx_i_DS[i] /= intdvxopen;
		//printf("%f %f\n", vx_i_DS[i], fx_i_DS[i]);//nitotal[0]);
		fprintf(output, "%f %f\n", vx_i_DS[i], fx_i_DS[i]);
	}

	fclose(output);
	fclose(outputyz);
	fclose(outputyzvy);
	fclose(outputyzvz);

	// Calculating correction due to non boltzmann electrons
	// AG: to the Bohm condition
	/*the gradient of the function */
	//double ne_p0=0.0;
	//ne_p0 = (ne_grid[1] - ne_grid[0]) / (phi[1] - phi[0]);
	//printf("ne_p0 = %f (%f)\n", ne_p0, (ne_grid[1] - ne_grid[0]) / (phi_grid[1] - phi_grid[0]));

	printf("flux calculated = %f\nflux from continuity = %f\n", *flux, 1.0/nitotal[0]);
	printf("in densfinorb: The density at x=0 obtained from the extracted distribution function is %f\nThe density at x=0 obtained from the ion density integral (more direct) is %f\n", intdvxopen, nitotal[0]);
	printf("Bohm integral = %f (obtained from the extracted distribution function at x=0)\nBohm integral = %f (obtained directly from the distribution function at infinity)\nThe Bohm condition is: Bohm integral = 2*dn_e/dphi/(Te*ni0)\n", intdvxopenBohm/intdvxopen, Bohm);//, (2.0*ne_p0)/(Te*nitotal[0])); 
	printf("in densfinorb: The flow at x=0 obtained from the extracted distribution function is %f\nThe flow at x=0 obtained from the distribution function at infinity is %f\nThe Bohm speed is %f\n", intdvxopenfluidBohm/intdvxopen, *flux, 1.0/sqrt(2.0)); 
	//printf("in densfinorb: kBohmsquared= %f (should be >0)\n", kBohmsquared);

}
else if (output_type ==2) {
	FILE *output;
	if ((output = fopen("OUTPUT/Umucutfile.txt", "w")) == NULL) {       
		printf("Cannot open %s\n", "outputfile.txt");
		exit(EXIT_FAILURE);
	}
	if (DEBUG == 0) printf("muopen Ucritf\n");
	for (j=0;j<dim_output;j++) {
		if (muopen[maxj] < vx_i_DS[j]) fx_i_DS[j] = sqrt(2.0)*Ucritf[maxj];
		else fx_i_DS[j] = sqrt(2.0*lin_interp(muopen, Ucritf, vx_i_DS[j], maxj, 986));
		//else fx_i_DS[j] = sqrt(2.0)*lin_interp(muopen, Ucritf, vx_i_DS[j], maxj, 986);
		fprintf(output, "%f %f\n", vx_i_DS[j], fx_i_DS[j]);
		//if (DEBUG == 1) printf("mu_crit = %f\tvpar_crit = %f\n", vx_i_DS[j], fx_i_DS[j]);
		if (DEBUG == 0) printf("mu_crit = %f\tU_crit = %f\n", vx_i_DS[j], 0.5*fx_i_DS[j]*fx_i_DS[j]);
	}
	fclose(output);
}

for (i=0; i<size_phigrid; i++) {
	n_grid[i] = nitotal[i];
	if (DEBUG == 1) printf("densfinorbspline.c: phi = %f\tnitotal = %f\n", phi_grid[i], nitotal[i]);
}

// If you love your variables (and your memory) set them free
free(xx);//
free(niclosed);//
free(niopen);//
free(jmclosed);//
free(jmopen);//
free(xifunction);//
free(xbar);//
free(chimin);//
free(chiMax);//
free(chimpp);//
free(crossed_min);//
free(crossed_max);//
free(twopimuprime);
free(openintegral);
free(openorbit);
free(upperlimit);
free(lowerlimit);
free(chiprimetop);
free(xtop);
free(itop);
free(imax);
free(imin);


for (int w = 0; w < sizexbar; w++)
{
	free(chi[w]);
	free(chiop[w]);
	free(Uperp[w]);
	free(mu[w]);
	free(upper[w]);//
	for (int s = 0; s < size_finegrid; s++)
	{
		free(vx[w][s]);//
	}
	free(vx[w]);//
}
free(chi);//
free(chiop);//
free(Uperp);//
free(mu);
free(upper);
free(vx);

free(phi);//
free(gg);//
free(ff);//
free(phip);//
free(nitotal);//

for (j=0; j<j_inf; j++) free(vxinf[j]);
free(muinf);
free(vxinf);
free(chiinf);

free(Ucritf);//
free(muopen);//
free(kdrop);//

if (charge > 0.0) {
	for (i=0; i<sizemumu; i++) free(distfunc_iprime[i]);
	free(distfunc_iprime);
}

clock_t end = clock(); // finds the end time of the computation
double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
printf("in densfinorb: module ran in %f seconds\n", jobtime);

return;
} // closes densfinorb function 
