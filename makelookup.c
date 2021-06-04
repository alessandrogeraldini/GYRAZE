// Author: Robbie Ewart 
// MODIFIED: 21/04/2021 by Alessandro Geraldini
/* This code calculates the density integral of a species when this is only a function of the potential at a given point */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mps.h"
#define TINY 1e-8

/* AG
   to use for ions in the Debye sheath, phi must be a list of POSITIVE numbers 
   equal to minus the potential relative to the Debye sheath entrance
   set number pointed at by vpar_cut_lookup to be equal to 0.0
   set size_mu = 0 while size_cut, size_vpar and mue_cut_lookup are redundant (can set to zero/NULL)
*/

void makelookup(double charge, double *phi_real, double *n_grid, int p_size, double *Phi_point, double **distfunc, double *vpar, double *mu, int size_vpar, int size_mu, double *mue_cut_lookup, double *vpar_cut_lookup, int size_cut) {
//defining des variables
int count = 0;
int vi, vi_1, p, len_F; //vi is a counting variable that will be saved for sums over velocity space, vi_1 is a special value in velocity space devoted to the first point, p is a counting variable that will be saved for counting over phi space, len_F saves the number of entries in the distribution;
double *phi;
double Phi=0.0, n_inf = 0.0, v_cut;
double v_max, v_s, v_min;//v_max is the maximum velocity the distribution function will go up to before effectively just reading 0 from there on, v_s is the separation in velocity space between ajacent points. v_cut is the velocity that would be require in order to just reach the plasma wall boundary, v_min is the minimum velocity (the entrance to the presheath) required to reach a given x
double sqrt_up, sqrt_lo, theta_up, theta_lo, sqrt_cut, theta_cut; //sqrt_up/lo reprresent the upper and lower limits of one strip integral, theta_up/lo are related to the hyperbolic arcsinh of some specific values (see document that hopefully exists)

double *F, *Fp, *Fpp, **ddistdvpar, **ddistdvpartwo; // the zeroth first and second derivatives of the distribution function respectively

double* in_err, *af_err, * n_res, w, w_up, w_lo; // for calculating errors in different parts of the code

double* n_pre; //n_pre[x] is an array that will contain the values that we expect the integral to give (for a maxwellian input)

int show_err;// change show_err to 1 or 0 if you do/don't want the error output to be printed to a file (note that the error output will only be meaningful for maxwellian inputs). Change p_option to control how the phi values are distributed with 1 being in an attempt to make n_grid roughly linear and 0 making phi linear;
int mu_ind;
double *nepart, mue_s = 0.05, vpar_cut;
FILE* fptr1;
show_err = 0;

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
		vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu_ind*mue_s, size_cut, 205);
		//printf("in makelookup.c: mu = %f and vpar_cut_DSE = %f\n", mu_ind*mue_s, vpar_cut);
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
			n_inf += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1])*mue_s;
		}
		//printf("count=%d\n", count);
	}
	printf("n_inf = %f\n", n_inf);


	// ELECTRON CURRENT
	for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
		nepart[mu_ind] = 0.0;
		vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu_ind*mue_s, size_cut, 205);
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
			Phi += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1 ])*mue_s;
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
				vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu_ind*mue_s, size_cut, 205);
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
					n_grid[p] += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1])*mue_s;
					//printf("n_grid %f\n", n_grid[p]);
				}
			}
		}
	}
}
else {
	for (vi = 0; vi < len_F - 1; vi++) {
		F[vi]     = distfunc[0][vi];
		Fp[vi]    = ddistdvpar[0][vi];
		Fpp[vi]   = ddistdvpartwo[0][vi];
		//printf("F = %f, Fp = %f, Fpp = %f\n", F[vi], Fp[vi], Fpp[vi]);
	}
	vpar_cut = *vpar_cut_lookup;
	if (phi[0] < 0.0) { // the case for electrons in MPS, not the case for ions in DS
		v_cut = sqrt(vpar_cut*vpar_cut + 2.0*(-phi[0]));
	}
	else v_cut = 0.0;



	// Normalization at infinity
	n_inf = 0.0;
	v_min = 0.0;

	if (phi[0] > 0.0) {
		n_inf = 0.0;	
		for (vi = 0; vi < len_F - 2; vi++)
		{
			n_inf += 0.5*(F[vi] + F[vi + 1]) * v_s;
		}
	}
	else {
	if (v_cut >= v_max) { //effectively no cut off
		v_min = 0.0;
		if ((int)floor(v_min / v_s) >= len_F - 1)
		{
			n_inf = 0;
			printf("v_min outside velocity grid");
		}
		else
		{
			if (fabs(v_min) < TINY)//special case where the integration becomes simple
			{
				for (vi = 0; vi < len_F - 2; vi++)
				{
					n_inf += (F[vi] + F[vi + 1]) * v_s;
				}
			}
			else //regualar case
			{
				vi_1 = (int)floor(v_min / v_s);
				sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0));
				theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0)) / v_min);

				n_inf += 2.0 * (((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
				n_inf += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));

				w = ((0.5 * pow((((vi_1 + 1) * v_s) - v_min), 2.0)) + (v_min * (((vi_1 + 1) * v_s) - v_min)));

				for (vi = vi_1 + 1; vi < len_F - 1; vi++)
				{
					sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0));
					sqrt_lo = sqrt(pow(((vi)* v_s), 2.0));
					theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0)) / (v_min));
					theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0)) / (v_min));

					n_inf += 2.0 * (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
					n_inf += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

					w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
					w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
				}
			}
		}
	}
	else {
		// for ions in DS always end up here
		//now instead there is a cutoff
			v_min = 0.0;
			vi_1 = 0;
			sqrt_up = ((vi_1 + 1) * v_s);

			n_inf += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
			n_inf += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));
			w = 0.5 * pow(v_s, 2.0);

			for (vi = vi_1 + 1; vi < len_F - 1; vi++)
			{
				sqrt_up = (vi + 1) * v_s;
				sqrt_lo = (vi * v_s);

				n_inf += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
				n_inf += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

				w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
				w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
			}

			for (vi = vi_1 + 1; vi < (int)floor(v_cut / v_s); vi++)
			{
				sqrt_up = (vi + 1) * v_s;
				sqrt_lo = (vi * v_s);

				n_inf += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
				n_inf += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

				w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
				w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
			}
			vi_1 = (int)floor(v_cut / v_s);
			sqrt_up = v_cut;
			sqrt_lo = vi_1 * v_s;

			n_inf += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
			n_inf += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));

			w_up = 0.5 * pow(v_cut, 2.0);
			w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
	}
	}

	printf("n_inf = %f\n", n_inf);

	// ELECTRON CURRENT
	// fill in


	if (phi[0] > 0.0) {
		for (p = 0;p < p_size; p++) {
			n_grid[p] = 0.0;
			for (vi = 0; vi < len_F - 2; vi++) {
				n_grid[p] += 0.5*(F[vi] + F[vi + 1]) * ( sqrt((vi+1)*v_s*(vi+1)*v_s + 2.0*phi[p]) - sqrt((vi)*v_s*(vi)*v_s + 2.0*phi[p]) );
			}
		}
	}
	else {
	if (v_cut >= v_max) { //effectively no cut off
		for (p = 0; p < p_size; p++) {
			n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] + (0.5 * pow(v_cut, 2.0)))));
			v_min = sqrt(-2.0 * phi[p]);
			//printf("v_min = %f\n", v_min);
			if ((int)floor(v_min / v_s) >= len_F - 1)
			{
				n_grid[p] = 0.0;
				printf("v_min outside velocity grid");
			}
			else
			{
				n_grid[p] = 0.0;
				if (fabs(phi[p]) < TINY) {
					n_grid[p] = n_inf;
				}
				else if (fabs(v_min) < TINY)//special case where the integration becomes simple
				{
					for (vi = 0; vi < len_F - 2; vi++)
					{
						n_grid[p] += (F[vi] + F[vi + 1]) * v_s;
					}
				}
				else //regualar case
				{
					vi_1 = (int)floor(v_min / v_s);
					sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p]));
					theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p])) / v_min);

					n_grid[p] += 2.0 * (((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
					n_grid[p] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));

					w = ((0.5 * pow((((vi_1 + 1) * v_s) - v_min), 2.0)) + (v_min * (((vi_1 + 1) * v_s) - v_min)));
					in_err[p] += 2.0 * ((((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w))));
					in_err[p] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));

					for (vi = vi_1 + 1; vi < len_F - 1; vi++)
					{
						sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
						sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
						theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
						theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

						n_grid[p] += 2.0 * (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
						n_grid[p] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
						w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
						af_err[p] += 2.0 * ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
						af_err[p] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
					}
				}

			}
			n_res[p] = n_grid[p] - n_pre[p];
			//printf("n_grid[%d/%d] = %f\n", p, p_size, n_grid[p]);
		}
	}
	else {
		// for ions in DS always end up here
		//now instead there is a cutoff
		for (p = 0;p < p_size; p++)
		{
			n_grid[p] = 0.0;
			if ((phi[p] + (0.5 * pow(v_cut, 2.0))) > 0)
			{
				n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] + (0.5 * pow(v_cut, 2.0)))));
			}
			else
			{
				n_pre[p] = 0.5 * exp(phi[p]);
			}
			//if (phi[0] < 0.0) {
			v_min = sqrt(-2.0 * phi[p]);
			//}
			//else { 
			//	v_min = 2.0*TINY;
			//}
			if ((int)floor(v_min / v_s) >= len_F - 1)
			{
				n_grid[p] = 0.0;
			}
			else
			{
				if (fabs(v_min) < TINY)
				{
					vi_1 = 0;
					sqrt_up = ((vi_1 + 1) * v_s);

					n_grid[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
					n_grid[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));
					w = 0.5 * pow(v_s, 2.0);
					af_err[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s)))) - erf(sqrt(w));
					af_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));

					for (vi = vi_1 + 1; vi < len_F - 1; vi++)
					{
						sqrt_up = (vi + 1) * v_s;
						sqrt_lo = (vi * v_s);

						n_grid[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
						n_grid[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
						w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
						af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
						af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
					}

					for (vi = vi_1 + 1; vi < (int)floor(v_cut / v_s); vi++)
					{
						sqrt_up = (vi + 1) * v_s;
						sqrt_lo = (vi * v_s);

						n_grid[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
						n_grid[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
						w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
						af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
						af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
					}
					vi_1 = (int)floor(v_cut / v_s);
					sqrt_up = v_cut;
					sqrt_lo = vi_1 * v_s;

					n_grid[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
					n_grid[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));

					w_up = 0.5 * pow(v_cut, 2.0);
					w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
					af_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
					af_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
				}
				else
				{
					//printf("am I always here for %d/%d\n", p, p_size);
					if ((int)floor(v_min / v_s) == (int)floor(v_cut / v_s))
					{
						vi_1 = (int)floor(v_min / v_s);
						sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p]));
						theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p])) / (v_min));

						n_grid[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
						n_grid[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

						w_up = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
						in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w_up))));
						in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

						sqrt_cut = sqrt(pow(v_cut, 2.0) + 2.0 * phi[p]);
						theta_cut = asinh(sqrt(pow(v_cut, 2.0) + 2.0 * phi[p]) / (v_min));

						n_grid[p] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * v_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
						n_grid[p] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (v_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));

						w_lo = 0.5 * (pow(v_cut, 2.0) - pow(v_min, 2.0));
						in_err[p] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * v_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
						in_err[p] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (v_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));

						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
						{
							sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
							sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
							theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
							theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

							n_grid[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							n_grid[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
							w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
							af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
							af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
						}
					}
					else //the most likely case
					{
						vi_1 = (int)floor(v_min / v_s);
						sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phi[p]);
						theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));

						n_grid[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
						n_grid[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

						w = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
						in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w))));
						in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
						{
							sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
							sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
							theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
							theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

							n_grid[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							n_grid[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
							w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
							af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
							af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
						}

						for (vi = vi_1 + 1; vi < (int)floor(v_cut / v_s); vi++)
						{
							sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
							sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
							theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
							theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

							n_grid[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							n_grid[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
							w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
							af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
							af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
						}

						vi_1 = (int)floor(v_cut / v_s);
						sqrt_up = sqrt(pow(v_cut, 2.0) + 2.0 * phi[p]);
						sqrt_lo = sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phi[p]);
						theta_up = asinh(sqrt(pow(v_cut, 2.0) + 2.0 * phi[p]) / (v_min));
						theta_lo = asinh(sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

						n_grid[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
						n_grid[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));

						w_up = 0.5 * (pow(v_cut, 2.0) - pow(v_min, 2.0));
						w_lo = 0.5 * (pow((vi_1 * v_s), 2.0) - pow(v_min, 2.0));
						in_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
						in_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));



					}
				}
			}
			n_res[p] = n_grid[p] - n_pre[p];
			//printf("n_grid[%d] = %f\n", p, n_grid[p]);
		}
	}
	}
}


Phi =  Phi / n_inf;
printf("Phi = %f\n", Phi);
printf("n_inf = %f\n", n_inf);
for (p = 0;p < p_size; p++)
{
	//n_grid[p] = n_grid[p] / n_grid[p_size - 1];
	n_grid[p] = n_grid[p] / n_inf;
	//if (n_grid[p] != n_grid[p]) n_grid[p] = 1.0;
	n_res[p] = n_grid[p] - n_pre[p];
	//printf("in makelookup.c: n[%d/%d] = %f, phi = %f\n", p, p_size, n_grid[p], phi[p]);
}

//only works if size_mu = 1 (or for the final value of mu)
if (show_err == 1)
{
	fptr1 = fopen("total_error.txt", "w");
	for (p = 0; p < p_size; p++)
	{
		fprintf(fptr1, "%.15f\n", n_res[p]);
	}
	fclose(fptr1);
	fptr1 = fopen("regular_point_error.txt", "w");
	for (p = 0; p < p_size; p++)
	{
		fprintf(fptr1, "%.15f\n", af_err[p]);
	}
	fclose(fptr1);
	fptr1 = fopen("initial_point_error.txt", "w");
	for (p = 0; p < p_size; p++)
	{
		fprintf(fptr1, "%.15f\n", in_err[p]);
	}
	fclose(fptr1);
	fptr1 = fopen("phi.txt", "w");
	for (p = 0; p < p_size; p++)
	{
		fprintf(fptr1, "%.15f\n", phi[p]);
	}
	fclose(fptr1);
	fptr1 = fopen("Edens.txt", "w");
	for (p = 0; p < p_size; p++)
	{
		fprintf(fptr1, "%.15f\n", in_err[p]);
	}
	fclose(fptr1);
}
printf("line 1163\n");
for (mu_ind = 0; mu_ind < size_mu; mu_ind++) {
	free(ddistdvpar[mu_ind]);
	free(ddistdvpartwo[mu_ind]);
}
free(ddistdvpar);
free(ddistdvpartwo);
free(F);
free(Fp);
free(Fpp);
free(in_err);
free(af_err);
free(n_res);
free(n_pre);
if (size_mu > 1) free(nepart);

*Phi_point = Phi;

printf("end of makelookup\n");
return;
}
