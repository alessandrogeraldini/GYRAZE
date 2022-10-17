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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "mps_renorm.h"


//double tophat(double x1, double x2, double x) {
//	double y;
//	if ((x1 <= x ) && (x <= x2 )) 
//		y = 1.0;
//	else
//		y = 0.0;
//	return y;
//}


// denszeroorb (previously makelookup) calculates the density integral of a species when this is only a function of the potential at a given point 

/* AG
   to use for ions in the Debye sheath, phi must be a list of POSITIVE numbers 
   equal to minus the potential relative to the Debye sheath entrance
   set number pointed at by vpar_cut_lookup to be equal to 0.0
   set size_mu = 0 while size_cut, size_vpar and mue_cut_lookup are redundant (can set to zero/NULL)
*/

//void denszeroorb(double charge, double TeovTs, double *phi_real, double *n_grid, int p_size, double *Phi_point, double **distfunc, double *vpar, double *mu, int size_vpar, int size_mu, double *vpar_cut_lookup, double gamma, double *x_grid) {
//	//define variables
//	int count = 0;
//	int vi, vi_1, p, len_F; //vi is a counting variable that will be saved for sums over velocity space, vi_1 is a special value in velocity space devoted to the first point, p is a counting variable that will be saved for counting over phi space, len_F saves the number of entries in the distribution;
//	double *phi, *vparacc, phip = -0.000000001 ;
//	double Phi=0.0, n_inf = 0.0;
//	double v_max, v_s, v_min;//v_max is the maximum velocity the distribution function will go up to before effectively just reading 0 from there on, v_s is the separation in velocity space between ajacent points. v_cut is the velocity that would be require in order to just reach the plasma wall boundary, v_min is the minimum velocity (the entrance to the presheath) required to reach a given x
//	double sqrt_up, sqrt_lo, theta_up, theta_lo, sqrt_cut, theta_cut; //sqrt_up/lo reprresent the upper and lower limits of one strip integral, theta_up/lo are related to the hyperbolic arcsinh of some specific values (see document that hopefully exists)
//
//	double *F, *Fp, *Fpp, **ddistdvpar, **ddistdvpartwo; // the zeroth first and second derivatives of the distribution function respectively
//
//	double* in_err, *af_err, * n_res, w, w_up, w_lo; // for calculating errors in different parts of the code
//
//	double* n_pre; //n_pre[x] is an array that will contain the values that we expect the integral to give (for a maxwellian input)
//
//	int show_err;// change show_err to 1 or 0 if you do/don't want the error output to be printed to a file (note that the error output will only be meaningful for maxwellian inputs). Change p_option to control how the phi values are distributed with 1 being in an attempt to make n_grid roughly linear and 0 making phi linear;
//	int mu_ind;
//	double *nepart, vpar_cut;
//	FILE* fptr1;
//	show_err = 0;
//
//	vparacc = malloc(size_vpar*sizeof(double));
//	phi = malloc(p_size*sizeof(double));
//	for (p=0; p<p_size; p++) {
//		phi[p] = -phi_real[p]*charge*TeovTs;
//	}
//
//
//	ddistdvpar = malloc(size_mu * sizeof(double));
//	ddistdvpartwo = malloc(size_mu * sizeof(double));
//	in_err = malloc(p_size * sizeof(double));
//	af_err = malloc(p_size * sizeof(double));
//	n_res = malloc(p_size * sizeof(double));
//	n_pre = malloc(p_size * sizeof(double));
//	len_F = size_vpar;
//	v_max = vpar[len_F - 1];
//	v_s = v_max / (len_F - 1);
//
//	F = malloc(len_F * sizeof(double));
//	Fp = malloc(len_F * sizeof(double));
//	Fpp = malloc(len_F * sizeof(double));
//	//Generate the first and second derivative of the distribution function
//	for (mu_ind = 0; mu_ind < size_mu; mu_ind++) {
//		ddistdvpar[mu_ind] = malloc(len_F * sizeof(double));
//		ddistdvpartwo[mu_ind] = malloc(len_F * sizeof(double));
//		for (vi = 0; vi < len_F; vi++)
//		{
//			if ((vi < len_F - 1) && (vi > 0))
//			{
//				ddistdvpar[mu_ind][vi] = (distfunc[mu_ind][vi + 1] - distfunc[mu_ind][vi - 1]) / (2.0 * v_s);
//				ddistdvpartwo[mu_ind][vi] = (distfunc[mu_ind][vi + 1] + distfunc[mu_ind][vi - 1] - (2.0 * distfunc[mu_ind][vi])) / (pow(v_s, 2.0));
//			}
//			else
//			{
//				if (vi < len_F - 1)
//				{
//					ddistdvpar[mu_ind][vi] = (distfunc[mu_ind][vi + 1] - distfunc[mu_ind][vi]) / (v_s);
//					ddistdvpartwo[mu_ind][vi] = (distfunc[mu_ind][vi + 2] + distfunc[mu_ind][vi] - (2.0 * distfunc[mu_ind][vi + 1])) / (pow(v_s, 2.0));
//				}
//				else
//				{
//					if (vi > 0)
//					{
//						ddistdvpar[mu_ind][vi] = (distfunc[mu_ind][vi] - distfunc[mu_ind][vi - 1]) / (v_s);
//						ddistdvpartwo[mu_ind][vi] = (distfunc[mu_ind][vi - 2] + distfunc[mu_ind][vi] - (2.0 * distfunc[mu_ind][vi - 1])) / (pow(v_s, 2.0));
//					}
//					else
//					{
//						printf("Fuck");
//					}
//				}
//
//			}
//		}
//	}
//
//	/* comment by AG
//	the integration has two options 
//	if size_mu is zero then the distribution function is 1D 
//	e.g. rho_e << lambda_D model of electrons in magnetic presheath x ~ rho_i (velocity coordinate = v_parallel) or ions in the Debye sheath x ~ lambda_D << rho_i (velocity coordinate = v_x)
//	if size_mu is non-zero then the distribution function is 2D but, differently from the integration in densfinorb.c, is still a local function of x
//	e.g. rho_e >~ lambda_D for electrons in the magnetic presheath x ~ rho_i
//	*/
//	if (size_mu > 1) {
//		printf("In 2D integration\n");
//		nepart = malloc(size_mu*sizeof(double));
//		// for finite electron gyroradius effects, we need to integrate in mu as well
//		// cutoff function for large rhoe
//
//		// Normalization at infinity
//		n_inf = 0.0;
//		v_min = 0.0;
//		for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
//			count = 0;
//			nepart[mu_ind] = 0.0;
//			vpar_cut = vpar_cut_lookup[mu_ind]; 
//			//lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
//			//printf("in makelookup.c: mu = %f and vpar_cut_DSE = %f\n", mu[mu_ind], vpar_cut);
//			vpar_cut = sqrt(vpar_cut*vpar_cut + 2.0*(-phi[0]) + TINY);
//			for (vi = 0; vi < len_F - 1; vi++) {
//				//F[vi]   = exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
//				//Fp[vi]   = -vi*v_s*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
//				//Fpp[vi]   = (vi*v_s*vi*v_s - 1.0)*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
//				F[vi]     = distfunc[mu_ind][vi];
//				Fp[vi]    = ddistdvpar[mu_ind][vi];
//				Fpp[vi]   = ddistdvpartwo[mu_ind][vi];
//			}
//			if (vpar_cut >= v_max) // effectively no cut off
//			{
//				for (vi = 0; vi < len_F - 2; vi++)
//					nepart[mu_ind] += (F[vi] + F[vi + 1]) * v_s;
//				count++;
//			}
//			else { // now there is a cut-off
//				vi_1 = 0;
//				sqrt_up = ((vi_1 + 1) * v_s);
//				nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
//				nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));
//
//				w = 0.5 * pow(v_s, 2.0);
//
//				for (vi = vi_1 + 1; vi < len_F - 1; vi++) {
//					sqrt_up = (vi + 1) * v_s;
//					sqrt_lo = (vi * v_s);
//
//					nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
//					nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//
//					w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//					w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
//				}
//
//				for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++) {
//					sqrt_up = (vi + 1) * v_s;
//					sqrt_lo = (vi * v_s);
//
//					nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
//					nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//
//					w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//					w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
//				}
//				vi_1 = (int)floor(vpar_cut / v_s);
//				sqrt_up = vpar_cut;
//				sqrt_lo = vi_1 * v_s;
//
//				nepart[mu_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
//				nepart[mu_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
//
//				w_up = 0.5 * pow(vpar_cut, 2.0);
//				w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
//			}
//			//printf("in makelookup.c: mu = %f and vpar_cut_DSE = %f, nepart[mu_ind] = %f\n", mu[mu_ind], vpar_cut, nepart[mu_ind]);
//			if (mu_ind != 0) {
//				n_inf += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1])*(mu[mu_ind] - mu[mu_ind-1]);
//			}
//			//printf("count=%d\n", count);
//		}
//		printf("n_inf = %f\n", n_inf);
//		
//	//else {
//		n_inf = 0.0;
//		v_min = sqrt(-2.0 * phip);
//		for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
//			nepart[mu_ind] = 0.0;
//			//vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
//			vpar_cut = vpar_cut_lookup[mu_ind] ; 
//			//vpar_cut = sqrt(vpar_cut*vpar_cut + TINY);
//			vpar_cut = sqrt(vpar_cut*vpar_cut + 2.0*(-phi[0]) + TINY);
//			//printf("vpar_cut = %f\n", vpar_cut);
//			for (vi = 0; vi < len_F - 1; vi++) {
//				//F[vi]   = exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
//				//Fp[vi]   = -vi*v_s*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
//				//Fpp[vi]   = (vi*v_s*vi*v_s - 1.0)*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
//				F[vi]     = distfunc[mu_ind][vi];
//				Fp[vi]    = ddistdvpar[mu_ind][vi];
//				Fpp[vi]   = ddistdvpartwo[mu_ind][vi];
//			}
//			if (vpar_cut >= v_max) // effectively no cut off
//			{
//				if ((int)floor(v_min / v_s) >= len_F - 1) {
//					//ne[p] = 0;
//					nepart[mu_ind] = 0.0;
//					printf("v_min outside velocity grid");
//				}
//				else {
//					if (fabs(v_min) < 2.0*TINY)//special case where the integration becomes simple
//					{
//						for (vi = 0; vi < len_F - 2; vi++)
//						{
//							nepart[mu_ind] += (F[vi] + F[vi + 1]) * v_s;
//						}
//					}
//					else //regular case
//					{
//
//
//						vi_1 = (int)floor(v_min / v_s);
//						sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phip));
//						theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phip)) / v_min);
//
//						nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//						nepart[mu_ind] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//
//						w = ((0.5 * pow((((vi_1 + 1) * v_s) - v_min), 2.0)) + (v_min * (((vi_1 + 1) * v_s) - v_min)));
//
//						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//						{
//							sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip);
//							sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip);
//							theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip) / (v_min));
//							theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip ) / (v_min));
//
//							nepart[mu_ind] += 2.0 * (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//							nepart[mu_ind] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//							w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//						}
//					}
//				}
//			//n_res[p] = ne[p] - n_pre[p];
//			}
//			else { // now there is a cut-off
//				v_min = sqrt(-2.0 * phip);
//				if ((int)floor(v_min / v_s) >= len_F - 1)
//				{
//					nepart[mu_ind] = 0.0;
//				}
//				else
//				{
//					if (fabs(v_min) < 1e-9)
//					{
//						vi_1 = 0;
//						sqrt_up = ((vi_1 + 1) * v_s);
//						nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
//						nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));
//
//						w = 0.5 * pow(v_s, 2.0);
//
//						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//						{
//							sqrt_up = (vi + 1) * v_s;
//							sqrt_lo = (vi * v_s);
//
//							nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
//							nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//
//							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//							w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
//						}
//
//						for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
//						{
//							sqrt_up = (vi + 1) * v_s;
//							sqrt_lo = (vi * v_s);
//
//							nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
//							nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//
//							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//							w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
//						}
//						vi_1 = (int)floor(vpar_cut / v_s);
//						sqrt_up = vpar_cut;
//						sqrt_lo = vi_1 * v_s;
//
//						nepart[mu_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
//						nepart[mu_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
//
//						w_up = 0.5 * pow(vpar_cut, 2.0);
//						w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
//					}
//					else
//					{
//						if ((int)floor(v_min / v_s) == (int)floor(vpar_cut / v_s))
//						{
//							vi_1 = (int)floor(v_min / v_s);
//							sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phip));
//							theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phip)) / (v_min));
//
//							nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//							nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));
//
//							w_up = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
//
//							sqrt_cut = sqrt(pow(vpar_cut, 2.0) + 2.0 * phip);
//							theta_cut = asinh(sqrt(pow(vpar_cut, 2.0) + 2.0 * phip) / (v_min));
//
//							nepart[mu_ind] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * vpar_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
//							nepart[mu_ind] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (vpar_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
//
//							w_lo = 0.5 * (pow(vpar_cut, 2.0) - pow(v_min, 2.0));
//
//							for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//							{
//								sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip);
//								sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip);
//								theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip) / (v_min));
//								theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip) / (v_min));
//
//								nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//								nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//								w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//								w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//							}
//						}
//						else //the most likely case
//						{
//							vi_1 = (int) (floor(v_min / v_s));
//							sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phip);
//							theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phip) / (v_min));
//							nepart[mu_ind] = 0.0;
//							
//
//							nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//							nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));
//
//							w = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
//
//							for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//							{
//								sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip);
//								sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip);
//								theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip) / (v_min));
//								theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip) / (v_min));
//
//								nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//								nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//								w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//								w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//							}
//
//							for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
//							{
//								sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip);
//								sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip);
//								theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phip) / (v_min));
//								theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phip) / (v_min));
//
//								nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//								nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//								w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//								w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//							}
//
//							vi_1 = (int)floor(vpar_cut / v_s);
//							sqrt_up = sqrt(pow(vpar_cut, 2.0) + 2.0 * phip);
//							sqrt_lo = sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phip);
//							theta_up = asinh(sqrt(pow(vpar_cut, 2.0) + 2.0 * phip) / (v_min));
//							theta_lo = asinh(sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phip) / (v_min));
//
//							nepart[mu_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//							nepart[mu_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));
//
//							w_up = 0.5 * (pow(vpar_cut, 2.0) - pow(v_min, 2.0));
//							w_lo = 0.5 * (pow((vi_1 * v_s), 2.0) - pow(v_min, 2.0));
//						}
//					}
//				//n_res[p] = ne[p] - n_pre[p];
//				}
//			}
//			if (mu_ind != 0) {
//				n_inf += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1])*(mu[mu_ind] - mu[mu_ind-1]);
//				//printf("n_grid %f\n", n_grid[p]);
//			}
//		}
//	//n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] -phi[0] )));
//	//printf("n_grid[%d] = %f, n_pre = %f\n", p, n_grid[p], n_pre[p]);
//	//}
//
//		// ELECTRON CURRENT
//		for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
//			nepart[mu_ind] = 0.0;
//			//vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
//			vpar_cut = vpar_cut_lookup[mu_ind] ; 
//			//vpar_cut = sqrt(vpar_cut*vpar_cut + TINY);
//			vpar_cut = sqrt(vpar_cut*vpar_cut + 2.0*(-phi[0]) + TINY);
//			for (vi = 0; vi < len_F - 1; vi++) {
//				//F[vi]   = vi*v_s*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
//				//Fp[vi]  = (1.0 - vi*v_s*vi*v_s)*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
//				//Fpp[vi] = (-3.0*vi*v_s + vi*v_s*vi*v_s*vi*v_s )*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
//				F[vi]   = vi*v_s*distfunc[mu_ind][vi];
//				Fp[vi]  = vi*v_s*ddistdvpar[mu_ind][vi] + distfunc[mu_ind][vi];
//				Fpp[vi] = vi*v_s*ddistdvpartwo[mu_ind][vi] + 2.0*ddistdvpar[mu_ind][vi];
//				//printf("F[%d]=%f (should be %f)\n", vi, F[vi], vi*v_s*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s));
//			}
//			//printf("v_max = %f\tvpar_cut = %f\n", v_max, vpar_cut);
//			if (vpar_cut >= v_max) // effectively no cut off
//				Phi = 0.0;
//			else { // now there is a cut-off
//				v_min = 0.0;
//				if ((int)floor(v_min / v_s) >= len_F - 1)
//				{
//					nepart[mu_ind] = 0.0;
//				}
//				else
//				{
//					vi_1 = -1;
//					for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//					{
//						sqrt_up = (vi + 1) * v_s;
//						sqrt_lo = (vi * v_s);
//
//						nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
//						nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//
//						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//						w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
//					}
//
//					for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
//					{
//						sqrt_up = (vi + 1) * v_s;
//						sqrt_lo = (vi * v_s);
//
//						nepart[mu_ind] -= (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
//						nepart[mu_ind] -= ( ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))) );
//
//						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//						w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
//					}
//					vi_1 = (int)floor(vpar_cut / v_s);
//					sqrt_up = vpar_cut;
//					sqrt_lo = vi_1 * v_s;
//
//					nepart[mu_ind] -= (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
//					nepart[mu_ind] -= ( ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))) );
//
//					w_up = 0.5 * pow(vpar_cut, 2.0);
//					w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
//				}
//			}
//			if (mu_ind != 0) 
//				Phi += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1 ])*(mu[mu_ind]-mu[mu_ind-1]);
//		}
//
//		// NOW CALCULATE THE DENSITY PROFILE
//		for (p = 0; p < p_size; p++)
//		{
//			//printf("phi[%d] = %f\n", p, phi[p]);
//			if (fabs(phi[p]) < TINY) {
//				n_grid[p] = n_inf;
//			}
//			else {
//				n_grid[p] = 0.0;
//				v_min = sqrt(-2.0 * phi[p]);
//				for (mu_ind=0; mu_ind<size_mu; mu_ind++) {
//					nepart[mu_ind] = 0.0;
//					//vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mu[mu_ind], size_cut, 205);
//					vpar_cut = vpar_cut_lookup[mu_ind] ; 
//					vpar_cut = sqrt(vpar_cut*vpar_cut +  (-2.0*phi[0]) + TINY);
//					//printf("vpar_cut = %f\n", vpar_cut);
//					for (vi = 0; vi < len_F - 1; vi++) {
//						//F[vi]   = exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
//						//Fp[vi]   = -vi*v_s*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
//						//Fpp[vi]   = (vi*v_s*vi*v_s - 1.0)*exp(- mu_ind*mue_s - 0.5*vi*vi*v_s*v_s);
//						F[vi]     = distfunc[mu_ind][vi];
//						Fp[vi]    = ddistdvpar[mu_ind][vi];
//						Fpp[vi]   = ddistdvpartwo[mu_ind][vi];
//					}
//					if (vpar_cut >= v_max) // effectively no cut off
//					{
//						if ((int)floor(v_min / v_s) >= len_F - 1)
//						{
//							//ne[p] = 0;
//							nepart[mu_ind] = 0.0;
//							printf("v_min outside velocity grid");
//						}
//						else
//						{
//							if (fabs(v_min) < 2.0*TINY)//special case where the integration becomes simple
//							{
//								for (vi = 0; vi < len_F - 2; vi++)
//								{
//									nepart[mu_ind] += (F[vi] + F[vi + 1]) * v_s;
//								}
//							}
//							else //regular case
//							{
//
//
//								vi_1 = (int)floor(v_min / v_s);
//								sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p]));
//								theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p])) / v_min);
//
//								nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//								nepart[mu_ind] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//
//								w = ((0.5 * pow((((vi_1 + 1) * v_s) - v_min), 2.0)) + (v_min * (((vi_1 + 1) * v_s) - v_min)));
//								in_err[p] += 2.0 * ((((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w))));
//								in_err[p] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//
//								for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//								{
//									sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
//									sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
//									theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//									theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//
//									nepart[mu_ind] += 2.0 * (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//									nepart[mu_ind] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//									w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//									w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//									af_err[p] += 2.0 * ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
//									af_err[p] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//								}
//							}
//						}
//						//n_res[p] = ne[p] - n_pre[p];
//					}
//					else { // now there is a cut-off
//						v_min = sqrt(-2.0 * phi[p]);
//						if ((int)floor(v_min / v_s) >= len_F - 1)
//						{
//							nepart[mu_ind] = 0.0;
//						}
//						else
//						{
//							if (fabs(v_min) < 1e-9)
//							{
//								vi_1 = 0;
//								sqrt_up = ((vi_1 + 1) * v_s);
//								nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
//								nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));
//
//								w = 0.5 * pow(v_s, 2.0);
//								af_err[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s)))) - erf(sqrt(w));
//								af_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));
//
//								for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//								{
//									sqrt_up = (vi + 1) * v_s;
//									sqrt_lo = (vi * v_s);
//
//									nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
//									nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//
//									w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//									w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
//									af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
//									af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//								}
//
//								for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
//								{
//									sqrt_up = (vi + 1) * v_s;
//									sqrt_lo = (vi * v_s);
//
//									nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
//									nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//
//									w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//									w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
//									af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
//									af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
//								}
//								vi_1 = (int)floor(vpar_cut / v_s);
//								sqrt_up = vpar_cut;
//								sqrt_lo = vi_1 * v_s;
//
//								nepart[mu_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
//								nepart[mu_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
//
//								w_up = 0.5 * pow(vpar_cut, 2.0);
//								w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
//								af_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
//								af_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
//							}
//							else
//							{
//								if ((int)floor(v_min / v_s) == (int)floor(vpar_cut / v_s))
//								{
//									vi_1 = (int)floor(v_min / v_s);
//									sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p]));
//									theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p])) / (v_min));
//
//									nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//									nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));
//
//									w_up = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
//									in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w_up))));
//									in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));
//
//									sqrt_cut = sqrt(pow(vpar_cut, 2.0) + 2.0 * phi[p]);
//									theta_cut = asinh(sqrt(pow(vpar_cut, 2.0) + 2.0 * phi[p]) / (v_min));
//
//									nepart[mu_ind] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * vpar_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
//									nepart[mu_ind] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (vpar_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
//
//									w_lo = 0.5 * (pow(vpar_cut, 2.0) - pow(v_min, 2.0));
//									in_err[p] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * vpar_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
//									in_err[p] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (vpar_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
//
//									for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//									{
//										sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
//										sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
//										theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//										theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//
//										nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//										nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//										w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//										w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//										af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
//										af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//									}
//								}
//								else //the most likely case
//								{
//									vi_1 = (int) (floor(v_min / v_s));
//									sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phi[p]);
//									theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//									nepart[mu_ind] = 0.0;
//									
//
//									nepart[mu_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
//									nepart[mu_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));
//
//									w = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
//									in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w))));
//									in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));
//
//									for (vi = vi_1 + 1; vi < len_F - 1; vi++)
//									{
//										sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
//										sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
//										theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//										theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//
//										nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//										nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//										w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//										w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//										af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
//										af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//									}
//
//									for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
//									{
//										sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
//										sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
//										theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//										theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//
//										nepart[mu_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//										nepart[mu_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//
//										w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
//										w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
//										af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
//										af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//									}
//
//									vi_1 = (int)floor(vpar_cut / v_s);
//									sqrt_up = sqrt(pow(vpar_cut, 2.0) + 2.0 * phi[p]);
//									sqrt_lo = sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phi[p]);
//									theta_up = asinh(sqrt(pow(vpar_cut, 2.0) + 2.0 * phi[p]) / (v_min));
//									theta_lo = asinh(sqrt(pow(((vi_1)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));
//
//									nepart[mu_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
//									nepart[mu_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));
//
//									w_up = 0.5 * (pow(vpar_cut, 2.0) - pow(v_min, 2.0));
//									w_lo = 0.5 * (pow((vi_1 * v_s), 2.0) - pow(v_min, 2.0));
//									in_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
//									in_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));
//								}
//							}
//						//n_res[p] = ne[p] - n_pre[p];
//						}
//					}
//					//if (nepart[mu_ind] != nepart[mu_ind]) printf("mu_ind = %d, p = %d\n", mu_ind, p);
//
//					if ( (gamma < SMALLGAMMA) && (gamma > TINY) ) nepart[mu_ind] *= smallrhoefrac(x_grid[p], mu[mu_ind]);
//					if (mu_ind != 0) {
//						n_grid[p] += 2.0*M_PI*0.5*(nepart[mu_ind] + nepart[mu_ind -1])*(mu[mu_ind] - mu[mu_ind-1]);
//						//printf("n_grid %f\n", n_grid[p]);
//					}
//				}
//			}
//			//n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] -phi[0] )));
//			//printf("n_grid[%d] = %f, n_pre = %f\n", p, n_grid[p], n_pre[p]);
//		}
//		free(nepart);
//	}
//	else {
//		printf("In 1D ion integration\n");
//		n_inf = 0.0;
//		for (vi = 0; vi < len_F - 1; vi++) {
//			F[vi] = distfunc[0][vi];
//			if (vi != 0)
//				n_inf += 0.5*(F[vi] + F[vi-1])*(vpar[vi] - vpar[vi-1]);
//		}
//		for (p=0; p<p_size;p++) {
//			n_grid[p] = 0.0;
//			//if (phi[p] < 0.0) n_grid[p] = 1.0;
//			//else {
//				for (vi = 0; vi < len_F - 1; vi++) {
//					vparacc[vi] = sqrt(vpar[vi]*vpar[vi] + 2.0*phi[p]);
//					if (vi != 0)
//						n_grid[p] += 0.5*(F[vi] + F[vi-1])*(vparacc[vi] - vparacc[vi-1]);
//				}
//			//}
//		}
//	}
//
//	Phi =  Phi / n_inf;
//	printf("Phi = %f\n", Phi);
//	printf("n_inf = %f\n", n_inf);
//	for (p = 0;p < p_size; p++)
//	{
//		//n_grid[p] = n_grid[p] / n_grid[p_size - 1];
//		n_grid[p] = n_grid[p] / n_inf;
//		//if (n_grid[p] != n_grid[p]) n_grid[p] = 1.0;
//		n_res[p] = n_grid[p] - n_pre[p];
//		if ( (DEBUG == 1) && (x_grid[p] < 10.0) )
//			printf("n_denszeroorb[%d/%d] = %f, x_grid = %f, phi = %f\n", p, p_size, n_grid[p], x_grid[p], phi[p]);
//	}
//
//	//only works if size_mu = 1 (or for the final value of mu)
//	if (show_err == 1)
//	{
//		fptr1 = fopen("OUTPUT/total_error.txt", "w");
//		for (p = 0; p < p_size; p++)
//		{
//			fprintf(fptr1, "%.15f\n", n_res[p]);
//		}
//		fclose(fptr1);
//		fptr1 = fopen("OUTPUT/regular_point_error.txt", "w");
//		for (p = 0; p < p_size; p++)
//		{
//			fprintf(fptr1, "%.15f\n", af_err[p]);
//		}
//		fclose(fptr1);
//		fptr1 = fopen("OUTPUT/initial_point_error.txt", "w");
//		for (p = 0; p < p_size; p++)
//		{
//			fprintf(fptr1, "%.15f\n", in_err[p]);
//		}
//		fclose(fptr1);
//		//fptr1 = fopen("phi.txt", "w");
//		//for (p = 0; p < p_size; p++)
//		//{
//		//	fprintf(fptr1, "%.15f\n", phi[p]);
//		//}
//		//fclose(fptr1);
//		fptr1 = fopen("OUTPUT/Edens.txt", "w");
//		for (p = 0; p < p_size; p++)
//		{
//			fprintf(fptr1, "%.15f\n", in_err[p]);
//		}
//		fclose(fptr1);
//	}
//
//	for (mu_ind = 0; mu_ind < size_mu; mu_ind++) {
//		free(ddistdvpar[mu_ind]);
//		free(ddistdvpartwo[mu_ind]);
//	}
//	free(ddistdvpar);
//	free(ddistdvpartwo);
//	free(vparacc);
//	free(F);
//	free(Fp);
//	free(Fpp);
//	free(in_err);
//	free(af_err);
//	free(n_res);
//	free(n_pre);
//	free(phi);
//
//	*Phi_point = Phi;
//
//	printf("end of denszeroorb\n");
//	return;
//}

void iondens2d(double Te, double alpha, int size_ygrid, int size_phigrid, int *size_ngrid, double **n_grid, double *y_grid, double *x_grid, double **phi_grid, double ***FF, double *mumu, double *UU, int sizemumu, int sizeUU, double grid_parameter, double *flux, int zoomfactor, double stopdens, double **vy_op, double **chiMax_op, double **dmudvy_op) {
	// declare variables
	clock_t begin = clock(); // Finds start time of computation
	double limit_rho = 8.0, *n_inf, yinf = 0.0;
	double Ucrit = 0.0, frac_reflected=0.0, frac = 0.0, vzcrit;
	double deltax, deltay, deltax_inf, deltaE = 0.1;
	double **chiinf, **muinf, ***vxinf;
	double *xx, **phi, **phi_x, **phi_y, ***chi;
	int s, w, ind, stop = 0, sizexbar, *maxj, j_inf, reflected = 0, *icrit;
	int p=0, i=0, ic=0, j=0, k=0, l=0;
	double **openorbit, openorbitantycal, ***mu, **muopen, **chiMopen, **xbaropen, **openorbitopen, **Ucritf, **xbar, *xbarcrit, *chiMcrit, **xifunction;
	int **jmclosed, **jmopen,  sizeU, size_finegrid;
	int **crossed_max, **crossed_min;
	int **lowerlimit, **upperlimit, ***upper, **imax, **imin;
	double ***Uperp, ****vx, **chiMax, **chimpp, **chimin;
	double oorbintgrd, oorbintgrdantycal, oorbintgrdBohm;
	double vz, U, dvz = 0.2, dvzopen = 0.2, dvx, dxbar, intdU=0.0, intdUopen=0.0, intdUopenBohm = 0.0;
	double intdUold=0.0, intdvx=0.0, intdvxold = 0.0, intdxbar=0.0, intdxbaropen=0.0, intdxbaropenBohm = 0.0, F, Fold=0.0, Fold_ref=0.0, Ucap, Bohm;
	double intdUopenflow = 0.0, intdUopenflowold = 0.0, intdxbaropenflow = 0.0, oorbintgrdflow = 0.0, oorbintgrdflowold = 0.0;
	double oorbintgrdold=0.0, oorbintgrdBohmold=0.0, Fopen=0.0, intdUopenold=0.0, intdUopenBohmold=0.0;
	double intdUantycal=0.0, intdvxantycal=0.0, vxnew=0.0, vxold = 0.0, Uperpnew = 0.0, **xtop, intdUopenantycal=0.0;
	double openorbitnew, chinew, munew = 0.0, vx0open;
	double oorbintgrdxbar, oorbintgrdsquare, intdUopenxbar = 0.0, intdUopensquare = 0.0, intdxbarxbar=0.0, intdxbarsquare=0.0;
	double oorbintgrdxbarold, oorbintgrdsquareold, intdUopenxbarold, intdUopensquareold;
	double xi, *gg, *ff;
	double du, fluxinf1old, fluxinfintgrdold, fluxinfintgrd, u, Chodura2, Chodura2old, Chodura1old, Chodura1, *Chodura;
	double *fluxinf, fluxinf1, densinf1, *densinf, densinf1old; 

/* 
 * Note: first index of FF and FFprime is y, second one is mu, third one is U-mu
 * xx is the array containing the position (distance from the wall) on the fine grid; phi contains the values of phi(x), extracted from a file. phi_x is phi prime, first derivative of phi. n; newphi is the new electrostatic potential guess; chi is the effective potential
 * n is the domain of the position x (the largest value of the index i, so e.g. x_grid[n] = L_2 in the paper); sizexbar is the domain of the position xbar (the largest value of the index j);
 * Uperp stores the possible values of Uperp associated with closed orbits, and so does vx; chiMax and chimin store the local maxima and minima of the effective potential maximum
 * openorbit is the Delta_M = 2*pi*dmu/dxbar; openorbitantycal is the analytical value of openorbit for a flat potential; mu is the array containing values of mu(xbar, Uperp), index j for values of xbar, k for values of Uperp; xbar is the grid of values used in the closed orbit integral; FF contains the distribution function, read from the file distfile.txt. UU and mumu contain the values of U and mu corresponding to the function FF (which is F(mu, U)); FFprime is the numerical first derivative of F with respect to U; xifunction is the function of x which defines the grid of values of xbar by finding a chi whose minimum lies exactly at each grid point x. 
 * crossed_max and crossed_min is non-zero when a minimum (or maximum) of chi is found for some xbar[j]
 * jmclosed represent minimum values of xbar above which we integrate open and closed orbit density integrals respectively (xbar_m,o and xbar_m in the paper); i is an index usually representing the positin x  j is an index usually representing the orbit position xbar; k is an index usually representing the energy Uperp (or velocity vx). It's always used in conjunction with j (and sometimes i); l is an index (used in for loops) usually representing the total energy (or velocity vz). It's used only in the DENSITY INTEGRALS part of the code; sizeU is the size of the integration range over U (or velocity vz). It is set later on in the code
 * lowerlimit represents the lower limit of k in the integrals over Uperp (or vx). It's needed because some of the earlies energies; (which are the largest because thy are values of chi stored after the maximum is found); may be so large that they are associated with very small values of the distribution function. This avoids integrating in an empty portion of phase space; upperlimit[j] represents the largest value of k (the smallest stored energy Uperp = chi_minimum) associated with some value of j; upper[j][i] represents the value of k associated with the smallest value of vx when integrating over Uperp. Going above upperlimit[j][i] makes Uperp < chi so velocities imaginary; imax/imin[j] stores the position of the maximum/minimum of the effective potential chi (It's x_M/x_m in the paper, which depends on xbar).
 * vz used in the density integral; U is the total energy, used in the density integral; dvz is the thickness of the vz grid used to take the integral over U (which is taken over vz in practice), dvzopen is the same for the open orbit piece; dvx is the thickness of the vx grid used to take the integral over Uperp ( which is taken over vx in practice). It must be evaluated because it depends on stored values of vx[j][i][k]; dxbar is the thickness of the xbar grid; intdU is the value of the integral over U in the closed orbit density integration process; intdUopen is the same as above, for the open orbit integral; intdUopenBohm same, for Bohm integral
 * oorbintgrd is the value of the integrand in the first open orbit integral (oorbintgrdantycal is the analytical result for flat potential); oorbintgrdBohm is the value of the integrand in the `Bohm' integral 
 * intdUold is a variable which stores the old intdU, so that the trapezium rule of integration can be applied (intdUold + intdU)*dvz; intdvx stores the integral over Uperp (hence over vx) in the closed orbit integral; intdxbar stores the value of the closed orbitintegral over xbar (which is the final result!), intdxbaropen does the same in the open orbit density integral; intdxbaropenBohm does the same for the Bohm integral; idealBohm is what the Bohm integral shoult be if Bohm condition is marginally satisfied; F is the value of the distribution function evaluated in the density integrals by interpolating FF, and Fold is the `old' needed to apply the trapezium rule; Fprime is the trilinearly interpolated value of FFprime, and Fprimeold is the same at the previous grid point (needed for trapezium rule); used in INTEGRALS OF DISTRIBUTION FUNCTION AT INFINITY; Ucap is the topmost total energy integrated to 
 * intdUantycal is the integral over U (or v_z) for a flat potential profile (phi =0) for some value of xbar and Uperp; intdvxantycal  is the integral over Uperp (or vx) for a flat potential profile for some value of xbar; vxnew is the value of vx at the 'new' grid point, used in the vx integral (taken using the trapezium rule); vxold is the value of vx at the 'old' grid point, used in the vx integral; Uperpnew is the value of Uperp (used in the closed orbit density integral); munew is the valye of mu (used in the closed orbit density integral); xtop is the top bounce point x_t of the last closed orbit; intdUopenantycal is the analytical value of the integral over U  in the open orbit density integral 
 * oorbintgrdxbar is an integral over the open orbit distribution function at x=0 which is needed to evaluate a coefficient that appears when; expanding quasineutrality near x=0. It was just for playing around and at the moment plays no role in the code; similarly with all other integrals here 
*/

	//make 2d phi grid which is same resolution or finer in x than input phi_grid

	// Make grids in g = sqrt(x) and the equidistant variable f
	gg = malloc(size_phigrid*sizeof(double)); 
	ff = malloc(size_phigrid*sizeof(double));
	for (i=0; i<size_phigrid; i++) {
		gg[i] = sqrt(x_grid[i]);
		ff[i] = pow(sqrt(grid_parameter)+gg[i], 2.0) - grid_parameter;
	}
	deltax = ff[1];
	deltay = y_grid[1] - y_grid[0];
	size_finegrid = zoomfactor*size_phigrid-zoomfactor;
	phi = malloc(size_ygrid*sizeof(double));
	for (p=0; p<size_ygrid; p++) 
		phi[p] = malloc(size_finegrid*sizeof(double)); // phi now has correct size
	printf("size_finegrid = %d\t it is %d x %d - %d\n", size_finegrid, zoomfactor, size_phigrid, zoomfactor);
	xx = malloc(size_finegrid*sizeof(double)); // Make finer x grid than input x_grid
	if (zoomfactor != 1) {
		for (p=0; p<size_ygrid;p++) {
			gsl_interp_accel *acc = gsl_interp_accel_alloc ();
			gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, size_phigrid);
			gsl_spline_init (spline, ff, phi_grid[p], size_phigrid);
			//gsl_spline_init (spline, gg, phi_grid, size_phigrid);
			for (i = 0; i < size_finegrid; i += 1) {
				xi = i*deltax/zoomfactor;
				if (i == 0) xi += 0.00001; //TINY;
				if (i == size_finegrid-1) xi -= 0.00001; //TINY;
				xx[i] = pow( pow(grid_parameter+xi, 0.5) - sqrt(grid_parameter), 2.0);
				//xx[i] = xi*xi;
				phi[p][i] = gsl_spline_eval (spline, xi, acc);
				phi[p][i] *= Te;
			}
			gsl_spline_free (spline);
			gsl_interp_accel_free (acc);
		}

	}
	else {
		for (p=0; p<size_ygrid;p++) {
			for (i = 0; i < size_finegrid; i += 1) {
				xi = i*deltax/zoomfactor;
				xx[i] = x_grid[i];
				phi[p][i] = phi_grid[p][i];
				phi[p][i] *= (Te);
			}
		}
	}
	printf("Te= %f\t\n", Te);
	for (p=0; p<size_ygrid;p++)
		printf("phi[%d][0] = %f\tphi_grid[%d][0] = %f\n", p, phi[p][0], p, phi_grid[p][0]);

	// Introduce a cap in energy (U, Uperp) high enough that we can safely assume F = 0
	Ucap = 18.0 + 2.0*Te;
	if (trilin_interp(yinf, 0.0, Ucap, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1) > 1e-9) printf("ERROR in iondens2d_renorm.c: increase Ucap please\n");
	else if (trilin_interp(yinf, Ucap, 0.0, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1) > 1e-9) printf("ERROR in iondens2d_renorm.c: increase Ucap please\n");
	else printf("Ucap has been checked to be large enough\n");

	/* We initialize all arrays that contain functions of position x with the correct size n */
	phi_x      = malloc(size_ygrid*sizeof(double)); 
	phi_y      = malloc(size_ygrid*sizeof(double)); 
	jmclosed   = malloc(size_ygrid*sizeof(int));
	jmopen     = malloc(size_ygrid*sizeof(int));
	xifunction = malloc(size_ygrid*sizeof(double));
	icrit      = malloc(size_ygrid*sizeof(int));
	/* FORM XBAR GRIDS 
	Take derivatives of phi and use them to obtain two grids for xbar, one to be used for closed orbits and one to be used for open orbits. */
	xbar = malloc(size_ygrid*sizeof(double));

	densinf = malloc(size_ygrid*sizeof(double));
	fluxinf = malloc(size_ygrid*sizeof(double));
	Chodura = malloc(size_ygrid*sizeof(double));
	// Evaluate integrals of distribution function at infinity
	fluxinfintgrd = fluxinf1 = densinf1 = Chodura2 = Chodura1 = 0.0;
	du = 0.01;
	printf("sizemumu = %d\tsizeUU = %d\n", sizemumu, sizeUU);
	for (p=0; p<size_ygrid; p++) {
		densinf[p] = Chodura[p] = fluxinf[p] = 0.0;
		for (i=0; i<sizemumu; i++) {
			fluxinf1old = fluxinf1;
			Chodura1old = Chodura1;
			Chodura1 = 0.0;
			fluxinf1 = 0.0;
			densinf1old = densinf1;
			densinf1 = 0.0;
			//Fprimeold = Fprime = 0.0;
			Fold = F = 0.0;
			fluxinfintgrdold = fluxinfintgrd = 0.0;
			for (j=0; j< sizeUU; j++) {
				u = sqrt(2.0*UU[j]);
				fluxinfintgrdold = fluxinfintgrd;
				Chodura2old = Chodura2;
				Fold = F;
				F = FF[p][i][j];
				fluxinfintgrd = F*u;
				if (j>0) Chodura2 = F/(u*u);
				else Chodura2 = 0.0;
				if (j != 0) {
					du = u - sqrt(2.0*UU[j-1]); 
					fluxinf1 += 0.5*(fluxinfintgrd + fluxinfintgrdold)*du;
					Chodura1 += 0.5*(Chodura2 + Chodura2old)*du;
					densinf1 += 0.5*(F + Fold)*du;
				}
			}
			if (i!=0) {
				fluxinf[p] += 0.5*(mumu[i]-mumu[i-1])*(fluxinf1 + fluxinf1old);
				Chodura[p] += 0.5*(mumu[i]-mumu[i-1])*(Chodura1 + Chodura1old);
				densinf[p] += 0.5*(mumu[i]-mumu[i-1])*(densinf1 + densinf1old);
			}
		}
		densinf[p] *= (4.0*M_PI);
		fluxinf[p] *= (4.0*M_PI);
		Chodura[p] *= (4.0*M_PI);
		flux[p] = fluxinf[p]/densinf[p];
		Chodura[p] /= densinf[p];
		printf("at y index = %d, densinf = %f\tfluxinf = %f\tChodura = %f\n", p, densinf[p], flux[p], Chodura[p]);

		xbar[p] = malloc(size_finegrid*sizeof(double));	
		phi_x[p] = malloc(size_finegrid*sizeof(double));
		phi_y[p] = malloc(size_finegrid*sizeof(double));
		for (i=0; i<size_finegrid; i++) {
			// Evaluate derivative of phi
			jmopen[p][i] = jmclosed[p][i] = 0; //*
			if (i == 0) 	
				phi_x[p][0] = (phi[p][1] - phi[p][0])/(xx[1]-xx[0]);
			else if (i == size_finegrid-1) 
				phi_x[p][i] = phi_x[p][i-1];
			else 
				phi_x[p][i] = ((xx[i] - xx[i-1])/(xx[i+1]- xx[i-1]))*(phi[p][i+1] - phi[p][i])/(xx[i+1] - xx[i]) + ((xx[i+1] - xx[i])/(xx[i+1]- xx[i-1]))*(phi[p][i] - phi[p][i-1])/(xx[i] - xx[i-1]); 

			if (p>1) 
				phi_y[p-1][i] = (phi[p][i] - phi[p-2][i])/(2*deltay);
			// xifunction is xbar corresponding to given position x being a stationary point
			xifunction[p][i] = xx[i] + phi_x[p][i];
			if (DEBUG == 1)
				printf("xifunction[%d] = %f\n", i, xifunction[p][i]); 
			if (i == 1) {	
				if  (xifunction[i] > xifunction[i-1]) { // immediately found that xi is increasing at x=0 telling us x_c = 0
					icrit[p] = i-1;
				}
			}
			else if (i > 1) {	
				if ( (xifunction[p][i] > xifunction[p][i-1]) && (xifunction[p][i-1] < xifunction[p][i-2] ) ) 
					icrit[p] = i-1;
				else if ( (xifunction[p][i] < xifunction[p][i-1]) && (xifunction[p][i-1]  > xifunction[p][i-2]) ) {
				// found maximum of xi: this only happens when phi has noise in second derivative
					printf("ERROR in iondens2d: too much noise in second derivative\n");
					printf("xifunction = %f\n", xifunction[p][i-1]);
					printf("i = %d\n", i);
				}
			} 
		}
		//printf("icrit = %d\n", icrit[p]);
	}
	//impose periodic boundary conditions in y
	for (i=0; i<size_finegrid; i++) {	
		phi_y[0][i] = (phi[1][i] - phi[size_ygrid-1][i])/(2*deltay);
		phi_y[size_ygrid-1][i] = (phi[0][i] - phi[size_ygrid-2][i])/(2*deltay);
	}

	chi         = (double***)  malloc(size_ygrid*sizeof(double*)); // chi(xbar, x) indices j and i 
	Uperp       = (double***)  malloc(size_ygrid*sizeof(double*)); 
	mu          = (double***)  malloc(size_ygrid*sizeof(double*)); 
	vx          = (double****) malloc(size_ygrid*sizeof(double**)); 
	upper       = (int***)     malloc(size_ygrid*sizeof(int*)); 
	chimin      = (double**)   malloc(size_ygrid*sizeof(double)); 
	chiMax      = (double**)   malloc(size_ygrid*sizeof(double));
	chimpp      = (double**)   malloc(size_ygrid*sizeof(double));
	crossed_min = (int**)      malloc(size_ygrid*sizeof(int)); 
	crossed_max = (int**)      malloc(size_ygrid*sizeof(int));
	openorbit   = (double**)   malloc(size_ygrid*sizeof(double));
	upperlimit  = (int**)      malloc(size_ygrid*sizeof(int));
	lowerlimit  = (int**)      malloc(size_ygrid*sizeof(int));
	xtop        = (double**)   malloc(size_ygrid*sizeof(double));
	imax        = (int**)      malloc(size_ygrid*sizeof(int));
	imin        = (int**)      malloc(size_ygrid*sizeof(int));
	n_inf       = (double*)    malloc(size_ygrid*sizeof(double));
	chiMcrit    = (double*)    malloc(size_ygrid*sizeof(double));
	xbarcrit    = (double*)    malloc(size_ygrid*sizeof(double));
	muopen      = (double**)   malloc(size_ygrid*sizeof(double));
	chiMopen      = (double**)   malloc(size_ygrid*sizeof(double));
	xbaropen      = (double**)   malloc(size_ygrid*sizeof(double));
	Ucritf      = (double**)   malloc(size_ygrid*sizeof(double));
	openorbitopen      = (double**)   malloc(size_ygrid*sizeof(double));
	vxinf      = (double***)   malloc(size_ygrid*sizeof(double));
	chiinf      = (double**)   malloc(size_ygrid*sizeof(double));
	muinf      = (double**)   malloc(size_ygrid*sizeof(double));

	FILE *fout; 
	if ((fout = fopen("OUTPUT/iondens2d_out.txt", "w")) == NULL) {	
		printf("Cannot open iondens2d_out.txt");
		exit(EXIT_FAILURE);
	}

	for (p=0; p<size_ygrid; p++) {
		j=0; // set counting index to zero
		for (k=icrit[p]+1; k<size_finegrid; k++) {
			xbar[p][j] = xifunction[p][k]; 
			j++;
		}
		xbarcrit[p] = xifunction[p][icrit[p]];
		chiMcrit[p] = 0.5*phi_x[p][icrit[p]]*phi_x[p][icrit[p]] + phi[p][icrit[p]];
		sizexbar = j;
		if (DEBUG == 1) for (j=0;j<sizexbar; j++) printf("xbar[0][%d] = %f\n", j, xbar[0][j]); 
		//printf("sizexbar = %d, size_finegrid (x) =%d\n", sizexbar, size_finegrid);
		//
		// Lots of array allocations now that size of xbar (vy + x)
		chi[p] = (double **) calloc(sizexbar,sizeof(double*)); // chi(xbar, x) indices j and i 
		Uperp[p] = (double**)  calloc(sizexbar,sizeof(double*)); 
		mu[p] = (double**)  calloc(sizexbar,sizeof(double*)); 
		vx[p] = (double***) calloc(sizexbar,sizeof(double**)); 
		upper[p] = (int**)     calloc(sizexbar,sizeof(int*)); 
		chimin[p] = (double*)   calloc(sizexbar,sizeof(double)); 
		chiMax[p] = (double*)   calloc(sizexbar,sizeof(double));
		chimpp[p] = (double*)   calloc(sizexbar,sizeof(double));
		crossed_min[p] = (int*)      calloc(sizexbar,sizeof(int)); 
		crossed_max[p] = (int*)      calloc(sizexbar,sizeof(int));
		openorbit[p] = (double*)   calloc(sizexbar,sizeof(double));
		upperlimit[p] = (int*)      calloc(sizexbar,sizeof(int));
		lowerlimit[p] = (int*)      calloc(sizexbar,sizeof(int));
		xtop[p] = (double*)   calloc(sizexbar,sizeof(double));
		imax[p] = (int*)      calloc(sizexbar,sizeof(int));
		imin[p] = (int*)      calloc(sizexbar,sizeof(int));

		/////////////////////////////////////////////////////
		/* CLOSED ORBIT ARRAY FILLING
		Set up the grid in xbar and also initialize all arrays that contain a different number at different values of xbar, indexed j */
		/* Loop below initializes all 2d arrays which are functions of xbar and x. It allocates the right amount of memory to arrays of pointers of size xbar. The result is a 2D array indexed j (size sizexbar) and i or k (size n, see above) */
		for (j=0;j<sizexbar;j++)  {
			imax[p][j] = imin[p][j] = -1;
			openorbit[p][j] = 0.0;
			crossed_min[p][j] = 0;
			crossed_max[p][j] = 0;
			xtop[p][j] = 0.0;
			chimin[p][j] = 0.0;
			chiMax[p][j] = 0.0;
			upperlimit[p][j] = -1;
			lowerlimit[p][j] = 0; 
			chi[p][j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing chi(xbar, x) indexes j and i 
			chi[p][j][0] =  0.5*pow((xx[0] - xbar[p][j]), 2.0) + phi[p][0];
			Uperp[p][j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing Uperp for closed orbits indexes j and k
			mu[p][j] = (double*)calloc(size_finegrid,sizeof(double)); // array containing adiab invariant mu(xbar, Uperp), indexes j and k
			upper[p][j] = (int*)calloc(size_finegrid,sizeof(int)); // array containing the index of energy Uperp corresponding to chi(x) 
			vx[p][j] = (double**)calloc(size_finegrid,sizeof(double*)); // see below
			/* The loop below initializes the 3d array containing the velocity of a particle at a given orbit position xbar, with particle position x and with energy Uperp, indexed i, j and k. */
			for (i = 0; i < size_finegrid; i++) 
				vx[p][j][i] = (double*)calloc(size_finegrid,sizeof(double)); // 4D array containing value of vx at different y, xbar, x and Uperp
		}
		/* This for loop fills in the arrays, calculating the integrals where necessary */
		for (i=1; i<size_finegrid; i++) {
			for (j=0; j<sizexbar; j++)
				chi[p][j][i] =  0.5*pow((xx[i] - xbar[p][j]), 2.0) + phi[p][i];
	//FINDING MAXIMA/MINIMA
	/* Below, we use the array elements to define arrays for the effective potential maxima and minima that exist for every xbar (index j). We search through the chi curve from x=0 to the largest value of x but we search through every chi curve first (so first scan in xbar at fixed x, then move to the next x). We create arrays of mu, Uperp and vx.
	Because I am comparing neighbouring values of x to find a maximum, I need to consider two distinct cases.The first one is treated in the if loop below, which the program should enter only if chi is decreasing at x=0. Implying that chi(x=0) is an effective potential maximum. The second one is the else if loop after that, which finds maxima of chi that are stationary points by comparing the value of the function before and after each point. Note that because we compare the function at a point with the function at points before and after, the point we consider at every iteration step in the index i is indexed (i-1), compared with (i-2) and i.  */
			if ( (i==1) && (chi[p][j][i] < chi[p][j][i-1]) ) {
				crossed_max[p][j] += 1;		
				imax[p][j] = 0;
				chiMax[p][j] = chi[p][j][0];
			} 
			else if ( (i > 1) && ((chi[p][j][i] < chi[p][j][i-1]) && (chi[p][j][i-1] > chi[p][j][i-2])) ) {
				crossed_max[p][j] += 1;
	// crossed_min and crossed_max let the code know whether we go up or down the effective potential well. We are starting to go down!
				if (crossed_max[p][j] > 1) {	
					printf("***WARNING*** There is more than one maximum!\n");
					printf("j is %d and maxima at %d and %d for first and second\n",  j, imax[p][j], i-1);
					for (k=imax[p][j]-1;k<i+1;k++) 
					printf("chi[%d][%d] = %f\n", j, k, chi[p][j][k]);
					crossed_max[p][j] = 1; 
					printf("exit code now\n");
					exit(-1);
				}
				imax[p][j] = i-1; // Store the index of the maximum 
				chiMax[p][j] = chi[p][j][i-1]; //Store chi Maximum itself, a function of xbar (index j)
			}
	/* We store the index corresponding to the position of a minimum for a given value of xbar.*/
			//watch out HERE
			else if ( (i > 1) && ((chi[p][j][i] > chi[p][j][i-1]) && (chi[p][j][i-1] < chi[p][j][i-2])) ) {
				crossed_min[p][j] += 1;
				//printf("i-1 = %d, imax[p][%d] = %d, crossed_min[p][%d] = %d, crossed_max[%d][%d] = %d\n", i-1, j,imax[p][j], j, crossed_min[p][j], p, j, crossed_max[p][j]);
				if (crossed_min[p][j] > 1) {	
					printf("***WARNING*** There is more than one minimum!\n");
					printf("j is %d and minima at %d and %d for first and second\n",  j, imin[p][j], i-1);
					for (k=imin[p][j]-1;k<i+1;k++) 
					printf("chi[%d][%d] = %f\n", j, k, chi[p][j][k]);
					crossed_min[p][j] = 1; 
				}
				imin[p][j] = i-1; // Store the index of the maximum 
				upperlimit[p][j] = imin[p][j] - imax[p][j];
				chimpp[p][j] = ( ( chi[p][j][i] - chi[p][j][i-1] ) / (xx[i] - xx[i-1]) - (chi[p][j][i-1] - chi[p][j][i-2]) / (xx[i-1] - xx[i-2]) ) *2.0/ (xx[i] - xx[i-2] ) ;  
				//printf("chimpp[%d] = %f\n", chimpp[j]);
				chimin[j] = chi[j][i-1];  //Store chi minimum itself, a function of xbar (index j)
				crossed_max[p][j] = 1;
			}
	/* We will now start going up the effective potential well! The temporary flag below is used because we still have to store the effective potential minimum as a possible value of Uperp. If we don't have this flag we miss the region near the minimum of chi. */
	/* FILLING IN ARRAYS */
	/* Once we cross a maximum, we start filling in arrays for mu, Uperp and vx related to the orbit we are in. As we go down the maximum we store the values of Uperp = chi(x) we encounter which will form a grid of (unevenly spaced) allowed values of Uperp. We also store the value of the small deltamu associated with every value of Uperp above the current value of chi(x), and add it to the previous values */                                                                               	
			if ( ( (crossed_max[p][j] == 1  && crossed_min[p][j] == 0) || (i-1 == imin[p][j]) )  ) {
				Uperp[p][j][i-1-imax[p][j]] = chi[p][j][i-1];
				if (Uperp[p][j][i-1-imax[p][j]] > Ucap && lowerlimit[p][j] != 0)
					lowerlimit[p][j] = i-1-imax[p][j]; 
				mu[p][j][i-1-imax[p][j]] = 0.0;
				upper[p][j][i-1] = i-1-imax[p][j];
	/* Note that the size of the dimension of the array with values of Uperp is set to n, which is larger than the size it will turn out to be. This is because in C there is no way to append elements to arrays as I go along, enough memory has to be given to the array from the start. n is the largest possible size the array could have. */
				for (k=0;k<=upper[p][j][i-1]; k++) 
				{	
					//printf("k=%d/%d\n", k, upper[j][i-1]);
					// replaced k with upper below
					if ( (upper[p][j][i-1] == 0) || (i-1 == 0) ) {
						vx[p][j][i-1][k] = sqrt((2.0*Uperp[p][j][k] - chi[p][j][i-1])); // equiv 0
						mu[p][j][k] += 0.0; 

					}
					else if ( (k == imin[p][j] - imax[p][j]) && (crossed_min[p][j] == 1) ) {
						vx[p][j][i-1][k] = 0.0;
						mu[p][j][k] = 0.0; 
					}
					else if ( (k == imin[p][j] - imax[p][j] - 1) && (crossed_min[p][j] == 1) ) {
						vx[p][j][i-1][k] = sqrt(2.0*(Uperp[p][j][k] - chi[p][j][i-1])); // equiv 0
						mu[p][j][k] = 0.5*pow(xx[imin[p][j]-1] - xx[imin[p][j]], 2.0)*pow(chimpp[p][j], 0.5);
						//printf("1..%f\n", sqrt(0.5*chimpp[j])*pow(xx[imin[j]-1] - xx[imin[j]], 2.0));
					}
					else if ( k == upper[p][j][i-1] - 1 )  {	
						vx[p][j][i-1][k] = sqrt(2.0*(Uperp[p][j][k] - chi[p][j][i-1])); // equiv 0
						mu[p][j][k] += (sqrt(2.0)/M_PI)*sqrt(chi[p][j][i-2]-chi[p][j][i-1])*(2.0/3.0)*(xx[i-1] - xx[i-2]);
					}
					else if (k == upper[p][j][i-1]) {
						vx[p][j][i-1][k] = 0.0;
					}
					else {
						vx[p][j][i-1][k] = sqrt(2.0*(Uperp[p][j][k] - chi[p][j][i-1]));
						mu[p][j][k] += (1.0/M_PI)*0.5*(vx[p][j][i-1][k] + vx[p][j][i-2][k])*(xx[i-1] - xx[i-2]); 
					}
					if (mu[p][j][k] != mu[p][j][k]) {
						printf("BEFORE: mu[%d][%d] is NAN\n", j, k); 
						exit(-1);
					}  
				}
			}
	/* Once we cross the minimum, we stop creating array elements with values of Uperp. However, we keep storing the value of vx associated with any given point x on an effective potential curve with xbar, with energy Uperp and using this value to finish performing the mu integral. This should happen as long the effective potential at the point under consideration is smaller than the effective potential maximum. */
			else if ( ( crossed_min[p][j] == 1 && crossed_max[p][j] == 1 && chi[p][j][i-1] < chiMax[p][j] && ( i-1 != imin[p][j] ) ) ) {
				for (k=0;k <= upperlimit[p][j] ;k++) {	
					if ( (chi[p][j][i-1] < Uperp[p][j][k]) && (chi[p][j][i-2] < Uperp[p][j][k]) ) {	
						vx[p][j][i-1][k] = sqrt(2.0*(Uperp[p][j][k] - chi[p][j][i-1]));
						mu[p][j][k] += (1.0/M_PI)*0.5*(vx[p][j][i-1][k] + vx[p][j][i-2][k])*(xx[i-1] - xx[i-2]); 
					}
					else if (Uperp[p][j][k] <= chi[p][j][i-1] && Uperp[p][j][k-1] > chi[p][j][i-1]) {
						upper[p][j][i-1] = k;
				//mu[j][k] += (2.0/M_PI)*(vx[j][i-2][k])*(xx[i-1] - xx[i-2])*(Uperp[j][k] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]);
						ind = 0;
						while (Uperp[p][j][k] < chi[p][j][i-2-ind]) 
							ind++;
						
						mu[p][j][k] += (sqrt(2.0)/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[p][j][k] - chi[p][j][i-2-ind], 1.5)/(chi[p][j][i-1-ind] - chi[p][j][i-2-ind]); // double check normalization
			
					}
					else if (Uperp[p][j][k] <= chi[p][j][i-1] && Uperp[p][j][k] > chi[p][j][i-2]) 
						mu[p][j][k] += (sqrt(2.0)/M_PI)*(2.0/3.0)*(xx[i-1] - xx[i-2])*pow(Uperp[p][j][k] - chi[p][j][i-2], 1.5)/(chi[p][j][i-1] - chi[p][j][i-2]); 
					if (mu[p][j][k] != mu[p][j][k]) {
						printf("mu[%d][%d] is NAN\n", j, k); 
						exit(-1);
					}  
				}
			}
	/* When the effective potential at the iteration point (i-1) under consideration becomes larger than the effective potential maximum, we finish performing the mu integral. We also store the position of the top of the orbit which has chi = chiMax, in order to perform the open orbit integral. If the loop below is accessed, a switch it turned off to signify the no more closed orbits can be present */
			//else if ( ( crossed_min[j] == 1 ) && ( crossed_max[j] == 1) && ( chi[j][i-1] > chiMax[j] - TINY) ) 	
			else if ( ( crossed_min[p][j] == 1 ) && ( crossed_max[p][j] == 1) && ( chi[p][j][i-1] > chiMax[p][j] ) ) {
				//itop[j] = i-2;
				xtop[p][j] = xx[i-2] + ((chiMax[p][j] - chi[p][j][i-2])/(chi[p][j][i-1] - chi[p][j][i-2]))*(xx[i-1] - xx[i-2]);
				for (k=0; k<upper[p][j][i-2]; k++) {
					ind = 0;
					//while (Uperp[j][0] < chi[j][i-2-ind]) 
					//	ind++; // find top bounce point index i-2-ind for a given xbar[j] and for Uperp[j][0] = chiMax
					//mu[j][k] += (sqrt(2.0)/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[j][0] - chi[j][i-2-ind], 1.5)/(chi[j][i-1-ind] - chi[j][i-2-ind]);  // CHECK NORMALIZATION
					while (Uperp[p][j][k] < chi[p][j][i-2-ind]) 
						ind++; // find top bounce point index i-2-ind, with i the smallest index such that chi[j][i-1] > chiMax[j] for a given xbar[j], for orbits with Uperp[j][k] > chi[j][i-2] 
					mu[p][j][k] += (sqrt(2.0)/M_PI)*(2.0/3.0)*(xx[i-1-ind] - xx[i-2-ind])*pow(Uperp[p][j][k] - chi[p][j][i-2-ind], 1.5)/(chi[p][j][i-1-ind] - chi[p][j][i-2-ind]);  // CHECK NORMALIZATION
					//printf("mu[%d][%d] = %f\n", j, k, mu[j][k]);
					//mu[j][k] += (1.0/M_PI)*(vx[j][i-2][k])*(xx[i-1] - xx[i-2])*(Uperp[j][k] - chi[j][i-2])/(chi[j][i-1] - chi[j][i-2]);
				}
				crossed_max[j] = 0; 
			}
			if (j!=0) {
				if ( ((chiMax[p][j-1] < chi[p][j-1][i-1] ) && (chiMax[p][j] > chi[p][j][i-1])) ) { // || (imax[p][j-1] == -1 && imax[p][j] != -1) )
					jmclosed[p][i-1] = j-1;
					if (i-1>icrit[p]) { 	
						jmopen[p][i-1] = j-1; 
						if (DEBUG == 1) {
							printf("i = %d, j = %d, icrit = %d\n", i-1, j-1, icrit[p]); 
							printf("jmopen[%d] = %d\n", i-1, jmopen[p][i-1]); 
						}
					} 
				}
			} 
			//printf("upper = %d, upperlimit = %d\n", upper[j][i-1], upperlimit[j]);
		} 

		upperlimit[p][sizexbar-1] = upperlimit[p][sizexbar-2] +1;
		printf("hello\n");
		maxj = malloc(size_ygrid*sizeof(int));
		//printf("xbar[maxj=%d] = %f\n", maxj, xbar[maxj]);	
		// OPEN ORBIT INTEGRAL
		/* Now we perform the open orbit integral. We use a change of variables which makes the integrand smooth at the top bounce point. The change of variables is to some var = sqrt(x_t - x) */
		maxj[p] = sizexbar +1 ; //*
		muopen[p] = malloc(maxj[p]*sizeof(double)); //*
		chiMopen[p] = malloc(maxj[p]*sizeof(double));//*	
		xbaropen[p] = malloc(maxj[p]*sizeof(double));//*
		openorbitopen[p] = malloc(maxj[p]*sizeof(double));
		Ucritf[p] = malloc(maxj[p]*sizeof(double));	

		xbaropen[p][0] = xbarcrit[p]; //*
		chiMopen[p][0] = chiMcrit[p]; //*
		muopen[p][0] = 0.0; //*
		openorbitopen[p][0] = 0.0; //*
		Ucritf[p][0] = chiMcrit[p]; 
		printf("sizexbar = %d\n", sizexbar);
		for (j=0;j<maxj[p]-1;j++) {
			//printf("lowerlimit[%d] = %d and upper[%d][0] = %d, upperlimit = %d\n", j, lowerlimit[j],j, upper[j][0], upperlimit[j]);
			mu[p][j][upperlimit[p][j]] = 0.0;
			chiMopen[p][j+1] = chiMax[p][j];
			muopen[p][j+1] = mu[p][j][0]; //*
			xbaropen[p][j+1] = xbar[p][j];
			Ucritf[p][j+1] = chiMax[p][j] - mu[p][j][0] ; //*
			//printf("mucrit , Ucrit = %f, %f\n", muopen[j+1], Ucritf[j+1]);
			if (DEBUG == 1)
				printf("jmopen[%d] = %d\n", j, jmopen[p][j]); 
			//if (j==0) twopimuprime[0] = 0.0; //*
			if (j==0) //*
				openorbit[p][j] = (2.0*M_PI) * ( ((xbar[p][j] - xbarcrit[p])/(xbar[p][j+1] - xbarcrit[p])) *(mu[p][j+1][0] - mu[p][j][0])/(xbar[p][j+1] - xbar[p][j]) + ((xbar[p][j+1] - xbar[p][j])/(xbar[p][j+1] - xbarcrit[p])) * (mu[p][j][0] - 0.0)/(xbar[p][j] - xbarcrit[p]) );	
			else if (j==maxj[p]-2) 
				openorbit[p][j] = openorbit[p][j-1];	
			else //* 
				openorbit[p][j] = (2.0*M_PI) * (((xbar[p][j] - xbar[p][j-1])/(xbar[p][j+1] - xbar[p][j-1])) *(mu[p][j+1][0] - mu[p][j][0])/(xbar[p][j+1] - xbar[p][j]) + ((xbar[p][j+1] - xbar[p][j])/(xbar[p][j+1] - xbar[p][j-1])) * (mu[p][j][0] - mu[p][j-1][0])/(xbar[p][j] - xbar[p][j-1]) );	
			openorbitopen[p][j+1] = openorbit[p][j];
			openorbitantycal = 4.0*M_PI*xbar[p][j];
			if (DEBUG == 1)
				printf("%f %f %f %f %f\n", xbar[p][j], mu[p][j][0], Uperp[p][j][0], openorbit[p][j], openorbitantycal); 
		}
		muopen[p][maxj[p]-1] = 5000.0;
		chiMopen[p][maxj[p]-1] = 5000.0; 
		xbaropen[p][maxj[p]-1] = 100.0;
		openorbitopen[p][maxj[p]-1] = 100.0;
		Ucritf[p][maxj[p]-1] = 0.0;
		if (DEBUG == 1) {
			printf("~~~~~The second element of FF is %f~~~~~\n",FF[0][0][1]);
			printf("~~~~~The second element of UU is %f~~~~~\n",UU[1]);
			printf("~~~~~The fourth element of mu is %f~~~~~\n",mumu[3]);
		}
		i=0;
		clock_t int1 = clock(); // Finds the time of the computation so far
		double inttime  = (double)(int1 - begin) / CLOCKS_PER_SEC;
		if (DEBUG == 1) 
			printf("in iondens2d: Array filling DONE: time is %f\n", inttime);

		/* DENSITY INTEGRALS 
		This part calculates the density integrals and outputs the result of the integration to a file fout and also the yz distribution function to three files one containing the distribution function the other two containing the velocity grid */
		////////////////////////////////////////////
		// calculate normalization at infinity
		deltax_inf = (xx[size_finegrid-1] - xx[size_finegrid-2]);
		j_inf = (int) sqrt(Ucap)/deltax_inf ;
		printf("j_inf = %d, deltax_inf = %f\n", j_inf, deltax_inf);
		muinf[p] = (double*)calloc(j_inf,sizeof(double*)); 
		vxinf[p] = (double**)calloc(j_inf,sizeof(double**)); 
		chiinf[p] = (double*)calloc(j_inf,sizeof(double*));
		for (j=j_inf-1; j>=0; j--) {	
			vxinf[p][j] = (double*)calloc(j_inf,sizeof(double*)); 
			chiinf[p][j] =  0.5*deltax_inf*j*deltax_inf*j;
			for (k=j; k<j_inf; k++) {
				vxinf[p][j][k] = sqrt(2.0*(chiinf[k] - chiinf[j]));
			}
			muinf[p][j] = chiinf[p][j];
		}

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
				munew = muinf[p][k]; 
				vxnew = vxinf[p][j][k];
				Ucrit = lin_interp(muopen[p], Ucritf[p], munew, maxj[p], 790);
				Uperpnew = chiinf[p][k];
				sizeU = (int) sqrt(Ucap - Uperpnew)/dvz;
				reflected = 1;
					for (l=0; l < sizeU; l++) {	
						if (l!=0) {	
							Fold = F;
							Fold_ref = F;
							vz = dvz*l;
							U = Uperpnew + 0.5*pow(vz, 2.0);
							if ( (U > munew) && (U - 0.5*vz*vz + 0.5*(vz-dvz)*(vz-dvz) < munew) ) {
								frac = (vz - sqrt(2.0*(munew - Uperpnew)))/dvz;
								Fold = trilin_interp(yinf, munew, 0.0, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
								F = trilin_interp(yinf, munew, U-munew, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							}
							else if (U > munew){
								frac = 1.0;
								F = trilin_interp(yinf, munew, U-munew, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
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
									vzcrit = sqrt(2.0*(Ucrit + munew - Uperpnew));
									Fold_ref = trilin_interp(yinf, munew, Ucrit, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
									frac_reflected = (vzcrit - (vz - dvz))/dvz;
								}
								//printf("frac_reflected = %f, l=%d, vzcrit = %f, vz - dvz = %f\n", frac_reflected, l, vzcrit, vz - dvz);
								//printf("U = %f, Ucrit = %f, munew = %f, Uperpnew = %f\n", U, Ucrit, munew, Uperpnew);
							}
							else frac_reflected = 0.0;
									
							intdU += frac*dvz*(F+Fold) + frac*frac_reflected*dvz*(F + Fold_ref);

							//if (reflected == 1) printf("particles are being reflected\n");
							//printf("frac = %f, frac_reflected = %f, F = %f, Fold = %f, Fold_ref = %f\n", frac, frac_reflected, F, Fold, Fold_ref);
						}
						else {	
							U = Uperpnew ;
							F = trilin_interp(yinf, munew, U-munew, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							intdU += 0.0;
						} 
					}

				intdUantycal = exp(-chiinf[p][k])*(1.0/(2.0*M_PI));// result with phi =0
				if (DEBUG == 1) printf("Analytical intdU is %f, numerical one is %f\n", intdUantycal, intdU);
				if (k!=j) {	
					dvx = vxnew - vxold; 
					intdvx += 4.0*0.5*dvx*(intdU+intdUold); // CHECK NORMALIZATION
				}
				intdvxantycal = (2.0/(2.0*sqrt(M_PI)))*exp(-deltax_inf*j*deltax_inf*j);//*erf(sqrt(xbar[j]*xbar[j]-(pos-xbar[j])*(pos-xbar[j])));
			}
			if (DEBUG == 1)	printf("intdvx is %f, analytical one is %f\n", intdvx, intdvxantycal); 
			//printf("intdvx is %f, while the analytical one is %f\n", intdvx, intdvxantycal);
			//intdvx = intdvxantycal;
			
			dxbar = deltax_inf;
			if (j!=0) intdxbar += (intdvx+intdvxold)*dxbar; // multiply by two because only half the xbars are contemplated here
		}
		n_inf[p] = intdxbar;
		if (DEBUG == 1) {	
			printf("n_inf = %f\n", n_inf[p]);
			if ( (n_inf[p]  != n_inf[p]) || (n_inf[p] < TINY) ) { 	
				printf("n_inf = %f\n", n_inf[p]);
				exit(-1);
			}
		}

		// density profile
		i=0; // start from x[0] = 0 i.e. the wall
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
				intdUopenflowold = intdUopenflow; intdUopenBohmold = intdUopenBohm;
				intdUopenxbarold = intdUopenxbar; intdUopensquareold = intdUopensquare;
				intdUopen = 0.0;
				intdUopenflow = 0.0; intdUopenBohm = 0.0;
				intdUopenxbar = 0.0; intdUopensquare = 0.0;
				if (j == jmopen[p][i]) {
					oorbintgrd = oorbintgrdold = 0.0;
					oorbintgrdflow = oorbintgrdflowold = 0.0; oorbintgrdBohm = oorbintgrdBohmold = 0.0;
					oorbintgrdxbar = oorbintgrdxbarold = 0.0; oorbintgrdsquare = oorbintgrdsquareold = 0.0;
					sizeU = (int) sqrt(2.0*(Ucap - chiMax[p][j]))/dvzopen;
					if (j==0) { //????????????
						munew = 0.0;
						openorbitnew = 0.0;
					}
					else {
						munew = mu[p][j+1][0] + ( (chiMax[p][j+1] - chi[p][j+1][i]) / (chiMax[p][j+1] - chi[p][j+1][i] + chi[p][j][i] - chiMax[p][j]) ) * (mu[p][j][0] - mu[p][j+1][0]);
						openorbitnew = openorbit[p][j+1] + ( (chiMax[p][j+1] - chi[p][j+1][i]) / (chiMax[p][j+1] - chi[p][j+1][i] + chi[p][j][i] - chiMax[p][j]) ) * (openorbit[p][j] - openorbit[p][j+1]);
					}
					for (l=0; l < sizeU; l++) {
						oorbintgrdold = oorbintgrd;
						oorbintgrdflowold = oorbintgrdflow; oorbintgrdBohmold = oorbintgrdBohm;
						oorbintgrdxbarold = oorbintgrdxbar; oorbintgrdsquareold = oorbintgrdsquare;
						vz = dvzopen*l;
						if (j!=0) { ///???????
							chinew = chi[p][j+1][i] + ( (chiMax[p][j+1] - chi[p][j+1][i]) / (chiMax[p][j+1] - chi[p][j+1][i] + chi[p][j][i] - chiMax[p][j]) ) * (chi[p][j][i] - chi[p][j+1][i]) ;
							vx0open = 0.0;
						}
						else {
							chinew = chi[p][j][i];
							if (chiMax[p][j] > chinew) 
								vx0open = sqrt(2.0*(chiMax[p][j] - chinew));
							else vx0open = 0.0;
						}
						U = chinew + 0.5*pow(vz, 2.0); 
						if ( (U > munew) && (U - 0.5*vz*vz + 0.5*(vz-dvzopen)*(vz-dvzopen) < munew) ) {
							frac = (vz - sqrt(2.0*(munew - chinew)))/dvzopen;
							Fopen = trilin_interp(yinf, munew, 0.0, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							oorbintgrdold = ( sqrt(vx0open*vx0open + alpha*sqrt(2.0*(munew - chinew))*openorbitnew) - vx0open )*Fopen;
							oorbintgrdflowold = 0.5*alpha*sqrt(2.0*(munew - chinew))*openorbitnew*Fopen;
							Fopen = trilin_interp(yinf, munew, U-munew, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							oorbintgrd = ( sqrt(vx0open*vx0open + alpha*vz*openorbitnew) - vx0open )*Fopen;
							oorbintgrdflow = 0.5*alpha*vz*openorbitnew*Fopen;
						}
						else {
							frac = 1.0;
							Fopen = trilin_interp(yinf, munew, U-munew, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							oorbintgrd = ( sqrt(vx0open*vx0open + alpha*vz*openorbitnew) - vx0open )*Fopen;
							oorbintgrdflow = 0.5*alpha*vz*openorbitnew*Fopen;
						}
						if (oorbintgrd != oorbintgrd) {	
							printf("vx0open = %f, openorbit[%d] = %f, imaginary oorbintgrd in initial piece of integral due to negative value of openorbit?\n", vx0open, j, openorbit[p][j]); 
							exit(-1);
						}
						if (i==0) 
							oorbintgrdantycal = sqrt(2.0*alpha*vz*M_PI*xbar[p][j])*(U-chiMax[p][j])*exp(-U)/pow(M_PI, 1.5);  // CHECK NORMALIZATION
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
					if (j == 0) { //*
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
				}
				else if (j > jmopen[p][i]) {
					oorbintgrd = oorbintgrdold = 0.0;
					oorbintgrdflow = oorbintgrdflowold = 0.0;
					oorbintgrdBohm = oorbintgrdBohmold = 0.0;
					oorbintgrdxbar = oorbintgrdxbarold = 0.0;
					oorbintgrdsquare = oorbintgrdsquareold = 0.0;
					sizeU = (int) sqrt(2.0*(Ucap - chiMax[p][j]))/dvzopen;
					for (l=0; l < sizeU; l++) {
						oorbintgrdold = oorbintgrd;
						oorbintgrdflowold = oorbintgrdflow;
						oorbintgrdBohmold = oorbintgrdBohm;
						oorbintgrdxbarold = oorbintgrdxbar;
						oorbintgrdsquareold = oorbintgrdsquare;
						vz = dvzopen*l;
						U = chiMax[p][j] + 0.5*pow(vz, 2.0);
						//vx0open = sqrt(TINY + chiMax[j] - chi[j][i]);
						vx0open = sqrt(2.0*(chiMax[p][j] - chi[p][j][i])); // SEARCH HERE TO FIND CURRENT POSITION OF NEW CHANGES TO NORMALIZATION
						if (vx0open != vx0open) {
							printf("HERE imaginary vx0open, j = %d, i is %d, chi[j][i] = %f, chiMax[j] = %f\n", j, i, chi[p][j][i], chiMax[p][j]); 
							exit(-1);
						}
						if ( (U >= mu[p][j][0]) && (U - 0.5*vz*vz + 0.5*(vz-dvzopen)*(vz-dvzopen) < mu[p][j][0]) ) {
							frac = (vz - sqrt(2.0*(mu[p][j][0] - chiMax[p][j])))/dvzopen;
							Fopen = trilin_interp(yinf, mu[p][j][0], 0.0, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							oorbintgrdold = sqrt(2.0*alpha*sqrt(2.0*(mu[p][j][0] - chiMax[p][j]))*openorbit[p][j])*Fopen;
							oorbintgrdflowold = alpha*sqrt(2.0*(mu[p][j][0] - chiMax[p][j]))*openorbit[p][j]*Fopen;
							oorbintgrdBohmold = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + 2.0*alpha*sqrt(2.0*(mu[p][j][0] - chiMax[p][j]))*openorbit[p][j]))*Fopen;
							oorbintgrdxbarold = xbar[p][j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + 2.0*alpha*sqrt(2.0*(mu[p][j][0] - chiMax[p][j]))*openorbit[p][j]))*Fopen;
							oorbintgrdsquareold = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + 2.0*alpha*sqrt(2.0*(mu[p][j][0] - chiMax[p][j]))*openorbit[p][j]), 1.5))*Fopen;
							Fopen = trilin_interp(yinf, mu[p][j][0], U-mu[p][j][0], FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							oorbintgrd = sqrt(2.0*alpha*vz*openorbit[p][j])*Fopen;
							oorbintgrdflow = alpha*vz*openorbit[p][j]*Fopen;
							oorbintgrdBohm = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbit[p][j]))*Fopen;
							oorbintgrdxbar = xbar[p][j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbit[p][j]))*Fopen;
							oorbintgrdsquare = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + 2.0*alpha*vz*openorbit[p][j]), 1.5))*Fopen;

						}
						else {
							frac = 1.0;
							Fopen = trilin_interp(yinf, mu[p][j][0], U-mu[p][j][0], FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
							//if (i ==0) printf("l = %d, Fopen = %f\n", l, Fopen);
							oorbintgrd = (sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbit[p][j]) - vx0open)*Fopen;
							oorbintgrdflow = alpha*vz*openorbit[p][j]*Fopen;
							oorbintgrdBohm = (1.0/vx0open - 1.0/sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbit[p][j]))*Fopen;
							oorbintgrdxbar = xbar[p][j]*(1.0/vx0open - 1.0/sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbit[p][j]))*Fopen;
							oorbintgrdsquare = (1.0/8.0)*(1.0/pow(vx0open, 3.0) - 1.0/pow((vx0open*vx0open + 2.0*alpha*vz*openorbit[p][j]), 1.5))*Fopen;
							//if (oorbintgrd != oorbintgrd) 
							if (frac != frac) 
								printf("AHA\n\n\n\n\n");
						}
						if (oorbintgrd != oorbintgrd) {
							printf("vx0open = %f, openorbit[%d] = %f, HERE imaginary oorbintgrd in initial piece of integral due to negative value of openorbit?\n", vx0open, j, openorbit[p][j]); 
							exit(-1);
						}
						if (i==0) {
							oorbintgrdantycal = sqrt(2.0*alpha*vz*M_PI*xbar[p][j])*(U-chiMax[p][j])*exp(-U)/pow(M_PI, 1.5); 
							if (DEBUG == 1) 
								printf("openorbit = %f (should be %f with flat potential profile)\n", oorbintgrd, oorbintgrdantycal);
							// check not working
						}
		// for a flat potential, oorbintgrd can be calculated analytically
						//oorbintgrd = oorbintgrdantycal;
		//printf("oorbintgrd is %f, analytical is %F\n", oorbintgrd, oorbintgrdantycal);
						if (l!=0) {
							intdUopen       += 2.0*frac*dvzopen*(oorbintgrd+oorbintgrdold);
							intdUopenflow   += 2.0*frac*dvzopen*(oorbintgrdflow+oorbintgrdflowold);
							intdUopenBohm   += 2.0*frac*dvzopen*(oorbintgrdBohm+oorbintgrdBohmold);
							intdUopenxbar   += 2.0*frac*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
							intdUopensquare += 2.0*frac*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
						}
						else {
							intdUopen       += 0.0;
							intdUopenflow   += 0.0;
							intdUopenBohm   += 0.0;
							intdUopenxbar   += 2.0*frac*dvzopen*(oorbintgrdxbar+oorbintgrdxbarold);
							intdUopensquare += 2.0*frac*dvzopen*(oorbintgrdsquare+oorbintgrdsquareold); 
						} 
					}
					if (intdUopen != intdUopen)  {
						printf("HERE, j is %d\n", j); 
						exit(-1);
					}
					if (i==0) {	
						intdUopenantycal = 2.0*0.919*sqrt(2.0*alpha*xbar[p][j])*exp(-xbar[p][j]*xbar[p][j])/M_PI; 
						if (DEBUG == 1) { 
							printf("xbar[p][%d] = %f, analytical is %f, numerical is %f\n", j, xbar[p][j], intdUopenantycal, intdUopen); printf("vx0open is %f (should be zero)\n", vx0open); 
						}
					} 
					dxbar = xbar[p][j] - xbar[p][j-1];
					if ( (j == jmopen[p][i]+1) &&  (j != 1) ) { // in this case dxbar is different //*
						dxbar = (xbar[p][j] - xbar[p][j-1])*(chiMax[p][j] - chi[p][j][i])/ (chiMax[p][j] - chi[p][j][i] + chi[p][j-1][i] - chiMax[p][j-1]);// open orbit density does not need to be so accurate at this point
						if (DEBUG == 1)
							printf("dxbar = %f\n", dxbar);
					}
					intdxbaropen += 0.5*(intdUopen+intdUopenold)*dxbar;
					intdxbaropenflow += 0.5*(intdUopenflow+intdUopenflowold)*dxbar;
					intdxbaropenBohm += 0.5*(intdUopenBohm+intdUopenBohmold)*dxbar;
					intdxbarxbar += 0.5*(intdUopenxbar+intdUopenxbarold)*dxbar;
					intdxbarsquare += 0.5*(intdUopensquare+intdUopensquareold)*dxbar; 
					if (intdxbaropen != intdxbaropen) {
						printf("PROBLEM HERE, j is %d\n", j); 
						exit(-1);
					} 
				}
				if ( j>=jmclosed[p][i] ) {	
					/* We have entered the closed orbit integral */
					intdvxold = intdvx;
					//intdvxflowold = intdvxflow;
					//intdvxflow = 0.0;
					intdvx = 0.0;
					if (j == jmclosed[p][i]) {
						intdvx = 0.0;
						//intdxbar += 0.0; 
					}
					else if (j > jmclosed[p][i]) {	
						if (DEBUG ==1) {
							printf("lowerlimit[%d] = %d\tupper[%d][%d] = %d\n", j, lowerlimit[p][j], j, i, upper[p][j][i]);
							printf("mu = %f\n", mu[p][j][lowerlimit[p][j]]);
							printf("imax[p][%d] = %d\n", j, imax[p][j]);
						}
						for (k=lowerlimit[p][j]; k<upper[p][j][i]+1; k++) {	
							vxold = vxnew;
							intdUold = intdU;
							intdU = 0.0;
							if (k == upper[p][j][i]) {
								Uperpnew = chi[p][j][i];
								if (lowerlimit[p][j] == upper[p][j][i] ) {
									//munew = 0.0;
									munew = mu[p][j][k]; 
									//printf("why would I ever be here? i = %d mu = %f\n\n\n", i, munew);
								}
								else if (k == upperlimit[p][j]) // correct but unnecessary as taken ino account below
									munew = 0.0;
								else {
									munew = ((chi[p][j][i] - Uperp[p][j][k])*mu[p][j][k-1] + (Uperp[p][j][k-1] - chi[p][j][i])*mu[p][j][k])/(Uperp[p][j][k-1] - Uperp[p][j][k]);
								}
								vxnew = 0.0;
								dvx = vxold - vxnew; 
							}
							else {	
								Uperpnew = Uperp[p][j][k];
								vxnew = vx[p][j][i][k];
								dvx = vxold - vxnew;
								munew = mu[p][j][k]; 
							}
							Ucrit = lin_interp(muopen[p], Ucritf[p], munew, maxj[p], 791);
							sizeU = (int) sqrt(2.0*(Ucap - Uperpnew))/dvz;
							reflected = 1;
							for (l=0; l < sizeU; l++) {	
								if (l!=0) {	
									Fold = F;
									Fold_ref = F;
									vz = dvz*l;
									U = Uperpnew + 0.5*vz*vz;
									if ( (U > munew) && (U - 0.5*vz*vz + 0.5*(vz-dvz)*(vz-dvz) < munew) ) {
										frac = (vz - sqrt(2.0*(munew - Uperpnew+TINY)))/dvz;
										Fold = trilin_interp(yinf, munew, 0.0, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
										F = trilin_interp(yinf, munew, U-munew, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
									}
									else if (U > munew){
										frac = 1.0;
										F = trilin_interp(yinf, munew, U-munew, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
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
											vzcrit = sqrt(2.0*(Ucrit + munew - Uperpnew));
											Fold_ref = trilin_interp(yinf, munew, Ucrit, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
											frac_reflected = (vzcrit - (vz - dvz))/dvz;
										}
									}
									else 
										frac_reflected = 0.0;
									intdU += frac*dvz*(F+Fold) + frac*frac_reflected*dvz*(F + Fold_ref);
								}
								else {	
									U = Uperpnew ;//+ 0.5*vz*vz;
									if (phi[0][0] < 0.0) reflected = 0.0;
									F = trilin_interp(yinf, munew, U-munew, FF, y_grid, mumu, UU, size_ygrid, sizemumu, sizeUU, -1, -1, -1); 
									intdU += 0.0;
								} 
							}
							intdUantycal = exp(-Uperpnew)*(1.0/(2.0*M_PI));// result with phi =0
							if (DEBUG == 1) 	
								printf("Analytical intdU is %f, numerical one is %f\n", intdUantycal, intdU);
							if (k==lowerlimit[p][j]) {	
								intdvx += 0.0; 
							}
							else {	
								intdvx += 4.0*0.5*dvx*(intdU+intdUold);
							}
							if (intdvx != intdvx) {	
								printf("intdvx is NAN, j=%d, i=%d\n", j, i); 
								exit(-1);
							} 
						}
						intdvxantycal = (2.0/(2.0*sqrt(M_PI)))*exp(-(xx[i]-xbar[p][j])*(xx[i]-xbar[p][j]))*erf(sqrt(xbar[p][j]*xbar[p][j]-(xx[i]-xbar[p][j])*(xx[i]-xbar[p][j])));
						if (DEBUG == 1) {
							printf("pos=%f, i=%d, j=%d, Uperp is %f, chi is %f\nvxnew and vxold are %f and %f and and vx is %f, dvx is %f\nintdvx is %f, analytical one is %f, upper is %d, upperlimit is %d\n", xx[i], i, j, Uperpnew, chi[p][j][i], vxnew, vxold, vx[p][j][i][k], dvx, intdvx, intdvxantycal, upper[p][j][i], upperlimit[p][j]); 
						}
						if (j==jmclosed[p][i]+1 && i!=0) {
							if (i==imax[p][j]) 
								dxbar = xbar[p][j] - xbar[p][j-1]; 
							else 
								dxbar = (xbar[p][j] - xbar[p][j-1])*(chiMax[p][j] - chi[p][j][i])/ (chiMax[p][j] - chi[p][j][i] + chi[p][j-1][i] - chiMax[p][j-1]);// open orbit density does not need to be so accurate at this point
							intdxbar += 0.5*(intdvx+intdvxold)*dxbar;
							if (intdxbar != intdxbar) {	
								printf("intdxbar is NAN, j=%d, i=%d\n", j, i);  
								exit(-1);
							} 
						}
						else
						{
							dxbar = xbar[p][j] - xbar[p][j-1];
							intdxbar += 0.5*(intdvx+intdvxold)*dxbar;
							if (intdxbar != intdxbar) {
								printf("intdxbar is NAN, j=%d, i=%d\n", j, i); 
								exit(-1);
							} 
						}
					} 
				} 
			}
			//niclosed = intdxbar; //niopen = intdxbaropen;
			n_grid[p][ic] = intdxbar + intdxbaropen;
			if (ic == 0) {
				Bohm = intdxbaropenBohm/n_grid[p][0];
				printf("flow velocity at x=0 = %f\n", (flux[p])/n_grid[p][0]);
				printf("flux evaluated at x=0 is %f\n", intdxbaropenflow/(n_inf[p]*alpha));
				printf("flux evaluated at x=0 is %f\n", *flux);
				//printf("intdxbarsquare = %f\tintdxbarxbar = %f\n", intdxbarsquare, intdxbarxbar);
				//kBohmsquared = intdxbarxbar/(intdxbarsquare - 0.5*n_grid[0]); 
			}
			if ( n_grid[p][ic] >= stopdens*n_inf[p]) {	
				printf("ic = %d/%d\n", ic, size_phigrid);
				stop = 1;
				*size_ngrid = ic; 
				printf("stopping because density larger than %f times the density at infinity\n", stopdens);
			}
			else if (x_grid[ic] > x_grid[size_phigrid-1] - limit_rho) {
				printf("ic = %d/%d\n", ic, size_phigrid);
				stop = 1;
				*size_ngrid = ic; 
				printf("WARNING in iondens2d.c: stopping for positive density\n");

			}
			if ( (DEBUG == 1) ) {	
				printf("%f, %f, %f is TOTAL, CLOSED and OPEN orbit density at position index %d, position %f, potential %f\n", n_grid[p][ic], intdxbar, intdxbaropen, ic, xx[i], phi_grid[p][ic]);
				if ( (n_grid[p][ic] != n_grid[p][ic]) || (n_grid[p][ic] < TINY) ) { 	
					printf("n_finorb[ic] = %f\n", n_grid[p][ic]);
					//exit(-1);
				} 
			} 
			fprintf(fout, "%f %f %f %f\n", xx[i], n_grid[p][ic]/n_inf[p], intdxbar/n_inf[p], intdxbaropen/n_inf[p]);
			ic += 1;
		} 
		fclose(fout);
		if (stop == 0) { 
			printf("ERROR: the density never reached stopdens*n_inf\n"); 
			exit(-1); 
		}
		for (ic = 0; ic < size_phigrid; ic++) {
			if (DEBUG == 1) 
				printf("before renormalizing n_finorb[%d] = %f, phi_grid[%d] = %f\n", ic, n_grid[p][ic], ic, phi_grid[p][ic]);
			if (ic < *size_ngrid) {
				n_grid[p][ic] /= n_inf[p];
				if (DEBUG ==1) {
					printf("x_grid[%d] = %f\tn_finorb[%d] = %f\tphi_grid[%d] = %f\n", ic, x_grid[ic], ic, n_grid[p][ic], ic, phi_grid[p][ic]);
				}
			}
			else n_grid[p][ic] = 0.0;
		}


		if (DEBUG == 1) { // change this to 0 to debug
			printf("mu vy  chiM dmudvy\n");
		}
		for (j=0; j<sizemumu; j++) {
			vy_op[p][j] = lin_interp(muopen[p], xbaropen[p], mumu[j], maxj[p], 2700);
			chiMax_op[p][j] = lin_interp(muopen[p], chiMopen[p], mumu[j], maxj[p], 2700);
			dmudvy_op[p][j] = lin_interp(muopen[p], openorbitopen[p], mumu[j], maxj[p], 2700);
			if (DEBUG == 1) 
				printf("%d/%d %f %f %f %f\n", j, sizemumu, mumu[j], vy_op[p][j], chiMax_op[p][j], dmudvy_op[p][j]);
		}

		// If you love your variables (and your memory) set them free // wise words Robbie
		free(chiMopen);//
		free(openorbitopen);//
		free(xbaropen);//
		free(Ucritf);//
		free(muopen);//
		free(xx);//
		free(jmclosed);//
		free(jmopen);//
		free(xifunction);//
		free(xbar);//
		free(chimin);//
		free(chiMax);//
		free(chimpp);//
		free(crossed_min);//
		free(crossed_max);//
		free(openorbit);//
		free(upperlimit);//
		free(lowerlimit);//
		free(xtop);//
		free(imax);//
		free(imin);//

		for (w = 0; w < sizexbar; w++) {
			free(chi[w]);//
			free(Uperp[w]);//
			free(mu[w]);//
			free(upper[w]);//
			for (s = 0; s < size_finegrid; s++)
				free(vx[w][s]);//
			free(vx[w]);//
		}
		free(chi);//
		free(Uperp);//
		free(mu);//
		free(upper);//
		free(vx);//

		free(phi);//
		free(gg);//
		free(ff);//
		free(phi_x);//

		for (j=0; j<j_inf; j++) {
			free(vxinf[j]);//
		}
		free(vxinf); //
		free(chiinf); //
		free(muinf);//

	}
	clock_t end = clock(); // finds the end time of the computation
	double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("in iondens2d: module ran in %f seconds\n", jobtime);
	return;
}
// closed densfinorb function 







// OLD COMMENTED OUT PART OF denszeroorb
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
