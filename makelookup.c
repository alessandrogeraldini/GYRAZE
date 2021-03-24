
// Author: Robbie Ewart 
// MODIFIED: 15/08/2019
/*This code generates a lookup table for values of phi given values of the electron density and should be run at the start of the program*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mps.h"

double FF(double beta) {
	double integrand, integral=0.0, beta_s = 0.005, betaval;
	int ind, num_beta= (int) (beta/beta_s);
	for (ind = 0; ind < num_beta+1; ind++) {
		betaval = beta*ind/num_beta;
		integrand = 0.5*(1.0-cos(2.0*betaval)) / ( M_PI - betaval + 0.5*sin(2.0*betaval) );
		integral += integrand*beta/num_beta;
	}
        //#returnvalue = beta**3/(3*np.pi) 
        return integral;
}

double FFp(double beta) {
	double integrand = 0.5*(1.0-cos(2.0*beta)) / ( M_PI - beta + 0.5*sin(2.0*beta) );
        return integrand;
}

double vparcut(double beta, double vcut, double phi_rel_DSE) {
	double vparcutval = sqrt((1.0-exp(-2.0*FF(beta)))*vcut*vcut/(sin(beta)*sin(beta)) + 2.0*phi_rel_DSE) ;
	return vparcutval;
}

double mucut(double beta, double vcut) {
	double mucutval = 0.5*exp(-2.0*FF(beta))*pow(vcut/sin(beta), 2.0);
	return mucutval;
}

void makelookup(double *phi,double *ne, int p_size,double v_cut, double v_cutDS, double *Phi_e_point)
{
	//defining des variables
	int i, vi, vi_1, p, len_F, cols = 0;//i is a counting variable, vi is a counting variable that will be saved for sums over velocity space, vi_1 is a special value in velocity space devoted to the first point, p is a counting variable that will be saved for counting over phi space, len_F saves the number of entries in the distribution;

	double Phi_e=0.0, mioverme = 1836.0;
	double gamma_e = 100000.0;
	double v_max, v_s, v_min;//v_max is the maximum velocity the distribution function will go up to before effectively just reading 0 from there on, v_s is the separation in velocity space between ajacent points. v_cut is the velocity that would be require in order to just reach the plasma wall boundary, v_min is the minimum velocity (the entrance to the presheath) required to reach a given x

	double sqrt_up, sqrt_lo, theta_up, theta_lo, sqrt_cut, theta_cut; //sqrt_up/lo reprresent the upper and lower limits of one strip integral, theta_up/lo are related to the hyperbolic arcsinh of some specific values (see document that hopefully exists)

	double* F, * Fp, * Fpp; // the zeroth first and second derivatives of the distribution function respectively

	double* in_err, * af_err, * n_res, w, w_up, w_lo; // for calculating errors in different parts of the code

	double* n_pre; //n_pre[x] is an array that will contain the values that we expect the integral to give (for a maxwellian input)

	int show_err, p_option;// change show_err to 1 or 0 if you do/don't want the error output to be printed to a file (note that the error output will only be meaningful for maxwellian inputs). Change p_option to control how the phi values are distributed with 1 being in an attempt to make ne roughly linear and 0 making phi linear;
	show_err = 0;
	p_option = 0;
	int upgraded = 1, mue_ind, mue_size = 140, size_cut=200;
	double *nepart, mue_s = 0.05, *mue_cut_lookup=malloc((1+size_cut)*sizeof(double)), *vpar_cut_lookup=malloc((1+size_cut)*sizeof(double)), vpar_cut;
	FILE* fptr1;

	char line_pot2[200], line_pot1[20];
	double* storevalspot1, * storevalspot2; //temporarily holds values to be moved around

	//getting some of the parameters
	fptr1 = fopen("dist_parameters.txt", "r");
	if (fptr1 == NULL)
	{
		printf("Crap");
	}
	i = 0;
	while (fgets(line_pot1, 20, fptr1) != NULL)
	{
		storevalspot1 = linetodata(line_pot1, strlen(line_pot1), &cols);
		if (i == 0)
		{
			v_max = *storevalspot1;
		}
		else
		{
			if (i == 1)
			{
				//v_cut = *storevalspot1;
			}
		}
		i++;
	}



	//Opening the distribution file to determine the the number of values of the distribution function
	fptr1 = fopen("dist_file.txt", "r");
	if (fptr1 == NULL)
	{
		printf("Crap");
	}

	//first reading through the file to find the values of the distribution function
	i = 0;
	while (fgets(line_pot2, 200, fptr1) != NULL)
	{
		i++;
	}
	len_F = i;
	F = malloc(len_F * sizeof(double));
	Fp = malloc(len_F * sizeof(double));
	Fpp = malloc(len_F * sizeof(double));
	in_err = malloc(p_size * sizeof(double));
	af_err = malloc(p_size * sizeof(double));
	n_res = malloc(p_size * sizeof(double));
	n_pre = malloc(p_size * sizeof(double));

	v_s = v_max / (len_F - 1);
	for (p = 0;p < p_size;p++)
	{
		in_err[p] = 0;
		af_err[p] = 0;
		n_res[p] = 0;
		n_pre[p] = 0;
		ne[p] = 0;
	}

	rewind(fptr1); // rewind through the files
	//reread throught the files now storing the values of the Distribution function
	i = 0;
	while (fgets(line_pot2, 200, fptr1) != NULL)
	{
		storevalspot2 = linetodata(line_pot2, strlen(line_pot2), &cols);
		F[i] = *storevalspot2;
		i++;
	}
	fclose(fptr1);

	//Generate the first and second derivative of the distribution function
	for (vi = 0; vi < len_F; vi++)
	{
		if ((vi < len_F - 1) && (vi > 0))
		{
//			printf("this part\n");
			Fp[vi] = (F[vi + 1] - F[vi - 1]) / (2.0 * v_s);
			Fpp[vi] = (F[vi + 1] + F[vi - 1] - (2.0 * F[vi])) / (pow(v_s, 2.0));
		}
		else
		{
			if (vi < len_F - 1)
			{
				Fp[vi] = (F[vi + 1] - F[vi]) / (v_s);
				Fpp[vi] = (F[vi + 2] + F[vi] - (2.0 * F[vi + 1])) / (pow(v_s, 2.0));
			}
			else
			{
				if (vi > 0)
				{
					Fp[vi] = (F[vi] - F[vi - 1]) / (v_s);
					Fpp[vi] = (F[vi - 2] + F[vi] - (2.0 * F[vi - 1])) / (pow(v_s, 2.0));
				}
				else
				{
					printf("Fuck");
				}
			}

		}


	}
	//Generating the test potential

	

	if (upgraded == 0) {
		if (p_option == 1)
		{
			for (p = 0; p < p_size; p++)
			{
				phi[p] = -(pow(v_cut, 2.0) / 2.0) + log(1 + (((double)p) * (exp(pow(v_cut, 2.0) / 2.0) - 1.0) / (((double)p_size) - 1.0)));
			}
		}
		else // we need a lookup table and so the potential in this module is just linear in the position index p for convenience
		{
			for (p = 0; p < p_size; p++)
			{
				phi[p] = -(pow(v_cut, 2.0) / 2.0) + ((pow(v_cut, 2.0) / 2.0) * ((double)p / ((double)p_size - 1.0)));
				//printf("%.15f\n", phi[p]);
			}
		}
	}

	//the integration

	if (upgraded == 1) {
	nepart = malloc(mue_size*sizeof(double));
	// for finite electron gyroradius effects, we need to integrate in mu as well
	printf("v_cut = %f, v_cutDS = %f, phiDSE = %f\n", v_cut, v_cutDS, phi[0]);
	//v_cutDS = sqrt(2.0*log(0.05) + log(1800/M_PI));
	// cutoff function for large rhoe
	vpar_cut_lookup[0] = 1e10;
	mue_cut_lookup[0] = 0.0;
	for (i=1; i< size_cut; i++) {
		vpar_cut_lookup[i] = vparcut(M_PI - i*M_PI/size_cut, v_cutDS, - phi[0]);
		mue_cut_lookup[i] = mucut(M_PI - i*M_PI/size_cut, v_cutDS);
		//printf("%f %f\n", vpar_cut_lookup[i], mue_cut_lookup[i]);
	}
	vpar_cut_lookup[size_cut] = sqrt(- 2.0*phi[0]) ;
	mue_cut_lookup[size_cut] = 1e10;
i=size_cut;
		//printf("%f %f\n", vpar_cut_lookup[i], mue_cut_lookup[i]);

	// ELECTRON CURRENT

	for (mue_ind=0; mue_ind<mue_size; mue_ind++) {
		nepart[mue_ind] = 0.0;
		if (gamma_e <= 1000.0) { //with finite gamma have a function to interpolate
			vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mue_ind*mue_s, size_cut, 205);
			//printf("\n");
			if (mue_ind == 0) vpar_cut = 1e10;
		}
		else { // with large gamma take the small gyrradius limit
			vpar_cut = v_cut;
		}

		for (vi = 0; vi < len_F - 1; vi++) {
			F[vi]   = vi*v_s*exp(- mue_ind*mue_s - 0.5*vi*vi*v_s*v_s);
			Fp[vi]  = (1.0 - vi*v_s*vi*v_s)*exp(- mue_ind*mue_s - 0.5*vi*vi*v_s*v_s);
			Fpp[vi] = (-3.0*vi*v_s + vi*v_s*vi*v_s*vi*v_s )*exp(- mue_ind*mue_s - 0.5*vi*vi*v_s*v_s);
		}
		//F = fullF[mue_ind];
		//Fp = fullFp[mue_ind];
		//Fpp = fullFpp[mue_ind];
		if (vpar_cut >= v_max) // effectively no cut off
			Phi_e = 0.0;
		else { // now there is a cut-off
			v_min = 0.0;
			if ((int)floor(v_min / v_s) >= len_F - 1)
			{
				nepart[mue_ind] = 0.0;
			}
			else
			{
				vi_1 = -1;
				for (vi = vi_1 + 1; vi < len_F - 1; vi++)
				{
					sqrt_up = (vi + 1) * v_s;
					sqrt_lo = (vi * v_s);

					nepart[mue_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
					nepart[mue_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

					w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
					w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
					af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
					af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
				}

				for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
				{
					sqrt_up = (vi + 1) * v_s;
					sqrt_lo = (vi * v_s);

					nepart[mue_ind] -= (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
					nepart[mue_ind] -= ( ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))) );

					w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
					w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
					af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
					af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
				}
				vi_1 = (int)floor(vpar_cut / v_s);
				sqrt_up = vpar_cut;
				sqrt_lo = vi_1 * v_s;

				nepart[mue_ind] -= (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
				nepart[mue_ind] -= ( ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))) );

				w_up = 0.5 * pow(vpar_cut, 2.0);
				w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
				af_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
				af_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
				//else
				//{
				//	if ((int)floor(v_min / v_s) == (int)floor(vpar_cut / v_s))
				//	{
				//		vi_1 = (int) floor(v_min / v_s);
				//		sqrt_up = (vi_1 + 1) * v_s ;
				//		theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) ) / (v_min));

				//		w_up = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
				//		in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * erf(sqrt(w_up))));
				//		in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

				//		sqrt_cut = vpar_cut ;
				//		theta_cut = asinh(vpar_cut / v_min);

				//		nepart[mue_ind] += (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * vpar_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
				//		nepart[mue_ind] += (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (vpar_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));

				//		w_lo = 0.5 * (pow(vpar_cut, 2.0) - pow(v_min, 2.0));
				//		in_err[p] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * vpar_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
				//		in_err[p] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (vpar_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));

				//		for (vi = vi_1 + 1; vi < len_F - 1; vi++)
				//		{
				//			sqrt_up = (vi + 1) * v_s ;
				//			sqrt_lo = vi * v_s;
				//			theta_up = asinh((vi + 1) * v_s  / v_min);
				//			theta_lo = asinh(vi * v_s / v_min);

				//			nepart[mue_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
				//			nepart[mue_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

				//			w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
				//			w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
				//			af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
				//			af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
				//		}
				//	}
				//	else //the most likely case
				//	{
				//		vi_1 = (int)floor(v_min / v_s);
				//		sqrt_up = (vi_1 + 1) * v_s;
				//		theta_up = asinh((vi_1 + 1) * v_s / v_min);

				//		//nepart[mue_ind] += (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
				//		//nepart[mue_ind] += (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

				//		w = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
				//		in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * erf(sqrt(w))));
				//		in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

				//		for (vi = vi_1 + 1; vi < len_F - 1; vi++)
				//		{
				//			sqrt_up = (vi + 1) * v_s;
				//			sqrt_lo = vi * v_s;
				//			theta_up = asinh((vi + 1) * v_s  / (v_min));
				//			theta_lo = asinh(vi * v_s / v_min);

				//			nepart[mue_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
				//			nepart[mue_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

				//			w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
				//			w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
				//			af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
				//			af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
				//		}

				//		for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
				//		{
				//			sqrt_up = (vi + 1) * v_s;
				//			sqrt_lo = vi * v_s;
				//			theta_up = asinh((vi + 1) * v_s / v_min);
				//			theta_lo = asinh(vi * v_s  / v_min);

				//			nepart[mue_ind] -= (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
				//			nepart[mue_ind] -= ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

				//			w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
				//			w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
				//			af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5  * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
				//			af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
				//		}

				//		vi_1 = (int)floor(vpar_cut / v_s);
				//		sqrt_up = vpar_cut;
				//		sqrt_lo = vi_1 * v_s;
				//		theta_up = asinh(vpar_cut  / v_min);
				//		theta_lo = asinh(vi_1 * v_s / v_min);

				//		nepart[mue_ind] -= (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
				//		nepart[mue_ind] -= ( ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))) );

				//		w_up = 0.5 * (pow(vpar_cut, 2.0) - pow(v_min, 2.0));
				//		w_lo = 0.5 * (pow((vi_1 * v_s), 2.0) - pow(v_min, 2.0));
				//		in_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
				//		in_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));
				//	}
				//}
			//n_res[p] = ne[p] - n_pre[p];
			}
		}
		if (mue_ind != 0) 
			Phi_e += 2.0*M_PI*0.5*(nepart[mue_ind] + nepart[mue_ind -1 ])*mue_s;
	}

	for (p = 0; p < p_size; p++)
	{
		ne[p] = 0.0;
		n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] + (0.5 * pow(v_cut, 2.0)))));
    		v_min = sqrt(-2.0 * phi[p]);
 		for (mue_ind=0; mue_ind<mue_size; mue_ind++) {
			nepart[mue_ind] = 0.0;
			vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mue_ind*mue_s, size_cut, 205);
			//if ( (v_cut*v_cut + 2.0*phi[p] < 1e-10) ) vpar_cut = v_cut;
			//else vpar_cut = lin_interp(mue_cut_lookup, vpar_cut_lookup, mue_ind*mue_s, size_cut, 205);
			//if (fabs(vpar_cut) < 1e-10) vpar_cut = 1e-10;
			if (mue_ind == 0) vpar_cut = 1e10;
			for (vi = 0; vi < len_F - 1; vi++) {
				F[vi]   = exp(- mue_ind*mue_s - 0.5*vi*vi*v_s*v_s);
				Fp[vi]   = -vi*v_s*exp(- mue_ind*mue_s - 0.5*vi*vi*v_s*v_s);
				Fpp[vi]   = (vi*v_s*vi*v_s - 1.0)*exp(- mue_ind*mue_s - 0.5*vi*vi*v_s*v_s);
			}
			//F = fullF[mue_ind];
			//Fp = fullFp[mue_ind];
			//Fpp = fullFpp[mue_ind];
			if (vpar_cut >= v_max) // effectively no cut off
			{
				if ((int)floor(v_min / v_s) >= len_F - 1)
				{
					//ne[p] = 0;
					nepart[mue_ind] = 0.0;
					printf("v_min outside velocity grid");
				}
				else
				{
					if (v_min == 0.0)//special case where the integration becomes simple
					{
						for (vi = 0; vi < len_F - 1; vi++)
						{
							//ne[p] += (F[vi] + F[vi + 1]) * v_s;
							nepart[mue_ind] += (F[vi] + F[vi + 1]) * v_s;
						}
					}
					else //regualar case
					{


						vi_1 = (int)floor(v_min / v_s);
						sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p]));
						theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p])) / v_min);

						nepart[mue_ind] += 2.0 * (((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
						nepart[mue_ind] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));

						w = ((0.5 * pow((((vi_1 + 1) * v_s) - v_min), 2.0)) + (v_min * (((vi_1 + 1) * v_s) - v_min)));
						in_err[p] += 2.0 * ((((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w))));
						in_err[p] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));

						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
						{
							sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
							sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
							theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
							theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

							nepart[mue_ind] += 2.0 * (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							nepart[mue_ind] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

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
				//if ((phi[p] + (0.5 * pow(vpar_cut, 2.0))) > 0)
				//{
				//	n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] + (0.5 * pow(vpar_cut, 2.0)))));
				//}
				//else
				//{
				//	n_pre[p] = 0.5 * exp(phi[p]);
				//}
				v_min = sqrt(-2.0 * phi[p]);
				if ((int)floor(v_min / v_s) >= len_F - 1)
				{
					nepart[mue_ind] = 0.0;
				}
				else
				{
					if (v_min == 0.0)
					{
						vi_1 = 0;
						sqrt_up = ((vi_1 + 1) * v_s);
						nepart[mue_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
						nepart[mue_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));

						w = 0.5 * pow(v_s, 2.0);
						af_err[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s)))) - erf(sqrt(w));
						af_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));

						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
						{
							sqrt_up = (vi + 1) * v_s;
							sqrt_lo = (vi * v_s);

							nepart[mue_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
							nepart[mue_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
							w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
							af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
							af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
						}

						for (vi = vi_1 + 1; vi < (int)floor(vpar_cut / v_s); vi++)
						{
							sqrt_up = (vi + 1) * v_s;
							sqrt_lo = (vi * v_s);

							nepart[mue_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
							nepart[mue_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

							w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
							w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
							af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
							af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
						}
						vi_1 = (int)floor(vpar_cut / v_s);
						sqrt_up = vpar_cut;
						sqrt_lo = vi_1 * v_s;

						nepart[mue_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
						nepart[mue_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));

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

							nepart[mue_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
							nepart[mue_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

							w_up = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
							in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w_up))));
							in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

							sqrt_cut = sqrt(pow(vpar_cut, 2.0) + 2.0 * phi[p]);
							theta_cut = asinh(sqrt(pow(vpar_cut, 2.0) + 2.0 * phi[p]) / (v_min));

							nepart[mue_ind] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * vpar_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
							nepart[mue_ind] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (vpar_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));

							w_lo = 0.5 * (pow(vpar_cut, 2.0) - pow(v_min, 2.0));
							in_err[p] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * vpar_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
							in_err[p] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (vpar_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));

							for (vi = vi_1 + 1; vi < len_F - 1; vi++)
							{
								sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
								sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
								theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
								theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

								nepart[mue_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
								nepart[mue_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

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

							nepart[mue_ind] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
							nepart[mue_ind] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

							w = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
							in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w))));
							in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

							for (vi = vi_1 + 1; vi < len_F - 1; vi++)
							{
								sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
								sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
								theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
								theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

								nepart[mue_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
								nepart[mue_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

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

								nepart[mue_ind] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
								nepart[mue_ind] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

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

							nepart[mue_ind] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							nepart[mue_ind] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));

							w_up = 0.5 * (pow(vpar_cut, 2.0) - pow(v_min, 2.0));
							w_lo = 0.5 * (pow((vi_1 * v_s), 2.0) - pow(v_min, 2.0));
							in_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
							in_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((vpar_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));
						}
					}
				//n_res[p] = ne[p] - n_pre[p];
				}
			}
			if (mue_ind != 0) 
				ne[p] += 2.0*M_PI*0.5*(nepart[mue_ind] + nepart[mue_ind -1])*mue_s;
		}
	}
	//END OF PART UNDER CONSTRUCTION
	}
	else {
	// PART BELOW SHOULD EVENTUALLY BE REPLACED (MORE LIKE UPGRADED) BY THE ABOVE
	if (v_cut >= v_max)//effectively no cut off
	{

		for (p = 0; p < p_size; p++)
		{
			n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] + (0.5 * pow(v_cut, 2.0)))));
			v_min = sqrt(-2.0 * phi[p]);
			if ((int)floor(v_min / v_s) >= len_F - 1)
			{
				ne[p] = 0;
				printf("v_min outside velocity grid");
			}
			else
			{
				if (v_min == 0.0)//special case where the integration becomes simple
				{
					for (vi = 0; vi < len_F - 1; vi++)
					{
						ne[p] += (F[vi] + F[vi + 1]) * v_s;
					}
				}
				else //regualar case
				{
					vi_1 = (int)floor(v_min / v_s);
					sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p]));
					theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p])) / v_min);

					ne[p] += 2.0 * (((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
					ne[p] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));

					w = ((0.5 * pow((((vi_1 + 1) * v_s) - v_min), 2.0)) + (v_min * (((vi_1 + 1) * v_s) - v_min)));
					in_err[p] += 2.0 * ((((F[vi_1 + 1] - (((vi_1 + 1) * v_s) * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w))));
					in_err[p] += 2.0 * ((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up) + ((1.0 / 6.0) * pow(sqrt_up, 3.0) * Fpp[vi_1 + 1]) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));

					for (vi = vi_1 + 1; vi < len_F - 1; vi++)
					{
						sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
						sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
						theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
						theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

						ne[p] += 2.0 * (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
						ne[p] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
						w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
						af_err[p] += 2.0 * ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
						af_err[p] += 2.0 * ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
					}
				}

			}
			n_res[p] = ne[p] - n_pre[p];
		}

	}
	else
	{
		//now instead there is a cutoff
		for (p = 0;p < p_size; p++)
		{
			if ((phi[p] + (0.5 * pow(v_cut, 2.0))) > 0)
			{
				n_pre[p] = 0.5 * exp(phi[p]) * (1 + erf(sqrt(phi[p] + (0.5 * pow(v_cut, 2.0)))));
			}
			else
			{
				n_pre[p] = 0.5 * exp(phi[p]);
			}
			v_min = sqrt(-2.0 * phi[p]);
			if ((int)floor(v_min / v_s) >= len_F - 1)
			{
				ne[p] = 0;
			}
			else
			{
				if (v_min == 0.0)
				{
					vi_1 = 0;
					sqrt_up = ((vi_1 + 1) * v_s);

					ne[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
					ne[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));

					w = 0.5 * pow(v_s, 2.0);
					af_err[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s)))) - erf(sqrt(w));
					af_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));

					for (vi = vi_1 + 1; vi < len_F - 1; vi++)
					{
						sqrt_up = (vi + 1) * v_s;
						sqrt_lo = (vi * v_s);

						ne[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
						ne[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
						w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
						af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
						af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
					}

					for (vi = vi_1 + 1; vi < (int)floor(v_cut / v_s); vi++)
					{
						sqrt_up = (vi + 1) * v_s;
						sqrt_lo = (vi * v_s);

						ne[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
						ne[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

						w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
						w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
						af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
						af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
					}
					vi_1 = (int)floor(v_cut / v_s);
					sqrt_up = v_cut;
					sqrt_lo = vi_1 * v_s;

					ne[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
					ne[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));

					w_up = 0.5 * pow(v_cut, 2.0);
					w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
					af_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
					af_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
				}
				else
				{
					if ((int)floor(v_min / v_s) == (int)floor(v_cut / v_s))
					{
						vi_1 = (int)floor(v_min / v_s);
						sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p]));
						theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + (2.0 * phi[p])) / (v_min));

						ne[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
						ne[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

						w_up = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
						in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w_up))));
						in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

						sqrt_cut = sqrt(pow(v_cut, 2.0) + 2.0 * phi[p]);
						theta_cut = asinh(sqrt(pow(v_cut, 2.0) + 2.0 * phi[p]) / (v_min));

						ne[p] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * v_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));
						ne[p] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (v_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));

						w_lo = 0.5 * (pow(v_cut, 2.0) - pow(v_min, 2.0));
						in_err[p] -= (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * (sqrt_up - sqrt_cut)) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) - (sqrt_cut * v_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
						in_err[p] -= (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * (sqrt_up - sqrt_cut))) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * (pow(sqrt_up, 3.0) - pow(sqrt_cut, 3.0))) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) - (v_cut * sqrt_cut) + (pow(v_min, 2.0) * (theta_up - theta_cut)))));

						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
						{
							sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
							sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
							theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
							theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

							ne[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							ne[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

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

						ne[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
						ne[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

						w = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
						in_err[p] += 2.0 * ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi[p]) * erf(sqrt(w))));
						in_err[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

						for (vi = vi_1 + 1; vi < len_F - 1; vi++)
						{
							sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]);
							sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]);
							theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi[p]) / (v_min));
							theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi[p]) / (v_min));

							ne[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							ne[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

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

							ne[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
							ne[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

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

						ne[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
						ne[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));

						w_up = 0.5 * (pow(v_cut, 2.0) - pow(v_min, 2.0));
						w_lo = 0.5 * (pow((vi_1 * v_s), 2.0) - pow(v_min, 2.0));
						in_err[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
						in_err[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))));



					}
				}
				n_res[p] = ne[p] - n_pre[p];
			}
		}
	}
	}


	Phi_e =  Phi_e / ne[p_size - 1];
	printf("Phi_e = %f\n", Phi_e);
	for (p = 0;p < p_size; p++)
	{
		//printf("ne[%d] = %f\n", p, ne[p]);
		ne[p] = ne[p] / ne[p_size - 1];
		n_res[p] = ne[p] - n_pre[p];
		//if (ne[p] < 0.0) 
		//printf("WARNING: The electron density is negative!\n ne[%d] = %f\n", p, ne[p]);
		//printf("n_pre[%d] = %f\n", p, n_pre[p]);
		//printf("ne[%d] = %f\n", p, ne[p]);
	}

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
	free(F);
	free(Fp);
	free(Fpp);
	free(in_err);
	free(af_err);
	free(n_res);
	free(n_pre);
	*Phi_e_point = Phi_e*sqrt(mioverme/2.0);
}
