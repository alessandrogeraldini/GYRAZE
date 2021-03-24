
// Author: Robbie Ewart 
// MODIFIED: 15/08/2019
/*This code is very similar to makelookup.c in that it generates values of electron density and is given values of phi however,
the makelookup.c code takes one cut off value and calculates the corresponding electron density for all phi below this. The code
on the other hand is taking all reasonable cut off values and calculating the corresponding electron density at that cut off*/

//This afords a couple of simplifications since at the cut off the second integral (of the returning electrons) evaluates to zero
// for all possible input functions 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mps.h"
//#include "Headerfiles/linetodata.h"
//#include "Headerfiles/cutlookup.h"


void cutlookup(double* phi_cut, double* ne_cut, int cut_points)
{
	//defining des variables
	int i, vi, vi_1, p, len_F, cols = 0;//i is a counting variable, vi is a counting variable that will be saved for sums over velocity space, vi_1 is a special value in velocity space devoted to the first point, p is a counting variable that will be saved for counting over phi space, len_F saves the number of entries in the distribution;

	double v_max, v_s, v_cut, v_min, v_top, *ne0;//v_max is the maximum velocity the distribution function will go up to before effectively just reading 0 from there on, v_s is the separation in velocity space between ajacent points. v_cut is the velocity that would be require in order to just reach the plasma wall boundary, v_min is the minimum velocity (the entrance to the presheath) required to reach a given x

	double sqrt_up, sqrt_lo, theta_up, theta_lo; //sqrt_up/lo reprresent the upper and lower limits of one strip integral, theta_up/lo are related to the hyperbolic arcsinh of some specific values (see document that hopefully exists)

	double* F, * Fp, * Fpp; // the zeroth first and second derivatives of the distribution function respectively

	double* in_err, * af_err, *af_err0, * n_res, w, w_up, w_lo; // for calculating errors in different parts of the code

	double* n_pre; //n_pre[x] is an array that will contain the values that we expect the integral to give (for a maxwellian input)

	int show_err;// change show_err to 1 or 0 if you do/don't want the error output to be printed to a file (note that the error output will only be meaningful for maxwellian inputs). Change p_option to control how the phi values are distributed with 1 being in an attempt to make ne roughly linear and 0 making phi linear;
	show_err = 1;
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
				v_cut = *storevalspot1;
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
	in_err = malloc(cut_points * sizeof(double));
	af_err = malloc(cut_points * sizeof(double));
	af_err0 = malloc(cut_points * sizeof(double));
	n_res = malloc(cut_points * sizeof(double));
	n_pre = malloc(cut_points * sizeof(double));
	ne0 = malloc(cut_points* sizeof(double));

	v_s = v_max / (len_F - 1);

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
	v_top = 4;
	for (p = 0;p < cut_points; p++)
	{
		phi_cut[p] = -(((double)p) / ((double)(cut_points - 1.0))) * pow(v_top, 2.0) / 2.0;
	}





	//the integration

	for (p = 0;p < cut_points; p++)
	{

		v_min = sqrt(-2.0 * phi_cut[p]);
		ne_cut[p] = 0;
		if ((int)floor(v_min / v_s) >= len_F - 1)
		{
			ne_cut[p] = 0;
		}
		else if (v_min == 0.0)
		{
			vi_1 = 0;
			sqrt_up = ((vi_1 + 1) * v_s);

			ne_cut[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
			ne_cut[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));

			w = 0.5 * pow(v_s, 2.0);
			af_err[p] += (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s)))) - 0.5 * erf(sqrt(w));
			af_err[p] += (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));

			for (vi = vi_1 + 1; vi < len_F - 1; vi++)
			{
				sqrt_up = (vi + 1) * v_s;
				sqrt_lo = (vi * v_s);

				ne_cut[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
				ne_cut[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

				w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
				w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
				af_err[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
				af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
			}
		}
		else //the most likely case
		{
			vi_1 = (int)floor(v_min / v_s);
			sqrt_up = sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phi_cut[p]);
			theta_up = asinh(sqrt(pow(((vi_1 + 1) * v_s), 2.0) + 2.0 * phi_cut[p]) / (v_min));

			ne_cut[p] += (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up))));
			ne_cut[p] += (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

			w = 0.5 * (pow(((vi_1 + 1) * v_s), 2.0) - pow(v_min, 2.0));
			in_err[p] += ((((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * ((sqrt_up * ((vi_1 + 1) * v_s)) + (pow(v_min, 2.0) * theta_up)))) - (0.5 * exp(phi_cut[p]) * erf(sqrt(w))));
			in_err[p] += (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up) + (pow(v_min, 2.0) * theta_up))));

			for (vi = vi_1 + 1; vi < len_F - 1; vi++)
			{
				sqrt_up = sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi_cut[p]);
				sqrt_lo = sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi_cut[p]);
				theta_up = asinh(sqrt(pow(((vi + 1) * v_s), 2.0) + 2.0 * phi_cut[p]) / (v_min));
				theta_lo = asinh(sqrt(pow(((vi)* v_s), 2.0) + 2.0 * phi_cut[p]) / (v_min));


				ne_cut[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
				ne_cut[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));

				w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
				w_lo = 0.5 * (pow(((vi)* v_s), 2.0) - pow(v_min, 2.0));
				af_err[p] += ((((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo))))) - (0.5 * exp(phi_cut[p]) * (erf(sqrt(w_up)) - erf(sqrt(w_lo)))));
				af_err[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo)) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)) + (pow(v_min, 2.0) * (theta_up - theta_lo)))));
			}
		}

		v_min = 0;
		v_cut = sqrt(-2.0 * phi_cut[p]);

		vi_1 = 0;
		sqrt_up = ((vi_1 + 1) * v_s);

		ne0[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s))));
		ne0[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));

		w = 0.5 * pow(v_s, 2.0);
		af_err0[p] += 2.0 * (((F[vi_1 + 1] - ((vi_1 + 1) * v_s * Fp[vi_1 + 1])) * sqrt_up) + (0.5 * Fp[vi_1 + 1] * (sqrt_up * ((vi_1 + 1) * v_s)))) - erf(sqrt(w));
		af_err0[p] += 2.0 * (((0.5 * (pow(((vi_1 + 1) * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1 + 1] * sqrt_up)) + ((1.0 / 6.0) * Fpp[vi_1 + 1] * pow(sqrt_up, 3.0)) - (0.5 * ((vi_1 + 1) * v_s) * Fpp[vi_1 + 1] * ((((vi_1 + 1) * v_s) * sqrt_up))));

		for (vi = vi_1 + 1; vi < len_F - 1; vi++)
		{
			sqrt_up = (vi + 1) * v_s;
			sqrt_lo = (vi * v_s);

			ne0[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
			ne0[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

			w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
			w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
			af_err0[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
			af_err0[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
		}

		for (vi = vi_1 + 1; vi < (int)floor(v_cut / v_s); vi++)
		{
			sqrt_up = (vi + 1) * v_s;
			sqrt_lo = (vi * v_s);

			ne0[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo)))));
			ne0[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));

			w_up = 0.5 * (pow(((vi + 1) * v_s), 2.0) - pow(v_min, 2.0));
			w_lo = 0.5 * (pow((vi * v_s), 2.0) - pow(v_min, 2.0));
			af_err0[p] += (((F[vi] - ((vi * v_s) * Fp[vi])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
			af_err0[p] += ((0.5 * (pow((vi * v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi * v_s) * Fpp[vi] * (((((vi + 1) * v_s) * sqrt_up) - ((vi * v_s) * sqrt_lo))));
		}
		vi_1 = (int)floor(v_cut / v_s);
		sqrt_up = v_cut;
		sqrt_lo = vi_1 * v_s;

		ne0[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo)))));
		ne0[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));

		w_up = 0.5 * pow(v_cut, 2.0);
		w_lo = 0.5 * pow((vi_1 * v_s), 2.0);
		af_err0[p] += (((F[vi_1] - ((vi_1 * v_s) * Fp[vi_1])) * (sqrt_up - sqrt_lo)) + (0.5 * Fp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))))) - (0.5 * (erf(sqrt(w_up)) - erf(sqrt(w_lo))));
		af_err0[p] += ((0.5 * (pow(((vi_1)* v_s), 2.0) + pow(v_min, 2.0)) * Fpp[vi_1] * (sqrt_up - sqrt_lo))) + ((1.0 / 6.0) * Fpp[vi_1] * (pow(sqrt_up, 3.0) - pow(sqrt_lo, 3.0))) - (0.5 * (vi_1 * v_s) * Fpp[vi_1] * (((v_cut * sqrt_up) - ((vi_1 * v_s) * sqrt_lo))));
	
		//if (p == 128)
		//{
		//	printf("%.10f\n", ne_cut[p]);
		//	printf("%.10f\n", ne0[p]);
		//}
	}


	for (p = 0;p < cut_points; p++)
	{
		ne_cut[p] = ne_cut[p] / ne0[p];
	}

	if (show_err == 1)
	{

	}
	free(F);
	free(Fp);
	free(Fpp);
	free(in_err);
	free(af_err);
	free(af_err0);
	free(n_res);
	free(n_pre);
	free(ne0);

}
