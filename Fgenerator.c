// Author: basically Alessandro but all the bugs are due to Robbie 


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "Headerfiles/linetodata.h"

#define M_PI acos(-1)

void Fgenerator(double phi_cut)
{
	double Te, alpha, Ts, maxE, dE, l, dEnarrow, U1, condition, flow, u, chodura;
	double normalisation, mu, U, FF, correction;
	char line[2000];
	double* storevals;
	int i = 0, j = 0 ,coldelectrons, col;

	correction = 1 + (exp(phi_cut) / ((1 + erf(sqrt(-phi_cut))) * sqrt(-M_PI * phi_cut))); //corrects the chodura condition for having a lower cutoff 

	FILE* input, *fp;
	if ((input = fopen("inputfile.txt", "r")) == NULL)
	{ // Check for presence of file
		printf("Cannot open %s\n", "inputfile.txt");
		exit(1);
	}
	while (fgets(line, 20, input) != NULL)
	{
		//string = malloc(strlen(line)*sizeof(char));
		//string = line;
		//printf("The string is %s\n", string);
		//storevals = linetodata(string, strlen(string), &col);
		storevals = linetodata(line, strlen(line), &col);

		if (i == 0)
		{
			alpha = *storevals;
		}
		else if (i == 1)
		{
			Te = *storevals;
		}
		i += 1;
	}
	fclose(input);
	printf("Te is %f\n", Te);
	U1 = 1.0;
	dE = 0.05;

	if (Te/correction > 0.999)
	{
		maxE = 18.0 + (4.0 * Te/correction);
		dE = 0.05;
		l = (int)floor(maxE / dE);
		coldelectrons = 0;
	}
	else if ((Te/correction > 0.00001) && (Te/correction < 0.999))
	{
		maxE = 18.0;
		dEnarrow = dE * Te/correction;
		l = (int)floor(maxE / dE) + (int)floor(U1 / dEnarrow) - (int)floor(U1 / dE);
		coldelectrons = 0;
	}
	else
	{
		coldelectrons = 1;
	}


	

	u = 0.0;
	condition = (2.0*correction) / Te;
	printf("condition is %f\n", condition);
	flow = 0.0; // this needs to be inputted but seems to have only been important in the prototyping phase???
	if (coldelectrons == 0)
	{
		if (Te/correction > 0.999)
		{
			chodura = 9999999.0;
			while ((chodura > condition)||(chodura<(condition -0.005)))
			{
				if (chodura > condition)
				{
					u += 0.0005;
				}
				else if (chodura < condition)
				{
					u -= 0.0005;
				}
				normalisation = ((1.0 + erf(u)) * (1.0 + (2.0 * pow(u, 2.0)))) + ((2.0 / sqrt(M_PI)) * u * exp(-pow(u, 2.0)));
				chodura = 2.0 * (1 + erf(u)) / normalisation;
				//flow = (0.5 * (1 + u * *2) * np.exp(-u * *2) + (1.5 + u * *2) * (np.sqrt(np.pi) * u / 2.0) * (1 + sp.erf(u))) * 4 / sqrt(np.pi) / normalization;
			}
		}
		else
		{
			chodura = 0.0;
			while ((chodura > condition)||(chodura<(condition -0.005)))
			{
				if (chodura > condition)
				{
					u -= 0.0005;
				}
				else if (chodura < condition)
				{
					u += 0.01;
				}
				if (u > 0.06)
				{
					normalisation = ((sqrt(M_PI * u)) - M_PI * exp(1.0 / u) * (1.0 - erf(1.0 / sqrt(u)))) / (2.0 * pow(u, 1.5));
					chodura = M_PI * exp(1.0 / u) * (1.0 - erf(1 / sqrt(u))) / (2.0 * sqrt(u) * normalisation);
				}
				else
				{
					normalisation = (sqrt(M_PI) / 4.0) - ((3.0 * sqrt(M_PI) / 8.0) * u) + ((15.0 * sqrt(M_PI) / 16.0) * (pow(u,2.0))) - (105.0 * sqrt(M_PI) / 32.0) * (pow(u,3.0)) + (945.0 * sqrt(M_PI) / 64.0) * (pow(u,4.0));
					chodura = (sqrt(M_PI) / normalisation) * (0.5 - (u / 4.0) + (3.0 * (pow(u,2.0)) / 8.0) - (15.0 * (pow(u,3.0)) / 16.0) + (105.0 * (pow(u,4.0)) / 32.0));
				}
            //printf("chodura is %f\n", chodura);
    			
			
			}
		}
	}
	printf("u is now %f\n", u);
	if ((fp = fopen("distfuncin.txt", "w")) == NULL)
	{
		printf("Unable to open distfuncin.txt");
		exit(1);
	}
	//fprintf(fp, " ");
	if (coldelectrons == 0)
	{
		if (Te/correction > 0.999) //In this case the distribution is a maxwellian tapered at the origin with velocity centered at the nessecary velocity to marginally satisfy the Bohm condition
		{
			for (i = 0; i < l; i++)
			{
				mu = ((double)i) * dE;
				for (j = 0; j < l; j++)
				{
					U = j * dE;
					if (U > mu)
					{
						FF = (1.0 / pow(M_PI, 1.5)) * (1.0 / normalisation) * (U - mu) * exp(-U + (2.0 * u * sqrt(U - mu)) - pow(u, 2.0));
						if (j == l - 1)
						{
							fprintf(fp, "%e\n", FF);
						}
						else
						{
							fprintf(fp, "%e ", FF);
						}
					}
					else
					{
						FF = 0.0;
						if (j == l - 1)
						{
							fprintf(fp, "%d\n", 0);
						}
						else
						{
							fprintf(fp, "%d ", 0);
						}
					}


				}
		
			}
			fclose(fp);
			if ((fp = fopen("Umufile.txt", "w")) == NULL)
			{
				printf("Unable to open Umufile.txt");
				exit(1);
			}
			for (i = 0;i < l; i++)
			{
				mu = ((double)i) * dE;
				if (i == l - 1)
				{
					fprintf(fp, "%f\n", mu);
				}
				else
				{
					fprintf(fp, "%f ", mu);
				}
			}
			for (j = 0;j < l; j++)
			{
				U = j * dE;
				if (i == l - 1)
				{
					fprintf(fp, "%f\n", U);
				}
				else
				{
					fprintf(fp, "%f ", U);
				}
			}
			fclose(fp);
		}
		else // Here the distribution is again a tapered maxwellian but instead it approaches a half maxwellian
		{
			for (i = 0;i < l;i++)
			{
				mu = ((double)i) * dE;
				for (j = 0;j < l;j++)
				{
					if (j < (int)floor(U1 / dEnarrow))
					{
						U = ((double)j) * dEnarrow;
					}
					else
					{
						U = (double)((j - ((int)floor(U1 / dEnarrow))) + ((int)floor(U1 / dE))) * dE;
					}
					if (U > mu)
					{
						FF = (1.0 / (4.0 * M_PI)) * (1.0 / normalisation) * ((U - mu) / (1.0 + (u * (U - mu)))) * exp(-U);
						if (j == l - 1)
						{
							fprintf(fp, "%e\n", FF);
						}
						else
						{
							fprintf(fp, "%e ", FF);
						}
					}
					else
					{
						FF = 0;
						if (j == l - 1)
						{
							fprintf(fp, "%d\n", 0);
						}
						else
						{
							fprintf(fp, "%d ", 0);
						}
					}
				}
			}
			fclose(fp);
			if ((fp = fopen("Umufile.txt", "w")) == NULL)
			{
				printf("Unable to open Umufile.txt");
				exit(1);
			}
			for (i = 0; i < l; i++)
			{
				mu = (double)i * dE;
				if (i == l - 1)
				{
					fprintf(fp, "%f\n", mu);
				}
				else
				{
					fprintf(fp, "%f ", mu);
				}
			}
			for (j = 0; j < l; j++)
			{
				if (j < (int)floor(U1 / dEnarrow))
				{
					U = ((double)j) * dEnarrow;
				}
				else
				{
					U = (double)((j - ((int)floor(U1 / dEnarrow))) + ((int)floor(U1 / dE))) * dE;
				}
				if (j == l - 1)
				{
					fprintf(fp, "%f\n", U);
				}
				else
				{
					fprintf(fp, "%f ", U);
				}
			}
			fclose(fp);
		}



	}





}
