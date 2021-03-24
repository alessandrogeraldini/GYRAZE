//Author: Robbie Ewart
//Date Created: 24/08/2019
/*generate_cspline generates a natural cubic spline that interpolates the n points (x[i],y[i]). Value of the spline at a
point x can be found with the function interpolatespline. The arguments of generate_cspline are as follows:
> struct spline *spl = a pointer to a spline structure you have already created (see spline_test for example). The
spline structure itself could be passed to such a function but with large structures it is better to pass pointers
> double *x = the x[i] array
> double *y = the y[i] array
> int n = the length of these arrays  
The function edit the spline provided giving it the members a,b,c,d. Calling the spline is done with call_spline
which takes the arguments
> struct spline *spl = again a pointer to the spline you created
> double *x = the x[i] array
> double given_x = the x value at which you wish to know the value of the spline
> int n = the size of the x[i] array
this will return the interpolated value of the spline.

(*basic facts about natural cubic splines*)
A spline interpolated n data points is made of n-1 piecewise cubic functions S_i(x), i E [1,n-1] such that:
> S_i(x[i]) = S_(i-1)(x[i]) = y[i] the functions agree on the boundary (knots)
> S_i'(x[i]) = S_(i-1)'(x[i]) the first derivative of the functions agree on the boundary
> S_i''(x[i]) = S_(i-1)''(x[i]) the second derivative of the functions agree on the boundary
> S_0''(x[0]) = S_(n-1)''x[n-1] = 0 the second derivative of the functions is zero at the end points

with our definition of the spline S_i is given by

S_i(x) = y[i] + b[i]*(x-x[i]) + c[i]*(x-x[i])^2 + d[i]*(x-x[i])^3
and this is the purpose of the arrays a,b,c,d

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mps.h"


//struct spline
//{
//	double* a;
//	double* b;
//	double* c;
//	double* d;
//
//};


void generate_cspline(struct spline *spl,double *x,double *y, int n)
{
	
	/* */
	int i;
	double* a, * b, * c, * d;
	a = malloc(n * sizeof(double));
	b = malloc((n - 1) * sizeof(double));
	c = malloc(n * sizeof(double));
	d = malloc((n - 1) * sizeof(double));
	double* h, * l, * u, * alpha, * z;
	h = malloc((n - 1) * sizeof(double));
	alpha = malloc((n - 1) * sizeof(double));
	l = malloc(n * sizeof(double));
	u = malloc(n * sizeof(double));
	z = malloc(n * sizeof(double));
	
	a[0] = y[0];
	a[n - 1] = y[n - 1];
	h[0] = x[1] - x[0];

	for (i = 1; i < n-1; i++)
	{
		a[i] = y[i];
		h[i] = x[i + 1] - x[i];
		alpha[i] = (3.0 * (y[i + 1] - y[i]) / (x[i + 1] - x[i])) - (3.0 * (y[i] - y[i - 1]) / (x[i] - x[i - 1]));
	}
	l[0] = 0;
	u[0] = 0;
	z[0] = 0;
	for (i = 1; i < n - 1; i++)
	{
		l[i] = (2 * (x[i + 1] - x[i])) - (h[i - 1] * u[i - 1]);
		u[i] = h[i] / l[i];
		z[i] = (alpha[i] - (h[i - 1] * z[i - 1])) / l[i];
	}
	l[n - 1] = 1;
	z[n - 1] = 0;
	c[n - 1] = 0;
	for (i = n - 2;i > -1; i--)
	{
		c[i] = z[i] - (u[i] * c[i + 1]);
		b[i] = ((a[i + 1] - a[i]) / h[i]) - (h[i] * (c[i + 1] + (2.0 * c[i])) / 3.0);
		d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
	}
	
	(*spl).a = malloc((n - 1) * (sizeof(double)));
	(*spl).b = malloc((n - 1) * (sizeof(double)));
	(*spl).c = malloc((n - 1) * (sizeof(double)));
	(*spl).d = malloc((n - 1) * (sizeof(double)));
	for (i = 0;i < n - 1;i++)
	{
		(*spl).a[i] = a[i];
		(*spl).b[i] = b[i];
		(*spl).c[i] = c[i];
		(*spl).d[i] = d[i];
	}
	
}

double call_spline(struct spline *spl,double *x,double given_x,int n,int line)
{
	
	double interp_y;
	int cont = 1; 
	int i;
	i = 0;
	
	while (cont == 1)
	{
		if (i == n-1)
		{
			printf("The x value exceeds the given x grid");
			printf("The x value was: %.10f\n", given_x);
			printf("Error on line %d",line);
			exit(EXIT_FAILURE);
		}
		if(given_x > x[i])
		{
			if (given_x <= x[i + 1])
			{
				interp_y = (*spl).a[i] + ((*spl).b[i]) * (given_x - x[i]) + ((*spl).c[i]) * pow(given_x - x[i], 2.0) + ((*spl).d[i]) * pow(given_x - x[i], 3.0);
				//printf("%.5f\n%.5f\n%.5f\n%.5f\n", (*spl).a[i], (*spl).a[i] + ((*spl).b[i]) * (given_x - x[i]), (*spl).a[i] + ((*spl).b[i]) * (given_x - x[i]) + ((*spl).c[i]) * pow(given_x - x[i], 2.0), (*spl).a[i] + ((*spl).b[i]) * (given_x - x[i]) + ((*spl).c[i]) * pow(given_x - x[i], 2.0) + ((*spl).d[i]) * pow(given_x - x[i], 3.0));
				//printf("%.5f",given_x - x[i]);
				//printf("%.5f", x[i + 1] - x[i]);
				cont = 0;
			}
		}
		i++;
	}
	return interp_y;
}


double lin_interp(double* x_grid, double* y_grid,double given_x ,int n, int line)
{
	double interp_y;
	int cont = 1;
	int i;
	i = 0;

	while (cont == 1)
	{
		if (i == n)
		{
			if (fabs(given_x - x_grid[0]) < 1e-10) {
				interp_y = y_grid[0];
				cont = 0;
			}
			else if (fabs(given_x - x_grid[n-1]) < 1e-10) {
				interp_y = y_grid[n-1];
				cont = 0;
			}
			else {
				printf("The x value exceeds the given x grid");
				printf("The x value was: %.10f\n", given_x);
				printf("Error on line %d", line);
				exit(EXIT_FAILURE);
			}
		}
		if (given_x > x_grid[i])
		{
			if (given_x <= x_grid[i + 1])
			{
				interp_y = y_grid[i] + (y_grid[i + 1] - y_grid[i]) * ((given_x - x_grid[i]) / (x_grid[i + 1] - x_grid[i]));
				cont = 0;
			}
		}
		i++;
	}
	return interp_y;

}
