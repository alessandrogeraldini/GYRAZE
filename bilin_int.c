#include <stdlib.h>
#include <stdio.h>
#include "Headerfiles/bilin_int.h"

double bilin_interp(double x, double y, double **FF, double *xx, double *yy, int cols, int rows, int guessi, int guessj)
{
	double answer, deltax, deltay, aminus, aplus, bminus, bplus;
	int i, j, ileft = 0, iright = cols-1, jdown = 0, jup = rows-1, out_of_bounds =0, A=0, B=0, C=0, D=0;
	//printf("Am I here?/n");
	if (guessi == -1)
	{
		guessi = (int) cols/2;
		//printf("%d is guessi\n", guessi);
	}
	i = guessi;
	while ( ( ( xx[i] - x > 0 ) || ( xx[i+1] - x < 0 ) ) && out_of_bounds != 1 )
	{
		A = (xx[i] - x > 0);
		B = (xx[i+1] - x < 0);	
		if (A==1)
		{
			if (i-ileft > 1)
			{
			iright = i;
			i = i - (i - ileft)/2;
			}
			else
			{
			iright = i;
			i = i - 1;
			}
		}
		else if (B==1)
		{
			if (iright - i > 1)
			{
				ileft = i;
				i = i + (iright - i)/2;
			}
			else
			{
				ileft = i;
				i = i +1;
			}
		}
		//printf("i is %d, A and B are %d and %d\n", i, A, B);
		if (i == -1 || i == cols)
		{
			out_of_bounds = 1;
			//printf("Argument is out of bounds on the right; extrapolation necessary instead of interpolation\n");
		}
	}

	if (guessj == -1)
	{
		guessj = (int) rows/2;
	}
	j = guessj;
	while ( ( ( yy[j] - y > 0 ) ||  ( yy[j+1] - y < 0 ) ) && out_of_bounds != 1 )
	{	
		C = (yy[j] - y > 0);
		D = (yy[j+1] - y < 0);
		//printf("j is %d and C and D are %d and %d\n", j, C, D);
		if (C==1)
		{
			if (j-jdown > 1)
			{
			jup = j;
			j = j - (j - jdown)/2;
			}
			else
			{
			jup = i;
			j = j - 1;
			}
		}
		else if (D==1)
		{
			if (jup - j > 1)
			{
				jdown = i;
				j = j + (jup - j)/2;
			}
			else
			{
				jdown = j;
				j = j +1;
			}
		}

		if (j == -1 || j == rows)
		{
			out_of_bounds = 1;
			//printf("Argument is out of bounds on the right; extrapolation necessary instead of interpolation\n");
		}
	//printf("jup is %d\n", jup);
	}
	//printf("out of loop\n");

if (out_of_bounds !=1)
{
	deltax = xx[i+1] - xx[i];
	aplus= (x - xx[i])/deltax;
	aminus = (xx[i+1] - x )/deltax;
	deltay = yy[j+1] - yy[j];
	bplus = (y - yy[j])/deltay;
	bminus = (yy[j+1] - y)/deltay;
	answer = bplus*aplus*FF[i+1][j+1] + bplus*aminus*FF[i][j+1] + bminus*aplus*FF[i+1][j] + bminus*aminus*FF[i][j];
	//printf("answer is %f\n", answer);
}
else
{
	answer = 0.0;
}

return answer;
}
