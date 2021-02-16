// The definition of a top hat function, returning 1.0 if the argument x of the function is in between the parameters x1 and x2
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Headerfiles/tophat.h"

double tophat(double x1, double x2, double x)
{
	double y;
	if ((x1 <= x ) && (x <= x2 ))
	{
		y = 1.0;
	}	
	else
	{
		y = 0.0;
	}
	return y;
}
