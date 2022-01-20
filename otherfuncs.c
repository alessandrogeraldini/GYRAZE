#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#define LEN 200


double lin_interp(double* x_grid, double* y_grid,double given_x ,int n, int line) {
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

double bilin_interp(double x, double y, double **FF, double *xx, double *yy, int cols, int rows, int guessi, int guessj) {
	double answer, deltax, deltay, aminus, aplus, bminus, bplus;
	int i, j, ileft = 0, iright = cols-1, jdown = 0, jup = rows-1, out_of_bounds =0, A=0, B=0, C=0, D=0;
	//for (i=0; i<cols; i++) printf("xx[%d] = %f\n", i, xx[i]);
	//for (i=0; i<rows; i++) printf("yy[%d] = %f\n", i, yy[i]);
	if (guessi == -1)
	{
		guessi = (int) cols/2;
		//printf("%d is guessi\n", guessi);
	}
	i = guessi;
	
	while ( ( ( xx[i] - x > 0 ) || ( xx[i+1] - x < 0 ) ) && out_of_bounds != 1 )
	{
		//printf("x = %f\nxx[i] = %f\nxx[i+1] = %f\n", x, xx[i], xx[i+1]);
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
			//printf("Argument is out of bounds horizontally; extrapolation necessary instead of interpolation\n");
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
			//printf("Argument is out of bounds; extrapolation necessary instead of interpolation\n");
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
	//printf("OUT OF BOUNDS\n");
}

return answer;
}


double *linetodata(char line[], int lenline, int *size) {
// takes a line from a data file and converts it to an array of numbers (in practice, to a double pointer)
// There's probably a better way to do it, but this works
	char testspace[2], numb[25];
	int i, lennumb, start_index, numb_index, *space_index, words=0, index;
	space_index = malloc(lenline * sizeof(int));
	double *numb_array, numb_double=0.0;
	/* Iterate over number of characters in the whole line
	 */
	for (i=0; i<lenline; i++)
	{
		testspace[0] = line[i];
		/* Compare each character to a white space. When a white space is found (or a 
		 * newline character), we break the string and extract the number from the 
		 * previous whitepace to the new one.
		 */
		if ( (strncmp(testspace, "\n", 1) == 0) || (strncmp(testspace, " ", 1) == 0) ) 
		{
			space_index[words] = i;
			words += 1;
		}
	}
	*size = words;
	numb_array = malloc(words*sizeof(double));
	for (index = 0;index < words; index++)
	{
		if (index == 0)
		{
			lennumb = space_index[index];
			start_index = 0;
		}
		else
		{
			lennumb = space_index[index] - space_index[index-1];
			start_index = space_index[index-1];
		}

		for (numb_index=0; numb_index < lennumb; numb_index++)
		{
			numb[numb_index] = line[start_index+numb_index];
			//printf("numb[%d] = %c\n", numb_index, numb[numb_index]);
		}
		for (numb_index=lennumb; numb_index < 25; numb_index++)
		{
			numb[numb_index] = ' ';
			//printf("numb[%d] = %c\n", numb_index, numb[numb_index]);
		}
		numb[25-1] = '\0';
		//printf("blblblablablablablabla %s\n", numb);
		//printf("%s\n", numb);
		sscanf(numb, "%lf", &numb_double);
		numb_array[index] = numb_double;
	}
	free(space_index);
	return numb_array;
}




double *linetodatanew(char *line, int *size) {
	int i, imin, lenline, words, started_word, save=0;
	char word[LEN];
	double *line_broken;
	/* Iterate over number of characters in the whole line
	 */
	if (*size<1)
	line_broken = malloc(LEN*sizeof(double));
	else 
	line_broken = malloc((*size)*sizeof(double));
	//printf("START of linetodata\n");
	lenline = strlen(line);
	imin = 0;
	words = 0;
	started_word = 0;
	for (i=0; i<lenline; i++)
	{ 
/* Compare each character to a white space (32) or newline (10) character. When a white space is found (or a 
 newline character), we break the string and extract the number from the 
 previous whitepace to the new one.
 */
		if ( (line[i] != 32) && (line[i] != 10) ) 
		{
			started_word = 1;
			//words += 1;
			if ( (isdigit(line[i]) != 0) && (words == 0) ) save = 1;
			if ( strncmp(line, "param", 5) == 0 ) save = 1;
		}
		else 
		{
			if (started_word == 1)
			{
				started_word = 0;
				words += 1;
			}
			imin = i+1;
		}
	}
	if (save == 1) { // CHANGE HERE
		//printf("saving line\n");
		imin = 0;
		words = 0;
		started_word = 0;
		for (i=0; i<lenline; i++)
		{
			/* Compare each character to a white space. When a white space is found (or a 
			 * newline character), we break the string and extract the number from the 
			 * previous whitepace to the new one.
			 */
			if ( (line[i] != 32) && (line[i] != 10) ) 
			{
				word[i-imin] = line[i];
				if (isdigit(word[i-imin]) != 0) {
					started_word = 1;
				}
			}
			else
			{
				if (started_word == 1)
				{	
					word[i-imin] = '\0';
					started_word = 0;
					sscanf(word, "%lf", (line_broken+words));
					words += 1;
					//line_broken += 1;
				}
				imin = i + 1;
			}
		}
	}
	*size = words;

	//printf("*line_broken   = %Le\n", *(line_broken)) ;
	//printf("*line_broken+1 = %Le\n", *(line_broken+1));
	//printf("*line_broken+2 = %Le\n", *(line_broken+2));
	//printf("*line_broken+3 = %Le\n", *(line_broken+3));
	//printf("%s\n", line);

	return line_broken;
}

