
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mps.h"
#define LEN 200

double lin_interp(double *xx, double *yy, double x, int n, int bla) {
	double answer, deltax, aminus, aplus;
	int i, ileft = 0, iright = n-1, out_of_bounds =0, guessi;
	guessi = n/2;
	i = guessi;
	//printf("x = %f\nxx[0] = %f\nxx[n-1] = %f\n", x, xx[0], xx[n-1]);
	if ( (x < xx[0]) || (x > xx[n-1]) ) out_of_bounds = 1;
	else out_of_bounds = 0;
	//printf("out_of_bounds = %d\n", out_of_bounds);
	while ( ( ( xx[i] - x > TINY ) || ( xx[i+1] - x < - TINY ) ) && out_of_bounds == 0 ) {
		//printf("x = %f\nxx[i] = %f\nxx[i+1] = %f\n", x, xx[i], xx[i+1]);
		if (xx[i] - x > TINY) {
			if (i-ileft > 1) {
				iright = i;
				i = i - (i - ileft)/2;
			}
			else {
				iright = i;
				i = i - 1;
			}
		}
		else if (xx[i+1] - x < - TINY) {
			if (iright - i > 1) {
				ileft = i;
				i = i + (iright - i)/2;
			}
			else
			{
				ileft = i;
				i = i +1;
			}
		}
		if (i == -1 || i == n) {
			out_of_bounds = 1;
			printf("Argument is out of bounds; extrapolation necessary instead of interpolation\n\n\n");
			exit(-1);
		}
	}

	if (out_of_bounds == 0) {
		deltax = xx[i+1] - xx[i];
		aplus= (x - xx[i])/deltax;
		aminus = (xx[i+1] - x )/deltax;
		answer = aplus*yy[i+1] + aminus*yy[i] ;
		//printf("answer is %f\n", answer);
		//if (answer < 0.0) {
		//	printf("whymain????, FF++, FF0+ FF00 FF+0 = %f %f %f %f\n\n", FF[i+1][j+1], FF[i][j+1], FF[i][j], FF[i+1][j]);
		//	printf(" mu UminmU = %f %f\txi xi+1 yi yi+1 = %f %f %f %f\n\n", x, y, xx[i], xx[i+1], yy[j], yy[j+1]);
		//	printf(" deltax aplus aminus bplus bminus  = %f %f %f %f %f\n\n", deltax, aplus, aminus, bplus, bminus);
		//	exit(-1);
		//}
	}
	else {
		if ( (x < xx[0]) && (x > xx[0] - (xx[1] - xx[0])) ) {
			answer = yy[0] - (xx[0] - x)*(yy[1] - yy[0])/(xx[1]-xx[0]);
			printf("WARNING in lin_interp: extrapolation carried out on the right. Make sure this does not happen often\n");
		}
		else if ( (x > xx[n-1]) && (x < xx[n-1] + (xx[n-1] - xx[n-2])) ) {
			answer = yy[0] + (xx[n-1] - x)*(yy[n-1] - yy[n-2])/(xx[n-1]-xx[n-2]); 
			printf("WARNING in lin_interp: extrapolation carried out on the right. Make sure this does not happen often\n");
		}
		else {
			printf("x = %f\nxx[0] = %f\nxx[n-1] = %f\n", x, xx[0], xx[n-1]);
			printf("ERROR in lin_interp function, %d\n", bla);
			exit(-1);
		}
		//answer = 0.0;
		//printf("OUT OF BOUNDS\n");
	}
	return answer;
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
	if ( (x < xx[0]) || (x > xx[cols-1]) ) out_of_bounds = 1;
	else out_of_bounds = 0;
	while ( ( ( xx[i] - x > 0 ) || ( xx[i+1] - x < 0 ) ) && out_of_bounds == 0 ) {
		//printf("x = %f\nxx[i] = %f\nxx[i+1] = %f\n", x, xx[i], xx[i+1]);
		A = (xx[i] - x > 0);
		B = (xx[i+1] - x < 0);	
		if (A==1)
		{
			if (i-ileft > 1) {
				iright = i;
				i = i - (i - ileft)/2;
			}
			else {
				iright = i;
				i = i - 1;
			}
		}
		else if (B==1)
		{
			if (iright - i > 1) {
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
			printf("Argument is out of bounds; extrapolation necessary instead of interpolation\n\n\n");
		}
	}

	if (guessj == -1) {
		guessj = (int) rows/2;
	}
	j = guessj;
	if ( (y < yy[0]) || (y > yy[rows-1]) ) out_of_bounds = 1;
	while ( ( ( yy[j] - y > 0 ) ||  ( yy[j+1] - y < 0 ) ) && out_of_bounds == 0 )
	{	
		C = (yy[j] - y > 0);
		D = (yy[j+1] - y < 0);
		//printf("j is %d and C and D are %d and %d\n", j, C, D);
		if (C==1) {
			if (j-jdown > 1) {
				jup = j;
				j = j - (j - jdown)/2;
			}
			else {
				jup = j;
				j = j - 1;
			}
		}
		else if (D==1) {
			if (jup - j > 1) {
				jdown = j;
				j = j + (jup - j)/2;
			}
			else {
				jdown = j;
				j = j +1;
			}
		}

		if (j == -1 || j == rows)
		{
			out_of_bounds = 1;
			printf("Argument is out of bounds; extrapolation necessary instead of interpolation\n\n\n");
		}
	//printf("jup is %d\n", jup);
	}
	//printf("out of loop\n");

	if (out_of_bounds == 0)
	{
		deltax = xx[i+1] - xx[i];
		aplus= (x - xx[i])/deltax;
		aminus = (xx[i+1] - x )/deltax;
		deltay = yy[j+1] - yy[j];
		bplus = (y - yy[j])/deltay;
		bminus = (yy[j+1] - y)/deltay;
		answer = bplus*aplus*FF[i+1][j+1] + bplus*aminus*FF[i][j+1] + bminus*aplus*FF[i+1][j] + bminus*aminus*FF[i][j];
		//printf("answer is %f\n", answer);
		if (answer < 0.0) {
			printf("distribution function value = %f\n", answer);
			printf("why is the interpolated distribution function so negative??, FF++, FF0+ FF00 FF+0 = %f %f %f %f\n\n", FF[i+1][j+1], FF[i][j+1], FF[i][j], FF[i+1][j]);
			printf(" mu UminmU = %f %f\txi xi+1 yi yi+1 = %f %f %f %f\n\n", x, y, xx[i], xx[i+1], yy[j], yy[j+1]);
			printf(" deltax aplus aminus bplus bminus  = %f %f %f %f %f\n\n", deltax, aplus, aminus, bplus, bminus);
			// this if loop is necessary if the input distribution function has small negative values due to numerical error in code simulating rest of plasma
			//if (fabs(answer) > 1e-10) 
			//	exit(-1);
			//else 
			//	answer = 0.0;

		}
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
	int i, lennumb, start_index, numb_index, *space_index, words=0, index, ignore_flag=0;
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
		if ( strncmp(testspace, "#", 1) == 0) {
			ignore_flag = 1;
			if (ignore_flag == 1) printf("YOLO\n");
		}
		if ( (strncmp(testspace, "\n", 1) == 0) || (strncmp(testspace, " ", 1) == 0) ) 
		{
			if (ignore_flag == 0) {
				space_index[words] = i;
				words += 1;
			}
			else ignore_flag = 0;
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

