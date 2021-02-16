MPS_it: MPS_iterate.o linetodata.o tophat.o iondens.o newguess.o cspline.o makelookup.o cutlookup.o Fgenerator.o bilin_int.o
	gcc -O2 -Wall -ffast-math -o MPS MPS_iterate.o linetodata.o tophat.o iondens.o newguess.o cspline.o makelookup.o cutlookup.o Fgenerator.o bilin_int.o
 
MPS_iterate.o: MPS_iterate.c 
	gcc -Wall -I/opt/local/include/gsl -c -o MPS_iterate.o MPS_iterate.c

bilin_int.o: bilin_int.c 
	gcc -Wall -I/opt/local/include/gsl -c -o bilin_int.o bilin_int.c

cspline.o: cspline.c
	gcc -Wall -I/opt/local/include/gsl -c -o cspline.o cspline.c

Fgenerator.o: Fgenerator.c
	gcc -Wall -I/opt/local/include/gsl -c -o Fgenerator.o Fgenerator.c

makelookup.o: makelookup.c
	gcc -Wall -I/opt/local/include/gsl -c -o makelookup.o makelookup.c

cutlookup.o: cutlookup.c
	gcc -Wall -I/opt/local/include/gsl -c -o cutlookup.o cutlookup.c

iondens.o: iondens.c
	gcc -Wall -I/opt/local/include/gsl -c -o iondens.o iondens.c

newguess.o: newguess.c
	gcc -Wall -I/opt/local/include/gsl -c -o newguess.o newguess.c

tophat.o: tophat.c
	gcc -Wall -I/opt/local/include/gsl -c -o tophat.o tophat.c

linetodata.o: linetodata.c
	gcc -Wall -I/opt/local/include/gsl -c -o linetodata.o linetodata.c
