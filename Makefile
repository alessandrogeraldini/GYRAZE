MPS: MPS.o denscalc.o potupdate.o otherfuncs.o
	gcc -Wall -lgsl -lgslcblas -lm -o MPS MPS.o denscalc.o potupdate.o otherfuncs.o 

MPS.o: MPS.c 
	gcc -Wall -I/usr/include/gsl -c -o MPS.o MPS.c

otherfuncs.o: otherfuncs.c 
	gcc -Wall -I/usr/include/gsl -c -o otherfuncs.o otherfuncs.c

denscalc.o: denscalc.c
	gcc -Wall -I/usr/include/gsl -c -o denscalc.o denscalc.c

potupdate.o: potupdate.c
	gcc -Wall -I/usr/include/gsl -c -o potupdate.o potupdate.c
