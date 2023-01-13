GIMS: GIMS.o denscalc.o potupdate.o otherfuncs.o
	gcc -Wall -g -lgsl -lgslcblas -lm -o GIMS GIMS.o denscalc.o potupdate.o otherfuncs.o 

2d: GIMS2d.o otherfuncs2d.o
	gcc -Wall -g -lgsl -lgslcblas -lm -o GIMS2d GIMS2d.o otherfuncs2d.o

GIMS.o: GIMS.c 
	gcc -Wall -c -o GIMS.o GIMS.c -g

otherfuncs.o: otherfuncs.c 
	gcc -Wall -g  -c -o otherfuncs.o otherfuncs.c

otherfuncs2d.o: otherfuncs2d.c 
	gcc -Wall -g  -c -o otherfuncs2d.o otherfuncs2d.c

denscalc.o: denscalc.c
	gcc -Wall -g -c -o denscalc.o denscalc.c

potupdate.o: potupdate.c
	gcc -Wall -g  -c -o potupdate.o potupdate.c

GIMS2d.o: GIMS2d.c 
	gcc -Wall -g  -c -o GIMS2d.o GIMS2d.c

clean:
	rm GIMS GIMS2d GIMS2d.o GIMS.o denscalc.o potupdate.o otherfuncs.o otherfuncs2d.o
