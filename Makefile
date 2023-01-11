MPS: MPS.o denscalc.o potupdate.o otherfuncs.o
	gcc -Wall -g -lgsl -lgslcblas -lm -o MPS MPS.o denscalc.o potupdate.o otherfuncs.o 

2d: MPS2d.o otherfuncs2d.o
	gcc -Wall -g -lgsl -lgslcblas -lm -o MPS2d MPS2d.o otherfuncs2d.o

MPS.o: MPS.c 
	gcc -Wall -c -o MPS.o MPS.c -g

otherfuncs.o: otherfuncs.c 
	gcc -Wall -g  -c -o otherfuncs.o otherfuncs.c

otherfuncs2d.o: otherfuncs2d.c 
	gcc -Wall -g  -c -o otherfuncs2d.o otherfuncs2d.c

denscalc.o: denscalc.c
	gcc -Wall -g -c -o denscalc.o denscalc.c

potupdate.o: potupdate.c
	gcc -Wall -g  -c -o potupdate.o potupdate.c

MPS2d.o: MPS2d.c 
	gcc -Wall -g  -c -o MPS2d.o MPS2d.c

clean:
	rm MPS MPS_renorm MPS2d MPS2d.o MPS.o MPS_renorm.o denscalc.o denscalc_renorm.o potupdate.o otherfuncs.o otherfuncs2d.o
