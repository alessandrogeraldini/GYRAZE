# set he compiler if CC and LD are unset
CC := gcc-15
LD?=${CC}

CFLAGS+= -Wall -Wextra -Wpedantic

# to build for release:
# $ make -B PROFILE=release
# to build for debug:
# $ make -B PROFILE=debug
PROFILE?=

ifeq ($(PROFILE), release)
CFLAGS+= -O3 -march=native -I/opt/homebrew/include
endif

ifeq ($(PROFILE), debug)
CFLAGS+= -O1 -g -ggdb -fsanitize=address -DDEBUG -fno-omit-frame-pointer -I/opt/homebrew/include
endif

LDFLAGS+= $(CFLAGS) $(OMPFLAG) -L/opt/homebrew/lib -lgsl -lgslcblas -lm

binaries=GYRAZE test_DS test_ion_DS

OMPFLAG?=-fopenmp

GYRAZE: GYRAZE.o denscalc.o densfinorb_par.o potupdate.o otherfuncs.o

test_DS: test_DS.o denscalc.o densfinorb_par.o otherfuncs.o

test_ion_DS: test_ion_DS.o denscalc.o densfinorb_par.o otherfuncs.o

densfinorb_par.o: densfinorb_par.c
	$(CC) $(CFLAGS) $(OMPFLAG) -c -o $@ $<

.phony: clean

all: clean $(binaries)

clean:
	rm -f $(binaries) *.o
