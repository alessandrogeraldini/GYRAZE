# set he compiler if CC and LD are unset
CC?=gcc
LD?=${CC}

CFLAGS+= -Wall -Wextra -Wpedantic

# to build for release:
# $ make -B PROFILE=release
# to build for debug:
# $ make -B PROFILE=debug
PROFILE?=

ifeq ($(PROFILE), release)
CFLAGS+= -O3 -march=native
endif

ifeq ($(PROFILE), debug)
CFLAGS+= -O1 -g -ggdb -fsanitize=address -DDEBUG -fno-omit-frame-pointer
endif

LDFLAGS+= $(CFLAGS) -lgsl -lgslcblas -lm

binaries=GYRAZE

GYRAZE: GYRAZE.o denscalc.o potupdate.o otherfuncs.o

.phony: clean

all: clean $(binaries)

clean:
	rm -f $(binaries) *.o
