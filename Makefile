# Platform-specific settings (CC, include/lib paths, ...) live in make/<name>.mk.
# The right one is auto-detected below; override explicitly with:
#   $ make SYSTEM=mac
#   $ make SYSTEM=stellar
#
# To add support for another system: copy make/generic.mk to make/<name>.mk,
# edit it, and either build with SYSTEM=<name> or add a detection rule below.

ifeq ($(origin SYSTEM), undefined)
  UNAME_S := $(shell uname -s)
  HOSTNAME := $(shell hostname -f 2>/dev/null || hostname)
  ifneq (,$(findstring stellar,$(HOSTNAME)))
    SYSTEM := stellar
  else ifeq ($(UNAME_S),Darwin)
    SYSTEM := mac
  else
    SYSTEM := generic
  endif
endif

$(info Building for SYSTEM=$(SYSTEM) (override with `make SYSTEM=<name>`; see make/*.mk))

include make/$(SYSTEM).mk

CFLAGS += -Wall -Wextra -Wpedantic

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

LDFLAGS+= $(CFLAGS) $(OMPFLAG) $(LIBFLAGS)

binaries=GYRAZE test_DS test_ion_DS test_jac_DS

GYRAZE: GYRAZE.o denscalc.o densfinorb_par.o potupdate.o otherfuncs.o

test_DS: test_DS.o denscalc.o densfinorb_par.o otherfuncs.o

test_ion_DS: test_ion_DS.o denscalc.o densfinorb_par.o otherfuncs.o

test_jac_DS: test_jac_DS.o denscalc.o densfinorb_par.o otherfuncs.o

# mps.h holds the compile-time switches; without this every object silently keeps the
# values it was built with when only the header changes.
GYRAZE.o denscalc.o densfinorb_par.o potupdate.o otherfuncs.o test_DS.o test_ion_DS.o test_jac_DS.o: mps.h

densfinorb_par.o: densfinorb_par.c mps.h
	$(CC) $(CFLAGS) $(OMPFLAG) -c -o $@ $<

.PHONY: clean

all: clean $(binaries)

clean:
	rm -f $(binaries) *.o
