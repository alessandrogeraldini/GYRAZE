# Fallback platform config, used when SYSTEM isn't recognized.
#
# Tries pkg-config to locate GSL; override CC/CFLAGS/LIBFLAGS on the command
# line if this doesn't work on your system, e.g.:
#   make CFLAGS=-I/path/to/gsl/include LIBFLAGS='-L/path/to/gsl/lib -lgsl -lgslcblas -lm'
#
# To add permanent support for a new system, copy this file to
# make/<name>.mk, edit it, and add a detection rule in the top-level
# Makefile (or just build with `make SYSTEM=<name>`).

CC ?= gcc
LD ?= $(CC)

CFLAGS += $(shell pkg-config --cflags gsl 2>/dev/null)
LIBFLAGS = $(shell pkg-config --libs gsl 2>/dev/null || echo -lgsl -lgslcblas -lm)

OMPFLAG ?= -fopenmp
