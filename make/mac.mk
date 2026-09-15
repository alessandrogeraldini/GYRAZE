# Platform config: macOS with Homebrew.
#
# Assumes GSL and an OpenMP-capable gcc are installed via Homebrew, e.g.:
#   brew install gsl gcc

CC := gcc-15
LD ?= $(CC)

HOMEBREW_PREFIX ?= /opt/homebrew

CFLAGS += -I$(HOMEBREW_PREFIX)/include
LIBFLAGS = -L$(HOMEBREW_PREFIX)/lib -lgsl -lgslcblas -lm

OMPFLAG ?= -fopenmp
