# Platform config: Princeton Stellar cluster.
#
# Before building, load the required module:
#   module load gsl/2.6
#
# The gsl module sets CPATH and LIBRARY_PATH, so no explicit -I/-L flags are
# needed here. Use the system default gcc (8.5.0) rather than gcc-toolset/13:
# the toolset module is missing libasan, so PROFILE=debug fails to link
# under it.

CC := gcc
LD ?= $(CC)

LIBFLAGS = -lgsl -lgslcblas -lm

OMPFLAG ?= -fopenmp
