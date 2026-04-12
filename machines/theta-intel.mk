# Machine-specific build settings for prj.
# Converted from machines_tmp/alcf/theta-intel.mk.

CC := cc
MACHINE_CFLAGS := -xMIC-AVX512 -diag-disable 161,3180,6843
MACHINE_LDLIBS := -L/usr/common/software/gsl/2.1/intel/lib -lgsl -lgslcblas
OMPFLAGS := -qopenmp
