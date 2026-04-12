# Machine-specific build settings for prj.
# Converted from machines_tmp/princeton/rhogar-gnu.mk.

CC := mpicc
MACHINE_CFLAGS := $(LOCAL_INCLUDE) -Wall -ftree-vectorize -fpeel-loops
HDF5_LIBS := -lhdf5
MACHINE_LDLIBS := $(LOCAL_LDFLAGS) -Wall -llapack -lblas -lgsl -lgslcblas -lmpi_mpifh -lmpi_usempif08
