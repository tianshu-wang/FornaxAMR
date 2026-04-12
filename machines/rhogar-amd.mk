# Machine-specific build settings for prj.
# Converted from machines_tmp/princeton/rhogar-amd.mk.

CC := mpicc
MACHINE_CFLAGS := $(LOCAL_INCLUDE) -fvectorize -floop-splitting
HDF5_LIBS := -lhdf5
MACHINE_LDLIBS := $(LOCAL_LDFLAGS) -llapack -lblas -lgsl -lgslcblas -lmpi_mpifh -lmpi_usempif08
