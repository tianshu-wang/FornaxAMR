# Machine-specific build settings for prj.
# Converted from machines_tmp/princeton/rhogar-intel.mk.

CC := /opt/intel/compilers_and_libraries_2020.1.217/linux/bin/intel64/icc
MACHINE_CFLAGS := $(LOCAL_INCLUDE) -xCORE-AVX512 -diag-disable=161,3180
HDF5_LIBS := -lhdf5
MACHINE_LDLIBS := $(LOCAL_LDFLAGS) -lgsl -mkl -lifcore -lmpi -ldl -lrt -lpthread
