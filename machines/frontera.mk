# Machine-specific build settings for prj.
# Converted from machines_tmp/tacc-xsede/frontera.mk.

CC := mpicc
MACHINE_CFLAGS := -xCORE-AVX512 -diag-disable 161,3180,6843 -I${TACC_GSL_INC}
HDF5_CFLAGS := -I${TACC_HDF5_INC}
HDF5_LIBS := -L${TACC_HDF5_LIB} -lhdf5
MACHINE_LDLIBS := -L${TACC_GSL_LIB} -lgsl -lifcore -mkl
OMPFLAGS := -qopenmp
