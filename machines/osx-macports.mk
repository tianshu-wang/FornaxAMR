# Machine-specific build settings for prj.
# Converted from machines_tmp/old/osx-macports.mk.

CC := mpicc
MACHINE_CFLAGS := -I/opt/local/include -Wall -ftree-vectorize -fpeel-loops
HDF5_LIBS := -lhdf5
MACHINE_LDLIBS := -L/opt/local/lib -Wall -lgsl -lgslcblas -llapack -lblas -lmpi -lpmpi
