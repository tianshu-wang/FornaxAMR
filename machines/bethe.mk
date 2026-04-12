# Machine-specific build settings for prj.
# Converted from machines_tmp/old/bethe.mk.

CC := h5pcc
MACHINE_CFLAGS := -ftree-vectorize -fpeel-loops
HDF5_LIBS := -lhdf5
MACHINE_LDLIBS := -llapack -lblas -lgsl -lgslcblas
