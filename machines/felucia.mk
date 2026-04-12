# Machine-specific build settings for prj.
# Converted from machines_tmp/old/felucia.mk.

CC := h5pcc
MACHINE_CFLAGS := -ftree-vectorize -fpeel-loops
HDF5_LIBS := -lhdf5
MACHINE_LDLIBS := -llapack -lblas -lgsl -lgslcblas
