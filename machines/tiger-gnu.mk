# Machine-specific build settings for prj.
# Converted from machines_tmp/princeton/tiger-gnu.mk.

CC := h5pcc
MACHINE_CFLAGS := -ftree-vectorize -fpeel-loops
HDF5_LIBS := -lhdf5
MACHINE_LDLIBS := -L/home/mskinner/TOOLS/gsl-2.1 -lgsl -lgslcblas -lblas -llapack
