# Machine-specific build settings for prj.
# Converted from machines_tmp/princeton/perseus-gnu.mk.

CC := h5pcc
MACHINE_CFLAGS := -ftree-vectorize -fpeel-loops
HDF5_LIBS := -lhdf5
MACHINE_LDLIBS := -lgsl -lgslcblas -lblas -llapack -lmpi_usempi -lmpi_mpifh -lmpi
