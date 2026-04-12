# Machine-specific build settings for prj.
# Converted from machines_tmp/old/hubert.mk.

CC := h5pcc
MACHINE_CFLAGS := -ftree-vectorize -ffast-math -fno-math-errno -fno-signed-zeros -fno-trapping-math -Wall
HDF5_LIBS := -lhdf5
MACHINE_LDLIBS := -llapack -lblas -lmpi_mpifh
