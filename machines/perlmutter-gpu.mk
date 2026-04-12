# Machine-specific build settings for prj.
# Converted from machines_tmp/nersc/perlmutter-gpu.mk.

CC := cc
MACHINE_CFLAGS := --diag_suppress=177,550 --display_error_number -Mautoinline
HDF5_LIBS := -lhdf5
MACHINE_LDLIBS := -lnvf -lmpi -v
OMPFLAGS := -mp=gpu -gpu=cc80,lineinfo,maxregcount:96 -Minfo=mp,accel
