# Machine-specific build settings for prj.
# Converted from machines_tmp/princeton/della-gpu.mk.

CC := mpicc
MACHINE_CFLAGS := --diag_suppress=177,550 --display_error_number
HDF5_LIBS := -L/usr/local/hdf5/nvhpc-22.5/openmpi-4.1.3/1.10.6/lib64 -lhdf5
MACHINE_LDLIBS := -lnvf -lmpi_mpifh -v
OMPFLAGS := -mp=gpu -gpu=cc80,lineinfo -Minfo=mp,accel
