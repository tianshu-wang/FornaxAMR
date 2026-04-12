# Machine-specific build settings for prj.
# Converted from machines_tmp/alcf/ALCF_ThetaGPU-NVIDIA_HPC_SDK.mk.

CC := mpicc
MACHINE_CFLAGS := --diag_suppress=177,550 --display_error_number
HDF5_LIBS := -L/lus/theta-fs0/software/thetagpu/hdf5/1.12.0/lib -lhdf5
MACHINE_LDLIBS := -lnvf -lmpi_mpifh -v
OMPFLAGS := -mp=gpu -gpu=cc80,lineinfo -Minfo=mp,accel
