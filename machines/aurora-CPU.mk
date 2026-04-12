# Machine-specific build settings for prj.
# Converted from machines_tmp/alcf/aurora-CPU.mk.

CC := mpicc
MACHINE_CFLAGS := -finline-functions -DHAVE_UNISTD_H
HDF5_LIBS := -L/soft/packaging/spack/e4s/23.05-2023.05.15.006/install/linux-sles15-x86_64/oneapi-2023.05.15.006.001/hdf5-1.14.2-mz7rnsojl74risazaqydz4vfyrmwelo5/lib -lhdf5
OMPFLAGS := -fp-model precise -mlong-double-64 -qmkl -fiopenmp -fopenmp-targets=spir64_gen="-mlong-double-64" -Xopenmp-target-backend=spir64_gen "-device 12.60.7" -D__STRICT_ANSI__
