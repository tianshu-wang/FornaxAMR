# Machine-specific build settings for prj.
# Converted from machines_tmp/tacc-xsede/frontera.mk.

CC := mpicc
# -xCORE-AVX512 already vectorizes sqrt/division via SVML, and Intel's default
# -fp-model fast=1 is errno-free (the Clang -fno-math-errno equivalent), so the
# sqrt-bound M1 radiation flux kernel is already optimized here. -fp-model fast
# is stated explicitly to lock in that intent. Do NOT use -fp-model fast=2 or
# -fast: the radiation closure relies on NaN/Inf and E<=0 guards.
MACHINE_CFLAGS := -xCORE-AVX512 -fp-model fast -diag-disable 161,3180,6843 -I${TACC_GSL_INC}
HDF5_CFLAGS := -I${TACC_HDF5_INC}
HDF5_LIBS := -L${TACC_HDF5_LIB} -lhdf5
MACHINE_LDLIBS := -L${TACC_GSL_LIB} -lgsl -lifcore -mkl
OMPFLAGS := -qopenmp
