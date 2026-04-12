# Machine-specific build settings for prj.
# Converted from machines_tmp/princeton/stellar-intel.mk.

CC := h5pcc
MACHINE_CFLAGS := -xCORE-AVX512 -diag-disable=161,3180
MACHINE_LDLIBS := -lgsl -mkl -lifcore -lmpi -ldl -lrt -lpthread
