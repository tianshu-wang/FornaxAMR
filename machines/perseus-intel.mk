# Machine-specific build settings for prj.
# Converted from machines_tmp/princeton/perseus-intel.mk.

CC := h5pcc
MACHINE_CFLAGS := -xCORE-AVX2
MACHINE_LDLIBS := -lgsl -mkl -lifcore -lmpi -lmpigi -ldl -lrt -lpthread
