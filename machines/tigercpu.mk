# Machine-specific build settings for prj.
# Converted from machines_tmp/princeton/tigercpu.mk.

CC := h5pcc
MACHINE_CFLAGS := -xCORE-AVX512 -qopt-zmm-usage=high -debug -diag-disable 161,3180,6843
MACHINE_LDLIBS := -lgsl -mkl -lifcore
OMPFLAGS := -qopenmp
