# Machine-specific build settings for prj.
# Converted from machines_tmp/princeton/tiger-intel.mk.

CC := h5pcc
MACHINE_CFLAGS := -diag-disable 161,3180,6843
MACHINE_LDLIBS := -lgsl -mkl
