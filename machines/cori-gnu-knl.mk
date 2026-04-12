# Machine-specific build settings for prj.
# Converted from machines_tmp/nersc/cori-gnu-knl.mk.

CC := cc
MACHINE_CFLAGS := -ffast-math -ftree-vectorize -funroll-loops -ftree-vectorizer-verbose=2 -march=knl
OMPFLAGS := -fopenmp
