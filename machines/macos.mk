# macOS-specific compiler and library settings for prj.
#
# Usage examples:
#   make -f Makefile -include machines/macos.mk
#   make CC=mpicc HDF5_CFLAGS=... HDF5_LIBS=...

CC := $(shell command -v mpicc >/dev/null 2>&1 && echo mpicc || echo cc)

MPI_CFLAGS := $(shell $(CC) --showme:compile 2>/dev/null)
MPI_LIBS := $(shell $(CC) --showme:link 2>/dev/null)

HOMEBREW_HDF5_MPI_PREFIX := $(shell brew --prefix hdf5-mpi 2>/dev/null)
HOMEBREW_HDF5_PREFIX := $(shell brew --prefix hdf5 2>/dev/null)

ifeq ($(strip $(HOMEBREW_HDF5_MPI_PREFIX)),)
HDF5_PREFIX := $(HOMEBREW_HDF5_PREFIX)
else
HDF5_PREFIX := $(HOMEBREW_HDF5_MPI_PREFIX)
endif

ifeq ($(strip $(HDF5_PREFIX)),)
HDF5_CFLAGS ?= $(shell pkg-config --cflags hdf5 2>/dev/null)
HDF5_LIBS ?= $(shell pkg-config --libs hdf5 2>/dev/null)
else
HDF5_CFLAGS ?= -I$(HDF5_PREFIX)/include
HDF5_LIBS ?= -L$(HDF5_PREFIX)/lib -lhdf5
endif
