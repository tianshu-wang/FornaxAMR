CC ?= $(shell command -v mpicc >/dev/null 2>&1 && echo mpicc || echo cc)
STD := -std=c99
WARN := -Wall -Wextra -pedantic

-include machine.mk
-include setup.mk

GRAVITY ?= 1
RADIATION ?= 0
MHD ?= 0
DUMP_SINGLE_PRECISION ?= 1
NRAD ?=
NEGROUP ?=
MPI_CFLAGS ?= $(shell $(CC) --showme:compile 2>/dev/null)
MPI_LIBS ?= $(shell $(CC) --showme:link 2>/dev/null)
HDF5_CFLAGS ?= $(shell pkg-config --cflags hdf5)
HDF5_LIBS ?= $(shell pkg-config --libs hdf5)
MACHINE_CPPFLAGS ?=
MACHINE_CFLAGS ?=
MACHINE_LDFLAGS ?=
MACHINE_LDLIBS ?=
OMPFLAGS ?=
NRAD_DEF := $(if $(NRAD),-DPRJ_NRAD=$(NRAD),)
NEGROUP_DEF := $(if $(NEGROUP),-DPRJ_NEGROUP=$(NEGROUP),)
CPPFLAGS := -Isrc -DPRJ_ENABLE_MPI -DPRJ_USE_GRAVITY=$(GRAVITY) -DPRJ_USE_RADIATION=$(RADIATION) -DPRJ_MHD=$(MHD) -DPRJ_DUMP_SINGLE_PRECISION=$(DUMP_SINGLE_PRECISION) $(NRAD_DEF) $(NEGROUP_DEF) $(MPI_CFLAGS) $(HDF5_CFLAGS) $(MACHINE_CPPFLAGS)
LDFLAGS := $(MACHINE_LDFLAGS)
LDLIBS := $(HDF5_LIBS) $(MPI_LIBS) $(MACHINE_LDLIBS)

ifeq ($(DEBUG),1)
CFLAGS := $(STD) $(WARN) -g -O0 -DPRJ_DEBUG $(MACHINE_CFLAGS) $(OMPFLAGS)
else
CFLAGS := $(STD) $(WARN) -O3 $(MACHINE_CFLAGS) $(OMPFLAGS)
endif

TARGET ?= prj
SRC_DIR := src
TEST_DIR := tests

SRCS := \
	$(SRC_DIR)/main.c \
	$(SRC_DIR)/prj_amr.c \
	$(SRC_DIR)/prj_boundary.c \
	$(SRC_DIR)/prj_eos.c \
	$(SRC_DIR)/prj_flux.c \
	$(SRC_DIR)/prj_gravity.c \
	$(SRC_DIR)/prj_io.c \
	$(SRC_DIR)/prj_mesh.c \
	$(SRC_DIR)/prj_mpi.c \
	$(SRC_DIR)/prj_radiation.c \
	$(SRC_DIR)/prj_rad3_opac.c \
	$(SRC_DIR)/prj_reconstruct.c \
	$(SRC_DIR)/prj_riemann.c \
	$(SRC_DIR)/prj_src.c \
	$(SRC_DIR)/prj_timeint.c \
	$(SRC_DIR)/prj_utils.c \
	problems/prj_problem_general.c \
	problems/prj_problem_cc.c \
	problems/prj_problem_ccsn.c \
	problems/prj_problem_sedov.c \
	problems/prj_problem_sedov_offcenter.c \
	problems/prj_problem_shock1d.c

OBJS := $(SRCS:.c=.o)
CORE_SRCS := $(filter-out $(SRC_DIR)/main.c,$(SRCS))
CORE_OBJS := $(CORE_SRCS:.c=.o)
TEST_SRCS := $(wildcard $(TEST_DIR)/*.c)
TEST_BINS := $(TEST_SRCS:.c=)

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LDLIBS) -o $@

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/prj.h $(SRC_DIR)/prj_defs.h $(SRC_DIR)/prj_types.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(TEST_DIR)/%: $(TEST_DIR)/%.c $(CORE_OBJS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(CORE_OBJS) $(LDFLAGS) $(LDLIBS) -o $@

test: $(TEST_BINS)
	@set -e; for test_bin in $(TEST_BINS); do ./$$test_bin; done

clean:
	rm -f $(TARGET) $(OBJS) $(CORE_OBJS) $(TEST_BINS)
