#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>

#include <hdf5.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#define PRJ_IO_METADATA_SIZE 639

static void prj_io_fail(const char *message);

static void prj_io_trim(char *text)
{
    char *start;
    char *end;

    if (text == 0) {
        return;
    }
    start = text;
    while (*start != '\0' && isspace((unsigned char)*start)) {
        start += 1;
    }
    if (start != text) {
        memmove(text, start, strlen(start) + 1);
    }
    end = text + strlen(text);
    while (end > text && isspace((unsigned char)end[-1])) {
        end -= 1;
    }
    *end = '\0';
}

static int prj_io_parse_bc(const char *value, int *bc_type)
{
    if (value == 0 || bc_type == 0) {
        return 1;
    }
    if (strcmp(value, "outflow") == 0) {
        *bc_type = PRJ_BC_OUTFLOW;
        return 0;
    }
    if (strcmp(value, "reflect") == 0) {
        *bc_type = PRJ_BC_REFLECT;
        return 0;
    }
    if (strcmp(value, "user") == 0) {
        *bc_type = PRJ_BC_USER;
        return 0;
    }
    return 1;
}

static int prj_io_parse_eos_kind(const char *value, int *eos_kind)
{
    if (value == 0 || eos_kind == 0) {
        return 1;
    }
    if (strcmp(value, "ideal") == 0) {
        *eos_kind = PRJ_EOS_KIND_IDEAL;
        return 0;
    }
    if (strcmp(value, "table") == 0 || strcmp(value, "tabulated") == 0) {
        *eos_kind = PRJ_EOS_KIND_TABLE;
        return 0;
    }
    return 1;
}

static int prj_io_parse_amr_estimator(const char *value, int *amr_estimator)
{
    if (value == 0 || amr_estimator == 0) {
        return 1;
    }
    if (strcmp(value, "lohner") == 0) {
        *amr_estimator = PRJ_AMR_ESTIMATOR_LOEHNER;
        return 0;
    }
    if (strcmp(value, "velocity") == 0) {
        *amr_estimator = PRJ_AMR_ESTIMATOR_VELOCITY;
        return 0;
    }
    if (strcmp(value, "pressure_scale_height") == 0) {
        *amr_estimator = PRJ_AMR_ESTIMATOR_PRESSURE_SCALE_HEIGHT;
        return 0;
    }
    if (strcmp(value, "density_jump") == 0) {
        *amr_estimator = PRJ_AMR_ESTIMATOR_DENSITY_JUMP;
        return 0;
    }
    return 1;
}

static int prj_io_parse_amr_slot_key(const char *key, const char *prefix, int *slot)
{
    size_t prefix_len;
    char *endptr;
    long idx;

    if (key == 0 || prefix == 0 || slot == 0) {
        return 0;
    }
    prefix_len = strlen(prefix);
    if (strncmp(key, prefix, prefix_len) != 0) {
        return 0;
    }
    if (key[prefix_len] == '\0') {
        *slot = 0;
        return 1;
    }
    idx = strtol(key + prefix_len, &endptr, 10);
    if (*endptr != '\0' || idx < 1 || idx > PRJ_AMR_N) {
        return 0;
    }
    *slot = (int)idx - 1;
    return 1;
}

static void prj_io_set_default_runtime(prj_sim *sim)
{
    if (sim == 0) {
        return;
    }

    sim->coord.x1min = -2.0;
    sim->coord.x1max = 2.0;
    sim->coord.x2min = -2.0;
    sim->coord.x2max = 2.0;
    sim->coord.x3min = -2.0;
    sim->coord.x3max = 2.0;
    sim->bc.bc_x1_inner = PRJ_BC_OUTFLOW;
    sim->bc.bc_x1_outer = PRJ_BC_OUTFLOW;
    sim->bc.bc_x2_inner = PRJ_BC_OUTFLOW;
    sim->bc.bc_x2_outer = PRJ_BC_OUTFLOW;
    sim->bc.bc_x3_inner = PRJ_BC_OUTFLOW;
    sim->bc.bc_x3_outer = PRJ_BC_OUTFLOW;
    sim->cfl = 0.8;
    sim->dt_factor = 1.2;
    sim->t_end = 0.1;
    sim->output_dt = -1.0;
    sim->restart_dt = -1.0;
    sim->max_steps = 100;
    sim->output_interval = -1;
    sim->restart_interval = -1;
    sim->amr_interval = -1;
    strncpy(sim->output_dir, "output/dump", sizeof(sim->output_dir) - 1);
    sim->output_dir[sizeof(sim->output_dir) - 1] = '\0';
    sim->progenitor_file[0] = '\0';
    strncpy(sim->problem_name, "general", sizeof(sim->problem_name) - 1);
    sim->problem_name[sizeof(sim->problem_name) - 1] = '\0';
    sim->restart_from_file = 0;
    sim->restart_file_name[0] = '\0';
    sim->mesh.root_nx[0] = 8;
    sim->mesh.root_nx[1] = 8;
    sim->mesh.root_nx[2] = 8;
    sim->mesh.max_level = 0;
    {
        int amr_idx;

        for (amr_idx = 0; amr_idx < PRJ_AMR_N; ++amr_idx) {
            sim->mesh.amr_refine_thresh[amr_idx] = 0.5;
            sim->mesh.amr_derefine_thresh[amr_idx] = 0.2;
            sim->mesh.amr_estimator[amr_idx] = PRJ_AMR_ESTIMATOR_LOEHNER;
            sim->mesh.amr_criterion_set[amr_idx] = 0;
        }
        sim->mesh.amr_estimator[0] = PRJ_AMR_ESTIMATOR_VELOCITY;
        sim->mesh.amr_criterion_set[0] = 1;
    }
    sim->mesh.amr_eps = 0.1;
    sim->mesh.use_amr_angle_resolution = 0;
    sim->mesh.amr_angle_resolution_limit = 0.0;
    sim->mesh.E_floor = -1.0;
    sim->eos.kind = PRJ_EOS_KIND_IDEAL;
    sim->eos.filename[0] = '\0';
    sim->rad.maxiter = 20;
    sim->rad.implicit_err_tol = 1.0e-6;
}

void prj_io_parser(prj_sim *sim, char *filename)
{
    FILE *fp;
    char line[1024];
    int lineno = 0;

    if (sim == 0) {
        prj_io_fail("prj_io_parser: sim is null");
    }

    prj_io_set_default_runtime(sim);
    if (filename == 0) {
        return;
    }

    fp = fopen(filename, "r");
    if (fp == 0) {
        prj_io_fail("prj_io_parser: failed to open parameter file");
    }

    while (fgets(line, sizeof(line), fp) != 0) {
        char *comment;
        char *eq;
        char *key;
        char *value;
        char *endptr;
        int amr_slot;

        lineno += 1;
        comment = strchr(line, '#');
        if (comment != 0) {
            *comment = '\0';
        }
        comment = strstr(line, "//");
        if (comment != 0) {
            *comment = '\0';
        }
        prj_io_trim(line);
        if (line[0] == '\0') {
            continue;
        }

        eq = strchr(line, '=');
        if (eq == 0) {
            fprintf(stderr, "prj_io_parser: invalid line %d in %s\n", lineno, filename);
            fclose(fp);
            exit(1);
        }
        *eq = '\0';
        key = line;
        value = eq + 1;
        prj_io_trim(key);
        prj_io_trim(value);

        if (strcmp(key, "x1min") == 0) {
            sim->coord.x1min = strtod(value, &endptr);
        } else if (strcmp(key, "x1max") == 0) {
            sim->coord.x1max = strtod(value, &endptr);
        } else if (strcmp(key, "x2min") == 0) {
            sim->coord.x2min = strtod(value, &endptr);
        } else if (strcmp(key, "x2max") == 0) {
            sim->coord.x2max = strtod(value, &endptr);
        } else if (strcmp(key, "x3min") == 0) {
            sim->coord.x3min = strtod(value, &endptr);
        } else if (strcmp(key, "x3max") == 0) {
            sim->coord.x3max = strtod(value, &endptr);
        } else if (strcmp(key, "cfl") == 0) {
            sim->cfl = strtod(value, &endptr);
        } else if (strcmp(key, "t_end") == 0) {
            sim->t_end = strtod(value, &endptr);
        } else if (strcmp(key, "dt_factor") == 0) {
            sim->dt_factor = strtod(value, &endptr);
        } else if (strcmp(key, "output_dt") == 0) {
            sim->output_dt = strtod(value, &endptr);
        } else if (strcmp(key, "restart_dt") == 0) {
            sim->restart_dt = strtod(value, &endptr);
        } else if (strcmp(key, "max_steps") == 0) {
            sim->max_steps = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "output_interval") == 0) {
            sim->output_interval = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "restart_interval") == 0) {
            sim->restart_interval = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "amr_interval") == 0) {
            sim->amr_interval = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "root_nx") == 0) {
            int root_n = (int)strtol(value, &endptr, 10);

            sim->mesh.root_nx[0] = root_n;
            sim->mesh.root_nx[1] = root_n;
            sim->mesh.root_nx[2] = root_n;
        } else if (strcmp(key, "root_nx1") == 0) {
            sim->mesh.root_nx[0] = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "root_nx2") == 0) {
            sim->mesh.root_nx[1] = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "root_nx3") == 0) {
            sim->mesh.root_nx[2] = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "max_level") == 0) {
            sim->mesh.max_level = (int)strtol(value, &endptr, 10);
        } else if (prj_io_parse_amr_slot_key(key, "amr_estimator", &amr_slot)) {
            if (prj_io_parse_amr_estimator(value, &sim->mesh.amr_estimator[amr_slot]) != 0) {
                endptr = value;
            } else {
                endptr = value + strlen(value);
                sim->mesh.amr_criterion_set[amr_slot] = 1;
            }
        } else if (prj_io_parse_amr_slot_key(key, "amr_refine_thresh", &amr_slot)) {
            sim->mesh.amr_refine_thresh[amr_slot] = strtod(value, &endptr);
            if (endptr != value) {
                sim->mesh.amr_criterion_set[amr_slot] = 1;
            }
        } else if (prj_io_parse_amr_slot_key(key, "amr_derefine_thresh", &amr_slot)) {
            sim->mesh.amr_derefine_thresh[amr_slot] = strtod(value, &endptr);
            if (endptr != value) {
                sim->mesh.amr_criterion_set[amr_slot] = 1;
            }
        } else if (strcmp(key, "use_amr_angle_resolution") == 0) {
            sim->mesh.use_amr_angle_resolution = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "amr_angle_resolution_limit") == 0) {
            sim->mesh.amr_angle_resolution_limit = strtod(value, &endptr);
        } else if (strcmp(key, "E_floor") == 0) {
            sim->mesh.E_floor = strtod(value, &endptr);
        } else if (strcmp(key, "output_dir") == 0) {
            strncpy(sim->output_dir, value, sizeof(sim->output_dir) - 1);
            sim->output_dir[sizeof(sim->output_dir) - 1] = '\0';
            endptr = value + strlen(value);
        } else if (strcmp(key, "progenitor_file") == 0) {
            strncpy(sim->progenitor_file, value, sizeof(sim->progenitor_file) - 1);
            sim->progenitor_file[sizeof(sim->progenitor_file) - 1] = '\0';
            endptr = value + strlen(value);
        } else if (strcmp(key, "problem") == 0) {
            strncpy(sim->problem_name, value, sizeof(sim->problem_name) - 1);
            sim->problem_name[sizeof(sim->problem_name) - 1] = '\0';
            endptr = value + strlen(value);
        } else if (strcmp(key, "restart_from_file") == 0) {
            sim->restart_from_file = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "restart_file_name") == 0) {
            strncpy(sim->restart_file_name, value, sizeof(sim->restart_file_name) - 1);
            sim->restart_file_name[sizeof(sim->restart_file_name) - 1] = '\0';
            endptr = value + strlen(value);
        } else if (strcmp(key, "eos_file") == 0) {
            strncpy(sim->eos.filename, value, sizeof(sim->eos.filename) - 1);
            sim->eos.filename[sizeof(sim->eos.filename) - 1] = '\0';
            sim->eos.kind = PRJ_EOS_KIND_TABLE;
            endptr = value + strlen(value);
        } else if (strcmp(key, "rad_maxiter") == 0) {
            sim->rad.maxiter = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "rad_implicit_err_tol") == 0) {
            sim->rad.implicit_err_tol = strtod(value, &endptr);
#if PRJ_NRAD > 0
        } else if (strcmp(key, "rad_table_param_file") == 0) {
            strncpy(sim->rad.table_param_file, value, sizeof(sim->rad.table_param_file) - 1);
            sim->rad.table_param_file[sizeof(sim->rad.table_param_file) - 1] = '\0';
            endptr = value + strlen(value);
        } else if (strcmp(key, "rad_table_file") == 0) {
            strncpy(sim->rad.table_file, value, sizeof(sim->rad.table_file) - 1);
            sim->rad.table_file[sizeof(sim->rad.table_file) - 1] = '\0';
            endptr = value + strlen(value);
        } else if (strncmp(key, "rad_emin_", 9) == 0) {
            char *kend;
            long idx = strtol(key + 9, &kend, 10);

            if (*kend != '\0' || idx < 0 || idx >= PRJ_NRAD) {
                fprintf(stderr, "prj_io_parser: bad rad_emin index '%s' in %s:%d\n", key, filename, lineno);
                fclose(fp);
                exit(1);
            }
            sim->rad.emin[idx] = strtod(value, &endptr);
        } else if (strncmp(key, "rad_emax_", 9) == 0) {
            char *kend;
            long idx = strtol(key + 9, &kend, 10);

            if (*kend != '\0' || idx < 0 || idx >= PRJ_NRAD) {
                fprintf(stderr, "prj_io_parser: bad rad_emax index '%s' in %s:%d\n", key, filename, lineno);
                fclose(fp);
                exit(1);
            }
            sim->rad.emax[idx] = strtod(value, &endptr);
#endif
        } else if (strcmp(key, "eos_type") == 0) {
            if (prj_io_parse_eos_kind(value, &sim->eos.kind) != 0) {
                endptr = value;
            } else {
                endptr = value + strlen(value);
            }
            endptr = value + strlen(value);
        } else if (strcmp(key, "bc_x1_inner") == 0) {
            if (prj_io_parse_bc(value, &sim->bc.bc_x1_inner) != 0) {
                endptr = value;
            } else {
                endptr = value + strlen(value);
            }
        } else if (strcmp(key, "bc_x1_outer") == 0) {
            if (prj_io_parse_bc(value, &sim->bc.bc_x1_outer) != 0) {
                endptr = value;
            } else {
                endptr = value + strlen(value);
            }
        } else if (strcmp(key, "bc_x2_inner") == 0) {
            if (prj_io_parse_bc(value, &sim->bc.bc_x2_inner) != 0) {
                endptr = value;
            } else {
                endptr = value + strlen(value);
            }
        } else if (strcmp(key, "bc_x2_outer") == 0) {
            if (prj_io_parse_bc(value, &sim->bc.bc_x2_outer) != 0) {
                endptr = value;
            } else {
                endptr = value + strlen(value);
            }
        } else if (strcmp(key, "bc_x3_inner") == 0) {
            if (prj_io_parse_bc(value, &sim->bc.bc_x3_inner) != 0) {
                endptr = value;
            } else {
                endptr = value + strlen(value);
            }
        } else if (strcmp(key, "bc_x3_outer") == 0) {
            if (prj_io_parse_bc(value, &sim->bc.bc_x3_outer) != 0) {
                endptr = value;
            } else {
                endptr = value + strlen(value);
            }
        } else {
            fprintf(stderr, "prj_io_parser: unknown key '%s' in %s:%d\n", key, filename, lineno);
            fclose(fp);
            exit(1);
        }

        prj_io_trim(endptr);
        if (*endptr != '\0') {
            fprintf(stderr, "prj_io_parser: invalid value for '%s' in %s:%d\n", key, filename, lineno);
            fclose(fp);
            exit(1);
        }
    }

    fclose(fp);
}

static void prj_io_dataset_name(int var, char *name, size_t size)
{
    if (var == PRJ_PRIM_RHO) {
        snprintf(name, size, "density");
    } else if (var == PRJ_PRIM_V1) {
        snprintf(name, size, "v1");
    } else if (var == PRJ_PRIM_V2) {
        snprintf(name, size, "v2");
    } else if (var == PRJ_PRIM_V3) {
        snprintf(name, size, "v3");
    } else if (var == PRJ_PRIM_EINT) {
        snprintf(name, size, "eint");
    } else if (var == PRJ_PRIM_YE) {
        snprintf(name, size, "ye");
    } else {
        int rad_idx = var - PRJ_NHYDRO;
        int field = rad_idx / (PRJ_NEGROUP * PRJ_RAD_GROUP_STRIDE);
        int local = rad_idx % (PRJ_NEGROUP * PRJ_RAD_GROUP_STRIDE);
        int group = local / PRJ_RAD_GROUP_STRIDE;
        int component = local % PRJ_RAD_GROUP_STRIDE;

        if (component == 0) {
            snprintf(name, size, "rad%d_g%02d_E", field, group);
        } else {
            snprintf(name, size, "rad%d_g%02d_F%d", field, group, component);
        }
    }
}

static void prj_io_fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}

static int prj_io_is_local_owner(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && (mpi == 0 || block->rank == mpi->rank);
}

static int prj_io_is_root_rank(void)
{
    prj_mpi *mpi = prj_mpi_current();

    return mpi == 0 || mpi->rank == 0;
}

static hid_t prj_io_create_file(const char *filename)
{
#if defined(PRJ_ENABLE_MPI)
    prj_mpi *mpi = prj_mpi_current();

    if (mpi != 0 && mpi->totrank > 1) {
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        hid_t file;

        H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
        file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        H5Pclose(fapl);
        return file;
    }
#endif
    return H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

static hid_t prj_io_open_file_readonly(const char *filename)
{
#if defined(PRJ_ENABLE_MPI)
    prj_mpi *mpi = prj_mpi_current();

    if (mpi != 0 && mpi->totrank > 1) {
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        hid_t file;

        H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
        file = H5Fopen(filename, H5F_ACC_RDONLY, fapl);
        H5Pclose(fapl);
        return file;
    }
#endif
    return H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
}

static hid_t prj_io_data_xfer_plist(void)
{
#if defined(PRJ_ENABLE_MPI)
    prj_mpi *mpi = prj_mpi_current();

    if (mpi != 0 && mpi->totrank > 1) {
        hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);

        H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
        return dxpl;
    }
#endif
    return H5P_DEFAULT;
}

static void prj_io_close_dxpl(hid_t dxpl)
{
    if (dxpl != H5P_DEFAULT) {
        H5Pclose(dxpl);
    }
}

static hid_t prj_io_dump_real_hdf5_type(void)
{
#if PRJ_DUMP_SINGLE_PRECISION
    return H5T_NATIVE_FLOAT;
#else
    return H5T_NATIVE_DOUBLE;
#endif
}

static void prj_io_write_attr_double(hid_t obj, const char *name, double value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(obj, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attr);
    H5Sclose(space);
}

static void prj_io_write_attr_int(hid_t obj, const char *name, int value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(obj, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr, H5T_NATIVE_INT, &value);
    H5Aclose(attr);
    H5Sclose(space);
}

static void prj_io_write_attr_int3(hid_t obj, const char *name, const int values[3])
{
    hsize_t dims[1] = {3};
    hid_t space = H5Screate_simple(1, dims, dims);
    hid_t attr = H5Acreate2(obj, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr, H5T_NATIVE_INT, values);
    H5Aclose(attr);
    H5Sclose(space);
}

static void prj_io_write_attr_double6(hid_t obj, const char *name, const prj_coord *coord)
{
    double values[6];
    hsize_t dims[1] = {6};
    hid_t space;
    hid_t attr;

    values[0] = coord->x1min;
    values[1] = coord->x1max;
    values[2] = coord->x2min;
    values[3] = coord->x2max;
    values[4] = coord->x3min;
    values[5] = coord->x3max;
    space = H5Screate_simple(1, dims, dims);
    attr = H5Acreate2(obj, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, values);
    H5Aclose(attr);
    H5Sclose(space);
}

static void prj_io_read_attr_double(hid_t obj, const char *name, double *value)
{
    hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);

    H5Aread(attr, H5T_NATIVE_DOUBLE, value);
    H5Aclose(attr);
}

static void prj_io_read_attr_int(hid_t obj, const char *name, int *value)
{
    hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);

    H5Aread(attr, H5T_NATIVE_INT, value);
    H5Aclose(attr);
}

static int prj_io_read_attr_int_optional(hid_t obj, const char *name, int default_value)
{
    int value = default_value;

    if (H5Aexists(obj, name) > 0) {
        prj_io_read_attr_int(obj, name, &value);
    }
    return value;
}

static double prj_io_read_attr_double_optional(hid_t obj, const char *name, double default_value)
{
    double value = default_value;

    if (H5Aexists(obj, name) > 0) {
        prj_io_read_attr_double(obj, name, &value);
    }
    return value;
}

static void prj_io_read_attr_int3(hid_t obj, const char *name, int values[3])
{
    hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);

    H5Aread(attr, H5T_NATIVE_INT, values);
    H5Aclose(attr);
}

static void prj_io_read_attr_double6(hid_t obj, const char *name, prj_coord *coord)
{
    double values[6];
    hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);

    H5Aread(attr, H5T_NATIVE_DOUBLE, values);
    H5Aclose(attr);
    coord->x1min = values[0];
    coord->x1max = values[1];
    coord->x2min = values[2];
    coord->x2max = values[3];
    coord->x3min = values[4];
    coord->x3max = values[5];
}

static void prj_io_fill_metadata(const prj_block *block, double *metadata_row)
{
    int n;
    int idx;

    for (n = 0; n < PRJ_IO_METADATA_SIZE; ++n) {
        metadata_row[n] = 0.0;
    }
    metadata_row[0] = (double)block->id;
    metadata_row[1] = 0.0;
    metadata_row[2] = (double)block->level;
    metadata_row[3] = (double)block->active;
    metadata_row[4] = block->xmin[0];
    metadata_row[5] = block->xmin[1];
    metadata_row[6] = block->xmin[2];
    metadata_row[7] = block->xmax[0];
    metadata_row[8] = block->xmax[1];
    metadata_row[9] = block->xmax[2];
    metadata_row[10] = block->dx[0];
    metadata_row[11] = block->dx[1];
    metadata_row[12] = block->dx[2];
    metadata_row[13] = (double)block->parent;
    metadata_row[14] = (double)block->base_block;
    for (n = 0; n < 8; ++n) {
        metadata_row[15 + n] = (double)block->children[n];
    }
    idx = 23;
    for (n = 0; n < 56; ++n) {
        int d;

        metadata_row[idx++] = (double)block->slot[n].id;
        metadata_row[idx++] = 0.0;
        for (d = 0; d < 3; ++d) {
            metadata_row[idx++] = block->slot[n].xmin[d];
        }
        for (d = 0; d < 3; ++d) {
            metadata_row[idx++] = block->slot[n].xmax[d];
        }
        for (d = 0; d < 3; ++d) {
            metadata_row[idx++] = block->slot[n].dx[d];
        }
    }
}

static void prj_io_unpack_metadata(prj_block *block, const double *metadata_row)
{
    int n;
    int idx;

    block->id = (int)metadata_row[0];
    block->rank = (int)metadata_row[1];
    block->level = (int)metadata_row[2];
    block->active = (int)metadata_row[3];
    block->xmin[0] = metadata_row[4];
    block->xmin[1] = metadata_row[5];
    block->xmin[2] = metadata_row[6];
    block->xmax[0] = metadata_row[7];
    block->xmax[1] = metadata_row[8];
    block->xmax[2] = metadata_row[9];
    block->dx[0] = metadata_row[10];
    block->dx[1] = metadata_row[11];
    block->dx[2] = metadata_row[12];
    block->parent = (int)metadata_row[13];
    block->base_block = (int)metadata_row[14];
    for (n = 0; n < 8; ++n) {
        block->children[n] = (int)metadata_row[15 + n];
    }
    idx = 23;
    for (n = 0; n < 56; ++n) {
        int d;

        block->slot[n].id = (int)metadata_row[idx++];
        block->slot[n].rank = (int)metadata_row[idx++];
        for (d = 0; d < 3; ++d) {
            block->slot[n].xmin[d] = metadata_row[idx++];
        }
        for (d = 0; d < 3; ++d) {
            block->slot[n].xmax[d] = metadata_row[idx++];
        }
        for (d = 0; d < 3; ++d) {
            block->slot[n].dx[d] = metadata_row[idx++];
        }
    }
}

void prj_io_write_restart(const prj_mesh *mesh, double time, int step, int dump_count,
    double next_output_time, double next_restart_time, double dt)
{
    char filename[64];
    hid_t file;
    hid_t space_data;
    hid_t space_meta;
    hid_t dset_data;
    hid_t dset_meta;
    hsize_t dims_data[3];
    hsize_t dims_meta[2];
    int bidx;

    snprintf(filename, sizeof(filename), "output/restart_%08d.h5", step);
    if (prj_io_is_root_rank()) {
        fprintf(stderr, "save restart file %s\n", filename);
    }
    dims_data[0] = (hsize_t)mesh->nblocks;
    dims_data[1] = (hsize_t)PRJ_NVAR_PRIM;
    dims_data[2] = (hsize_t)(PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE);
    dims_meta[0] = (hsize_t)mesh->nblocks;
    dims_meta[1] = (hsize_t)PRJ_IO_METADATA_SIZE;
    file = prj_io_create_file(filename);
    prj_io_write_attr_double(file, "time", time);
    prj_io_write_attr_int(file, "step", step);
    prj_io_write_attr_int(file, "dump_count", dump_count);
    prj_io_write_attr_double(file, "next_output_time", next_output_time);
    prj_io_write_attr_double(file, "next_restart_time", next_restart_time);
    prj_io_write_attr_double(file, "dt", dt);
    prj_io_write_attr_int(file, "nblocks", mesh->nblocks);
    prj_io_write_attr_int(file, "nvar_prim", PRJ_NVAR_PRIM);
    prj_io_write_attr_int(file, "block_size", PRJ_BLOCK_SIZE);
    prj_io_write_attr_int(file, "max_level", mesh->max_level);
    prj_io_write_attr_int3(file, "root_nx", mesh->root_nx);
    prj_io_write_attr_double6(file, "coord", &mesh->coord);

    space_data = H5Screate_simple(3, dims_data, dims_data);
    dset_data = H5Dcreate2(file, "Data", H5T_NATIVE_DOUBLE, space_data, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    space_meta = H5Screate_simple(2, dims_meta, dims_meta);
    dset_meta = H5Dcreate2(file, "MetaData", H5T_NATIVE_DOUBLE, space_meta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];

        if (prj_io_is_root_rank()) {
            hsize_t start_meta[2] = {(hsize_t)bidx, 0};
            hsize_t count_meta[2] = {1, (hsize_t)PRJ_IO_METADATA_SIZE};
            hid_t mem_meta = H5Screate_simple(2, count_meta, count_meta);
            hid_t file_meta = H5Dget_space(dset_meta);
            double metadata_row[PRJ_IO_METADATA_SIZE];
            hid_t dxpl = prj_io_data_xfer_plist();

            prj_io_fill_metadata(block, metadata_row);
            H5Sselect_hyperslab(file_meta, H5S_SELECT_SET, start_meta, 0, count_meta, 0);
            H5Dwrite(dset_meta, H5T_NATIVE_DOUBLE, mem_meta, file_meta, dxpl, metadata_row);
            prj_io_close_dxpl(dxpl);
            H5Sclose(file_meta);
            H5Sclose(mem_meta);
        }
        if (block->active == 1 && block->W != 0) {
            hsize_t start_data[3] = {(hsize_t)bidx, 0, 0};
            hsize_t count_data[3] = {1, (hsize_t)PRJ_NVAR_PRIM, (hsize_t)(PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE)};
            hid_t mem_data = H5Screate_simple(3, count_data, count_data);
            hid_t file_data = H5Dget_space(dset_data);
            double *buffer = (double *)calloc((size_t)PRJ_NVAR_PRIM * (size_t)(PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE), sizeof(*buffer));
            hid_t dxpl = prj_io_data_xfer_plist();
            int v;
            int i;
            int j;
            int k;

            if (buffer == 0) {
                prj_io_fail("prj_io_write_restart: allocation failed");
            }
            for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                    for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                        for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                            size_t cell = (size_t)i * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE + (size_t)j * PRJ_BLOCK_SIZE + (size_t)k;
                            buffer[(size_t)v * count_data[2] + cell] = block->W[VIDX(v, i, j, k)];
                        }
                    }
                }
            }
            H5Sselect_hyperslab(file_data, H5S_SELECT_SET, start_data, 0, count_data, 0);
            H5Dwrite(dset_data, H5T_NATIVE_DOUBLE, mem_data, file_data, dxpl, buffer);
            prj_io_close_dxpl(dxpl);
            free(buffer);
            H5Sclose(file_data);
            H5Sclose(mem_data);
        }
    }
    H5Dclose(dset_data);
    H5Sclose(space_data);
    H5Dclose(dset_meta);
    H5Sclose(space_meta);
    H5Fclose(file);
}

void prj_io_read_restart(prj_mesh *mesh, const prj_eos *eos, const char *filename,
    double *time, int *step, int *dump_count, double *next_output_time, double *next_restart_time,
    double *dt)
{
    hid_t file;
    hid_t dset_data;
    hid_t dset_meta;
    int nblocks;
    int block_size;
    int nvar_prim;
    int root_nx[3];
    int max_level;
    prj_coord coord;
    double *metadata;
    int bidx;
    prj_mpi *mpi = prj_mpi_current();
    prj_bc bc = {
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW,
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW,
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW
    };

    file = prj_io_open_file_readonly(filename);
    prj_io_read_attr_double(file, "time", time);
    prj_io_read_attr_int(file, "step", step);
    if (dump_count != 0) {
        *dump_count = prj_io_read_attr_int_optional(file, "dump_count", 0);
    }
    if (next_output_time != 0) {
        *next_output_time = prj_io_read_attr_double_optional(file, "next_output_time", -1.0);
    }
    if (next_restart_time != 0) {
        *next_restart_time = prj_io_read_attr_double_optional(file, "next_restart_time", -1.0);
    }
    if (dt != 0) {
        *dt = prj_io_read_attr_double_optional(file, "dt", 0.0);
    }
    prj_io_read_attr_int(file, "nblocks", &nblocks);
    prj_io_read_attr_int(file, "nvar_prim", &nvar_prim);
    prj_io_read_attr_int(file, "block_size", &block_size);
    prj_io_read_attr_int(file, "max_level", &max_level);
    prj_io_read_attr_int3(file, "root_nx", root_nx);
    prj_io_read_attr_double6(file, "coord", &coord);
    if (nvar_prim != PRJ_NVAR_PRIM || block_size != PRJ_BLOCK_SIZE) {
        prj_io_fail("prj_io_read_restart: incompatible restart metadata");
    }

    if (prj_mesh_init(mesh, root_nx[0], root_nx[1], root_nx[2], max_level, &coord) != 0) {
        prj_io_fail("prj_io_read_restart: mesh init failed");
    }
    mesh->nblocks = nblocks;

    metadata = (double *)calloc((size_t)nblocks * PRJ_IO_METADATA_SIZE, sizeof(*metadata));
    if (metadata == 0) {
        prj_io_fail("prj_io_read_restart: allocation failed");
    }

    dset_meta = H5Dopen2(file, "MetaData", H5P_DEFAULT);
    H5Dread(dset_meta, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata);
    H5Dclose(dset_meta);

    for (bidx = 0; bidx < nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        const double *meta_row = &metadata[(size_t)bidx * PRJ_IO_METADATA_SIZE];

        prj_io_unpack_metadata(block, meta_row);
        prj_block_setup_geometry(block, &coord);
    }

    prj_mpi_decompose(mesh);
    if (mpi != 0) {
        prj_mpi_prepare(mesh, mpi);
    }

    dset_data = H5Dopen2(file, "Data", H5P_DEFAULT);
    for (bidx = 0; bidx < nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];

        if (block->id >= 0 && block->active == 1 && prj_io_is_local_owner(block) && block->W != 0) {
            hsize_t start_data[3] = {(hsize_t)bidx, 0, 0};
            hsize_t count_data[3] = {1, (hsize_t)PRJ_NVAR_PRIM, (hsize_t)(PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE)};
            hid_t mem_data = H5Screate_simple(3, count_data, count_data);
            hid_t file_data;
            double *buffer = (double *)calloc((size_t)PRJ_NVAR_PRIM * (size_t)(PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE), sizeof(*buffer));
            hid_t dxpl = prj_io_data_xfer_plist();
            int v;
            int i;
            int j;
            int k;

            if (buffer == 0) {
                prj_io_fail("prj_io_read_restart: allocation failed");
            }
            file_data = H5Dget_space(dset_data);
            H5Sselect_hyperslab(file_data, H5S_SELECT_SET, start_data, 0, count_data, 0);
            H5Dread(dset_data, H5T_NATIVE_DOUBLE, mem_data, file_data, dxpl, buffer);
            prj_io_close_dxpl(dxpl);
            H5Sclose(file_data);
            for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                    for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                        for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                            size_t cell = (size_t)i * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE + (size_t)j * PRJ_BLOCK_SIZE + (size_t)k;
                            block->W[VIDX(v, i, j, k)] = buffer[(size_t)v * count_data[2] + cell];
                        }
                    }
                }
            }
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        double Wcell[PRJ_NVAR_PRIM];
                        double Ucell[PRJ_NVAR_CONS];

                        for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                            Wcell[v] = block->W[VIDX(v, i, j, k)];
                        }
                        prj_eos_prim2cons((prj_eos *)eos, Wcell, Ucell);
                        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                            block->U[VIDX(v, i, j, k)] = Ucell[v];
                            block->dUdt[VIDX(v, i, j, k)] = 0.0;
                            block->flux[0][VIDX(v, i, j, k)] = 0.0;
                            block->flux[1][VIDX(v, i, j, k)] = 0.0;
                            block->flux[2][VIDX(v, i, j, k)] = 0.0;
                        }
                    }
                }
            }
            free(buffer);
            H5Sclose(mem_data);
        }
    }
    H5Dclose(dset_data);
    H5Fclose(file);
    prj_mesh_update_max_active_level(mesh);
    if (prj_io_is_root_rank()) {
        fprintf(stderr, "read restart file %s with %d blocks\n", filename, prj_mesh_count_active(mesh));
    }
    prj_eos_fill_active_cells(mesh, (prj_eos *)eos, 1);
    prj_boundary_fill_ghosts(mesh, &bc, 1);
    prj_eos_fill_mesh(mesh, (prj_eos *)eos, 1);
    free(metadata);
}

void prj_io_write_dump(const prj_mesh *mesh, const char *basename, int dump_index, int step, double time)
{
    char filename[256];
    hid_t file;
    hid_t dset_level;
    hid_t dset_coord;
    hid_t dset_var[PRJ_NVAR_PRIM];
    hid_t space_level;
    hid_t space_coord;
    hid_t space_var;
    hsize_t dims_level[1];
    hsize_t dims_coord[2];
    hsize_t dims_var[4];
    int local_active = 0;
    int nactive = 0;
    int bidx;
    int active_idx = 0;
    hid_t dump_real_type = prj_io_dump_real_hdf5_type();

    snprintf(filename, sizeof(filename), "%s_%05d.h5", basename, dump_index);
    if (prj_io_is_root_rank()) {
        fprintf(stderr, "save dump file %s\n", filename);
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        if (mesh->blocks[bidx].id >= 0 && mesh->blocks[bidx].active == 1 && prj_io_is_local_owner(&mesh->blocks[bidx])) {
            local_active += 1;
        }
    }
    nactive = (int)prj_mpi_global_sum((double)local_active);
    dims_level[0] = (hsize_t)nactive;
    dims_coord[0] = (hsize_t)nactive;
    dims_coord[1] = 3;
    dims_var[0] = (hsize_t)nactive;
    dims_var[1] = PRJ_BLOCK_SIZE;
    dims_var[2] = PRJ_BLOCK_SIZE;
    dims_var[3] = PRJ_BLOCK_SIZE;
    file = prj_io_create_file(filename);
    prj_io_write_attr_int(file, "dump_index", dump_index);
    prj_io_write_attr_int(file, "step", step);
    prj_io_write_attr_double(file, "time", time);
    prj_io_write_attr_int3(file, "root_nx", mesh->root_nx);
    prj_io_write_attr_double6(file, "coord", &mesh->coord);
    space_level = H5Screate_simple(1, dims_level, dims_level);
    dset_level = H5Dcreate2(file, "level", H5T_NATIVE_INT, space_level, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    space_coord = H5Screate_simple(2, dims_coord, dims_coord);
    dset_coord = H5Dcreate2(file, "coordinate", H5T_NATIVE_DOUBLE, space_coord, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    space_var = H5Screate_simple(4, dims_var, dims_var);
    for (bidx = 0; bidx < PRJ_NVAR_PRIM; ++bidx) {
        char name[32];

        prj_io_dataset_name(bidx, name, sizeof(name));
        dset_var[bidx] = H5Dcreate2(file, name, dump_real_type, space_var, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    {
        int offset = 0;
#if defined(PRJ_ENABLE_MPI)
        prj_mpi *mpi = prj_mpi_current();

        if (mpi != 0 && mpi->totrank > 1) {
            MPI_Exscan(&local_active, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (mpi->rank == 0) {
                offset = 0;
            }
        }
#endif
        active_idx = offset;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int var;

        if (block->id < 0 || block->active != 1 || !prj_io_is_local_owner(block) || block->W == 0) {
            continue;
        }
        {
            hsize_t start_level[1] = {(hsize_t)active_idx};
            hsize_t count_level[1] = {1};
            hid_t mem_level = H5Screate_simple(1, count_level, count_level);
            hid_t file_level = H5Dget_space(dset_level);
            int level_value = block->level;
            hid_t dxpl = prj_io_data_xfer_plist();

            H5Sselect_hyperslab(file_level, H5S_SELECT_SET, start_level, 0, count_level, 0);
            H5Dwrite(dset_level, H5T_NATIVE_INT, mem_level, file_level, dxpl, &level_value);
            prj_io_close_dxpl(dxpl);
            H5Sclose(file_level);
            H5Sclose(mem_level);
        }
        {
            hsize_t start_coord[2] = {(hsize_t)active_idx, 0};
            hsize_t count_coord[2] = {1, 3};
            hid_t mem_coord = H5Screate_simple(2, count_coord, count_coord);
            hid_t file_coord = H5Dget_space(dset_coord);
            double coord_value[3];
            hid_t dxpl = prj_io_data_xfer_plist();

            coord_value[0] = block->xmin[0];
            coord_value[1] = block->xmin[1];
            coord_value[2] = block->xmin[2];
            H5Sselect_hyperslab(file_coord, H5S_SELECT_SET, start_coord, 0, count_coord, 0);
            H5Dwrite(dset_coord, H5T_NATIVE_DOUBLE, mem_coord, file_coord, dxpl, coord_value);
            prj_io_close_dxpl(dxpl);
            H5Sclose(file_coord);
            H5Sclose(mem_coord);
        }
        for (var = 0; var < PRJ_NVAR_PRIM; ++var) {
            hsize_t start_var[4] = {(hsize_t)active_idx, 0, 0, 0};
            hsize_t count_var[4] = {1, PRJ_BLOCK_SIZE, PRJ_BLOCK_SIZE, PRJ_BLOCK_SIZE};
            hid_t mem_var = H5Screate_simple(4, count_var, count_var);
            hid_t file_var = H5Dget_space(dset_var[var]);
            hid_t dxpl = prj_io_data_xfer_plist();
            int i;
            int j;
            int k;
            size_t ncells = (size_t)PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE;

            #if PRJ_DUMP_SINGLE_PRECISION
            float *buffer = (float *)calloc(ncells, sizeof(*buffer));

            if (buffer == 0) {
                prj_io_fail("prj_io_write_dump: allocation failed");
            }
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        size_t cell = (size_t)i * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE + (size_t)j * PRJ_BLOCK_SIZE + (size_t)k;
                        buffer[cell] = (float)block->W[VIDX(var, i, j, k)];
                    }
                }
            }
            H5Sselect_hyperslab(file_var, H5S_SELECT_SET, start_var, 0, count_var, 0);
            H5Dwrite(dset_var[var], H5T_NATIVE_FLOAT, mem_var, file_var, dxpl, buffer);
            #else
            double *buffer = (double *)calloc(ncells, sizeof(*buffer));

            if (buffer == 0) {
                prj_io_fail("prj_io_write_dump: allocation failed");
            }
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        size_t cell = (size_t)i * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE + (size_t)j * PRJ_BLOCK_SIZE + (size_t)k;
                        buffer[cell] = block->W[VIDX(var, i, j, k)];
                    }
                }
            }
            H5Sselect_hyperslab(file_var, H5S_SELECT_SET, start_var, 0, count_var, 0);
            H5Dwrite(dset_var[var], H5T_NATIVE_DOUBLE, mem_var, file_var, dxpl, buffer);
            #endif
            prj_io_close_dxpl(dxpl);
            free(buffer);
            H5Sclose(file_var);
            H5Sclose(mem_var);
        }
        active_idx += 1;
    }
    H5Dclose(dset_level);
    H5Sclose(space_level);
    H5Dclose(dset_coord);
    H5Sclose(space_coord);
    for (bidx = 0; bidx < PRJ_NVAR_PRIM; ++bidx) {
        H5Dclose(dset_var[bidx]);
    }
    H5Sclose(space_var);
    {
        const prj_grav_mono *grav_mono = prj_gravity_active_monopole();

        if (grav_mono != 0 && grav_mono->nbins > 0 && grav_mono->accel != 0 && grav_mono->lapse != 0) {
            hsize_t dims_accel[1] = {(hsize_t)grav_mono->nbins};
            hsize_t dims_lapse[1] = {(hsize_t)grav_mono->nbins + 1};
            hid_t space_accel = H5Screate_simple(1, dims_accel, dims_accel);
            hid_t space_lapse = H5Screate_simple(1, dims_lapse, dims_lapse);
            hid_t dset_accel = H5Dcreate2(file, "accel", H5T_NATIVE_DOUBLE, space_accel,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            hid_t dset_lapse = H5Dcreate2(file, "lapse", H5T_NATIVE_DOUBLE, space_lapse,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            if (prj_io_is_root_rank()) {
                hid_t dxpl_accel = prj_io_data_xfer_plist();
                hid_t dxpl_lapse = prj_io_data_xfer_plist();

                H5Dwrite(dset_accel, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dxpl_accel, grav_mono->accel);
                H5Dwrite(dset_lapse, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dxpl_lapse, grav_mono->lapse);
                prj_io_close_dxpl(dxpl_accel);
                prj_io_close_dxpl(dxpl_lapse);
            }
            H5Dclose(dset_accel);
            H5Dclose(dset_lapse);
            H5Sclose(space_accel);
            H5Sclose(space_lapse);
        }
    }
    H5Fclose(file);
}
