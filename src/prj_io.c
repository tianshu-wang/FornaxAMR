#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <dirent.h>

#include <hdf5.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#define PRJ_IO_METADATA_SIZE 639
#define PRJ_IO_DUMP_NAME_SIZE 32

#if PRJ_DUMP_SINGLE_PRECISION
typedef float prj_io_dump_real;
#else
typedef double prj_io_dump_real;
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
    if (strcmp(value, "fractional_jump") == 0) {
        *amr_estimator = PRJ_AMR_ESTIMATOR_FRACTIONAL_JUMP;
        return 0;
    }
    return 1;
}

static int prj_io_parse_lohner_var(const char *value, int *lohner_var)
{
    if (value == 0 || lohner_var == 0) {
        return 1;
    }
    if (strcmp(value, "density") == 0) {
        *lohner_var = PRJ_LOHNER_VAR_DENSITY;
        return 0;
    }
    if (strcmp(value, "pressure") == 0) {
        *lohner_var = PRJ_LOHNER_VAR_PRESSURE;
        return 0;
    }
    if (strcmp(value, "temperature") == 0) {
        *lohner_var = PRJ_LOHNER_VAR_TEMPERATURE;
        return 0;
    }
    return 1;
}

static int prj_io_parse_fractional_jump_var(const char *value, int *jump_var)
{
    if (value == 0 || jump_var == 0) {
        return 1;
    }
    if (strcmp(value, "density") == 0) {
        *jump_var = PRJ_FRACTIONAL_JUMP_VAR_DENSITY;
        return 0;
    }
    if (strcmp(value, "pressure") == 0) {
        *jump_var = PRJ_FRACTIONAL_JUMP_VAR_PRESSURE;
        return 0;
    }
    return 1;
}

#if PRJ_MHD
static int prj_io_parse_mhd_init_type(const char *value, int *mhd_init_type)
{
    if (value == 0 || mhd_init_type == 0) {
        return 1;
    }
    if (strcmp(value, "uniform") == 0) {
        *mhd_init_type = PRJ_MHD_INIT_UNIFORM;
        return 0;
    }
    if (strcmp(value, "dipole") == 0 ||
        strcmp(value, "dipole_core") == 0 ||
        strcmp(value, "dipole_uniform_core") == 0) {
        *mhd_init_type = PRJ_MHD_INIT_DIPOLE_CORE;
        return 0;
    }
    return 1;
}
#endif

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
    sim->x_com_err_tol = 0.5;
    sim->grav.use_multipole_gravity = 1;
    sim->t_end = 0.1;
    sim->output_dt = -1.0;
    sim->restart_dt = -1.0;
    sim->max_steps = 100;
    sim->output_interval = -1;
    sim->restart_interval = -1;
    sim->timer_interval = 100;
    sim->amr_interval = -1;
    sim->progenitor_file[0] = '\0';
    sim->perturbation_gaussian_norm = 0.0;
    sim->perturbation_seed = 0ULL;
    strncpy(sim->problem_name, "general", sizeof(sim->problem_name) - 1);
    sim->problem_name[sizeof(sim->problem_name) - 1] = '\0';
    sim->restart_from_file = 0;
    sim->restart_from_latest = 0;
    sim->restart_file_name[0] = '\0';
    sim->mesh.root_nx[0] = 8;
    sim->mesh.root_nx[1] = 8;
    sim->mesh.root_nx[2] = 8;
    sim->mesh.x_com[0] = 0.0;
    sim->mesh.x_com[1] = 0.0;
    sim->mesh.x_com[2] = 0.0;
    sim->mesh.max_level = 0;
    sim->mesh.min_dx = 0.0;
    sim->mesh.min_allowable_cell_size = 0.0;
    sim->mesh.max_blocks = 65536;
    {
        int amr_idx;

        for (amr_idx = 0; amr_idx < PRJ_AMR_N; ++amr_idx) {
            sim->mesh.amr_refine_thresh[amr_idx] = 0.5;
            sim->mesh.amr_derefine_thresh[amr_idx] = 0.2;
            sim->mesh.amr_estimator[amr_idx] = PRJ_AMR_ESTIMATOR_LOEHNER;
            sim->mesh.amr_lohner_var[amr_idx] = PRJ_LOHNER_VAR_PRESSURE;
            sim->mesh.amr_fractional_jump_var[amr_idx] = PRJ_FRACTIONAL_JUMP_VAR_PRESSURE;
            sim->mesh.amr_lohner_eps[amr_idx] = 0.1;
            sim->mesh.amr_criterion_set[amr_idx] = 0;
        }
        sim->mesh.amr_estimator[0] = PRJ_AMR_ESTIMATOR_VELOCITY;
        sim->mesh.amr_criterion_set[0] = 1;
    }
    sim->mesh.use_amr_angular_resolution_limit = 0;
    sim->mesh.use_BJ = 0;
    sim->mesh.amr_init_scale_factor = 0.5;
    sim->mesh.E_floor = -1.0;
    sim->eos.kind = PRJ_EOS_KIND_IDEAL;
    sim->eos.filename[0] = '\0';
    sim->eos.E_injected = 0.0;
    sim->rad.maxiter = 20;
    sim->rad.implicit_err_tol = 1.0e-6;
#if PRJ_NRAD > 0
    sim->rad.kom_epsilon = 0.1;
    sim->rad.kom_delta = 0.01;
    sim->rad.kom_dtmin = 1.0e-20;
    sim->rad.kom_rhocut = 1.0e13;
    sim->rad.min_inel_density = 1.0e8;
    {
        double kom_c = 299792458.0;
        double kom_e_unit = 1e6 * 1.60217662e-19;
        double kom_t_unit = 1.0545718e-34 / kom_e_unit;
        double kom_l_unit = kom_t_unit * kom_c;
        double kom_G2 = 1.327817e-22;
        double kom_m = 939.0;
        double kom_prot = (-0.5 * (1.0 - 4.0 * 0.23122)) * (-0.5 * (1.0 - 4.0 * 0.23122))
            + 5.0 * (-1.2723 / 2.0) * (-1.2723 / 2.0);
        double kom_neut = (-0.5) * (-0.5) + 5.0 * (1.2723 / 2.0) * (1.2723 / 2.0);
        sim->rad.kom_nucinel_const = 2.0 * kom_G2 * 1e3 * kom_c * kom_c
            * kom_l_unit * kom_l_unit * kom_l_unit
            / (kom_m * kom_e_unit * 3.0 * M_PI * kom_m * kom_t_unit);
        sim->rad.kom_nucinel_prot = kom_prot;
        sim->rad.kom_nucinel_neut = kom_neut;
    }
    {
        int nu_i;
        for (nu_i = 0; nu_i < PRJ_NRAD; ++nu_i) {
            sim->rad.kom_Ecut[nu_i] = 300.0;
        }
    }
#endif
#if PRJ_MHD
    sim->mhd_init_type = PRJ_MHD_INIT_UNIFORM;
    sim->mhd_B_norm = 0.0;
    sim->mhd_B_scale = 1.0;
#endif
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
        } else if (strcmp(key, "x_com_err_tol") == 0) {
            sim->x_com_err_tol = strtod(value, &endptr);
        } else if (strcmp(key, "use_multipole_gravity") == 0) {
            sim->grav.use_multipole_gravity = (int)strtol(value, &endptr, 10);
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
        } else if (strcmp(key, "timer_interval") == 0) {
            sim->timer_interval = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "amr_interval") == 0) {
            sim->amr_interval = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "perturbation_gaussian_norm") == 0) {
            sim->perturbation_gaussian_norm = strtod(value, &endptr);
        } else if (strcmp(key, "perturbation_seed") == 0) {
            sim->perturbation_seed = strtoull(value, &endptr, 10);
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
        } else if (strcmp(key, "max_blocks") == 0) {
            sim->mesh.max_blocks = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "min_dx") == 0) {
            sim->mesh.min_dx = strtod(value, &endptr);
        } else if (prj_io_parse_amr_slot_key(key, "amr_estimator", &amr_slot)) {
            if (prj_io_parse_amr_estimator(value, &sim->mesh.amr_estimator[amr_slot]) != 0) {
                endptr = value;
            } else {
                endptr = value + strlen(value);
                sim->mesh.amr_criterion_set[amr_slot] = 1;
            }
        } else if (prj_io_parse_amr_slot_key(key, "amr_lohner_var", &amr_slot)) {
            if (prj_io_parse_lohner_var(value, &sim->mesh.amr_lohner_var[amr_slot]) != 0) {
                endptr = value;
            } else {
                endptr = value + strlen(value);
            }
        } else if (prj_io_parse_amr_slot_key(key, "amr_fractional_jump_var", &amr_slot)) {
            if (prj_io_parse_fractional_jump_var(value, &sim->mesh.amr_fractional_jump_var[amr_slot]) != 0) {
                endptr = value;
            } else {
                endptr = value + strlen(value);
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
        } else if (prj_io_parse_amr_slot_key(key, "amr_lohner_eps", &amr_slot)) {
            sim->mesh.amr_lohner_eps[amr_slot] = strtod(value, &endptr);
        } else if (strcmp(key, "use_amr_angular_resolution_limit") == 0) {
            sim->mesh.use_amr_angular_resolution_limit = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "use_BJ") == 0 || strcmp(key, "use_bj") == 0) {
            sim->mesh.use_BJ = (int)strtol(value, &endptr, 10);
        } else if (strcmp(key, "amr_init_scale_factor") == 0) {
            sim->mesh.amr_init_scale_factor = strtod(value, &endptr);
        } else if (strcmp(key, "E_floor") == 0) {
            sim->mesh.E_floor = strtod(value, &endptr);
#if PRJ_MHD
        } else if (strcmp(key, "B_norm") == 0 || strcmp(key, "mhd_B_norm") == 0) {
            sim->mhd_B_norm = strtod(value, &endptr);
        } else if (strcmp(key, "B_scale") == 0 || strcmp(key, "mhd_B_scale") == 0) {
            sim->mhd_B_scale = strtod(value, &endptr);
        } else if (strcmp(key, "B_init") == 0 ||
            strcmp(key, "mhd_init") == 0 ||
            strcmp(key, "mhd_init_type") == 0) {
            if (prj_io_parse_mhd_init_type(value, &sim->mhd_init_type) != 0) {
                endptr = value;
            } else {
                endptr = value + strlen(value);
            }
#endif
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
        } else if (strcmp(key, "restart_from_latest") == 0) {
            sim->restart_from_latest = (int)strtol(value, &endptr, 10);
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
        } else if (strcmp(key, "eleinel_table_dir") == 0) {
            strncpy(sim->rad.eleinel_table_dir, value, sizeof(sim->rad.eleinel_table_dir) - 1);
            sim->rad.eleinel_table_dir[sizeof(sim->rad.eleinel_table_dir) - 1] = '\0';
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
        } else if (strcmp(key, "kom_epsilon") == 0) {
            sim->rad.kom_epsilon = strtod(value, &endptr);
        } else if (strcmp(key, "kom_delta") == 0) {
            sim->rad.kom_delta = strtod(value, &endptr);
        } else if (strcmp(key, "kom_dtmin") == 0) {
            sim->rad.kom_dtmin = strtod(value, &endptr);
        } else if (strcmp(key, "kom_rhocut") == 0) {
            sim->rad.kom_rhocut = strtod(value, &endptr);
        } else if (strcmp(key, "min_inel_density") == 0) {
            sim->rad.min_inel_density = strtod(value, &endptr);
        } else if (strncmp(key, "kom_Ecut_", 9) == 0) {
            char *kend;
            long idx = strtol(key + 9, &kend, 10);

            if (*kend != '\0' || idx < 0 || idx >= PRJ_NRAD) {
                fprintf(stderr, "prj_io_parser: bad kom_Ecut index '%s' in %s:%d\n", key, filename, lineno);
                fclose(fp);
                exit(1);
            }
            sim->rad.kom_Ecut[idx] = strtod(value, &endptr);
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

static void prj_io_fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}

static int prj_io_is_local_owner(const prj_mpi *mpi, const prj_block *block)
{
    return block != 0 && block->id >= 0 && (mpi == 0 || block->rank == mpi->rank);
}

static int prj_io_is_root_rank(const prj_mpi *mpi)
{
    return mpi == 0 || mpi->rank == 0;
}

static int prj_io_restart_write_data_block(const prj_block *block)
{
    return block != 0 && block->active == 1 && block->W != 0;
}

static int prj_io_restart_read_data_block(const prj_mpi *mpi, const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 &&
        prj_io_is_local_owner(mpi, block) && block->W != 0;
}

#if PRJ_MHD
static int prj_io_restart_write_bf_block(const prj_block *block)
{
    return block != 0 && block->active == 1 &&
        block->Bf[0] != 0 && block->Bf[1] != 0 && block->Bf[2] != 0;
}
#endif

static hid_t prj_io_create_file(const prj_mpi *mpi, const char *filename)
{
#if defined(PRJ_ENABLE_MPI)
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

static hid_t prj_io_open_file_readonly(const prj_mpi *mpi, const char *filename)
{
#if defined(PRJ_ENABLE_MPI)
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

static hid_t prj_io_data_xfer_plist(const prj_mpi *mpi)
{
#if defined(PRJ_ENABLE_MPI)
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

static void prj_io_write_hyperslab(hid_t dset, const prj_mpi *mpi, hid_t mem_type, int ndims,
    const hsize_t *start, const hsize_t *count, const void *buffer)
{
    hid_t mem_space = H5Screate_simple(ndims, count, count);
    hid_t file_space = H5Dget_space(dset);
    hid_t dxpl = prj_io_data_xfer_plist(mpi);

    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, 0, count, 0);
    H5Dwrite(dset, mem_type, mem_space, file_space, dxpl, buffer);
    prj_io_close_dxpl(dxpl);
    H5Sclose(file_space);
    H5Sclose(mem_space);
}

static void prj_io_read_hyperslab(hid_t dset, const prj_mpi *mpi, hid_t mem_type, int ndims,
    const hsize_t *start, const hsize_t *count, void *buffer)
{
    hid_t mem_space = H5Screate_simple(ndims, count, count);
    hid_t file_space = H5Dget_space(dset);
    hid_t dxpl = prj_io_data_xfer_plist(mpi);

    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, 0, count, 0);
    H5Dread(dset, mem_type, mem_space, file_space, dxpl, buffer);
    prj_io_close_dxpl(dxpl);
    H5Sclose(file_space);
    H5Sclose(mem_space);
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

static void prj_io_write_attr_double3(hid_t obj, const char *name, const double values[3])
{
    hsize_t dims[1] = {3};
    hid_t space = H5Screate_simple(1, dims, dims);
    hid_t attr = H5Acreate2(obj, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr, H5T_NATIVE_DOUBLE, values);
    H5Aclose(attr);
    H5Sclose(space);
}

static void prj_io_read_attr_double3(hid_t obj, const char *name, double values[3])
{
    hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);

    H5Aread(attr, H5T_NATIVE_DOUBLE, values);
    H5Aclose(attr);
}

static void prj_io_read_attr_double3_optional(hid_t obj, const char *name, double values[3],
    const double defaults[3])
{
    if (H5Aexists(obj, name) > 0) {
        prj_io_read_attr_double3(obj, name, values);
        return;
    }
    values[0] = defaults[0];
    values[1] = defaults[1];
    values[2] = defaults[2];
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
    metadata_row[14] = 0.0;
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

int prj_io_find_latest_restart(const char *dir, char *out_filename, size_t out_size, int *out_id)
{
    DIR *dh;
    struct dirent *de;
    int max_id = -1;

    if (dir == 0 || out_filename == 0 || out_size == 0) {
        return 1;
    }
    dh = opendir(dir);
    if (dh == 0) {
        return 1;
    }
    while ((de = readdir(dh)) != 0) {
        int id;
        char extra;
        if (sscanf(de->d_name, "restart_%d.h5%c", &id, &extra) == 1 && id >= 0) {
            if (id > max_id) {
                max_id = id;
            }
        }
    }
    closedir(dh);
    if (max_id < 0) {
        return 1;
    }
    if ((size_t)snprintf(out_filename, out_size, "%s/restart_%08d.h5", dir, max_id) >= out_size) {
        return 1;
    }
    if (out_id != 0) {
        *out_id = max_id;
    }
    return 0;
}

void prj_io_write_restart(const prj_mesh *mesh, const prj_mpi *mpi, double time, int step, int dump_count,
    double last_output_time, double last_restart_time, double dt)
{
    char filename[64];
    hid_t file;
    hid_t space_data;
    hid_t space_meta;
    hid_t dset_data;
    hid_t dset_meta;
#if PRJ_MHD
    hid_t space_bf;
    hid_t dset_bf;
#endif
    hsize_t dims_data[3];
    hsize_t dims_meta[2];
#if PRJ_MHD
    hsize_t dims_bf[3];
#endif
    int bidx;

    snprintf(filename, sizeof(filename), "output/restart_%08d.h5", step);
    if (prj_io_is_root_rank(mpi)) {
        fprintf(stderr, "save restart file %s\n", filename);
    }
    dims_data[0] = (hsize_t)mesh->nblocks;
    dims_data[1] = (hsize_t)PRJ_NVAR_PRIM;
    dims_data[2] = (hsize_t)(PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE);
    dims_meta[0] = (hsize_t)mesh->nblocks;
    dims_meta[1] = (hsize_t)PRJ_IO_METADATA_SIZE;
#if PRJ_MHD
    dims_bf[0] = (hsize_t)mesh->nblocks;
    dims_bf[1] = 3;
    dims_bf[2] = (hsize_t)PRJ_BLOCK_NFACES;
#endif
    file = prj_io_create_file(mpi, filename);
    prj_io_write_attr_double(file, "time", time);
    prj_io_write_attr_int(file, "step", step);
    prj_io_write_attr_int(file, "dump_count", dump_count);
    prj_io_write_attr_double(file, "last_output_time", last_output_time);
    prj_io_write_attr_double(file, "last_restart_time", last_restart_time);
    prj_io_write_attr_double(file, "dt", dt);
    prj_io_write_attr_int(file, "nblocks", mesh->nblocks);
    prj_io_write_attr_int(file, "nvar_prim", PRJ_NVAR_PRIM);
    prj_io_write_attr_int(file, "block_size", PRJ_BLOCK_SIZE);
    prj_io_write_attr_int(file, "max_level", mesh->max_level);
    prj_io_write_attr_double(file, "min_dx", mesh->min_dx);
    prj_io_write_attr_int3(file, "root_nx", mesh->root_nx);
    prj_io_write_attr_double6(file, "coord", &mesh->coord);
    {
        double x_com_values[3] = {0.0, 0.0, 0.0};

        if (mesh != 0) {
            x_com_values[0] = mesh->x_com[0];
            x_com_values[1] = mesh->x_com[1];
            x_com_values[2] = mesh->x_com[2];
        }
        prj_io_write_attr_double3(file, "x_com", x_com_values);
    }

    space_data = H5Screate_simple(3, dims_data, dims_data);
    dset_data = H5Dcreate2(file, "Data", H5T_NATIVE_DOUBLE, space_data, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    space_meta = H5Screate_simple(2, dims_meta, dims_meta);
    dset_meta = H5Dcreate2(file, "MetaData", H5T_NATIVE_DOUBLE, space_meta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#if PRJ_MHD
    space_bf = H5Screate_simple(3, dims_bf, dims_bf);
    dset_bf = H5Dcreate2(file, "Bf", H5T_NATIVE_DOUBLE, space_bf, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
    if (prj_io_is_root_rank(mpi) && mesh->nblocks > 0) {
        hsize_t start_meta[2] = {0, 0};
        hsize_t count_meta[2] = {(hsize_t)mesh->nblocks, (hsize_t)PRJ_IO_METADATA_SIZE};
        double *metadata = (double *)calloc((size_t)mesh->nblocks * PRJ_IO_METADATA_SIZE, sizeof(*metadata));

        if (metadata == 0) {
            prj_io_fail("prj_io_write_restart: metadata allocation failed");
        }
        for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
            prj_io_fill_metadata(&mesh->blocks[bidx], &metadata[(size_t)bidx * PRJ_IO_METADATA_SIZE]);
        }
        prj_io_write_hyperslab(dset_meta, mpi, H5T_NATIVE_DOUBLE, 2, start_meta, count_meta, metadata);
        free(metadata);
    }
    bidx = 0;
    while (bidx < mesh->nblocks) {
        int run_start;
        int run_len;
        hsize_t start_data[3];
        hsize_t count_data[3];
        double *buffer;
        int ridx;
        int v;
        int i;
        int j;
        int k;
        size_t ncells = (size_t)PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE;

        while (bidx < mesh->nblocks && !prj_io_restart_write_data_block(&mesh->blocks[bidx])) {
            bidx += 1;
        }
        run_start = bidx;
        while (bidx < mesh->nblocks && prj_io_restart_write_data_block(&mesh->blocks[bidx])) {
            bidx += 1;
        }
        run_len = bidx - run_start;
        if (run_len == 0) {
            continue;
        }
        start_data[0] = (hsize_t)run_start;
        start_data[1] = 0;
        start_data[2] = 0;
        count_data[0] = (hsize_t)run_len;
        count_data[1] = (hsize_t)PRJ_NVAR_PRIM;
        count_data[2] = (hsize_t)ncells;
        buffer = (double *)calloc((size_t)run_len * (size_t)PRJ_NVAR_PRIM * ncells, sizeof(*buffer));
        if (buffer == 0) {
            prj_io_fail("prj_io_write_restart: allocation failed");
        }
        for (ridx = 0; ridx < run_len; ++ridx) {
            const prj_block *block = &mesh->blocks[run_start + ridx];

            for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                double wscale = 1.0;
#if PRJ_NRAD > 0
                /* Radiation vars (v >= PRJ_NHYDRO) live in RAD_SCALE*erg units
                   internally; write them in physical erg so restart files stay
                   bit-identical to the pre-scaling format. */
                if (v >= PRJ_NHYDRO) {
                    wscale = RAD_SCALE;
                }
#endif
                for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                    for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                        for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                            size_t cell = (size_t)i * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE + (size_t)j * PRJ_BLOCK_SIZE + (size_t)k;
                            size_t offset = ((size_t)ridx * (size_t)PRJ_NVAR_PRIM + (size_t)v) * ncells + cell;

                            buffer[offset] = block->W[WIDX(v, i, j, k)] * wscale;
                        }
                    }
                }
            }
        }
        prj_io_write_hyperslab(dset_data, mpi, H5T_NATIVE_DOUBLE, 3, start_data, count_data, buffer);
        free(buffer);
    }
#if PRJ_MHD
    bidx = 0;
    while (bidx < mesh->nblocks) {
        int run_start;
        int run_len;
        hsize_t start_bf[3];
        hsize_t count_bf[3];
        double *buffer;
        int ridx;
        int d;
        int n;

        while (bidx < mesh->nblocks && !prj_io_restart_write_bf_block(&mesh->blocks[bidx])) {
            bidx += 1;
        }
        run_start = bidx;
        while (bidx < mesh->nblocks && prj_io_restart_write_bf_block(&mesh->blocks[bidx])) {
            bidx += 1;
        }
        run_len = bidx - run_start;
        if (run_len == 0) {
            continue;
        }
        start_bf[0] = (hsize_t)run_start;
        start_bf[1] = 0;
        start_bf[2] = 0;
        count_bf[0] = (hsize_t)run_len;
        count_bf[1] = 3;
        count_bf[2] = (hsize_t)PRJ_BLOCK_NFACES;
        buffer = (double *)calloc((size_t)run_len * 3U * (size_t)PRJ_BLOCK_NFACES, sizeof(*buffer));
        if (buffer == 0) {
            prj_io_fail("prj_io_write_restart: Bf allocation failed");
        }
        for (ridx = 0; ridx < run_len; ++ridx) {
            const prj_block *block = &mesh->blocks[run_start + ridx];

            for (d = 0; d < 3; ++d) {
                for (n = 0; n < PRJ_BLOCK_NFACES; ++n) {
                    size_t offset = ((size_t)ridx * 3U + (size_t)d) * (size_t)PRJ_BLOCK_NFACES + (size_t)n;

                    buffer[offset] = block->Bf[d][n];
                }
            }
        }
        prj_io_write_hyperslab(dset_bf, mpi, H5T_NATIVE_DOUBLE, 3, start_bf, count_bf, buffer);
        free(buffer);
    }
#endif
#if PRJ_MHD
    H5Dclose(dset_bf);
    H5Sclose(space_bf);
#endif
    H5Dclose(dset_data);
    H5Sclose(space_data);
    H5Dclose(dset_meta);
    H5Sclose(space_meta);
    H5Fclose(file);
}

void prj_io_read_restart(prj_mesh *mesh, const prj_eos *eos, prj_mpi *mpi, const char *filename,
    double *time, int *step, int *dump_count, double *last_output_time, double *last_restart_time,
    double *dt)
{
    hid_t file;
    hid_t dset_data;
    hid_t dset_meta;
#if PRJ_MHD
    hid_t dset_bf;
#endif
    int nblocks;
    int block_size;
    int nvar_prim;
    int root_nx[3];
    int max_level;
    double min_dx;
    prj_coord coord;
    double *metadata;
    int bidx;
    prj_bc bc = {
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW,
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW,
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW
    };

    file = prj_io_open_file_readonly(mpi, filename);
    prj_io_read_attr_double(file, "time", time);
    prj_io_read_attr_int(file, "step", step);
    if (dump_count != 0) {
        *dump_count = prj_io_read_attr_int_optional(file, "dump_count", 0);
    }
    if (last_output_time != 0) {
        *last_output_time = prj_io_read_attr_double_optional(file, "last_output_time", -1.0);
    }
    if (last_restart_time != 0) {
        *last_restart_time = prj_io_read_attr_double_optional(file, "last_restart_time", -1.0);
    }
    if (dt != 0) {
        *dt = prj_io_read_attr_double_optional(file, "dt", 0.0);
    }
    prj_io_read_attr_int(file, "nblocks", &nblocks);
    prj_io_read_attr_int(file, "nvar_prim", &nvar_prim);
    prj_io_read_attr_int(file, "block_size", &block_size);
    prj_io_read_attr_int(file, "max_level", &max_level);
    min_dx = prj_io_read_attr_double_optional(file, "min_dx", 0.0);
    prj_io_read_attr_int3(file, "root_nx", root_nx);
    prj_io_read_attr_double6(file, "coord", &coord);
    if (nvar_prim != PRJ_NVAR_PRIM || block_size != PRJ_BLOCK_SIZE) {
        prj_io_fail("prj_io_read_restart: incompatible restart metadata");
    }

    if (prj_mesh_init(mesh, root_nx[0], root_nx[1], root_nx[2], max_level, &coord) != 0) {
        prj_io_fail("prj_io_read_restart: mesh init failed");
    }
    {
        double defaults[3] = {0.0, 0.0, 0.0};

        prj_io_read_attr_double3_optional(file, "x_com", mesh->x_com, defaults);
    }
    mesh->min_dx = min_dx;
    prj_mesh_update_min_allowable_cell_size(mesh);
    if (nblocks > mesh->nblocks_max) {
        fprintf(stderr,
            "prj_io_read_restart: restart has nblocks=%d > max_blocks=%d. "
            "Increase max_blocks in the param file.\n",
            nblocks, mesh->nblocks_max);
        prj_io_fail("prj_io_read_restart: block capacity exceeded");
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
        prj_block_update_can_refine(block, mesh);
    }

    /* The metadata staging buffer (nblocks * PRJ_IO_METADATA_SIZE doubles,
       allocated on every rank) is only needed to unpack the block tree above.
       Free it here, before prj_mpi_prepare()/assign_block_storage and the Data
       read allocate the bulk per-block cell storage; otherwise its ~5 KB/block
       sits on every rank through the peak-memory phase and can OOM a restart
       that runs on the same rank count as the writer (which only allocated this
       buffer on the root rank). */
    free(metadata);
    metadata = 0;

    for (bidx = 0; bidx < nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (block->id < 0) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            int nid = block->slot[n].id;

            if (nid >= 0 && nid < mesh->nblocks && mesh->blocks[nid].id >= 0) {
                prj_neighbor_compute_geometry(block, &mesh->blocks[nid], &block->slot[n]);
            }
        }
    }
    if (prj_mesh_rebuild_morton_lookup(mesh) != 0) {
        prj_io_fail("prj_io_read_restart: Morton lookup rebuild failed");
    }

    prj_mpi_decompose(mesh, mpi);
    if (mpi != 0) {
        prj_mpi_prepare(mesh, mpi);
    }
    prj_mesh_update_r_com(mesh);

    dset_data = H5Dopen2(file, "Data", H5P_DEFAULT);
#if PRJ_MHD
    if (H5Lexists(file, "Bf", H5P_DEFAULT) <= 0) {
        prj_io_fail("prj_io_read_restart: MHD restart is missing Bf dataset");
    }
    dset_bf = H5Dopen2(file, "Bf", H5P_DEFAULT);
#endif
    bidx = 0;
    while (bidx < nblocks) {
        int run_start;
        int run_len;
        hsize_t start_data[3];
        hsize_t count_data[3];
        double *buffer;
        int ridx;
        int v;
        int i;
        int j;
        int k;
        size_t ncells = (size_t)PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE;

        while (bidx < nblocks && !prj_io_restart_read_data_block(mpi, &mesh->blocks[bidx])) {
            bidx += 1;
        }
        run_start = bidx;
        while (bidx < nblocks && prj_io_restart_read_data_block(mpi, &mesh->blocks[bidx])) {
            bidx += 1;
        }
        run_len = bidx - run_start;
        if (run_len == 0) {
            continue;
        }
        start_data[0] = (hsize_t)run_start;
        start_data[1] = 0;
        start_data[2] = 0;
        count_data[0] = (hsize_t)run_len;
        count_data[1] = (hsize_t)PRJ_NVAR_PRIM;
        count_data[2] = (hsize_t)ncells;
        buffer = (double *)calloc((size_t)run_len * (size_t)PRJ_NVAR_PRIM * ncells, sizeof(*buffer));
        if (buffer == 0) {
            prj_io_fail("prj_io_read_restart: allocation failed");
        }
        prj_io_read_hyperslab(dset_data, mpi, H5T_NATIVE_DOUBLE, 3, start_data, count_data, buffer);
        for (ridx = 0; ridx < run_len; ++ridx) {
            prj_block *block = &mesh->blocks[run_start + ridx];

            for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                double wscale = 1.0;
#if PRJ_NRAD > 0
                /* On-disk radiation vars are in physical erg; convert back to
                   the internal RAD_SCALE*erg units on read. */
                if (v >= PRJ_NHYDRO) {
                    wscale = 1.0 / RAD_SCALE;
                }
#endif
                for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                    for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                        for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                            size_t cell = (size_t)i * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE + (size_t)j * PRJ_BLOCK_SIZE + (size_t)k;
                            size_t offset = ((size_t)ridx * (size_t)PRJ_NVAR_PRIM + (size_t)v) * ncells + cell;

                            block->W[WIDX(v, i, j, k)] = buffer[offset] * wscale;
                        }
                    }
                }
            }
        }
        free(buffer);
#if PRJ_MHD
        {
            hsize_t start_bf[3] = {(hsize_t)run_start, 0, 0};
            hsize_t count_bf[3] = {(hsize_t)run_len, 3, (hsize_t)PRJ_BLOCK_NFACES};
            double *bf_buffer = (double *)calloc((size_t)run_len * 3U * (size_t)PRJ_BLOCK_NFACES, sizeof(*bf_buffer));
            int d;
            int n;

            if (bf_buffer == 0) {
                prj_io_fail("prj_io_read_restart: Bf allocation failed");
            }
            for (ridx = 0; ridx < run_len; ++ridx) {
                prj_block *block = &mesh->blocks[run_start + ridx];

                if (block->Bf[0] == 0 || block->Bf[1] == 0 || block->Bf[2] == 0) {
                    prj_io_fail("prj_io_read_restart: missing MHD block storage");
                }
            }
            prj_io_read_hyperslab(dset_bf, mpi, H5T_NATIVE_DOUBLE, 3, start_bf, count_bf, bf_buffer);
            for (ridx = 0; ridx < run_len; ++ridx) {
                prj_block *block = &mesh->blocks[run_start + ridx];

                for (d = 0; d < 3; ++d) {
                    for (n = 0; n < PRJ_BLOCK_NFACES; ++n) {
                        size_t offset = ((size_t)ridx * 3U + (size_t)d) * (size_t)PRJ_BLOCK_NFACES + (size_t)n;
                        double value = bf_buffer[offset];
                        double *bf1 = prj_block_bf_stage(block, d, 1);

                        block->Bf[d][n] = value;
                        bf1[n] = value;
                    }
                }
                prj_mhd_bf2bc_all((prj_eos *)eos, block, 0);
            }
            free(bf_buffer);
        }
#endif
        for (ridx = 0; ridx < run_len; ++ridx) {
            prj_block *block = &mesh->blocks[run_start + ridx];

            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        double Wcell[PRJ_NVAR_PRIM];
                        double Ucell[PRJ_NVAR_CONS];

                        for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                            Wcell[v] = block->W[WIDX(v, i, j, k)];
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
        }
    }
#if PRJ_MHD
    H5Dclose(dset_bf);
#endif
    H5Dclose(dset_data);
    H5Fclose(file);
    prj_mesh_update_max_active_level(mesh);
    if (prj_io_is_root_rank(mpi)) {
        fprintf(stderr, "read restart file %s with %d blocks\n", filename, prj_mesh_count_active(mesh));
    }
    prj_eos_fill_active_cells(mesh, (prj_eos *)eos, mpi, 1, PRJ_EOS_CTX_MAIN);
    prj_boundary_fill_ghosts(mesh, mpi, &bc, 1);
    prj_eos_fill_mesh(mesh, (prj_eos *)eos, mpi, 1, PRJ_EOS_CTX_MAIN);
}

static int prj_io_dump_write_eos_block(const prj_block *block)
{
    return block != 0 && block->active == 1 && block->eosvar != 0;
}

/* The dump now mirrors the restart on-disk layout: a small fixed number of fat
   datasets written as contiguous per-rank runs, instead of ~18 thin per-variable
   datasets. This collapses the collective HDF5 metadata work (dataset creates +
   attributes), which dominated dump time at high MPI-rank/block counts.

   Datasets, all indexed by global block id (dim 0 = mesh->nblocks):
     MetaData : double[nblocks][PRJ_IO_METADATA_SIZE]   (geometry/AMR tree, root)
     Data     : real[nblocks][PRJ_NVAR_PRIM][BS^3]       (all primitive vars)
     Eos      : real[nblocks][PRJ_NVAR_EOSVAR][BS^3]      (pressure,temp,gamma)
     Bf       : double[nblocks][3][NFACES]                (MHD face fields)
   Cell data uses the dump real type (single precision by default). Radiation
   prims are written in physical erg (* RAD_SCALE) to match the restart format;
   B is written unscaled (no sqrt(4*pi)) like the restart W. Only active blocks
   that this rank owns (W/eosvar/Bf allocated) are filled; inactive rows are left
   at the HDF5 fill value, exactly as prj_io_write_restart does. */
void prj_io_write_dump(const prj_mesh *mesh, const prj_grav *grav, const prj_mpi *mpi,
    int dump_index, int step, double time)
{
    char filename[256];
    hid_t file;
    hid_t space_data;
    hid_t dset_data;
    hid_t space_eos;
    hid_t dset_eos;
    hid_t space_meta;
    hid_t dset_meta;
#if PRJ_MHD && PRJ_MHD_DEBUG
    hid_t space_bf;
    hid_t dset_bf;
#endif
    hsize_t dims_data[3];
    hsize_t dims_eos[3];
    hsize_t dims_meta[2];
#if PRJ_MHD && PRJ_MHD_DEBUG
    hsize_t dims_bf[3];
#endif
    hid_t dump_real_type = prj_io_dump_real_hdf5_type();
    size_t ncells = (size_t)PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE;
    int bidx;

    snprintf(filename, sizeof(filename), "output/dump_%05d.h5", dump_index);
    if (prj_io_is_root_rank(mpi)) {
        fprintf(stderr, "save dump file %s\n", filename);
    }

    dims_data[0] = (hsize_t)mesh->nblocks;
    dims_data[1] = (hsize_t)PRJ_NVAR_PRIM;
    dims_data[2] = (hsize_t)ncells;
    dims_eos[0] = (hsize_t)mesh->nblocks;
    dims_eos[1] = (hsize_t)PRJ_NVAR_EOSVAR;
    dims_eos[2] = (hsize_t)ncells;
    dims_meta[0] = (hsize_t)mesh->nblocks;
    dims_meta[1] = (hsize_t)PRJ_IO_METADATA_SIZE;
#if PRJ_MHD && PRJ_MHD_DEBUG
    dims_bf[0] = (hsize_t)mesh->nblocks;
    dims_bf[1] = 3;
    dims_bf[2] = (hsize_t)PRJ_BLOCK_NFACES;
#endif

    file = prj_io_create_file(mpi, filename);
    prj_io_write_attr_int(file, "dump_format_version", 3);
    prj_io_write_attr_int(file, "dump_index", dump_index);
    prj_io_write_attr_int(file, "step", step);
    prj_io_write_attr_double(file, "time", time);
    prj_io_write_attr_int(file, "nblocks", mesh->nblocks);
    prj_io_write_attr_int(file, "nvar_prim", PRJ_NVAR_PRIM);
    prj_io_write_attr_int(file, "nvar_eos", PRJ_NVAR_EOSVAR);
    prj_io_write_attr_int(file, "block_size", PRJ_BLOCK_SIZE);
    prj_io_write_attr_int(file, "max_level", mesh->max_level);
    prj_io_write_attr_double(file, "min_dx", mesh->min_dx);
    prj_io_write_attr_int3(file, "root_nx", mesh->root_nx);
    prj_io_write_attr_double6(file, "coord", &mesh->coord);
    {
        double x_com_values[3] = {0.0, 0.0, 0.0};

        if (mesh != 0) {
            x_com_values[0] = mesh->x_com[0];
            x_com_values[1] = mesh->x_com[1];
            x_com_values[2] = mesh->x_com[2];
        }
        prj_io_write_attr_double3(file, "x_com", x_com_values);
    }

    space_data = H5Screate_simple(3, dims_data, dims_data);
    dset_data = H5Dcreate2(file, "Data", dump_real_type, space_data, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    space_eos = H5Screate_simple(3, dims_eos, dims_eos);
    dset_eos = H5Dcreate2(file, "Eos", dump_real_type, space_eos, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    space_meta = H5Screate_simple(2, dims_meta, dims_meta);
    dset_meta = H5Dcreate2(file, "MetaData", H5T_NATIVE_DOUBLE, space_meta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#if PRJ_MHD && PRJ_MHD_DEBUG
    space_bf = H5Screate_simple(3, dims_bf, dims_bf);
    dset_bf = H5Dcreate2(file, "Bf", dump_real_type, space_bf, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

    if (prj_io_is_root_rank(mpi) && mesh->nblocks > 0) {
        hsize_t start_meta[2] = {0, 0};
        hsize_t count_meta[2] = {(hsize_t)mesh->nblocks, (hsize_t)PRJ_IO_METADATA_SIZE};
        double *metadata = (double *)calloc((size_t)mesh->nblocks * PRJ_IO_METADATA_SIZE, sizeof(*metadata));

        if (metadata == 0) {
            prj_io_fail("prj_io_write_dump: metadata allocation failed");
        }
        for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
            prj_io_fill_metadata(&mesh->blocks[bidx], &metadata[(size_t)bidx * PRJ_IO_METADATA_SIZE]);
        }
        prj_io_write_hyperslab(dset_meta, mpi, H5T_NATIVE_DOUBLE, 2, start_meta, count_meta, metadata);
        free(metadata);
    }

    bidx = 0;
    while (bidx < mesh->nblocks) {
        int run_start;
        int run_len;
        hsize_t start_data[3];
        hsize_t count_data[3];
        prj_io_dump_real *buffer;
        int ridx;
        int v;
        int i;
        int j;
        int k;

        while (bidx < mesh->nblocks && !prj_io_restart_write_data_block(&mesh->blocks[bidx])) {
            bidx += 1;
        }
        run_start = bidx;
        while (bidx < mesh->nblocks && prj_io_restart_write_data_block(&mesh->blocks[bidx])) {
            bidx += 1;
        }
        run_len = bidx - run_start;
        if (run_len == 0) {
            continue;
        }
        start_data[0] = (hsize_t)run_start;
        start_data[1] = 0;
        start_data[2] = 0;
        count_data[0] = (hsize_t)run_len;
        count_data[1] = (hsize_t)PRJ_NVAR_PRIM;
        count_data[2] = (hsize_t)ncells;
        buffer = (prj_io_dump_real *)calloc((size_t)run_len * (size_t)PRJ_NVAR_PRIM * ncells, sizeof(*buffer));
        if (buffer == 0) {
            prj_io_fail("prj_io_write_dump: data allocation failed");
        }
        for (ridx = 0; ridx < run_len; ++ridx) {
            const prj_block *block = &mesh->blocks[run_start + ridx];

            for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                double wscale = 1.0;
#if PRJ_NRAD > 0
                if (v >= PRJ_NHYDRO) {
                    wscale = RAD_SCALE;
                }
#endif
                for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                    for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                        for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                            size_t cell = (size_t)i * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE + (size_t)j * PRJ_BLOCK_SIZE + (size_t)k;
                            size_t offset = ((size_t)ridx * (size_t)PRJ_NVAR_PRIM + (size_t)v) * ncells + cell;

                            buffer[offset] = (prj_io_dump_real)(block->W[WIDX(v, i, j, k)] * wscale);
                        }
                    }
                }
            }
        }
        prj_io_write_hyperslab(dset_data, mpi, dump_real_type, 3, start_data, count_data, buffer);
        free(buffer);
    }

    bidx = 0;
    while (bidx < mesh->nblocks) {
        int run_start;
        int run_len;
        hsize_t start_eos[3];
        hsize_t count_eos[3];
        prj_io_dump_real *buffer;
        int ridx;
        int v;
        int i;
        int j;
        int k;

        while (bidx < mesh->nblocks && !prj_io_dump_write_eos_block(&mesh->blocks[bidx])) {
            bidx += 1;
        }
        run_start = bidx;
        while (bidx < mesh->nblocks && prj_io_dump_write_eos_block(&mesh->blocks[bidx])) {
            bidx += 1;
        }
        run_len = bidx - run_start;
        if (run_len == 0) {
            continue;
        }
        start_eos[0] = (hsize_t)run_start;
        start_eos[1] = 0;
        start_eos[2] = 0;
        count_eos[0] = (hsize_t)run_len;
        count_eos[1] = (hsize_t)PRJ_NVAR_EOSVAR;
        count_eos[2] = (hsize_t)ncells;
        buffer = (prj_io_dump_real *)calloc((size_t)run_len * (size_t)PRJ_NVAR_EOSVAR * ncells, sizeof(*buffer));
        if (buffer == 0) {
            prj_io_fail("prj_io_write_dump: eos allocation failed");
        }
        for (ridx = 0; ridx < run_len; ++ridx) {
            const prj_block *block = &mesh->blocks[run_start + ridx];

            for (v = 0; v < PRJ_NVAR_EOSVAR; ++v) {
                for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                    for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                        for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                            size_t cell = (size_t)i * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE + (size_t)j * PRJ_BLOCK_SIZE + (size_t)k;
                            size_t offset = ((size_t)ridx * (size_t)PRJ_NVAR_EOSVAR + (size_t)v) * ncells + cell;

                            buffer[offset] = (prj_io_dump_real)block->eosvar[EIDX(v, i, j, k)];
                        }
                    }
                }
            }
        }
        prj_io_write_hyperslab(dset_eos, mpi, dump_real_type, 3, start_eos, count_eos, buffer);
        free(buffer);
    }

#if PRJ_MHD && PRJ_MHD_DEBUG
    bidx = 0;
    while (bidx < mesh->nblocks) {
        int run_start;
        int run_len;
        hsize_t start_bf[3];
        hsize_t count_bf[3];
        prj_io_dump_real *buffer;
        int ridx;
        int d;
        int n;

        while (bidx < mesh->nblocks && !prj_io_restart_write_bf_block(&mesh->blocks[bidx])) {
            bidx += 1;
        }
        run_start = bidx;
        while (bidx < mesh->nblocks && prj_io_restart_write_bf_block(&mesh->blocks[bidx])) {
            bidx += 1;
        }
        run_len = bidx - run_start;
        if (run_len == 0) {
            continue;
        }
        start_bf[0] = (hsize_t)run_start;
        start_bf[1] = 0;
        start_bf[2] = 0;
        count_bf[0] = (hsize_t)run_len;
        count_bf[1] = 3;
        count_bf[2] = (hsize_t)PRJ_BLOCK_NFACES;
        buffer = (prj_io_dump_real *)calloc((size_t)run_len * 3U * (size_t)PRJ_BLOCK_NFACES, sizeof(*buffer));
        if (buffer == 0) {
            prj_io_fail("prj_io_write_dump: Bf allocation failed");
        }
        for (ridx = 0; ridx < run_len; ++ridx) {
            const prj_block *block = &mesh->blocks[run_start + ridx];

            for (d = 0; d < 3; ++d) {
                for (n = 0; n < PRJ_BLOCK_NFACES; ++n) {
                    size_t offset = ((size_t)ridx * 3U + (size_t)d) * (size_t)PRJ_BLOCK_NFACES + (size_t)n;

                    buffer[offset] = (prj_io_dump_real)block->Bf[d][n];
                }
            }
        }
        prj_io_write_hyperslab(dset_bf, mpi, dump_real_type, 3, start_bf, count_bf, buffer);
        free(buffer);
    }
#endif

    {
        if (grav != 0 && grav->nbins > 0 && grav->phi != 0 && grav->lapse != 0) {
            hsize_t dims_phi[1] = {(hsize_t)grav->nbins + 1};
            hsize_t dims_lapse[1] = {(hsize_t)grav->nbins + 1};
            hid_t space_phi = H5Screate_simple(1, dims_phi, dims_phi);
            hid_t space_lapse = H5Screate_simple(1, dims_lapse, dims_lapse);
            hid_t dset_phi = H5Dcreate2(file, "phi", H5T_NATIVE_DOUBLE, space_phi,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            hid_t dset_lapse = H5Dcreate2(file, "lapse", H5T_NATIVE_DOUBLE, space_lapse,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            if (prj_io_is_root_rank(mpi)) {
                hid_t dxpl_phi = prj_io_data_xfer_plist(mpi);
                hid_t dxpl_lapse = prj_io_data_xfer_plist(mpi);

                H5Dwrite(dset_phi, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dxpl_phi, grav->phi);
                H5Dwrite(dset_lapse, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dxpl_lapse, grav->lapse);
                prj_io_close_dxpl(dxpl_phi);
                prj_io_close_dxpl(dxpl_lapse);
            }
            H5Dclose(dset_phi);
            H5Dclose(dset_lapse);
            H5Sclose(space_phi);
            H5Sclose(space_lapse);
        }
    }
#if PRJ_MHD && PRJ_MHD_DEBUG
    H5Dclose(dset_bf);
    H5Sclose(space_bf);
#endif
    H5Dclose(dset_meta);
    H5Sclose(space_meta);
    H5Dclose(dset_eos);
    H5Sclose(space_eos);
    H5Dclose(dset_data);
    H5Sclose(space_data);
    H5Fclose(file);
}
