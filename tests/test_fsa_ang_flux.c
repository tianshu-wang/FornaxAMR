#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

#if PRJ_USE_RADIATION_FSA
static void die(const char *msg)
{
    fprintf(stderr, "test_fsa_ang_flux: %s\n", msg);
    exit(1);
}

static void *xcalloc(size_t n, size_t size)
{
    void *ptr = calloc(n, size);

    if (ptr == 0) {
        die("allocation failed");
    }
    return ptr;
}

static void setup_block(prj_block *block)
{
    int d;

    memset(block, 0, sizeof(*block));
    block->dx[0] = 1.0;
    block->dx[1] = 1.0;
    block->dx[2] = 1.0;
    block->W_mhd = (double *)xcalloc((size_t)PRJ_NVAR_MHD_PRIM *
        (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_BLOCK_NCELLS, sizeof(*block->W_mhd));
    block->W_rad = (double *)xcalloc((size_t)PRJ_NVAR_RAD_PRIM *
        (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_BLOCK_NCELLS, sizeof(*block->W_rad));
#if PRJ_USE_RADIAL_FRAME_FSA
    block->ang_geom_fsa = (double *)xcalloc(3U * (size_t)PRJ_NARC *
        (size_t)PRJ_BLOCK_NCELLS, sizeof(*block->ang_geom_fsa));
#endif
    for (d = 0; d < 3; ++d) {
        block->grav[d] = (double *)xcalloc(PRJ_BLOCK_NCELLS, sizeof(*block->grav[d]));
    }
}

static void free_block(prj_block *block)
{
    int d;

    free(block->W_mhd);
    free(block->W_rad);
    block->W_mhd = 0;
    block->W_rad = 0;
#if PRJ_USE_RADIAL_FRAME_FSA
    free(block->ang_geom_fsa);
    block->ang_geom_fsa = 0;
#endif
    for (d = 0; d < 3; ++d) {
        free(block->grav[d]);
        block->grav[d] = 0;
    }
}

static double pattern_intensity(int field, int group, int angle, double scale)
{
    double val = 1.0 + 0.03 * (double)(field + 1) +
        0.005 * (double)(group + 1) + 0.0007 * (double)(angle + 1);

    return scale * val;
}

static void fill_state(const prj_rad *rad, prj_block *block,
    double u[PRJ_NVAR_CONS], double scale)
{
    int field;
    int group;
    int angle;

    memset(u, 0, (size_t)PRJ_NVAR_CONS * sizeof(*u));
    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int v = PRJ_CONS_RAD_I(field, group, angle);
                double I = pattern_intensity(field, group, angle, scale);
                double J = rad->solid_angle[angle] * I;

                block->W_rad[WIDX(PRJ_RAD_PRIM_I(field, group, angle), 0, 0, 0)] = J;
                u[v] = J;
            }
        }
    }
}

static double integrated_total(const double u[PRJ_NVAR_CONS])
{
    double total = 0.0;
    int field;
    int group;
    int angle;

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                total += u[PRJ_CONS_RAD_I(field, group, angle)];
            }
        }
    }
    return total;
}

static void set_gravity(prj_block *block, double ax, double ay, double az)
{
    int idx = IDX(0, 0, 0);

    block->grav[0][idx] = ax;
    block->grav[1][idx] = ay;
    block->grav[2][idx] = az;
}

static void check_zero_speed_noop(const prj_rad *rad, prj_block *block)
{
    double u[PRJ_NVAR_CONS];
    double before[PRJ_NVAR_CONS];
    int field;
    int group;
    int angle;

    set_gravity(block, 0.0, 0.0, 0.0);
    fill_state(rad, block, u, 1.0);
    memcpy(before, u, sizeof(u));

    prj_rad_ang_flux_apply(rad, block, block->W_rad, u, 0, 0, 0, 1.0, 0.25);

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int v = PRJ_CONS_RAD_I(field, group, angle);

                if (u[v] != before[v]) {
                    die("zero angular speed changed angular-cell energy");
                }
            }
        }
    }
}

static void check_conservation(const prj_rad *rad, prj_block *block)
{
    double u[PRJ_NVAR_CONS];
    double total_before;
    double total_after;
    double max_delta = 0.0;
    int field;
    int group;
    int angle;

    set_gravity(block, 0.6 * PRJ_CLIGHT, -0.25 * PRJ_CLIGHT, 0.4 * PRJ_CLIGHT);
    fill_state(rad, block, u, 1.0);
    total_before = integrated_total(u);

    prj_rad_ang_flux_apply(rad, block, block->W_rad, u, 0, 0, 0, 1.0, 1.0e-4);

    total_after = integrated_total(u);
    if (fabs(total_after - total_before) >
        1.0e-12 * fmax(fabs(total_before), 1.0)) {
        die("angular flux did not conserve angular-cell-integrated energy");
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int v = PRJ_CONS_RAD_I(field, group, angle);
                double ref = rad->solid_angle[angle] *
                    pattern_intensity(field, group, angle, 1.0);
                double delta = fabs(u[v] - ref);

                if (delta > max_delta) {
                    max_delta = delta;
                }
            }
        }
    }
    if (max_delta <= 0.0) {
        die("nonzero angular speed produced no update");
    }
}

static void check_positivity_limiter(const prj_rad *rad, prj_block *block)
{
    double u[PRJ_NVAR_CONS];
    double total_before;
    double total_after;
    double max_delta = 0.0;
    int field;
    int group;
    int angle;

    set_gravity(block, 80.0 * PRJ_CLIGHT, -55.0 * PRJ_CLIGHT, 35.0 * PRJ_CLIGHT);
    fill_state(rad, block, u, 1.0e-3);
    total_before = integrated_total(u);

    prj_rad_ang_flux_apply(rad, block, block->W_rad, u, 0, 0, 0, 1.0, 1.0);

    total_after = integrated_total(u);
    if (fabs(total_after - total_before) >
        1.0e-10 * fmax(fabs(total_before), 1.0)) {
        die("limited angular flux did not conserve angular-cell-integrated energy");
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                int v = PRJ_CONS_RAD_I(field, group, angle);
                double ref = rad->solid_angle[angle] *
                    pattern_intensity(field, group, angle, 1.0e-3);
                double delta = fabs(u[v] - ref);

                if (u[v] < -1.0e-14) {
                    die("angular flux limiter left a negative angular-cell energy");
                }
                if (delta > max_delta) {
                    max_delta = delta;
                }
            }
        }
    }
    if (max_delta <= 0.0) {
        die("large-step angular flux produced no update");
    }
}

static void check_geometry_flux_direction(const prj_rad *rad, prj_block *block)
{
#if PRJ_USE_RADIAL_FRAME_FSA
    double u[PRJ_NVAR_CONS];
    double before[PRJ_NVAR_CONS];
    double total_before;
    double total_after;
    const int arc = 0;
    int c0 = rad->arc_neighbor[2 * arc];
    int c1 = rad->arc_neighbor[2 * arc + 1];
    double speed = 0.35;
    int d;

    set_gravity(block, 0.0, 0.0, 0.0);
    fill_state(rad, block, u, 1.0);
    memcpy(before, u, sizeof(u));
    total_before = integrated_total(u);

    for (d = 0; d < 3; ++d) {
        block->ang_geom_fsa[PRJ_FSA_ANG_GEOM_IDX(arc, d, 0, 0, 0)] =
            -speed * rad->arc_vec[3 * arc + d];
    }

    prj_rad_ang_flux_apply(rad, block, block->W_rad, u, 0, 0, 0, 1.0, 1.0e-4);

    total_after = integrated_total(u);
    if (fabs(total_after - total_before) >
        1.0e-12 * fmax(fabs(total_before), 1.0)) {
        die("geometry angular flux did not conserve angular-cell-integrated energy");
    }
    if (u[PRJ_CONS_RAD_I(0, 0, c0)] <= before[PRJ_CONS_RAD_I(0, 0, c0)] ||
        u[PRJ_CONS_RAD_I(0, 0, c1)] >= before[PRJ_CONS_RAD_I(0, 0, c1)]) {
        die("geometry angular flux used the wrong donor direction");
    }
#else
    (void)rad;
    (void)block;
#endif
}

int main(void)
{
    prj_rad rad;
    prj_block block;

    memset(&rad, 0, sizeof(rad));
    prj_rad_fsa_calculate_directions(&rad);
    setup_block(&block);

    check_zero_speed_noop(&rad, &block);
    check_conservation(&rad, &block);
    check_positivity_limiter(&rad, &block);
    check_geometry_flux_direction(&rad, &block);

    free_block(&block);
    prj_rad_fsa_free_geometry(&rad);

    printf("test_fsa_ang_flux: ok (N_ANGLE_LEV=%d, cells=%d, arcs=%d)\n",
        PRJ_N_ANGLE_LEV, PRJ_NANGLE, PRJ_NARC);
    return 0;
}
#else
int main(void)
{
    printf("test_fsa_ang_flux: skipped (RADIATION_FSA=0)\n");
    return 0;
}
#endif
