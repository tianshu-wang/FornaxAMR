#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

static void die(const char *msg)
{
    fprintf(stderr, "test_z4c: %s\n", msg);
    exit(1);
}

static void check_enum_names(void)
{
    static const char *const expected[PRJ_NZ4C] = {
        "z4c_chi",
        "z4c_gxx", "z4c_gxy", "z4c_gxz", "z4c_gyy", "z4c_gyz", "z4c_gzz",
        "z4c_Khat",
        "z4c_Axx", "z4c_Axy", "z4c_Axz", "z4c_Ayy", "z4c_Ayz", "z4c_Azz",
        "z4c_Gamx", "z4c_Gamy", "z4c_Gamz",
        "z4c_Theta",
        "z4c_alpha",
        "z4c_betax", "z4c_betay", "z4c_betaz"
    };
    static const char *const expected_tmunu[PRJ_NTMUNU] = {
        "z4c_tmunu_Sxx", "z4c_tmunu_Sxy", "z4c_tmunu_Sxz",
        "z4c_tmunu_Syy", "z4c_tmunu_Syz", "z4c_tmunu_Szz",
        "z4c_tmunu_E", "z4c_tmunu_Sx", "z4c_tmunu_Sy", "z4c_tmunu_Sz"
    };
    int v;

    if (PRJ_NZ4C != 22) {
        die("unexpected PRJ_NZ4C");
    }
    for (v = 0; v < PRJ_NZ4C; ++v) {
        if (strcmp(prj_z4c_var_name(v), expected[v]) != 0) {
            die("unexpected Z4c variable name/order");
        }
    }
    if (PRJ_NTMUNU != 10) {
        die("unexpected PRJ_NTMUNU");
    }
    for (v = 0; v < PRJ_NTMUNU; ++v) {
        if (strcmp(prj_z4c_tmunu_var_name(v), expected_tmunu[v]) != 0) {
            die("unexpected Tmunu variable name/order");
        }
    }
}

#if PRJ_DYNAMIC_GR
static void assert_close(const char *name, double got, double expected, double tol)
{
    if (!isfinite(got) || fabs(got - expected) > tol) {
        fprintf(stderr, "test_z4c: %s got %.17e expected %.17e tol %.3e\n",
            name, got, expected, tol);
        exit(1);
    }
}

static void assert_close_rel(const char *name, double got, double expected, double rel)
{
    double scale = expected != 0.0 ? fabs(expected) : 1.0;
    double tol = rel * scale;

    if (tol < 1.0e-36) {
        tol = 1.0e-36;
    }
    assert_close(name, got, expected, tol);
}

static double geo_factor(void)
{
    double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    return PRJ_GNEWT / (c2 * c2);
}

static void disable_z4c_damping(prj_mesh *mesh)
{
    mesh->z4c_params.damp_kappa1_inv_cm = 0.0;
    mesh->z4c_params.shift_eta_inv_cm = 0.0;
    mesh->z4c_params.ssl_damping_amp_inv_cm = 0.0;
}

static double flat_gauge_denom_inv_cm(void)
{
    return 3.0 * sqrt(2.0) / 3.0e10;
}

static void set_constant_diagonal_z4c_metric(prj_block *block,
    double chi, double gxx, double gyy, double gzz)
{
    double *z;
    int i;
    int j;
    int k;

    if (block == 0) {
        die("missing block for metric setup");
    }
    z = prj_block_z4c_stage(block, 0);
    if (z == 0) {
        die("missing z4c storage for metric setup");
    }
    for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
        for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
            for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)] = chi;
                z[Z4CIDX(PRJ_Z4C_GXX, i, j, k)] = gxx;
                z[Z4CIDX(PRJ_Z4C_GXY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GXZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GYY, i, j, k)] = gyy;
                z[Z4CIDX(PRJ_Z4C_GYZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GZZ, i, j, k)] = gzz;
                z[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AXX, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AXY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AXZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AYY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AYZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_AZZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GAMX, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GAMY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_GAMZ, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_THETA, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_ALPHA, i, j, k)] = 1.0;
                z[Z4CIDX(PRJ_Z4C_BETAX, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_BETAY, i, j, k)] = 0.0;
                z[Z4CIDX(PRJ_Z4C_BETAZ, i, j, k)] = 0.0;
            }
        }
    }
}

static void init_one_block_mesh(prj_mesh *mesh, prj_coord *coord)
{
    memset(mesh, 0, sizeof(*mesh));
    prj_z4c_init_params(&mesh->z4c_params);
    mesh->use_full_dynamic_gr = 0;
    coord->x1min = 0.0;
    coord->x1max = (double)PRJ_BLOCK_SIZE * 3.0e10;
    coord->x2min = 0.0;
    coord->x2max = (double)PRJ_BLOCK_SIZE * 3.0e10;
    coord->x3min = 0.0;
    coord->x3max = (double)PRJ_BLOCK_SIZE * 3.0e10;
    if (prj_mesh_init(mesh, 1, 1, 1, 0, coord, 0) != 0) {
        die("mesh init failed");
    }
    if (mesh->blocks[0].W_mhd != 0) {
        prj_fill(mesh->blocks[0].W_mhd, (size_t)PRJ_NVAR_MHD_PRIM *
            (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_BLOCK_NCELLS, 0.0);
    }
    if (mesh->blocks[0].W_rad != 0) {
        prj_fill(mesh->blocks[0].W_rad, (size_t)PRJ_NVAR_RAD_PRIM *
            (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_BLOCK_NCELLS, 0.0);
    }
    if (mesh->blocks[0].eosvar != 0) {
        prj_fill(mesh->blocks[0].eosvar, (size_t)PRJ_NVAR_EOSVAR *
            (size_t)PRJ_BLOCK_NCELLS, 0.0);
    }
}

static void check_flat_state(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    double *z;
    double *rhs;
    double dt;
    int i, j, k, v;

    init_one_block_mesh(&mesh, &coord);
    block = &mesh.blocks[0];
    if (block->z4c == 0 || block->z4c_rhs == 0) {
        die("Z4c block storage missing");
    }

    prj_z4c_init_mesh_flat(&mesh, 0);
    disable_z4c_damping(&mesh);
    dt = prj_z4c_calc_dt_seconds(&mesh, 0, 0.5);
    assert_close("dt_z4c", dt,
        0.5 * 3.0e10 / (3.0 * sqrt(2.0) * PRJ_CLIGHT), 1.0e-12);

    z = prj_block_z4c_stage(block, 0);
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                assert_close("chi", z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)], 1.0, 0.0);
                assert_close("gxx", z[Z4CIDX(PRJ_Z4C_GXX, i, j, k)], 1.0, 0.0);
                assert_close("gyy", z[Z4CIDX(PRJ_Z4C_GYY, i, j, k)], 1.0, 0.0);
                assert_close("gzz", z[Z4CIDX(PRJ_Z4C_GZZ, i, j, k)], 1.0, 0.0);
                assert_close("alpha", z[Z4CIDX(PRJ_Z4C_ALPHA, i, j, k)], 1.0, 0.0);
                for (v = 0; v < PRJ_NZ4C; ++v) {
                    if (v != PRJ_Z4C_CHI && v != PRJ_Z4C_GXX &&
                        v != PRJ_Z4C_GYY && v != PRJ_Z4C_GZZ &&
                        v != PRJ_Z4C_ALPHA) {
                        assert_close("flat zero", z[Z4CIDX(v, i, j, k)], 0.0, 0.0);
                    }
                }
            }
        }
    }

    prj_z4c_compute_rhs(&mesh, 0, 0, 0, 0, 0.0);
    rhs = prj_block_z4c_rhs_stage(block, 0);
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                for (v = 0; v < PRJ_NZ4C; ++v) {
                    assert_close("flat rhs", rhs[Z4CIDX(v, i, j, k)], 0.0, 1.0e-30);
                }
            }
        }
    }
    prj_mesh_destroy(&mesh);
}

static void check_dt_includes_damping_terms(void)
{
    prj_mesh mesh;
    prj_coord coord;
    double cfl = 0.5;
    double gauge_denom = flat_gauge_denom_inv_cm();
    double dt;

    init_one_block_mesh(&mesh, &coord);
    prj_z4c_init_mesh_flat(&mesh, 0);

    disable_z4c_damping(&mesh);
    mesh.z4c_params.damp_kappa1_inv_cm = 1.0e-5;
    mesh.z4c_params.damp_kappa2 = 0.0;
    dt = prj_z4c_calc_dt_seconds(&mesh, 0, cfl);
    assert_close_rel("dt_z4c kappa damping", dt,
        cfl / (PRJ_CLIGHT * (gauge_denom + 2.0 * mesh.z4c_params.damp_kappa1_inv_cm)),
        1.0e-12);

    disable_z4c_damping(&mesh);
    mesh.z4c_params.shift_eta_inv_cm = 7.0e-6;
    dt = prj_z4c_calc_dt_seconds(&mesh, 0, cfl);
    assert_close_rel("dt_z4c shift damping", dt,
        cfl / (PRJ_CLIGHT * (gauge_denom + mesh.z4c_params.shift_eta_inv_cm)),
        1.0e-12);

    disable_z4c_damping(&mesh);
    mesh.z4c_params.slow_start_lapse = 1;
    mesh.z4c_params.ssl_damping_amp_inv_cm = 3.0e-6;
    mesh.z4c_params.ssl_damping_index = 1;
    mesh.time_seconds = 0.0;
    dt = prj_z4c_calc_dt_seconds(&mesh, 0, cfl);
    assert_close_rel("dt_z4c lapse damping", dt,
        cfl / (PRJ_CLIGHT * (gauge_denom + mesh.z4c_params.ssl_damping_amp_inv_cm)),
        1.0e-12);

    prj_mesh_destroy(&mesh);
}

static double z4c_poly_axis(double x)
{
    return (((((((0.071 * x - 0.043) * x + 0.019) * x - 0.031) * x +
        0.053) * x - 0.079) * x + 0.097) * x - 0.113);
}

static double z4c_poly_pattern(int var, double i, double j, double k, int stage)
{
    double x = 0.125 * i;
    double y = 0.125 * j;
    double z = 0.125 * k;
    double cross = x * y * z + 0.25 * x * x * y * y * z -
        0.125 * y * z * z * z;

    return 100.0 * (double)(stage + 1) + 10.0 * (double)(var + 1) +
        (1.0 + 0.01 * (double)(var + 1)) *
        (z4c_poly_axis(x) + 0.75 * z4c_poly_axis(y + 0.125) -
        0.5 * z4c_poly_axis(z - 0.25)) + 0.0625 * cross;
}

static void fill_z4c_polynomial(prj_block *block)
{
    int stage, var, i, j, k;

    for (stage = 0; stage < PRJ_BLOCK_NSTAGES; ++stage) {
        double *z = prj_block_z4c_stage(block, stage);

        for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
            for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
                for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                    for (var = 0; var < PRJ_NZ4C; ++var) {
                        z[Z4CIDX(var, i, j, k)] =
                            z4c_poly_pattern(var, (double)i, (double)j, (double)k, stage);
                    }
                }
            }
        }
    }
}

static void init_allocated_test_block(prj_block *block, int id)
{
    int n;

    memset(block, 0, sizeof(*block));
    block->id = id;
    block->active = 1;
    block->rank = 0;
    block->parent = -1;
    for (n = 0; n < 8; ++n) {
        block->children[n] = -1;
    }
    for (n = 0; n < 56; ++n) {
        block->slot[n].id = -1;
        block->slot[n].rank = 0;
        block->slot[n].rel_level = 0;
        block->slot[n].type = PRJ_NEIGHBOR_NONE;
    }
    if (prj_block_alloc_data(block) != 0) {
        die("test block allocation failed");
    }
}

static void check_z4c_aux_cleared(const prj_block *block)
{
    size_t n;

    for (n = 0; n < (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_NZ4C *
        (size_t)PRJ_BLOCK_NCELLS_Z4C; ++n) {
        if (block->z4c_rhs[n] != 0.0) {
            die("Z4c AMR transfer did not clear rhs");
        }
    }
}

static void check_z4c_amr_transfer(void)
{
    prj_block parent;
    prj_block child;
    prj_block children_storage[8];
    const prj_block *children[8];
    const int oct = 5;
    const int restrict_var = PRJ_Z4C_THETA;
    const int restrict_stage = 1;
    double *z_child;
    double expected;
    int i, j, k, var, stage;
    int o;

    init_allocated_test_block(&parent, 100);
    init_allocated_test_block(&child, 101);
    fill_z4c_polynomial(&parent);
    prj_fill(child.z4c_rhs, (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_NZ4C *
        (size_t)PRJ_BLOCK_NCELLS_Z4C, 42.0);
    prj_z4c_amr_prolongate_child(&parent, &child, oct);
    for (stage = 0; stage < PRJ_BLOCK_NSTAGES; ++stage) {
        z_child = prj_block_z4c_stage(&child, stage);
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            int gi = ((oct & 1) ? PRJ_BLOCK_SIZE : 0) + i;
            double x = (double)(gi / 2) + ((gi & 1) ? 0.25 : -0.25);

            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                int gj = ((oct & 2) ? PRJ_BLOCK_SIZE : 0) + j;
                double y = (double)(gj / 2) + ((gj & 1) ? 0.25 : -0.25);

                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    int gk = ((oct & 4) ? PRJ_BLOCK_SIZE : 0) + k;
                    double z = (double)(gk / 2) + ((gk & 1) ? 0.25 : -0.25);

                    for (var = 0; var < PRJ_NZ4C; ++var) {
                        expected = z4c_poly_pattern(var, x, y, z, stage);
                        assert_close("z4c prolong degree-7",
                            z_child[Z4CIDX(var, i, j, k)], expected, 1.0e-10);
                    }
                }
            }
        }
    }
    check_z4c_aux_cleared(&child);

    prj_fill(parent.z4c_rhs, (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_NZ4C *
        (size_t)PRJ_BLOCK_NCELLS_Z4C, 44.0);
    for (o = 0; o < 8; ++o) {
        int st, v, ii, jj, kk;

        init_allocated_test_block(&children_storage[o], 200 + o);
        children[o] = &children_storage[o];
        for (st = 0; st < PRJ_BLOCK_NSTAGES; ++st) {
            double *z = prj_block_z4c_stage(&children_storage[o], st);

            for (ii = -PRJ_NGHOST_Z4C; ii < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++ii) {
                for (jj = -PRJ_NGHOST_Z4C; jj < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++jj) {
                    for (kk = -PRJ_NGHOST_Z4C; kk < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++kk) {
                        for (v = 0; v < PRJ_NZ4C; ++v) {
                            z[Z4CIDX(v, ii, jj, kk)] =
                                1000.0 + 100.0 * (double)o + 10.0 * (double)st +
                                0.01 * (double)v;
                        }
                    }
                }
            }
        }
    }
    prj_z4c_amr_restrict_parent(children, &parent);
    for (o = 0; o < 8; ++o) {
        int ioff = (o & 1) ? PRJ_BLOCK_SIZE / 2 : 0;
        int joff = (o & 2) ? PRJ_BLOCK_SIZE / 2 : 0;
        int koff = (o & 4) ? PRJ_BLOCK_SIZE / 2 : 0;
        double *z_parent = prj_block_z4c_stage(&parent, restrict_stage);

        expected = 1000.0 + 100.0 * (double)o + 10.0 * (double)restrict_stage +
            0.01 * (double)restrict_var;
        assert_close("z4c restrict",
            z_parent[Z4CIDX(restrict_var, ioff + 1, joff + 1, koff + 1)],
            expected, 1.0e-12);
    }
    check_z4c_aux_cleared(&parent);
    for (o = 0; o < 8; ++o) {
        prj_block_free_data(&children_storage[o]);
    }
    prj_block_free_data(&child);
    prj_block_free_data(&parent);
}

static int z4c_test_fine_index_odd(int index)
{
    return (index % 2) != 0;
}

static double z4c_slot_parent_coord(const prj_neighbor *slot, int axis, int offset)
{
    int fine_index = offset + slot->recv_loc_start_z4c[axis];
    int parent_index = offset / 2 + slot->send_loc_start_z4c[axis];

    return (double)parent_index + (z4c_test_fine_index_odd(fine_index) ? 0.25 : -0.25);
}

static void set_test_block_box(prj_block *block, int level,
    double x0, double x1, double y0, double y1, double z0, double z1)
{
    block->level = level;
    block->xmin[0] = x0;
    block->xmax[0] = x1;
    block->xmin[1] = y0;
    block->xmax[1] = y1;
    block->xmin[2] = z0;
    block->xmax[2] = z1;
    block->dx[0] = (x1 - x0) / (double)PRJ_BLOCK_SIZE;
    block->dx[1] = (y1 - y0) / (double)PRJ_BLOCK_SIZE;
    block->dx[2] = (z1 - z0) / (double)PRJ_BLOCK_SIZE;
}

static void init_z4c_two_block_mesh(prj_mesh *mesh, int fine_on_high_side)
{
    prj_block *coarse;
    prj_block *fine;
    prj_neighbor *slot;
    int d;

    memset(mesh, 0, sizeof(*mesh));
    prj_z4c_init_params(&mesh->z4c_params);
    mesh->nblocks = 2;
    mesh->nblocks_max = 2;
    mesh->max_blocks = 2;
    mesh->max_level = 1;
    mesh->blocks = (prj_block *)calloc(2U, sizeof(*mesh->blocks));
    if (mesh->blocks == 0) {
        die("two-block mesh allocation failed");
    }

    coarse = &mesh->blocks[0];
    fine = &mesh->blocks[1];
    init_allocated_test_block(coarse, 0);
    init_allocated_test_block(fine, 1);
    set_test_block_box(coarse, 0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    if (fine_on_high_side) {
        set_test_block_box(fine, 1, 1.0, 1.5, 0.0, 0.5, 0.0, 0.5);
    } else {
        set_test_block_box(fine, 1, -0.5, 0.0, 0.0, 0.5, 0.0, 0.5);
    }

    slot = &coarse->slot[0];
    prj_neighbor_compute_geometry(coarse, fine, slot);
    slot->id = 1;
    slot->rank = 0;
    for (d = 0; d < 3; ++d) {
        slot->xmin[d] = fine->xmin[d];
        slot->xmax[d] = fine->xmax[d];
        slot->dx[d] = fine->dx[d];
    }
    if (slot->type != PRJ_NEIGHBOR_FACE || slot->rel_level != 1) {
        die("unexpected two-block Z4c neighbor geometry");
    }
}

static void destroy_z4c_two_block_mesh(prj_mesh *mesh)
{
    int b;

    if (mesh->blocks != 0) {
        for (b = 0; b < mesh->nblocks; ++b) {
            prj_block_free_data(&mesh->blocks[b]);
        }
        free(mesh->blocks);
    }
    mesh->blocks = 0;
    mesh->nblocks = 0;
}

static void check_z4c_one_prolong_ghost_side(int fine_on_high_side)
{
    prj_mesh mesh;
    prj_block *coarse;
    prj_block *fine;
    const prj_neighbor *slot;
    const double sentinel = -9.87654321e90;
    int stage = PRJ_BLOCK_NSTAGES > 1 ? 1 : 0;
    int ni, nj, nk;
    int i, j, k, var;
    int checked = 0;

    init_z4c_two_block_mesh(&mesh, fine_on_high_side);
    coarse = &mesh.blocks[0];
    fine = &mesh.blocks[1];
    slot = &coarse->slot[0];
    fill_z4c_polynomial(coarse);
    prj_fill(prj_block_z4c_stage(fine, stage),
        (size_t)PRJ_NZ4C * (size_t)PRJ_BLOCK_NCELLS_Z4C, sentinel);
    prj_z4c_fill_ghosts(&mesh, 0, 0, stage);

    if (fine_on_high_side && slot->recv_loc_start_z4c[0] >= 0) {
        die("expected negative fine ghost receive range");
    }
    if (!fine_on_high_side && slot->recv_loc_start_z4c[0] != PRJ_BLOCK_SIZE) {
        die("expected high-side fine ghost receive range");
    }

    ni = slot->recv_loc_end_z4c[0] - slot->recv_loc_start_z4c[0];
    nj = slot->recv_loc_end_z4c[1] - slot->recv_loc_start_z4c[1];
    nk = slot->recv_loc_end_z4c[2] - slot->recv_loc_start_z4c[2];
    for (i = 0; i < ni; ++i) {
        int di = i + slot->recv_loc_start_z4c[0];
        double x = z4c_slot_parent_coord(slot, 0, i);

        for (j = 0; j < nj; ++j) {
            int dj = j + slot->recv_loc_start_z4c[1];
            double y = z4c_slot_parent_coord(slot, 1, j);

            for (k = 0; k < nk; ++k) {
                int dk = k + slot->recv_loc_start_z4c[2];
                double z = z4c_slot_parent_coord(slot, 2, k);
                double *dst;

                if (di < -PRJ_NGHOST_Z4C || di >= PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C ||
                    dj < -PRJ_NGHOST_Z4C || dj >= PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C ||
                    dk < -PRJ_NGHOST_Z4C || dk >= PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C) {
                    continue;
                }
                dst = prj_block_z4c_stage(fine, stage);
                for (var = 0; var < PRJ_NZ4C; ++var) {
                    double expected = z4c_poly_pattern(var, x, y, z, stage);

                    assert_close("z4c prolong ghost degree-7",
                        dst[Z4CIDX(var, di, dj, dk)], expected, 1.0e-10);
                    ++checked;
                }
            }
        }
    }
    if (checked == 0) {
        die("Z4c prolong ghost test checked no cells");
    }
    destroy_z4c_two_block_mesh(&mesh);
}

static void check_z4c_amr_prolong_ghost_fill(void)
{
    check_z4c_one_prolong_ghost_side(1);
    check_z4c_one_prolong_ghost_side(0);
}

static void set_z4c_chi_step(prj_block *block, double jump)
{
    double *z = prj_block_z4c_stage(block, 0);
    int i;
    int j;
    int k;

    for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
        for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
            for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                z[Z4CIDX(PRJ_Z4C_CHI, i, j, k)] = 1.0 + (i >= PRJ_BLOCK_SIZE / 2 ? jump : 0.0);
            }
        }
    }
}

static void configure_z4c_dchi_amr(prj_mesh *mesh, double refine_thresh)
{
    int n;

    for (n = 0; n < PRJ_AMR_N; ++n) {
        mesh->amr_criterion_set[n] = 0;
        mesh->amr_refine_thresh[n] = 0.0;
        mesh->amr_derefine_thresh[n] = 0.0;
    }
    mesh->amr_estimator[0] = PRJ_AMR_ESTIMATOR_Z4C_DCHI;
    mesh->amr_refine_thresh[0] = refine_thresh;
    mesh->amr_derefine_thresh[0] = 0.0;
    mesh->amr_criterion_set[0] = 1;
}

static void check_z4c_dchi_amr_tag(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;

    init_one_block_mesh(&mesh, &coord);
    mesh.max_level = 1;
    block = &mesh.blocks[0];
    prj_z4c_init_mesh_flat(&mesh, 0);
    configure_z4c_dchi_amr(&mesh, 0.01);
    if (prj_amr_criteria_need_eosvar(&mesh) != 0) {
        die("z4c_dchi AMR should not require EOS variables");
    }

    prj_amr_tag(&mesh, 0, 0);
    if (block->refine_flag != 0) {
        die("constant chi should not refine");
    }

    set_z4c_chi_step(block, 0.2);
    configure_z4c_dchi_amr(&mesh, 0.05);
    prj_amr_tag(&mesh, 0, 0);
    if (block->refine_flag != 1) {
        die("z4c_dchi jump should refine");
    }

    configure_z4c_dchi_amr(&mesh, 0.5);
    prj_amr_tag(&mesh, 0, 0);
    if (block->refine_flag != 0) {
        die("z4c_dchi should respect high refine threshold");
    }

    prj_mesh_destroy(&mesh);
}

static void set_linear_sommerfeld_var(prj_block *block, int var, const double coeff[4])
{
    double *z = prj_block_z4c_stage(block, 0);
    int i, j, k;

    for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
        for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
            for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                double x = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                double y = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                double zc = block->xmin[2] + ((double)k + 0.5) * block->dx[2];

                z[Z4CIDX(var, i, j, k)] =
                    coeff[0] + coeff[1] * x + coeff[2] * y + coeff[3] * zc;
            }
        }
    }
}

static double expected_sommerfeld(const prj_block *block, int i, int j, int k,
    const double coeff[4], double speed)
{
    double x = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
    double y = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
    double z = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
    double r = sqrt(x * x + y * y + z * z);
    double value = coeff[0] + coeff[1] * x + coeff[2] * y + coeff[3] * z;
    double adv = (x * coeff[1] + y * coeff[2] + z * coeff[3]) / r;

    return -speed * value / r - speed * adv;
}

static prj_bc make_uniform_bc(int type)
{
    prj_bc bc;

    bc.bc_x1_inner = type;
    bc.bc_x1_outer = type;
    bc.bc_x2_inner = type;
    bc.bc_x2_outer = type;
    bc.bc_x3_inner = type;
    bc.bc_x3_outer = type;
    return bc;
}

static void check_z4c_sommerfeld_rhs(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    double *rhs;
    prj_bc bc;
    const int i = 0, j = 2, k = 3;
    const double theta[4] = {2.0e-6, 1.0e-24, -2.0e-24, 3.0e-24};
    const double khat[4] = {-3.0e-6, 2.0e-24, 1.0e-24, -1.0e-24};
    const double gamx[4] = {4.0e-6, -1.0e-24, 4.0e-24, 2.0e-24};
    const double axx[4] = {5.0e-6, 3.0e-24, -1.0e-24, 2.0e-24};
    const double sentinel = 9.0e99;

    init_one_block_mesh(&mesh, &coord);
    block = &mesh.blocks[0];
    prj_z4c_init_mesh_flat(&mesh, 0);
    set_linear_sommerfeld_var(block, PRJ_Z4C_THETA, theta);
    set_linear_sommerfeld_var(block, PRJ_Z4C_KHAT, khat);
    set_linear_sommerfeld_var(block, PRJ_Z4C_GAMX, gamx);
    set_linear_sommerfeld_var(block, PRJ_Z4C_AXX, axx);
    rhs = prj_block_z4c_rhs_stage(block, 0);

    bc = make_uniform_bc(PRJ_BC_REFLECT);
    prj_fill(rhs, (size_t)PRJ_NZ4C * (size_t)PRJ_BLOCK_NCELLS_Z4C, sentinel);
    prj_z4c_apply_sommerfeld_rhs(&mesh, 0, &bc, 0, 0);
    assert_close("sommerfeld reflect",
        rhs[Z4CIDX(PRJ_Z4C_THETA, i, j, k)], sentinel, 0.0);

    bc = make_uniform_bc(PRJ_BC_REFLECT);
    bc.bc_x1_inner = PRJ_BC_OUTFLOW;
    prj_fill(rhs, (size_t)PRJ_NZ4C * (size_t)PRJ_BLOCK_NCELLS_Z4C, sentinel);
    prj_z4c_apply_sommerfeld_rhs(&mesh, 0, &bc, 0, 0);
    assert_close_rel("sommerfeld theta",
        rhs[Z4CIDX(PRJ_Z4C_THETA, i, j, k)],
        expected_sommerfeld(block, i, j, k, theta, 1.0), 1.0e-12);
    assert_close_rel("sommerfeld khat",
        rhs[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)],
        expected_sommerfeld(block, i, j, k, khat, sqrt(2.0)), 1.0e-12);
    assert_close_rel("sommerfeld gamx",
        rhs[Z4CIDX(PRJ_Z4C_GAMX, i, j, k)],
        expected_sommerfeld(block, i, j, k, gamx, 1.0), 1.0e-12);
    assert_close_rel("sommerfeld axx",
        rhs[Z4CIDX(PRJ_Z4C_AXX, i, j, k)],
        expected_sommerfeld(block, i, j, k, axx, 1.0), 1.0e-12);
    assert_close("sommerfeld leaves alpha",
        rhs[Z4CIDX(PRJ_Z4C_ALPHA, i, j, k)], sentinel, 0.0);
    assert_close("sommerfeld leaves interior",
        rhs[Z4CIDX(PRJ_Z4C_THETA, 1, j, k)], sentinel, 0.0);

    bc.bc_x1_inner = PRJ_BC_USER;
    mesh.z4c_params.user_Sbc = 0;
    prj_fill(rhs, (size_t)PRJ_NZ4C * (size_t)PRJ_BLOCK_NCELLS_Z4C, sentinel);
    prj_z4c_apply_sommerfeld_rhs(&mesh, 0, &bc, 0, 0);
    assert_close("sommerfeld user gated",
        rhs[Z4CIDX(PRJ_Z4C_THETA, i, j, k)], sentinel, 0.0);
    mesh.z4c_params.user_Sbc = 1;
    prj_fill(rhs, (size_t)PRJ_NZ4C * (size_t)PRJ_BLOCK_NCELLS_Z4C, sentinel);
    prj_z4c_apply_sommerfeld_rhs(&mesh, 0, &bc, 0, 0);
    assert_close_rel("sommerfeld user enabled",
        rhs[Z4CIDX(PRJ_Z4C_THETA, i, j, k)],
        expected_sommerfeld(block, i, j, k, theta, 1.0), 1.0e-12);
    prj_mesh_destroy(&mesh);
}

static void check_rhs_hydro_matter_projection(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    double *W;
    double *rhs;
    const int i = 1, j = 2, k = 3;
    const double pi = acos(-1.0);
    const double rho = 2.0;
    const double v1 = 3.0;
    const double v2 = 4.0;
    const double v3 = 5.0;
    const double eint = 7.0;
    const double pressure = 11.0;
    double fac = geo_factor();
    double B1 = 0.0;
    double B2 = 0.0;
    double B3 = 0.0;
    double Bmag2 = 0.0;
    double E;
    double Sx;
    double Sxx;
    double Sxy;
    double Syy;
    double Szz;
    double S;

    init_one_block_mesh(&mesh, &coord);
    block = &mesh.blocks[0];
    prj_z4c_init_mesh_flat(&mesh, 0);
    disable_z4c_damping(&mesh);
    W = prj_block_mhd_stage(block, 0);
    W[WIDX(PRJ_PRIM_RHO, i, j, k)] = rho;
    W[WIDX(PRJ_PRIM_V1, i, j, k)] = v1;
    W[WIDX(PRJ_PRIM_V2, i, j, k)] = v2;
    W[WIDX(PRJ_PRIM_V3, i, j, k)] = v3;
    W[WIDX(PRJ_PRIM_EINT, i, j, k)] = eint;
    block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)] = pressure;
#if PRJ_MHD
    B1 = 13.0;
    B2 = 17.0;
    B3 = 19.0;
    W[WIDX(PRJ_PRIM_B1, i, j, k)] = B1;
    W[WIDX(PRJ_PRIM_B2, i, j, k)] = B2;
    W[WIDX(PRJ_PRIM_B3, i, j, k)] = B3;
    Bmag2 = B1 * B1 + B2 * B2 + B3 * B3;
#endif

    E = fac * (rho * PRJ_CLIGHT * PRJ_CLIGHT + rho * eint +
        0.5 * rho * (v1 * v1 + v2 * v2 + v3 * v3) + 0.5 * Bmag2);
    Sx = fac * rho * v1 * PRJ_CLIGHT;
    Sxx = fac * (rho * v1 * v1 + pressure + 0.5 * Bmag2 - B1 * B1);
    Sxy = fac * (rho * v1 * v2 - B1 * B2);
    Syy = fac * (rho * v2 * v2 + pressure + 0.5 * Bmag2 - B2 * B2);
    Szz = fac * (rho * v3 * v3 + pressure + 0.5 * Bmag2 - B3 * B3);
    S = Sxx + Syy + Szz;

    prj_z4c_compute_rhs(&mesh, 0, 0, 0, 0, 0.0);
    rhs = prj_block_z4c_rhs_stage(block, 0);
    assert_close_rel("rhs hydro Khat", rhs[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)],
        4.0 * pi * (S + E), 1.0e-10);
    assert_close_rel("rhs hydro Theta", rhs[Z4CIDX(PRJ_Z4C_THETA, i, j, k)],
        -8.0 * pi * E, 1.0e-10);
    assert_close_rel("rhs hydro Gamx", rhs[Z4CIDX(PRJ_Z4C_GAMX, i, j, k)],
        -16.0 * pi * Sx, 1.0e-10);
    assert_close_rel("rhs hydro Axx", rhs[Z4CIDX(PRJ_Z4C_AXX, i, j, k)],
        -8.0 * pi * (Sxx - S / 3.0), 1.0e-10);
    assert_close_rel("rhs hydro Axy", rhs[Z4CIDX(PRJ_Z4C_AXY, i, j, k)],
        -8.0 * pi * Sxy, 1.0e-10);
    prj_mesh_destroy(&mesh);
}

#if !PRJ_MHD
static void check_rhs_hydro_matter_uses_z4c_metric(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    double *W;
    double *rhs;
    const int i = 1, j = 2, k = 3;
    const double pi = acos(-1.0);
    const double chi = 0.25;
    const double gxx = 2.0;
    const double gyy = 0.5;
    const double gzz = 1.0;
    const double gamma_xx = gxx / chi;
    const double gamma_yy = gyy / chi;
    const double gamma_zz = gzz / chi;
    const double gu_xx = 1.0 / gxx;
    const double gamma_inv_xx = chi / gxx;
    const double gamma_inv_yy = chi / gyy;
    const double gamma_inv_zz = chi / gzz;
    const double rho = 2.0;
    const double v1 = 1.0e9;
    const double eint = 7.0;
    const double pressure = 11.0;
    const double beta1 = v1 / PRJ_CLIGHT;
    const double beta_cov1 = gamma_xx * beta1;
    const double beta2 = gamma_xx * beta1 * beta1;
    const double wlor2 = 1.0 / (1.0 - beta2);
    const double rhoh = rho * PRJ_CLIGHT * PRJ_CLIGHT + rho * eint + pressure;
    double fac = geo_factor();
    double E;
    double Sx;
    double Sxx;
    double Syy;
    double Szz;
    double S;

    init_one_block_mesh(&mesh, &coord);
    mesh.use_full_dynamic_gr = 1;
    block = &mesh.blocks[0];
    prj_z4c_init_mesh_flat(&mesh, 0);
    disable_z4c_damping(&mesh);
    set_constant_diagonal_z4c_metric(block, chi, gxx, gyy, gzz);
    W = prj_block_mhd_stage(block, 0);
    W[WIDX(PRJ_PRIM_RHO, i, j, k)] = rho;
    W[WIDX(PRJ_PRIM_V1, i, j, k)] = v1;
    W[WIDX(PRJ_PRIM_V2, i, j, k)] = 0.0;
    W[WIDX(PRJ_PRIM_V3, i, j, k)] = 0.0;
    W[WIDX(PRJ_PRIM_EINT, i, j, k)] = eint;
    block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)] = pressure;

    E = fac * (rhoh * wlor2 - pressure);
    Sx = fac * rhoh * wlor2 * beta_cov1;
    Sxx = fac * (rhoh * wlor2 * beta_cov1 * beta_cov1 +
        pressure * gamma_xx);
    Syy = fac * pressure * gamma_yy;
    Szz = fac * pressure * gamma_zz;
    S = gamma_inv_xx * Sxx + gamma_inv_yy * Syy + gamma_inv_zz * Szz;

    prj_z4c_compute_rhs(&mesh, 0, 0, 0, 0, 0.0);
    rhs = prj_block_z4c_rhs_stage(block, 0);
    assert_close_rel("rhs curved hydro Khat", rhs[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)],
        4.0 * pi * (S + E), 1.0e-10);
    assert_close_rel("rhs curved hydro Theta", rhs[Z4CIDX(PRJ_Z4C_THETA, i, j, k)],
        -8.0 * pi * E, 1.0e-10);
    assert_close_rel("rhs curved hydro Gamx", rhs[Z4CIDX(PRJ_Z4C_GAMX, i, j, k)],
        -16.0 * pi * gu_xx * Sx, 1.0e-10);
    assert_close_rel("rhs curved hydro Axx", rhs[Z4CIDX(PRJ_Z4C_AXX, i, j, k)],
        -8.0 * pi * (chi * Sxx - S * gxx / 3.0), 1.0e-10);
    assert_close_rel("rhs curved hydro Ayy", rhs[Z4CIDX(PRJ_Z4C_AYY, i, j, k)],
        -8.0 * pi * (chi * Syy - S * gyy / 3.0), 1.0e-10);
    prj_mesh_destroy(&mesh);
}
#endif

#if PRJ_MHD
static void check_rhs_grmhd_matter_uses_z4c_metric(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    double *W;
    double *rhs;
    const int i = 1, j = 2, k = 3;
    const double pi = acos(-1.0);
    const double chi = 0.25;
    const double gxx = 2.0;
    const double gyy = 0.5;
    const double gzz = 1.0;
    const double gamma_xx = gxx / chi;
    const double gamma_yy = gyy / chi;
    const double gamma_zz = gzz / chi;
    const double gu_xx = 1.0 / gxx;
    const double gu_yy = 1.0 / gyy;
    const double gamma_inv_xx = chi / gxx;
    const double gamma_inv_yy = chi / gyy;
    const double gamma_inv_zz = chi / gzz;
    const double rho = 2.0;
    const double v1 = 0.02 * PRJ_CLIGHT;
    const double v2 = -0.015 * PRJ_CLIGHT;
    const double v3 = 0.01 * PRJ_CLIGHT;
    const double eint = 7.0;
    const double pressure = 11.0;
    const double B1 = 13.0;
    const double B2 = 17.0;
    const double B3 = -19.0;
    const double beta1 = v1 / PRJ_CLIGHT;
    const double beta2v = v2 / PRJ_CLIGHT;
    const double beta3 = v3 / PRJ_CLIGHT;
    const double beta_cov1 = gamma_xx * beta1;
    const double beta_cov2 = gamma_yy * beta2v;
    const double beta_cov3 = gamma_zz * beta3;
    const double beta_sq = beta_cov1 * beta1 + beta_cov2 * beta2v +
        beta_cov3 * beta3;
    const double wlor2 = 1.0 / (1.0 - beta_sq);
    const double Bcov1 = gamma_xx * B1;
    const double Bcov2 = gamma_yy * B2;
    const double Bcov3 = gamma_zz * B3;
    const double Bsq = Bcov1 * B1 + Bcov2 * B2 + Bcov3 * B3;
    const double Bbeta = Bcov1 * beta1 + Bcov2 * beta2v + Bcov3 * beta3;
    const double mag_pressure = 0.5 * (Bsq / wlor2 + Bbeta * Bbeta);
    const double ptot = pressure + mag_pressure;
    const double rhoh = rho * PRJ_CLIGHT * PRJ_CLIGHT + rho * eint + pressure;
    double fac = geo_factor();
    double E;
    double Sx;
    double Sy;
    double Sxx;
    double Sxy;
    double Syy;
    double Szz;
    double S;

    init_one_block_mesh(&mesh, &coord);
    mesh.use_full_dynamic_gr = 1;
    block = &mesh.blocks[0];
    prj_z4c_init_mesh_flat(&mesh, 0);
    disable_z4c_damping(&mesh);
    set_constant_diagonal_z4c_metric(block, chi, gxx, gyy, gzz);
    W = prj_block_mhd_stage(block, 0);
    W[WIDX(PRJ_PRIM_RHO, i, j, k)] = rho;
    W[WIDX(PRJ_PRIM_V1, i, j, k)] = v1;
    W[WIDX(PRJ_PRIM_V2, i, j, k)] = v2;
    W[WIDX(PRJ_PRIM_V3, i, j, k)] = v3;
    W[WIDX(PRJ_PRIM_EINT, i, j, k)] = eint;
    W[WIDX(PRJ_PRIM_B1, i, j, k)] = B1;
    W[WIDX(PRJ_PRIM_B2, i, j, k)] = B2;
    W[WIDX(PRJ_PRIM_B3, i, j, k)] = B3;
    block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)] = pressure;

    E = fac * (rhoh * wlor2 - pressure + Bsq - mag_pressure);
    Sx = fac * ((rhoh + Bsq) * wlor2 * beta_cov1 - Bbeta * Bcov1);
    Sy = fac * ((rhoh + Bsq) * wlor2 * beta_cov2 - Bbeta * Bcov2);
    Sxx = fac * ((rhoh + Bsq) * wlor2 * beta_cov1 * beta_cov1 +
        ptot * gamma_xx - Bcov1 * Bcov1 / wlor2 -
        2.0 * Bbeta * Bcov1 * beta_cov1);
    Sxy = fac * ((rhoh + Bsq) * wlor2 * beta_cov1 * beta_cov2 -
        Bcov1 * Bcov2 / wlor2 -
        Bbeta * (Bcov1 * beta_cov2 + Bcov2 * beta_cov1));
    Syy = fac * ((rhoh + Bsq) * wlor2 * beta_cov2 * beta_cov2 +
        ptot * gamma_yy - Bcov2 * Bcov2 / wlor2 -
        2.0 * Bbeta * Bcov2 * beta_cov2);
    Szz = fac * ((rhoh + Bsq) * wlor2 * beta_cov3 * beta_cov3 +
        ptot * gamma_zz - Bcov3 * Bcov3 / wlor2 -
        2.0 * Bbeta * Bcov3 * beta_cov3);
    S = gamma_inv_xx * Sxx + gamma_inv_yy * Syy + gamma_inv_zz * Szz;

    prj_z4c_compute_rhs(&mesh, 0, 0, 0, 0, 0.0);
    rhs = prj_block_z4c_rhs_stage(block, 0);
    assert_close_rel("rhs curved GRMHD Khat", rhs[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)],
        4.0 * pi * (S + E), 1.0e-10);
    assert_close_rel("rhs curved GRMHD Theta", rhs[Z4CIDX(PRJ_Z4C_THETA, i, j, k)],
        -8.0 * pi * E, 1.0e-10);
    assert_close_rel("rhs curved GRMHD Gamx", rhs[Z4CIDX(PRJ_Z4C_GAMX, i, j, k)],
        -16.0 * pi * gu_xx * Sx, 1.0e-10);
    assert_close_rel("rhs curved GRMHD Gamy", rhs[Z4CIDX(PRJ_Z4C_GAMY, i, j, k)],
        -16.0 * pi * gu_yy * Sy, 1.0e-10);
    assert_close_rel("rhs curved GRMHD Axx", rhs[Z4CIDX(PRJ_Z4C_AXX, i, j, k)],
        -8.0 * pi * (chi * Sxx - S * gxx / 3.0), 1.0e-10);
    assert_close_rel("rhs curved GRMHD Axy", rhs[Z4CIDX(PRJ_Z4C_AXY, i, j, k)],
        -8.0 * pi * chi * Sxy, 1.0e-10);
    assert_close_rel("rhs curved GRMHD Ayy", rhs[Z4CIDX(PRJ_Z4C_AYY, i, j, k)],
        -8.0 * pi * (chi * Syy - S * gyy / 3.0), 1.0e-10);
    prj_mesh_destroy(&mesh);
}
#endif

static void check_rhs_hydro_matter_uses_minkowski_without_full_gr(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    double *W;
    double *rhs;
    const int i = 1, j = 2, k = 3;
    const double pi = acos(-1.0);
    const double chi = 0.25;
    const double gxx = 2.0;
    const double gyy = 0.5;
    const double gzz = 1.0;
    const double gu_xx = 1.0 / gxx;
    const double gu_yy = 1.0 / gyy;
    const double gu_zz = 1.0 / gzz;
    const double rho = 2.0;
    const double v1 = 1.0e9;
    const double eint = 7.0;
    const double pressure = 11.0;
    const double beta1 = v1 / PRJ_CLIGHT;
    const double beta2 = beta1 * beta1;
    const double wlor2 = 1.0 / (1.0 - beta2);
    const double rhoh = rho * PRJ_CLIGHT * PRJ_CLIGHT + rho * eint + pressure;
    double fac = geo_factor();
    double E;
    double Sx;
    double Sxx;
    double Syy;
    double Szz;
    double S;

    init_one_block_mesh(&mesh, &coord);
    mesh.use_full_dynamic_gr = 0;
    block = &mesh.blocks[0];
    prj_z4c_init_mesh_flat(&mesh, 0);
    disable_z4c_damping(&mesh);
    set_constant_diagonal_z4c_metric(block, chi, gxx, gyy, gzz);
    W = prj_block_mhd_stage(block, 0);
    W[WIDX(PRJ_PRIM_RHO, i, j, k)] = rho;
    W[WIDX(PRJ_PRIM_V1, i, j, k)] = v1;
    W[WIDX(PRJ_PRIM_V2, i, j, k)] = 0.0;
    W[WIDX(PRJ_PRIM_V3, i, j, k)] = 0.0;
    W[WIDX(PRJ_PRIM_EINT, i, j, k)] = eint;
    block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)] = pressure;

    E = fac * (rhoh * wlor2 - pressure);
    Sx = fac * rhoh * wlor2 * beta1;
    Sxx = fac * (rhoh * wlor2 * beta1 * beta1 + pressure);
    Syy = fac * pressure;
    Szz = fac * pressure;
    S = chi * (gu_xx * Sxx + gu_yy * Syy + gu_zz * Szz);

    prj_z4c_compute_rhs(&mesh, 0, 0, 0, 0, 0.0);
    rhs = prj_block_z4c_rhs_stage(block, 0);
    assert_close_rel("rhs coevolved hydro Khat", rhs[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)],
        4.0 * pi * (S + E), 1.0e-10);
    assert_close_rel("rhs coevolved hydro Theta", rhs[Z4CIDX(PRJ_Z4C_THETA, i, j, k)],
        -8.0 * pi * E, 1.0e-10);
    assert_close_rel("rhs coevolved hydro Gamx", rhs[Z4CIDX(PRJ_Z4C_GAMX, i, j, k)],
        -16.0 * pi * gu_xx * Sx, 1.0e-10);
    assert_close_rel("rhs coevolved hydro Axx", rhs[Z4CIDX(PRJ_Z4C_AXX, i, j, k)],
        -8.0 * pi * (chi * Sxx - S * gxx / 3.0), 1.0e-10);
    assert_close_rel("rhs coevolved hydro Ayy", rhs[Z4CIDX(PRJ_Z4C_AYY, i, j, k)],
        -8.0 * pi * (chi * Syy - S * gyy / 3.0), 1.0e-10);
    prj_mesh_destroy(&mesh);
}

#if PRJ_USE_RADIATION_M1
static double test_m1_chi_exact(double f)
{
    return (3.0 + 4.0 * f * f) / (5.0 + 2.0 * sqrt(4.0 - 3.0 * f * f));
}

static void check_rhs_m1_matter_projection(void)
{
    prj_mesh mesh;
    prj_coord coord;
    prj_block *block;
    prj_rad rad;
    double *W_rad;
    double *rhs;
    const int i = 0, j = 0, k = 0;
    const double pi = acos(-1.0);
    const double Erad = 2.0e20;
    const double f = 0.5;
    const double Frad = f * PRJ_CLIGHT * Erad;
    const double chi = test_m1_chi_exact(f);
    const double a_c = 0.5 * (1.0 - chi);
    const double b_c = 0.5 * (3.0 * chi - 1.0);
    double fac = geo_factor();
    double E = fac * Erad;
    double Sx = fac * Frad / PRJ_CLIGHT;
    double Sxx = fac * Erad * (a_c + b_c);
    double Syy = fac * Erad * a_c;
    double Szz = Syy;
    double S = Sxx + Syy + Szz;
    int n;

    init_one_block_mesh(&mesh, &coord);
    block = &mesh.blocks[0];
    memset(&rad, 0, sizeof(rad));
    for (n = 0; n <= NCLOSURE; ++n) {
        rad.chi[n] = test_m1_chi_exact((double)n / (double)NCLOSURE);
    }
    prj_z4c_init_mesh_flat(&mesh, 0);
    disable_z4c_damping(&mesh);
    W_rad = prj_block_rad_stage(block, 0);
    W_rad[WIDX(PRJ_RAD_PRIM_E(0, 0), i, j, k)] = Erad / RAD_SCALE;
    W_rad[WIDX(PRJ_RAD_PRIM_F1(0, 0), i, j, k)] = Frad / RAD_SCALE;
    prj_z4c_compute_rhs(&mesh, 0, &rad, 0, 0, 0.0);
    rhs = prj_block_z4c_rhs_stage(block, 0);
    assert_close_rel("rhs m1 Khat", rhs[Z4CIDX(PRJ_Z4C_KHAT, i, j, k)],
        4.0 * pi * (S + E), 1.0e-10);
    assert_close_rel("rhs m1 Theta", rhs[Z4CIDX(PRJ_Z4C_THETA, i, j, k)],
        -8.0 * pi * E, 1.0e-10);
    assert_close_rel("rhs m1 Gamx", rhs[Z4CIDX(PRJ_Z4C_GAMX, i, j, k)],
        -16.0 * pi * Sx, 1.0e-10);
    assert_close_rel("rhs m1 Axx", rhs[Z4CIDX(PRJ_Z4C_AXX, i, j, k)],
        -8.0 * pi * (Sxx - S / 3.0), 1.0e-10);
    assert_close_rel("rhs m1 Ayy", rhs[Z4CIDX(PRJ_Z4C_AYY, i, j, k)],
        -8.0 * pi * (Syy - S / 3.0), 1.0e-10);
    prj_mesh_destroy(&mesh);
}
#endif

#endif

int main(int argc, char **argv)
{
#if defined(PRJ_ENABLE_MPI)
    MPI_Init(&argc, &argv);
#else
    (void)argc;
    (void)argv;
#endif
    check_enum_names();
#if PRJ_DYNAMIC_GR
    check_flat_state();
    check_dt_includes_damping_terms();
    check_z4c_amr_transfer();
    check_z4c_amr_prolong_ghost_fill();
    check_z4c_dchi_amr_tag();
    check_z4c_sommerfeld_rhs();
    check_rhs_hydro_matter_projection();
    check_rhs_hydro_matter_uses_minkowski_without_full_gr();
#if !PRJ_MHD
    check_rhs_hydro_matter_uses_z4c_metric();
#else
    check_rhs_grmhd_matter_uses_z4c_metric();
#endif
#if PRJ_USE_RADIATION_M1
    check_rhs_m1_matter_projection();
#endif
#endif
    printf("test_z4c: ok\n");
#if defined(PRJ_ENABLE_MPI)
    MPI_Finalize();
#endif
    return 0;
}
