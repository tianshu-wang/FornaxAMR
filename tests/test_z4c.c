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

    dt = prj_z4c_calc_dt_seconds(&mesh, 0, 0.5);
    assert_close("dt_z4c", dt, 0.5, 1.0e-3);

    prj_z4c_init_mesh_flat(&mesh, 0);
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

static double z4c_linear_pattern(int var, double i, double j, double k, int stage)
{
    return 100.0 * (double)(stage + 1) + 10.0 * (double)(var + 1) +
        0.125 * i + 0.25 * j + 0.5 * k;
}

static void fill_z4c_linear(prj_block *block)
{
    int stage, var, i, j, k;

    for (stage = 0; stage < PRJ_BLOCK_NSTAGES; ++stage) {
        double *z = prj_block_z4c_stage(block, stage);

        for (i = -PRJ_NGHOST_Z4C; i < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++i) {
            for (j = -PRJ_NGHOST_Z4C; j < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++j) {
                for (k = -PRJ_NGHOST_Z4C; k < PRJ_BLOCK_SIZE + PRJ_NGHOST_Z4C; ++k) {
                    for (var = 0; var < PRJ_NZ4C; ++var) {
                        z[Z4CIDX(var, i, j, k)] =
                            z4c_linear_pattern(var, (double)i, (double)j, (double)k, stage);
                    }
                }
            }
        }
    }
}

static void init_allocated_test_block(prj_block *block, int id)
{
    memset(block, 0, sizeof(*block));
    block->id = id;
    block->active = 1;
    block->rank = 0;
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
    const int i = 3, j = 4, k = 5;
    const int var = PRJ_Z4C_THETA;
    const int stage = 1;
    double *z_child;
    double expected;
    int gi, gj, gk;
    int pi, pj, pk;
    int o;

    init_allocated_test_block(&parent, 100);
    init_allocated_test_block(&child, 101);
    fill_z4c_linear(&parent);
    prj_fill(child.z4c_rhs, (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_NZ4C *
        (size_t)PRJ_BLOCK_NCELLS_Z4C, 42.0);
    prj_z4c_amr_prolongate_child(&parent, &child, oct);
    z_child = prj_block_z4c_stage(&child, stage);
    gi = ((oct & 1) ? PRJ_BLOCK_SIZE : 0) + i;
    gj = ((oct & 2) ? PRJ_BLOCK_SIZE : 0) + j;
    gk = ((oct & 4) ? PRJ_BLOCK_SIZE : 0) + k;
    pi = gi / 2;
    pj = gj / 2;
    pk = gk / 2;
    expected = z4c_linear_pattern(var,
        (double)pi + ((gi & 1) ? 0.25 : -0.25),
        (double)pj + ((gj & 1) ? 0.25 : -0.25),
        (double)pk + ((gk & 1) ? 0.25 : -0.25), stage);
    assert_close("z4c prolong", z_child[Z4CIDX(var, i, j, k)], expected, 1.0e-12);
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
        double *z_parent = prj_block_z4c_stage(&parent, stage);

        expected = 1000.0 + 100.0 * (double)o + 10.0 * (double)stage +
            0.01 * (double)var;
        assert_close("z4c restrict",
            z_parent[Z4CIDX(var, ioff + 1, joff + 1, koff + 1)], expected, 1.0e-12);
    }
    check_z4c_aux_cleared(&parent);
    for (o = 0; o < 8; ++o) {
        prj_block_free_data(&children_storage[o]);
    }
    prj_block_free_data(&child);
    prj_block_free_data(&parent);
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
    check_z4c_amr_transfer();
    check_z4c_sommerfeld_rhs();
    check_rhs_hydro_matter_projection();
    check_rhs_hydro_matter_uses_minkowski_without_full_gr();
#if !PRJ_MHD
    check_rhs_hydro_matter_uses_z4c_metric();
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
