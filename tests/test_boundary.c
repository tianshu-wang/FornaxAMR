#include "prj.h"

static void prj_set_block_primitive(prj_block *block)
{
    int i;
    int j;
    int k;
    double eint = 1.5 * (double)(block->id + 1);

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                block->W[VIDX(PRJ_PRIM_RHO, i, j, k)] = 1.0;
                block->W[VIDX(PRJ_PRIM_V1, i, j, k)] = 0.0;
                block->W[VIDX(PRJ_PRIM_V2, i, j, k)] = 0.0;
                block->W[VIDX(PRJ_PRIM_V3, i, j, k)] = 0.0;
                block->W[VIDX(PRJ_PRIM_EINT, i, j, k)] = eint;
                block->W[VIDX(PRJ_PRIM_YE, i, j, k)] = 0.1;
                block->W1[VIDX(PRJ_PRIM_RHO, i, j, k)] = 0.0;
                block->W1[VIDX(PRJ_PRIM_V1, i, j, k)] = 0.0;
                block->W1[VIDX(PRJ_PRIM_V2, i, j, k)] = 0.0;
                block->W1[VIDX(PRJ_PRIM_V3, i, j, k)] = 0.0;
                block->W1[VIDX(PRJ_PRIM_EINT, i, j, k)] = 0.0;
                block->W1[VIDX(PRJ_PRIM_YE, i, j, k)] = 0.0;
            }
        }
    }
}

static int prj_almost_equal(double a, double b)
{
    double diff = a - b;

    if (diff < 0.0) {
        diff = -diff;
    }
    return diff < 1.0e-10;
}

static void prj_set_reflective_state(prj_block *block)
{
    int i;
    int j;
    int k;

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                block->W[VIDX(PRJ_PRIM_RHO, i, j, k)] = 1.0 + 0.1 * (double)(i + PRJ_NGHOST);
                block->W[VIDX(PRJ_PRIM_V1, i, j, k)] = 10.0 + (double)i;
                block->W[VIDX(PRJ_PRIM_V2, i, j, k)] = 20.0 + (double)j;
                block->W[VIDX(PRJ_PRIM_V3, i, j, k)] = 30.0 + (double)k;
                block->W[VIDX(PRJ_PRIM_EINT, i, j, k)] = 40.0 + (double)(i + j + k);
                block->W[VIDX(PRJ_PRIM_YE, i, j, k)] = 0.25;
            }
        }
    }
}

static int prj_point_inside_strict(const prj_block *block, double x1, double x2, double x3)
{
    return x1 > block->xmin[0] && x1 < block->xmax[0] &&
        x2 > block->xmin[1] && x2 < block->xmax[1] &&
        x3 > block->xmin[2] && x3 < block->xmax[2];
}

static int prj_check_neighbor_ghosts(prj_mesh *mesh)
{
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *src = &mesh->blocks[bidx];
        int n;

        if (src->id < 0 || src->active != 1) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            int nid = src->slot[n].id;

            if (nid >= 0 && mesh->blocks[nid].id >= 0 && mesh->blocks[nid].active == 1) {
                prj_block *dst = &mesh->blocks[nid];
                int i;
                int j;
                int k;

                for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
                    for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                        for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                            double x1;
                            double x2;
                            double x3;

                            if (i >= 0 && i < PRJ_BLOCK_SIZE &&
                                j >= 0 && j < PRJ_BLOCK_SIZE &&
                                k >= 0 && k < PRJ_BLOCK_SIZE) {
                                continue;
                            }
                            x1 = dst->xmin[0] + ((double)i + 0.5) * dst->dx[0];
                            x2 = dst->xmin[1] + ((double)j + 0.5) * dst->dx[1];
                            x3 = dst->xmin[2] + ((double)k + 0.5) * dst->dx[2];
                            if (prj_point_inside_strict(src, x1, x2, x3)) {
                                double expected = 1.5 * (double)(src->id + 1);

                                if (!prj_almost_equal(dst->W[VIDX(PRJ_PRIM_EINT, i, j, k)], expected)) {
                                    return 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}

static int prj_check_physical_outflow_edges(const prj_mesh *mesh)
{
    const prj_block *b = &mesh->blocks[0];

    if (!prj_almost_equal(b->W[VIDX(PRJ_PRIM_EINT, -1, -1, -1)], b->W[VIDX(PRJ_PRIM_EINT, 0, 0, 0)])) {
        return 1;
    }
    if (!prj_almost_equal(b->W[VIDX(PRJ_PRIM_EINT, -2, 3, -2)], b->W[VIDX(PRJ_PRIM_EINT, 1, 3, 1)])) {
        return 2;
    }
    if (!prj_almost_equal(b->W[VIDX(PRJ_PRIM_EINT, 4, -1, -2)], b->W[VIDX(PRJ_PRIM_EINT, 4, 0, 1)])) {
        return 3;
    }
    return 0;
}

static int prj_check_reflective_edges(prj_mesh *mesh)
{
    prj_block *b = &mesh->blocks[0];

    if (!prj_almost_equal(b->W[VIDX(PRJ_PRIM_V1, -1, -1, -1)], -b->W[VIDX(PRJ_PRIM_V1, 0, 0, 0)])) {
        return 1;
    }
    if (!prj_almost_equal(b->W[VIDX(PRJ_PRIM_V2, -1, -1, -1)], -b->W[VIDX(PRJ_PRIM_V2, 0, 0, 0)])) {
        return 2;
    }
    if (!prj_almost_equal(b->W[VIDX(PRJ_PRIM_V3, -1, -1, -1)], -b->W[VIDX(PRJ_PRIM_V3, 0, 0, 0)])) {
        return 3;
    }
    if (!prj_almost_equal(b->W[VIDX(PRJ_PRIM_RHO, -1, -1, -1)], b->W[VIDX(PRJ_PRIM_RHO, 0, 0, 0)])) {
        return 4;
    }
    return 0;
}

static int prj_check_get_prim_cases(prj_mesh *mesh)
{
    prj_block *block = &mesh->blocks[0];
    double w[PRJ_NVAR_PRIM];
    int i;
    int j;
    int k;

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                block->W[VIDX(PRJ_PRIM_RHO, i, j, k)] = 1.0 + 2.0 * (double)i + 3.0 * (double)j + 5.0 * (double)k;
                block->W[VIDX(PRJ_PRIM_V1, i, j, k)] = -2.0 + (double)i - (double)j + 0.5 * (double)k;
                block->W[VIDX(PRJ_PRIM_V2, i, j, k)] = 0.25 * (double)i + 0.5 * (double)j + 0.75 * (double)k;
                block->W[VIDX(PRJ_PRIM_V3, i, j, k)] = 7.0 - (double)i + 2.0 * (double)j - 3.0 * (double)k;
                block->W[VIDX(PRJ_PRIM_EINT, i, j, k)] = 11.0 + 4.0 * (double)i - 2.0 * (double)j + (double)k;
                block->W[VIDX(PRJ_PRIM_YE, i, j, k)] = 0.1 + 0.01 * (double)(i + j + k);
            }
        }
    }

    prj_boundary_get_prim(
        block, 1,
        block->xmin[0] + 3.5 * block->dx[0],
        block->xmin[1] + 4.5 * block->dx[1],
        block->xmin[2] + 5.5 * block->dx[2],
        w);
    if (!prj_almost_equal(w[PRJ_PRIM_EINT], block->W[VIDX(PRJ_PRIM_EINT, 3, 4, 5)])) {
        return 1;
    }

    prj_boundary_get_prim(
        block, 1,
        block->xmin[0] + 4.0 * block->dx[0],
        block->xmin[1] + 5.0 * block->dx[1],
        block->xmin[2] + 6.0 * block->dx[2],
        w);
    if (!prj_almost_equal(
            w[PRJ_PRIM_EINT],
            0.125 * (
                block->W[VIDX(PRJ_PRIM_EINT, 3, 4, 5)] +
                block->W[VIDX(PRJ_PRIM_EINT, 4, 4, 5)] +
                block->W[VIDX(PRJ_PRIM_EINT, 3, 5, 5)] +
                block->W[VIDX(PRJ_PRIM_EINT, 4, 5, 5)] +
                block->W[VIDX(PRJ_PRIM_EINT, 3, 4, 6)] +
                block->W[VIDX(PRJ_PRIM_EINT, 4, 4, 6)] +
                block->W[VIDX(PRJ_PRIM_EINT, 3, 5, 6)] +
                block->W[VIDX(PRJ_PRIM_EINT, 4, 5, 6)]))) {
        return 2;
    }

    prj_boundary_get_prim(
        block, 1,
        block->xmin[0] + 3.75 * block->dx[0],
        block->xmin[1] + 4.5 * block->dx[1],
        block->xmin[2] + 5.5 * block->dx[2],
        w);
    if (!prj_almost_equal(w[PRJ_PRIM_EINT], 11.0 + 4.0 * 3.25 - 2.0 * 4.0 + 5.0)) {
        return 3;
    }

    return 0;
}

int main(void)
{
    prj_coord coord = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    prj_bc bc = {
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW,
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW,
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW
    };
    prj_mesh mesh;
    int i;

    if (prj_mesh_init(&mesh, 2, 2, 2, 1, &coord) != 0) {
        return 1;
    }
    for (i = 0; i < mesh.nblocks; ++i) {
        if (mesh.blocks[i].id >= 0 && mesh.blocks[i].active == 1) {
            prj_set_block_primitive(&mesh.blocks[i]);
        }
    }
    prj_boundary_fill_ghosts(&mesh, &bc, 1);
    if (prj_check_neighbor_ghosts(&mesh) != 0) {
        prj_mesh_destroy(&mesh);
        return 2;
    }
    if (prj_check_physical_outflow_edges(&mesh) != 0) {
        prj_mesh_destroy(&mesh);
        return 3;
    }
    if (prj_check_get_prim_cases(&mesh) != 0) {
        prj_mesh_destroy(&mesh);
        return 4;
    }

    prj_amr_refine_block(&mesh, 7);
    for (i = 0; i < mesh.nblocks; ++i) {
        if (mesh.blocks[i].id >= 0 && mesh.blocks[i].active == 1) {
            prj_set_block_primitive(&mesh.blocks[i]);
        }
    }
    prj_boundary_fill_ghosts(&mesh, &bc, 1);
    if (prj_check_neighbor_ghosts(&mesh) != 0) {
        prj_mesh_destroy(&mesh);
        return 5;
    }

    {
        prj_bc reflect_bc = {
            PRJ_BC_REFLECT, PRJ_BC_OUTFLOW,
            PRJ_BC_REFLECT, PRJ_BC_OUTFLOW,
            PRJ_BC_REFLECT, PRJ_BC_OUTFLOW
        };

        prj_set_reflective_state(&mesh.blocks[0]);
        prj_boundary_fill_ghosts(&mesh, &reflect_bc, 1);
        if (prj_check_reflective_edges(&mesh) != 0) {
            prj_mesh_destroy(&mesh);
            return 6;
        }
    }

    prj_mesh_destroy(&mesh);
    return 0;
}
