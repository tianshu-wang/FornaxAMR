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
                            if (x1 > src->xmin[0] && x1 < src->xmax[0] &&
                                x2 > src->xmin[1] && x2 < src->xmax[1] &&
                                x3 > src->xmin[2] && x3 < src->xmax[2]) {
                                double expected = 1.5 * (double)(src->id + 1);

                                if (dst->W[VIDX(PRJ_PRIM_EINT, i, j, k)] != expected) {
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
    if (mesh.blocks[0].W[VIDX(PRJ_PRIM_EINT, -1, -1, -1)] != 1.5) {
        prj_mesh_destroy(&mesh);
        return 3;
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
        return 4;
    }

    prj_mesh_destroy(&mesh);
    return 0;
}
