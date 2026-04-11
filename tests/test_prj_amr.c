#include <stddef.h>

#include "prj.h"

static int prj_almost_equal(double a, double b)
{
    double diff = a - b;

    if (diff < 0.0) {
        diff = -diff;
    }
    return diff < 1.0e-10;
}

static int prj_point_in_block(const prj_block *b, double x, double y, double z)
{
    return b->xmin[0] <= x && x <= b->xmax[0] &&
        b->xmin[1] <= y && y <= b->xmax[1] &&
        b->xmin[2] <= z && z <= b->xmax[2];
}

static double prj_total_mass(const prj_mesh *mesh)
{
    double total = 0.0;
    int bidx;

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *b = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (b->id < 0 || b->active != 1) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    total += b->U[VIDX(PRJ_CONS_RHO, i, j, k)] * b->vol;
                }
            }
        }
    }
    return total;
}

static int prj_has_back_neighbor(const prj_mesh *mesh, int src_id, int dst_id)
{
    int n;
    const prj_block *dst = &mesh->blocks[dst_id];

    for (n = 0; n < 56; ++n) {
        if (dst->slot[n].id == src_id) {
            return 1;
        }
    }
    return 0;
}

static void prj_fill_linear_block(prj_block *b)
{
    int i;
    int j;
    int k;

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                double x = b->xmin[0] + ((double)i + 0.5) * b->dx[0];
                double y = b->xmin[1] + ((double)j + 0.5) * b->dx[1];
                double z = b->xmin[2] + ((double)k + 0.5) * b->dx[2];

                b->U[VIDX(PRJ_CONS_RHO, i, j, k)] = 1.0 + 2.0 * x - 3.0 * y + 0.5 * z;
                b->U[VIDX(PRJ_CONS_ETOT, i, j, k)] = 2.0 + x + y + z;
                b->U[VIDX(PRJ_CONS_MOM1, i, j, k)] = 0.0;
                b->U[VIDX(PRJ_CONS_MOM2, i, j, k)] = 0.0;
                b->U[VIDX(PRJ_CONS_MOM3, i, j, k)] = 0.0;
                b->U[VIDX(PRJ_CONS_YE, i, j, k)] = 0.1 * b->U[VIDX(PRJ_CONS_RHO, i, j, k)];
            }
        }
    }
}

int main(void)
{
    prj_coord coord = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    prj_eos eos;
    prj_mesh mesh;
    prj_mesh adapt_mesh;
    prj_block *parent;
    prj_block *children[8];
    prj_block roundtrip;
    double mass_before;
    double mass_after_refine;
    double mass_after_coarsen;
    int oct;
    int active_children = 0;
    int found_refined = 0;

    eos.filename[0] = '\0';
    roundtrip.id = -1;
    roundtrip.W = 0;
    roundtrip.W1 = 0;
    roundtrip.U = 0;
    roundtrip.dUdt = 0;
    roundtrip.flux[0] = 0;
    roundtrip.flux[1] = 0;
    roundtrip.flux[2] = 0;
    if (prj_mesh_init(&mesh, 2, 2, 2, 2, &coord) != 0) {
        return 1;
    }
    mesh.amr_refine_thresh = 0.05;
    mesh.amr_derefine_thresh = 0.0;

    parent = prj_mesh_get_block(&mesh, 7);
    if (parent == 0 || !prj_point_in_block(parent, 0.5, 0.5, 0.5)) {
        prj_mesh_destroy(&mesh);
        return 2;
    }

    for (oct = 0; oct < mesh.nblocks; ++oct) {
        prj_block *b = &mesh.blocks[oct];
        int i;
        int j;
        int k;

        if (b->id < 0 || b->active != 1) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double x = b->xmin[0] + ((double)i + 0.5) * b->dx[0];
                    double rho = x < 0.5 ? 1.0 : 4.0;

                    b->U[VIDX(PRJ_CONS_RHO, i, j, k)] = rho;
                    b->U[VIDX(PRJ_CONS_ETOT, i, j, k)] = rho;
                    b->U[VIDX(PRJ_CONS_MOM1, i, j, k)] = 0.0;
                    b->U[VIDX(PRJ_CONS_MOM2, i, j, k)] = 0.0;
                    b->U[VIDX(PRJ_CONS_MOM3, i, j, k)] = 0.0;
                    b->U[VIDX(PRJ_CONS_YE, i, j, k)] = 0.1 * rho;
                }
            }
        }
    }

    mass_before = prj_total_mass(&mesh);
    prj_amr_refine_block(&mesh, 7);
    parent = prj_mesh_get_block(&mesh, 7);
    if (parent->active != 0) {
        prj_mesh_destroy(&mesh);
        return 3;
    }
    for (oct = 0; oct < 8; ++oct) {
        children[oct] = prj_mesh_get_block(&mesh, parent->children[oct]);
        if (children[oct] == 0 || children[oct]->parent != 7) {
            prj_mesh_destroy(&mesh);
            return 4;
        }
        if (!prj_almost_equal(children[oct]->xmax[0] - children[oct]->xmin[0], 0.5) ||
            !prj_almost_equal(children[oct]->xmax[1] - children[oct]->xmin[1], 0.5) ||
            !prj_almost_equal(children[oct]->xmax[2] - children[oct]->xmin[2], 0.5)) {
            prj_mesh_destroy(&mesh);
            return 5;
        }
        {
            int ncount = 0;
            int n;

            for (n = 0; n < 56; ++n) {
                if (children[oct]->slot[n].id >= 0) {
                    ncount += 1;
                }
            }
            if (ncount < 7 || ncount > 14) {
                prj_mesh_destroy(&mesh);
                return 6;
            }
        }
        active_children += 1;
    }
    if (active_children != 8) {
        prj_mesh_destroy(&mesh);
        return 7;
    }
    mass_after_refine = prj_total_mass(&mesh);
    if (!prj_almost_equal(mass_before, mass_after_refine)) {
        prj_mesh_destroy(&mesh);
        return 8;
    }

    if (prj_block_alloc_data(&roundtrip) != 0) {
        prj_mesh_destroy(&mesh);
        return 9;
    }
    roundtrip.xmin[0] = parent->xmin[0];
    roundtrip.xmax[0] = parent->xmax[0];
    roundtrip.xmin[1] = parent->xmin[1];
    roundtrip.xmax[1] = parent->xmax[1];
    roundtrip.xmin[2] = parent->xmin[2];
    roundtrip.xmax[2] = parent->xmax[2];
    roundtrip.dx[0] = parent->dx[0];
    roundtrip.dx[1] = parent->dx[1];
    roundtrip.dx[2] = parent->dx[2];
    prj_fill_linear_block(parent);
    for (oct = 0; oct < 8; ++oct) {
        prj_amr_prolongate(parent, children[oct], oct);
    }
    {
        const prj_block *child_ptrs[8];
        int i;
        int j;
        int k;

        for (oct = 0; oct < 8; ++oct) {
            child_ptrs[oct] = children[oct];
        }
        prj_amr_restrict(child_ptrs, &roundtrip);
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    if (!prj_almost_equal(roundtrip.U[VIDX(PRJ_CONS_RHO, i, j, k)],
                            parent->U[VIDX(PRJ_CONS_RHO, i, j, k)])) {
                        prj_block_free_data(&roundtrip);
                        prj_mesh_destroy(&mesh);
                        return 10;
                    }
                }
            }
        }
    }
    prj_block_free_data(&roundtrip);

    for (oct = 0; oct < mesh.nblocks; ++oct) {
        prj_block *b = &mesh.blocks[oct];
        int i;
        int j;
        int k;

        if (b->id < 0 || b->active != 1) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double x = b->xmin[0] + ((double)i + 0.5) * b->dx[0];
                    double rho = x < 0.5 ? 1.0 : 4.0;

                    b->U[VIDX(PRJ_CONS_RHO, i, j, k)] = rho;
                    b->U[VIDX(PRJ_CONS_ETOT, i, j, k)] = rho;
                    b->U[VIDX(PRJ_CONS_MOM1, i, j, k)] = 0.0;
                    b->U[VIDX(PRJ_CONS_MOM2, i, j, k)] = 0.0;
                    b->U[VIDX(PRJ_CONS_MOM3, i, j, k)] = 0.0;
                    b->U[VIDX(PRJ_CONS_YE, i, j, k)] = 0.1 * rho;
                }
            }
        }
    }

    prj_amr_coarsen_block(&mesh, 7);
    mass_after_coarsen = prj_total_mass(&mesh);
    if (!prj_almost_equal(mass_after_refine, mass_after_coarsen)) {
        prj_mesh_destroy(&mesh);
        return 11;
    }

    prj_mesh_destroy(&mesh);

    if (prj_mesh_init(&adapt_mesh, 2, 2, 2, 2, &coord) != 0) {
        return 12;
    }
    adapt_mesh.amr_refine_thresh = 0.05;
    adapt_mesh.amr_derefine_thresh = 0.0;

    parent = prj_mesh_get_block(&adapt_mesh, 7);
    prj_amr_refine_block(&adapt_mesh, 7);
    parent = prj_mesh_get_block(&adapt_mesh, 7);
    children[0] = prj_mesh_get_block(&adapt_mesh, parent->children[0]);
    children[0]->refine_flag = 1;
    found_refined = 0;
    prj_amr_adapt(&adapt_mesh, &eos);

    for (oct = 0; oct < adapt_mesh.nblocks; ++oct) {
        prj_block *b = &adapt_mesh.blocks[oct];
        int n;

        if (b->id < 0 || b->active != 1) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            int nid = b->slot[n].id;

            if (nid >= 0 && adapt_mesh.blocks[nid].id >= 0 && adapt_mesh.blocks[nid].active == 1) {
                int dlevel = b->level - adapt_mesh.blocks[nid].level;

                if (dlevel < 0) {
                    dlevel = -dlevel;
                }
                if (dlevel > 1) {
                    prj_mesh_destroy(&adapt_mesh);
                    return 13;
                }
            }
        }
        if (b->level == 2) {
            found_refined = 1;
        }
    }
    if (found_refined == 0) {
        prj_mesh_destroy(&adapt_mesh);
        return 14;
    }

    for (oct = 0; oct < adapt_mesh.nblocks; ++oct) {
        prj_block *b = &adapt_mesh.blocks[oct];
        int n;

        if (b->id < 0 || b->active != 1) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            if (b->slot[n].id >= 0 && !prj_has_back_neighbor(&adapt_mesh, b->id, b->slot[n].id)) {
                prj_mesh_destroy(&adapt_mesh);
                return 15;
            }
        }
    }

    for (oct = 0; oct < adapt_mesh.nblocks; ++oct) {
        prj_block *b = &adapt_mesh.blocks[oct];
        int i;
        int j;
        int k;

        if (b->id < 0 || b->active != 1) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double x = b->xmin[0] + ((double)i + 0.5) * b->dx[0];
                    double rho = x < 0.5 ? 1.0 : 4.0;

                    b->U[VIDX(PRJ_CONS_RHO, i, j, k)] = rho;
                    b->U[VIDX(PRJ_CONS_ETOT, i, j, k)] = rho;
                    b->U[VIDX(PRJ_CONS_MOM1, i, j, k)] = 0.0;
                    b->U[VIDX(PRJ_CONS_MOM2, i, j, k)] = 0.0;
                    b->U[VIDX(PRJ_CONS_MOM3, i, j, k)] = 0.0;
                    b->U[VIDX(PRJ_CONS_YE, i, j, k)] = 0.1 * rho;
                }
            }
        }
    }
    prj_amr_adapt(&adapt_mesh, &eos);

    {
        int refined_blocks = 0;

        for (oct = 0; oct < adapt_mesh.nblocks; ++oct) {
            if (adapt_mesh.blocks[oct].id >= 0 && adapt_mesh.blocks[oct].active == 1 &&
                adapt_mesh.blocks[oct].xmin[0] < 0.5 && adapt_mesh.blocks[oct].xmax[0] > 0.0) {
                refined_blocks += adapt_mesh.blocks[oct].level > 0;
            }
        }
        if (refined_blocks == 0) {
            prj_mesh_destroy(&adapt_mesh);
            return 16;
        }
    }

    prj_mesh_destroy(&adapt_mesh);
    return 0;
}
