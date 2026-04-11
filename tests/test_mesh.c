#include <stddef.h>

#include "prj_defs.h"
#include "prj_types.h"
#include "prj_mesh.h"

static int prj_neighbor_slot_index(int ox, int oy, int oz)
{
    int index;

    if (ox == 0 && oy == 0 && oz == 0) {
        return -1;
    }

    index = (ox + 1) * 9 + (oy + 1) * 3 + (oz + 1);
    if (index > 13) {
        index -= 1;
    }
    return index;
}

static int prj_almost_equal(double a, double b)
{
    double diff = a - b;

    if (diff < 0.0) {
        diff = -diff;
    }
    return diff < 1.0e-12;
}

int main(void)
{
    prj_coord coord = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    prj_mesh mesh;
    int i;
    int j;
    int k;
    size_t prim_count = (size_t)PRJ_NVAR_PRIM * (size_t)PRJ_BLOCK_NCELLS;
    size_t cons_count = (size_t)PRJ_NVAR_CONS * (size_t)PRJ_BLOCK_NCELLS;

    if (prj_mesh_init(&mesh, 2, 2, 2, 0, &coord) != 0) {
        return 1;
    }

    if (mesh.nblocks != 8) {
        prj_mesh_destroy(&mesh);
        return 2;
    }
    if (prj_mesh_count_active(&mesh) != 8) {
        prj_mesh_destroy(&mesh);
        return 3;
    }

    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            for (k = 0; k < 2; ++k) {
                int id = (i * 2 + j) * 2 + k;
                prj_block *b = prj_mesh_get_block(&mesh, id);
                int ox;
                int oy;
                int oz;

                if (b == 0) {
                    prj_mesh_destroy(&mesh);
                    return 4;
                }

                if (!prj_almost_equal(b->xmin[0], -1.0 + (double)i) ||
                    !prj_almost_equal(b->xmax[0], (double)i) ||
                    !prj_almost_equal(b->xmin[1], -1.0 + (double)j) ||
                    !prj_almost_equal(b->xmax[1], (double)j) ||
                    !prj_almost_equal(b->xmin[2], -1.0 + (double)k) ||
                    !prj_almost_equal(b->xmax[2], (double)k)) {
                    prj_mesh_destroy(&mesh);
                    return 5;
                }

                for (ox = -1; ox <= 1; ++ox) {
                    for (oy = -1; oy <= 1; ++oy) {
                        for (oz = -1; oz <= 1; ++oz) {
                            int slot_index = prj_neighbor_slot_index(ox, oy, oz);
                            int ni = i + ox;
                            int nj = j + oy;
                            int nk = k + oz;

                            if (slot_index < 0) {
                                continue;
                            }

                            if (ni < 0 || ni >= 2 || nj < 0 || nj >= 2 || nk < 0 || nk >= 2) {
                                if (b->slot[slot_index].id != -1) {
                                    prj_mesh_destroy(&mesh);
                                    return 6;
                                }
                                continue;
                            }

                            if (b->slot[slot_index].id != (ni * 2 + nj) * 2 + nk) {
                                prj_mesh_destroy(&mesh);
                                return 7;
                            }

                            {
                                prj_block *neighbor = prj_mesh_get_block(&mesh, b->slot[slot_index].id);
                                int opposite = prj_neighbor_slot_index(-ox, -oy, -oz);

                                if (neighbor == 0 || neighbor->slot[opposite].id != b->id) {
                                    prj_mesh_destroy(&mesh);
                                    return 8;
                                }
                            }
                        }
                    }
                }

                if ((size_t)(b->W1 - b->W) != prim_count) {
                    prj_mesh_destroy(&mesh);
                    return 9;
                }
                if ((size_t)(b->U - b->W1) != prim_count) {
                    prj_mesh_destroy(&mesh);
                    return 10;
                }
                if ((size_t)(b->dUdt - b->U) != cons_count) {
                    prj_mesh_destroy(&mesh);
                    return 11;
                }
                if ((size_t)(b->flux[0] - b->dUdt) != cons_count) {
                    prj_mesh_destroy(&mesh);
                    return 12;
                }
                if ((size_t)(b->flux[1] - b->flux[0]) != cons_count) {
                    prj_mesh_destroy(&mesh);
                    return 13;
                }
                if ((size_t)(b->flux[2] - b->flux[1]) != cons_count) {
                    prj_mesh_destroy(&mesh);
                    return 14;
                }
            }
        }
    }

    prj_mesh_destroy(&mesh);
    return 0;
}
