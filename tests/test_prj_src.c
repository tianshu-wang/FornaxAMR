#include "prj.h"

int main(void)
{
    prj_coord coord = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    prj_eos eos;
    prj_mesh mesh;
    prj_block *block;
    int n;

    eos.filename[0] = '\0';

    if (prj_mesh_init(&mesh, 1, 1, 1, 0, &coord) != 0) {
        return 1;
    }
    block = &mesh.blocks[0];
    for (n = 0; n < PRJ_NVAR_PRIM * PRJ_BLOCK_NCELLS; ++n) {
        block->W[n] = (double)(n % 7);
    }
    for (n = 0; n < PRJ_NVAR_CONS * PRJ_BLOCK_NCELLS; ++n) {
        block->U[n] = (double)(n % 5);
    }
    for (n = 0; n < PRJ_NVAR_CONS * PRJ_BLOCK_NCELLS; ++n) {
        block->dUdt[n] = 3.5;
    }

    prj_src_update(&mesh, &eos, block->W, block->U, block->dUdt);

    for (n = 0; n < PRJ_NVAR_CONS * PRJ_BLOCK_NCELLS; ++n) {
        if (block->dUdt[n] != 3.5) {
            prj_mesh_destroy(&mesh);
            return 2;
        }
    }

    prj_mesh_destroy(&mesh);
    return 0;
}
