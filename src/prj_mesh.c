#include <stddef.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

static int prj_mesh_block_is_local_active(const prj_mpi *mpi, const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static int prj_neighbor_slot_index(int ox, int oy, int oz)
{
    int index;

    if (ox < -1 || ox > 1 || oy < -1 || oy > 1 || oz < -1 || oz > 1) {
        return -1;
    }
    if (ox == 0 && oy == 0 && oz == 0) {
        return -1;
    }

    index = (ox + 1) * 9 + (oy + 1) * 3 + (oz + 1);
    if (index > 13) {
        index -= 1;
    }
    return index;
}

static double prj_neighbor_abs(double x)
{
    return x < 0.0 ? -x : x;
}

void prj_neighbor_compute_geometry(const prj_block *a, const prj_block *b, prj_neighbor *slot)
{
    double tol=0;
    int axisrel[3];
    int touching;
    int d;

    if (a == 0 || b == 0 || slot == 0) {
        return;
    }

    touching = 0;
    for (d = 0; d < 3; ++d) {
        tol = 1.0e-2*PRJ_MIN(a->dx[d],b->dx[d]);
        if (prj_neighbor_abs(a->xmax[d] - b->xmin[d]) < tol) {
            axisrel[d] = 1;
            touching += 1;
        } else if (prj_neighbor_abs(b->xmax[d] - a->xmin[d]) < tol) {
            axisrel[d] = -1;
            touching += 1;
        } else {
            axisrel[d] = 0;
        }
    }

    if (touching == 1) {
        slot->type = PRJ_NEIGHBOR_FACE;
    } else if (touching == 2) {
        slot->type = PRJ_NEIGHBOR_EDGE;
    } else if (touching == 3) {
        slot->type = PRJ_NEIGHBOR_CORNER;
    } else {
        slot->type = PRJ_NEIGHBOR_NONE;
    }

    slot->rel_level = b->level - a->level;

    for (d = 0; d < 3; ++d) {
        if (slot->rel_level==0){
            // Same level
            if (axisrel[d]==1) {
                slot->recv_loc_start[d] = -PRJ_NGHOST;
                slot->recv_loc_end[d] = 0;
                slot->send_loc_start[d] = PRJ_BLOCK_SIZE-PRJ_NGHOST;
                slot->send_loc_end[d] = PRJ_BLOCK_SIZE;
            } else if (axisrel[d]==-1) {
                slot->recv_loc_start[d] = PRJ_BLOCK_SIZE;
                slot->recv_loc_end[d] = PRJ_BLOCK_SIZE + PRJ_NGHOST;
                slot->send_loc_start[d] = 0;
                slot->send_loc_end[d] = PRJ_NGHOST;
            } else {
                slot->recv_loc_start[d] = 0;
                slot->recv_loc_end[d] = PRJ_BLOCK_SIZE;
                slot->send_loc_start[d] = 0;
                slot->send_loc_end[d] = PRJ_BLOCK_SIZE;
            }
        } else if (slot->rel_level==-1) {
            // Neighbor is coarser
            if (axisrel[d]==1) {
                slot->recv_loc_start[d] = -PRJ_NGHOST;
                slot->recv_loc_end[d] = 0;
                slot->send_loc_start[d] = PRJ_BLOCK_SIZE-2*PRJ_NGHOST;
                slot->send_loc_end[d] = PRJ_BLOCK_SIZE;
            } else if (axisrel[d]==-1) {
                slot->recv_loc_start[d] = PRJ_BLOCK_SIZE;
                slot->recv_loc_end[d] = PRJ_BLOCK_SIZE + PRJ_NGHOST;
                slot->send_loc_start[d] = 0;
                slot->send_loc_end[d] = 2*PRJ_NGHOST;
            } else {
                if ((b->xmin[d]+b->xmax[d]) > (a->xmin[d]+a->xmax[d])) {
                    slot->recv_loc_start[d] = 0;
                    slot->recv_loc_end[d] = PRJ_BLOCK_SIZE/2;
                } else {
                    slot->recv_loc_start[d] = PRJ_BLOCK_SIZE/2;
                    slot->recv_loc_end[d] = PRJ_BLOCK_SIZE;
                }
                slot->send_loc_start[d] = 0;
                slot->send_loc_end[d] = PRJ_BLOCK_SIZE;
            }
        } else if (slot->rel_level==1) {
            // Neighbor is finer
            if (axisrel[d]==1) {
                slot->recv_loc_start[d] = -PRJ_NGHOST;
                slot->recv_loc_end[d] = 0;
                slot->send_loc_start[d] = PRJ_BLOCK_SIZE-PRJ_NGHOST/2;
                slot->send_loc_end[d] = PRJ_BLOCK_SIZE;
            } else if (axisrel[d]==-1) {
                slot->recv_loc_start[d] = PRJ_BLOCK_SIZE;
                slot->recv_loc_end[d] = PRJ_BLOCK_SIZE + PRJ_NGHOST;
                slot->send_loc_start[d] = 0;
                slot->send_loc_end[d] = PRJ_NGHOST/2;
            } else {
                if ((b->xmin[d]+b->xmax[d]) > (a->xmin[d]+a->xmax[d])) {
                    slot->recv_loc_start[d] = -PRJ_NGHOST;
                    slot->recv_loc_end[d] = PRJ_BLOCK_SIZE;
                    slot->send_loc_start[d] = PRJ_BLOCK_SIZE/2-PRJ_NGHOST/2;
                    slot->send_loc_end[d] = PRJ_BLOCK_SIZE;
                } else {
                    slot->recv_loc_start[d] = 0;
                    slot->recv_loc_end[d] = PRJ_BLOCK_SIZE+PRJ_NGHOST;
                    slot->send_loc_start[d] = 0;
                    slot->send_loc_end[d] = PRJ_BLOCK_SIZE/2+PRJ_NGHOST/2;
                }
            }
        } else {
            fprintf(stderr,"Neighbor geometry error!\n");
            exit(1);
        }
    }

    /* Radiation uses a narrower ghost band (PRJ_NGHOST_RAD <= PRJ_NGHOST).
     * Clip the recv box to the rad zone and map the result back to a send box
     * using the same rel_level mapping the send_loc_* assignments above use. */
    for (d = 0; d < 3; ++d) {
        int rs = slot->recv_loc_start[d];
        int re = slot->recv_loc_end[d];
        int rs_rad = rs > -PRJ_NGHOST_RAD ? rs : -PRJ_NGHOST_RAD;
        int re_rad = re < PRJ_BLOCK_SIZE + PRJ_NGHOST_RAD ? re : PRJ_BLOCK_SIZE + PRJ_NGHOST_RAD;
        int base;

        if (rs_rad > re_rad) {
            rs_rad = re_rad;
        }
        slot->recv_loc_start_rad[d] = rs_rad;
        slot->recv_loc_end_rad[d] = re_rad;

        base = slot->send_loc_start[d];
        if (slot->rel_level == 0) {
            int shift = base - rs;
            slot->send_loc_start_rad[d] = rs_rad + shift;
            slot->send_loc_end_rad[d] = re_rad + shift;
        } else if (slot->rel_level == -1) {
            slot->send_loc_start_rad[d] = base + 2 * (rs_rad - rs);
            slot->send_loc_end_rad[d] = base + 2 * (re_rad - rs);
        } else {
            slot->send_loc_start_rad[d] = base + (rs_rad - rs) / 2;
            slot->send_loc_end_rad[d] = base + (re_rad - rs) / 2;
        }
    }
}

static void prj_neighbor_clear_derived(prj_neighbor *slot)
{
    int d;

    if (slot == 0) {
        return;
    }
    slot->rel_level = 0;
    slot->type = PRJ_NEIGHBOR_NONE;
    for (d = 0; d < 3; ++d) {
        slot->send_loc_start[d] = 0;
        slot->send_loc_end[d] = 0;
        slot->recv_loc_start[d] = 0;
        slot->recv_loc_end[d] = 0;
        slot->send_loc_start_rad[d] = 0;
        slot->send_loc_end_rad[d] = 0;
        slot->recv_loc_start_rad[d] = 0;
        slot->recv_loc_end_rad[d] = 0;
    }
}

static void prj_block_init_empty(prj_block *b)
{
    int n;

    b->id = -1;
    b->rank = 0;
    b->level = 0;
    b->active = 1;
    b->refine_flag = 0;
    b->can_refine = 1;
    b->W = 0;
    b->W1 = 0;
    b->eosvar = 0;
    b->cell_derived_done = 0;
    b->U = 0;
    b->dUdt = 0;
    b->flux[0] = 0;
    b->flux[1] = 0;
    b->flux[2] = 0;
    b->v_riemann[0] = 0;
    b->v_riemann[1] = 0;
    b->v_riemann[2] = 0;
    b->kappa_cell = 0;
    b->sigma_cell = 0;
    b->lapse = 0;
    for (n = 0; n < 3; ++n) {
        b->grav[n] = 0;
    }
    b->r_com = 0;
    b->Ylm = 0;
#if PRJ_MHD
    for (n = 0; n < 3; ++n) {
        b->face_fidelity[n] = 0;
        b->edge_fidelity[n] = 0;
        b->Bf[n] = 0;
        b->Bf1[n] = 0;
        b->Bv1[n] = 0;
        b->Bv2[n] = 0;
        b->emf[n] = 0;
    }
#endif
    b->ridx = 0;
    b->fr = 0;
    b->vol = 0.0;
    b->area[0] = 0.0;
    b->area[1] = 0.0;
    b->area[2] = 0.0;
    b->parent = -1;

    for (n = 0; n < 8; ++n) {
        b->children[n] = -1;
    }
    for (n = 0; n < 56; ++n) {
        b->slot[n].id = -1;
        b->slot[n].rank = 0;
        b->slot[n].xmin[0] = 0.0;
        b->slot[n].xmin[1] = 0.0;
        b->slot[n].xmin[2] = 0.0;
        b->slot[n].xmax[0] = 0.0;
        b->slot[n].xmax[1] = 0.0;
        b->slot[n].xmax[2] = 0.0;
        b->slot[n].dx[0] = 0.0;
        b->slot[n].dx[1] = 0.0;
        b->slot[n].dx[2] = 0.0;
        prj_neighbor_clear_derived(&b->slot[n]);
    }
}

int prj_block_alloc_data(prj_block *b)
{
    size_t prim_count;
    size_t eosvar_count;
    size_t cons_count;
    size_t total_count;
    double *base;
    int *cell_derived_done;
#if PRJ_MHD
    int *face_fidelity[3] = {0, 0, 0};
    int *edge_fidelity[3] = {0, 0, 0};
#endif
    int *ridx;
    double *fr;

    if (b == 0) {
        return 1;
    }

    prj_block_free_data(b);

    prim_count = (size_t)PRJ_NVAR_PRIM * (size_t)PRJ_BLOCK_NCELLS;
    eosvar_count = (size_t)PRJ_NVAR_EOSVAR * (size_t)PRJ_BLOCK_NCELLS;
    cons_count = (size_t)PRJ_NVAR_CONS * (size_t)PRJ_BLOCK_NCELLS;
    total_count = 2U * prim_count + eosvar_count + 5U * cons_count + 9U * (size_t)PRJ_BLOCK_NCELLS;
    total_count += 5U * (size_t)PRJ_BLOCK_NCELLS;
    total_count += (size_t)(LMAX*LMAX) * (size_t)PRJ_BLOCK_NCELLS;
#if PRJ_MHD
    total_count += 6U * (size_t)PRJ_BLOCK_NFACES + 6U * (size_t)PRJ_BLOCK_NCELLS + 3U * (size_t)PRJ_BLOCK_NEDGES;
#endif
#if PRJ_NRAD > 0
    total_count += 2U * (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP * (size_t)PRJ_BLOCK_NCELLS;
#endif

    base = (double *)malloc(total_count * sizeof(*base));
    if (base == 0) {
        return 2;
    }
    cell_derived_done = (int *)calloc((size_t)PRJ_BLOCK_NCELLS, sizeof(*cell_derived_done));
#if PRJ_MHD
    for (int d = 0; d < 3; ++d) {
        face_fidelity[d] = (int *)calloc((size_t)PRJ_BLOCK_NFACES, sizeof(*face_fidelity[d]));
        edge_fidelity[d] = (int *)calloc((size_t)PRJ_BLOCK_NEDGES, sizeof(*edge_fidelity[d]));
    }
#endif
    ridx = (int *)malloc((size_t)PRJ_BLOCK_NCELLS * sizeof(*ridx));
    fr = (double *)malloc((size_t)PRJ_BLOCK_NCELLS * sizeof(*fr));
    if (cell_derived_done == 0 ||
#if PRJ_MHD
        face_fidelity[0] == 0 || face_fidelity[1] == 0 || face_fidelity[2] == 0 ||
        edge_fidelity[0] == 0 || edge_fidelity[1] == 0 || edge_fidelity[2] == 0 ||
#endif
        ridx == 0 || fr == 0) {
        free(fr);
        free(ridx);
#if PRJ_MHD
        for (int d = 0; d < 3; ++d) {
            free(edge_fidelity[d]);
            free(face_fidelity[d]);
        }
#endif
        free(cell_derived_done);
        free(base);
        return 2;
    }

    b->W = base;
    base += prim_count;
    b->W1 = base;
    base += prim_count;
    b->eosvar = base;
    base += eosvar_count;
    b->cell_derived_done = cell_derived_done;
    b->U = base;
    base += cons_count;
    b->dUdt = base;
    base += cons_count;
    b->flux[0] = base;
    base += cons_count;
    b->flux[1] = base;
    base += cons_count;
    b->flux[2] = base;
    base += cons_count;
    b->v_riemann[0] = base;
    base += 3U * (size_t)PRJ_BLOCK_NCELLS;
    b->v_riemann[1] = base;
    base += 3U * (size_t)PRJ_BLOCK_NCELLS;
    b->v_riemann[2] = base;
    base += 3U * (size_t)PRJ_BLOCK_NCELLS;
    b->lapse = base;
    base += (size_t)PRJ_BLOCK_NCELLS;
    for (int d = 0; d < 3; ++d) {
        b->grav[d] = base;
        base += (size_t)PRJ_BLOCK_NCELLS;
    }
    b->r_com = base;
    base += (size_t)PRJ_BLOCK_NCELLS;
    b->Ylm = base;
    base += (size_t)(LMAX * LMAX) * (size_t)PRJ_BLOCK_NCELLS;
#if PRJ_MHD
    for (int d = 0; d < 3; ++d) {
        b->face_fidelity[d] = face_fidelity[d];
        b->edge_fidelity[d] = edge_fidelity[d];
        b->Bf[d] = base;
        base += (size_t)PRJ_BLOCK_NFACES;
    }
    for (int d = 0; d < 3; ++d) {
        b->Bf1[d] = base;
        base += (size_t)PRJ_BLOCK_NFACES;
    }
    for (int d = 0; d < 3; ++d) {
        b->Bv1[d] = base;
        base += (size_t)PRJ_BLOCK_NCELLS;
    }
    for (int d = 0; d < 3; ++d) {
        b->Bv2[d] = base;
        base += (size_t)PRJ_BLOCK_NCELLS;
    }
    for (int d = 0; d < 3; ++d) {
        b->emf[d] = base;
        base += (size_t)PRJ_BLOCK_NEDGES;
    }
#endif
#if PRJ_NRAD > 0
    b->kappa_cell = base;
    base += (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP * (size_t)PRJ_BLOCK_NCELLS;
    b->sigma_cell = base;
    base += (size_t)PRJ_NRAD * (size_t)PRJ_NEGROUP * (size_t)PRJ_BLOCK_NCELLS;
#endif
    b->ridx = ridx;
    b->fr = fr;
    prj_fill(b->lapse, (size_t)PRJ_BLOCK_NCELLS, 1.0);
    for (int d = 0; d < 3; ++d) {
        prj_fill(b->grav[d], (size_t)PRJ_BLOCK_NCELLS, 0.0);
    }
    prj_fill(b->r_com, (size_t)PRJ_BLOCK_NCELLS, 0.0);
    prj_fill(b->Ylm, (size_t)(LMAX * LMAX) * (size_t)PRJ_BLOCK_NCELLS, 0.0);
    return 0;
}

void prj_block_free_data(prj_block *b)
{
    if (b == 0) {
        return;
    }

    free(b->W);
    free(b->cell_derived_done);
#if PRJ_MHD
    for (int d = 0; d < 3; ++d) {
        free(b->face_fidelity[d]);
        free(b->edge_fidelity[d]);
    }
#endif
    free(b->ridx);
    free(b->fr);
    b->W = 0;
    b->W1 = 0;
    b->eosvar = 0;
    b->cell_derived_done = 0;
    b->U = 0;
    b->dUdt = 0;
    b->flux[0] = 0;
    b->flux[1] = 0;
    b->flux[2] = 0;
    b->v_riemann[0] = 0;
    b->v_riemann[1] = 0;
    b->v_riemann[2] = 0;
    b->kappa_cell = 0;
    b->sigma_cell = 0;
    b->lapse = 0;
    for (int d = 0; d < 3; ++d) {
        b->grav[d] = 0;
    }
    b->r_com = 0;
    b->Ylm = 0;
#if PRJ_MHD
    for (int d = 0; d < 3; ++d) {
        b->face_fidelity[d] = 0;
        b->edge_fidelity[d] = 0;
        b->Bf[d] = 0;
        b->Bf1[d] = 0;
        b->Bv1[d] = 0;
        b->Bv2[d] = 0;
        b->emf[d] = 0;
    }
#endif
    b->ridx = 0;
    b->fr = 0;
}

void prj_block_setup_geometry(prj_block *b, const prj_coord *coord)
{
    (void)coord;

    if (b == 0) {
        return;
    }

    b->vol = b->dx[0] * b->dx[1] * b->dx[2];
    b->area[0] = b->dx[1] * b->dx[2];
    b->area[1] = b->dx[0] * b->dx[2];
    b->area[2] = b->dx[0] * b->dx[1];
}

static double prj_block_max_cell_dx_over_rcom(const prj_block *b, const prj_mesh *mesh)
{
    double x_com[3] = {0.0, 0.0, 0.0};
    double dx_max;
    double max_ratio = 0.0;
    int n;

    if (b == 0) {
        return 0.0;
    }
    if (mesh != 0) {
        x_com[0] = mesh->x_com[0];
        x_com[1] = mesh->x_com[1];
        x_com[2] = mesh->x_com[2];
    }

    dx_max = b->dx[0];
    if (b->dx[1] > dx_max) {
        dx_max = b->dx[1];
    }
    if (b->dx[2] > dx_max) {
        dx_max = b->dx[2];
    }

    for (n = 0; n < 8; ++n) {
        double center[3];
        double dx1;
        double dx2;
        double dx3;
        double r;
        int d;

        for (d = 0; d < 3; ++d) {
            if ((n & (1 << d)) != 0) {
                center[d] = b->xmax[d] - 0.5 * b->dx[d];
            } else {
                center[d] = b->xmin[d] + 0.5 * b->dx[d];
            }
        }
        dx1 = center[0] - x_com[0];
        dx2 = center[1] - x_com[1];
        dx3 = center[2] - x_com[2];
        r = sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
        if (r <= 0.0) {
            return HUGE_VAL;
        }
        {
            double ratio = dx_max / r;

            if (ratio > max_ratio) {
                max_ratio = ratio;
            }
        }
    }
    return max_ratio;
}

void prj_block_update_can_refine(prj_block *b, const prj_mesh *mesh)
{
    double ratio;
    double xc[3];
    double x_com[3] = {0.0, 0.0, 0.0};
    double dx1;
    double dx2;
    double dx3;
    double r_com;
    double limit;

    if (b == 0) {
        return;
    }
    b->can_refine = 1;
    if (mesh == 0 || mesh->use_amr_angular_resolution_limit == 0) {
        return;
    }

    if (mesh != 0) {
        x_com[0] = mesh->x_com[0];
        x_com[1] = mesh->x_com[1];
        x_com[2] = mesh->x_com[2];
    }

    xc[0] = 0.5 * (b->xmin[0] + b->xmax[0]);
    xc[1] = 0.5 * (b->xmin[1] + b->xmax[1]);
    xc[2] = 0.5 * (b->xmin[2] + b->xmax[2]);
    dx1 = xc[0] - x_com[0];
    dx2 = xc[1] - x_com[1];
    dx3 = xc[2] - x_com[2];
    r_com = sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
    limit = prj_amr_angular_resolution_limit(xc[0], xc[1], xc[2], r_com);
    if (limit <= 0.0) {
        return;
    }

    ratio = prj_block_max_cell_dx_over_rcom(b, mesh);
    if (ratio < limit) {
        b->can_refine = 0;
    }
}

int prj_block_cache_index(int i, int j, int k)
{
    return IDX(i, j, k);
}

void prj_mesh_update_block_r_com(prj_block *block, const prj_mesh *mesh)
{
    /* 3-point Gauss-Legendre quadrature on [-1, 1]: nodes +-sqrt(3/5), 0 and
       weights 5/9, 8/9, 5/9. The weights are pre-divided by 2 so that, mapped
       onto a cell, the per-axis weights sum to 1 (volume average). */
    static const double gq_node[3] = {
        -0.77459666924148337704, 0.0, 0.77459666924148337704
    };
    static const double gq_wnorm[3] = {5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0};
    double x_com[3] = {0.0, 0.0, 0.0};
    int i;
    int j;
    int k;

    if (block == 0 || block->r_com == 0) {
        return;
    }
    if (block->id < 0 || block->dx[0] <= 0.0 || block->dx[1] <= 0.0 || block->dx[2] <= 0.0) {
        prj_fill(block->r_com, (size_t)PRJ_BLOCK_NCELLS, 0.0);
        return;
    }
    if (mesh != 0) {
        x_com[0] = mesh->x_com[0];
        x_com[1] = mesh->x_com[1];
        x_com[2] = mesh->x_com[2];
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                double xc = block->xmin[0] + ((double)i + 0.5) * block->dx[0] - x_com[0];
                double yc = block->xmin[1] + ((double)j + 0.5) * block->dx[1] - x_com[1];
                double zc = block->xmin[2] + ((double)k + 0.5) * block->dx[2] - x_com[2];
                double hx = 0.5 * block->dx[0];
                double hy = 0.5 * block->dx[1];
                double hz = 0.5 * block->dx[2];
                int cache_idx = prj_block_cache_index(i, j, k);
                double r_avg = 0.0;
                int a;
                int b;
                int c;

                /* Volume-averaged radius (1/V) * integral of |x - x_com| dV,
                   evaluated by 3x3x3 tensor-product Gauss quadrature. */
                for (a = 0; a < 3; ++a) {
                    double dx1 = xc + hx * gq_node[a];

                    for (b = 0; b < 3; ++b) {
                        double dx2 = yc + hy * gq_node[b];
                        double wab = gq_wnorm[a] * gq_wnorm[b];

                        for (c = 0; c < 3; ++c) {
                            double dx3 = zc + hz * gq_node[c];

                            r_avg += wab * gq_wnorm[c] *
                                sqrt(dx1 * dx1 + dx2 * dx2 + dx3 * dx3);
                        }
                    }
                }
                block->r_com[cache_idx] = r_avg;
            }
        }
    }
}

void prj_mesh_update_r_com(prj_mesh *mesh)
{
    int bidx;

    if (mesh == 0) {
        return;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_mesh_update_block_r_com(&mesh->blocks[bidx], mesh);
    }
}

prj_block *prj_mesh_get_block(prj_mesh *mesh, int id)
{
    if (mesh == 0 || id < 0 || id >= mesh->nblocks) {
        return 0;
    }

    return &mesh->blocks[id];
}

int prj_mesh_count_active(const prj_mesh *mesh)
{
    int i;
    int count;

    if (mesh == 0) {
        return 0;
    }

    count = 0;
    for (i = 0; i < mesh->nblocks; ++i) {
        if (mesh->blocks[i].active == 1) {
            count += 1;
        }
    }
    return count;
}

static double prj_mesh_block_cell_size(const prj_block *block)
{
    double cell_size;

    if (block == 0) {
        return 0.0;
    }
    cell_size = block->dx[0];
    if (block->dx[1] > cell_size) {
        cell_size = block->dx[1];
    }
    if (block->dx[2] > cell_size) {
        cell_size = block->dx[2];
    }
    return cell_size;
}

double prj_mesh_min_cell_size(const prj_mesh *mesh)
{
    double min_cell_size = 1.0e99;
    int i;
    int found = 0;

    if (mesh == 0) {
        return 0.0;
    }
    for (i = 0; i < mesh->nblocks; ++i) {
        if (mesh->blocks[i].active == 1) {
            double cell_size = prj_mesh_block_cell_size(&mesh->blocks[i]);

            if (cell_size < min_cell_size) {
                min_cell_size = cell_size;
            }
            found = 1;
        }
    }
    return found != 0 ? min_cell_size : 0.0;
}

int prj_mesh_update_center_of_mass(prj_mesh *mesh, const prj_mpi *mpi, double x_com_err_tol)
{
    double local[4] = {0.0, 0.0, 0.0, 0.0};
    double global[4] = {0.0, 0.0, 0.0, 0.0};
    double x_com_new[3];
    double dx0;
    double dx1;
    double dx2;
    double distance;
    double min_cell;
    double threshold;
    int bidx;
    int d;

    if (mesh == 0) {
        return 0;
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (!prj_mesh_block_is_local_active(mpi, block) || block->W == 0) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    double rho = block->W[VIDX(PRJ_PRIM_RHO, i, j, k)];
                    double dm = rho * block->vol;
                    double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                    double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                    double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];

                    local[0] += dm;
                    local[1] += dm * x1;
                    local[2] += dm * x2;
                    local[3] += dm * x3;
                }
            }
        }
    }

#if defined(PRJ_ENABLE_MPI)
    if (mpi != 0 && mpi->totrank > 1) {
        MPI_Allreduce(local, global, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    } else {
        for (d = 0; d < 4; ++d) {
            global[d] = local[d];
        }
    }
#else
    for (d = 0; d < 4; ++d) {
        global[d] = local[d];
    }
#endif

    if (global[0] > 0.0) {
        x_com_new[0] = global[1] / global[0];
        x_com_new[1] = global[2] / global[0];
        x_com_new[2] = global[3] / global[0];
    } else {
        for (d = 0; d < 3; ++d) {
            x_com_new[d] = mesh->x_com[d];
        }
    }

    dx0 = x_com_new[0] - mesh->x_com[0];
    dx1 = x_com_new[1] - mesh->x_com[1];
    dx2 = x_com_new[2] - mesh->x_com[2];
    distance = sqrt(dx0 * dx0 + dx1 * dx1 + dx2 * dx2);
    min_cell = mesh->min_allowable_cell_size;
    threshold = x_com_err_tol * min_cell;
    if (threshold < 0.0) {
        threshold = 0.0;
    }

    if (min_cell <= 0.0 || distance <= threshold) {
        return 0;
    }

    for (d = 0; d < 3; ++d) {
        mesh->x_com[d] = x_com_new[d];
    }
    prj_mesh_update_r_com(mesh);
    return 1;
}

void prj_mesh_update_min_allowable_cell_size(prj_mesh *mesh)
{
    double dx[3];
    double root_cell_size;
    int level = 0;

    if (mesh == 0 ||
        mesh->root_nx[0] <= 0 || mesh->root_nx[1] <= 0 || mesh->root_nx[2] <= 0) {
        return;
    }

    dx[0] = fabs(mesh->coord.x1max - mesh->coord.x1min) /
        ((double)mesh->root_nx[0] * (double)PRJ_BLOCK_SIZE);
    dx[1] = fabs(mesh->coord.x2max - mesh->coord.x2min) /
        ((double)mesh->root_nx[1] * (double)PRJ_BLOCK_SIZE);
    dx[2] = fabs(mesh->coord.x3max - mesh->coord.x3min) /
        ((double)mesh->root_nx[2] * (double)PRJ_BLOCK_SIZE);
    root_cell_size = dx[0];
    if (dx[1] > root_cell_size) {
        root_cell_size = dx[1];
    }
    if (dx[2] > root_cell_size) {
        root_cell_size = dx[2];
    }
    if (root_cell_size <= 0.0) {
        mesh->min_allowable_cell_size = 0.0;
        return;
    }

    if (mesh->max_level >= 0) {
        level = mesh->max_level;
    }
    if (mesh->min_dx > 0.0) {
        int min_dx_level = 0;
        double cell_size = root_cell_size;

        /* Match AMR's cap: a parent whose cell size is still above min_dx
         * may refine once more, so the finest allowed child can fall below it. */
        while (cell_size > mesh->min_dx) {
            cell_size *= 0.5;
            ++min_dx_level;
        }
        if (mesh->max_level < 0 || min_dx_level < level) {
            level = min_dx_level;
        }
    }
    if (level < 0) {
        level = 0;
    }
    mesh->min_allowable_cell_size = ldexp(root_cell_size, -level);
}

void prj_mesh_update_max_active_level(prj_mesh *mesh)
{
    int i;
    int local_max;

    if (mesh == 0) {
        return;
    }

    local_max = -1;
    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block *b = &mesh->blocks[i];

        if (b->active == 0) {
            continue;
        }
        if (b->level > local_max) {
            local_max = b->level;
        }
    }

#if defined(PRJ_ENABLE_MPI)
    MPI_Allreduce(&local_max, &mesh->max_active_level, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
    mesh->max_active_level = local_max;
#endif
}

static void prj_mesh_mark_cell_derived_range(prj_block *block,
    const prj_neighbor *slot, int value)
{
    int i;
    int j;
    int k;

    if (block == 0 || slot == 0 || block->cell_derived_done == 0) {
        return;
    }
    for (i = slot->recv_loc_start[0]; i < slot->recv_loc_end[0]; ++i) {
        for (j = slot->recv_loc_start[1]; j < slot->recv_loc_end[1]; ++j) {
            for (k = slot->recv_loc_start[2]; k < slot->recv_loc_end[2]; ++k) {
                block->cell_derived_done[IDX(i, j, k)] = value;
            }
        }
    }
}

void prj_mesh_update_cell_derived_mask(prj_mesh *mesh)
{
    int bidx;

    if (mesh == 0) {
        return;
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1 || block->cell_derived_done == 0) {
            continue;
        }
        memset(block->cell_derived_done, 0,
            (size_t)PRJ_BLOCK_NCELLS * sizeof(*block->cell_derived_done));
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    block->cell_derived_done[IDX(i, j, k)] = 1;
                }
            }
        }
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (block->id < 0 || block->active != 1) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &block->slot[n];
            int nid = slot->id;
            prj_block *neighbor;

            if (slot->rel_level != 0 || nid < 0 || nid >= mesh->nblocks ||
                mesh->blocks[nid].active != 1) {
                continue;
            }
            neighbor = &mesh->blocks[nid];
            prj_mesh_mark_cell_derived_range(neighbor, slot, 1);
        }
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        int n;

        if (block->id < 0 || block->active != 1) {
            continue;
        }
        for (n = 0; n < 56; ++n) {
            const prj_neighbor *slot = &block->slot[n];
            int nid = slot->id;
            prj_block *neighbor;

            if (slot->rel_level == 0 || nid < 0 || nid >= mesh->nblocks ||
                mesh->blocks[nid].active != 1) {
                continue;
            }
            neighbor = &mesh->blocks[nid];
            prj_mesh_mark_cell_derived_range(neighbor, slot, 0);
        }
    }
}

static int prj_mesh_is_active_block(const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1;
}

static unsigned long long prj_mesh_morton_part(unsigned int value)
{
    unsigned long long x = (unsigned long long)(value & 0x1fffffU);

    x = (x | (x << 32)) & 0x1f00000000ffffULL;
    x = (x | (x << 16)) & 0x1f0000ff0000ffULL;
    x = (x | (x << 8)) & 0x100f00f00f00f00fULL;
    x = (x | (x << 4)) & 0x10c30c30c30c30c3ULL;
    x = (x | (x << 2)) & 0x1249249249249249ULL;
    return x;
}

static unsigned long long prj_mesh_morton_key(int level, int ix, int iy, int iz)
{
    unsigned long long key;

    key = ((unsigned long long)(unsigned int)level) * 0x9e3779b97f4a7c15ULL;
    key ^= prj_mesh_morton_part((unsigned int)ix) << 2;
    key ^= prj_mesh_morton_part((unsigned int)iy) << 1;
    key ^= prj_mesh_morton_part((unsigned int)iz);
    return key;
}

static unsigned long long prj_mesh_morton_hash(unsigned long long key)
{
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key;
}

static int prj_mesh_level_axis_count(const prj_mesh *mesh, int level, int axis)
{
    double count;

    if (mesh == 0 || level < 0 || axis < 0 || axis >= 3 ||
        mesh->root_nx[axis] <= 0) {
        return 0;
    }
    count = ldexp((double)mesh->root_nx[axis], level);
    if (!isfinite(count) || count <= 0.0 || count > (double)INT_MAX) {
        return 0;
    }
    return (int)(count + 0.5);
}

static int prj_mesh_block_logical_coords(const prj_mesh *mesh,
    const prj_block *block, int coord[3])
{
    double coord_min[3];
    double coord_max[3];
    int axis;

    if (mesh == 0 || block == 0 || coord == 0 || block->level < 0) {
        return 0;
    }
    coord_min[0] = mesh->coord.x1min;
    coord_min[1] = mesh->coord.x2min;
    coord_min[2] = mesh->coord.x3min;
    coord_max[0] = mesh->coord.x1max;
    coord_max[1] = mesh->coord.x2max;
    coord_max[2] = mesh->coord.x3max;

    for (axis = 0; axis < 3; ++axis) {
        double extent = coord_max[axis] - coord_min[axis];
        int level_count = prj_mesh_level_axis_count(mesh, block->level, axis);
        double scale;
        double scaled;
        int c;

        if (extent <= 0.0 || level_count <= 0) {
            return 0;
        }
        scale = (double)level_count / extent;
        scaled = (block->xmin[axis] - coord_min[axis]) * scale;
        if (!isfinite(scaled) || scaled < -0.5 || scaled > (double)INT_MAX) {
            return 0;
        }
        c = (int)(scaled + 0.5);
        if (c < 0 || c >= level_count) {
            return 0;
        }
        coord[axis] = c;
    }
    return 1;
}

static int prj_mesh_morton_lookup_capacity(int count)
{
    int capacity = 16;
    int target;

    if (count <= 0) {
        return 0;
    }
    if (count > INT_MAX / 4) {
        return 0;
    }
    target = count * 4;
    while (capacity < target) {
        if (capacity > INT_MAX / 2) {
            return 0;
        }
        capacity *= 2;
    }
    return capacity;
}

static int prj_mesh_morton_lookup_insert(prj_mesh *mesh, const prj_block *block)
{
    int coord[3];
    unsigned long long key;
    unsigned long long hash;
    unsigned int mask;
    int probe;

    if (mesh == 0 || block == 0 || mesh->morton_lookup == 0 ||
        mesh->morton_lookup_capacity <= 0) {
        return 1;
    }
    if (!prj_mesh_block_logical_coords(mesh, block, coord)) {
        return 1;
    }

    key = prj_mesh_morton_key(block->level, coord[0], coord[1], coord[2]);
    hash = prj_mesh_morton_hash(key);
    mask = (unsigned int)(mesh->morton_lookup_capacity - 1);
    for (probe = 0; probe < mesh->morton_lookup_capacity; ++probe) {
        unsigned int idx = (unsigned int)(hash + (unsigned long long)probe) & mask;
        prj_morton_lookup_entry *entry = &mesh->morton_lookup[idx];

        if (entry->occupied == 0 ||
            (entry->level == block->level &&
             entry->coord[0] == coord[0] &&
             entry->coord[1] == coord[1] &&
             entry->coord[2] == coord[2])) {
            entry->occupied = 1;
            entry->level = block->level;
            entry->coord[0] = coord[0];
            entry->coord[1] = coord[1];
            entry->coord[2] = coord[2];
            entry->id = block->id;
            entry->key = key;
            return 0;
        }
    }
    return 1;
}

int prj_mesh_rebuild_morton_lookup(prj_mesh *mesh)
{
    int active_count;
    int capacity;
    int i;

    if (mesh == 0) {
        return 1;
    }

    active_count = 0;
    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_mesh_is_active_block(&mesh->blocks[i])) {
            active_count += 1;
        }
    }
    mesh->morton_lookup_count = active_count;
    if (active_count <= 0) {
        return 0;
    }

    capacity = prj_mesh_morton_lookup_capacity(active_count);
    if (capacity <= 0) {
        mesh->morton_lookup_count = 0;
        return 1;
    }
    if (mesh->morton_lookup_capacity != capacity) {
        prj_morton_lookup_entry *lookup =
            (prj_morton_lookup_entry *)calloc((size_t)capacity, sizeof(*lookup));

        if (lookup == 0) {
            mesh->morton_lookup_count = 0;
            return 1;
        }
        free(mesh->morton_lookup);
        mesh->morton_lookup = lookup;
        mesh->morton_lookup_capacity = capacity;
    } else {
        memset(mesh->morton_lookup, 0,
            (size_t)mesh->morton_lookup_capacity * sizeof(*mesh->morton_lookup));
    }

    for (i = 0; i < mesh->nblocks; ++i) {
        if (prj_mesh_is_active_block(&mesh->blocks[i]) &&
            prj_mesh_morton_lookup_insert(mesh, &mesh->blocks[i]) != 0) {
            mesh->morton_lookup_count = 0;
            return 1;
        }
    }
    return 0;
}

int prj_mesh_morton_lookup_block(const prj_mesh *mesh, int level, int ix, int iy, int iz)
{
    unsigned long long key;
    unsigned long long hash;
    unsigned int mask;
    int probe;

    if (mesh == 0 || mesh->morton_lookup == 0 ||
        mesh->morton_lookup_capacity <= 0 || mesh->morton_lookup_count <= 0 ||
        level < 0 || ix < 0 || iy < 0 || iz < 0) {
        return -1;
    }

    key = prj_mesh_morton_key(level, ix, iy, iz);
    hash = prj_mesh_morton_hash(key);
    mask = (unsigned int)(mesh->morton_lookup_capacity - 1);
    for (probe = 0; probe < mesh->morton_lookup_capacity; ++probe) {
        unsigned int idx = (unsigned int)(hash + (unsigned long long)probe) & mask;
        const prj_morton_lookup_entry *entry = &mesh->morton_lookup[idx];

        if (entry->occupied == 0) {
            return -1;
        }
        if (entry->level == level &&
            entry->coord[0] == ix &&
            entry->coord[1] == iy &&
            entry->coord[2] == iz) {
            return entry->id;
        }
    }
    return -1;
}

int prj_mesh_init(prj_mesh *mesh, int root_nx1, int root_nx2, int root_nx3, int max_level, const prj_coord *coord)
{
    int i;
    int j;
    int k;
    int id;
    int nroot;
    int capacity;
    double block_dx[3];
    double saved_amr_refine_thresh[PRJ_AMR_N];
    double saved_amr_derefine_thresh[PRJ_AMR_N];
    int saved_amr_estimator[PRJ_AMR_N];
    int saved_amr_lohner_var[PRJ_AMR_N];
    int saved_amr_fractional_jump_var[PRJ_AMR_N];
    int saved_amr_criterion_set[PRJ_AMR_N];
    double saved_amr_lohner_eps[PRJ_AMR_N];
    int saved_use_amr_angular_resolution_limit;
    int saved_use_BJ;
    double saved_amr_init_scale_factor;
    double saved_E_floor;
    double saved_min_dx;
    int saved_max_blocks;
    int amr_idx;

    if (mesh == 0 || coord == 0 || root_nx1 <= 0 || root_nx2 <= 0 || root_nx3 <= 0) {
        return 1;
    }

    for (amr_idx = 0; amr_idx < PRJ_AMR_N; ++amr_idx) {
        saved_amr_refine_thresh[amr_idx] = mesh->amr_refine_thresh[amr_idx];
        saved_amr_derefine_thresh[amr_idx] = mesh->amr_derefine_thresh[amr_idx];
        saved_amr_estimator[amr_idx] = mesh->amr_estimator[amr_idx];
        saved_amr_lohner_var[amr_idx] = mesh->amr_lohner_var[amr_idx];
        saved_amr_fractional_jump_var[amr_idx] = mesh->amr_fractional_jump_var[amr_idx];
        saved_amr_lohner_eps[amr_idx] = mesh->amr_lohner_eps[amr_idx];
        saved_amr_criterion_set[amr_idx] = mesh->amr_criterion_set[amr_idx];
    }
    saved_use_amr_angular_resolution_limit = mesh->use_amr_angular_resolution_limit;
    saved_use_BJ = mesh->use_BJ;
    saved_amr_init_scale_factor = mesh->amr_init_scale_factor;
    saved_E_floor = mesh->E_floor;
    saved_min_dx = mesh->min_dx;
    saved_max_blocks = mesh->max_blocks;

    mesh->nblocks = 0;
    mesh->nblocks_max = 0;
    mesh->max_level = max_level;
    mesh->min_dx = saved_min_dx;
    mesh->max_active_level = -1;
    mesh->root_nx[0] = root_nx1;
    mesh->root_nx[1] = root_nx2;
    mesh->root_nx[2] = root_nx3;
    mesh->x_com[0] = 0.0;
    mesh->x_com[1] = 0.0;
    mesh->x_com[2] = 0.0;
    mesh->coord = *coord;
    prj_mesh_update_min_allowable_cell_size(mesh);
    for (amr_idx = 0; amr_idx < PRJ_AMR_N; ++amr_idx) {
        mesh->amr_refine_thresh[amr_idx] = saved_amr_refine_thresh[amr_idx];
        mesh->amr_derefine_thresh[amr_idx] = saved_amr_derefine_thresh[amr_idx];
        mesh->amr_estimator[amr_idx] = saved_amr_estimator[amr_idx];
        mesh->amr_lohner_var[amr_idx] = saved_amr_lohner_var[amr_idx];
        mesh->amr_fractional_jump_var[amr_idx] = saved_amr_fractional_jump_var[amr_idx];
        mesh->amr_lohner_eps[amr_idx] = saved_amr_lohner_eps[amr_idx];
        mesh->amr_criterion_set[amr_idx] = saved_amr_criterion_set[amr_idx];
    }
    mesh->use_amr_angular_resolution_limit = saved_use_amr_angular_resolution_limit;
    mesh->use_BJ = saved_use_BJ;
    mesh->amr_init_scale_factor = saved_amr_init_scale_factor;
    mesh->E_floor = saved_E_floor;
    mesh->max_blocks = saved_max_blocks;
    mesh->amr_init_refine_fn = 0;
    mesh->amr_init_refine_userdata = 0;
    mesh->blocks = 0;
    free(mesh->morton_lookup);
    mesh->morton_lookup = 0;
    mesh->morton_lookup_capacity = 0;
    mesh->morton_lookup_count = 0;

    nroot = root_nx1 * root_nx2 * root_nx3;
    capacity = mesh->max_blocks > 0 ? mesh->max_blocks : 65536;
    if (capacity < nroot) {
        capacity = nroot;
    }
    mesh->max_blocks = capacity;

    mesh->blocks = (prj_block *)malloc((size_t)capacity * sizeof(*mesh->blocks));
    if (mesh->blocks == 0) {
        return 2;
    }
    mesh->nblocks_max = capacity;

    for (i = 0; i < capacity; ++i) {
        prj_block_init_empty(&mesh->blocks[i]);
    }

    block_dx[0] = (coord->x1max - coord->x1min) / (double)root_nx1;
    block_dx[1] = (coord->x2max - coord->x2min) / (double)root_nx2;
    block_dx[2] = (coord->x3max - coord->x3min) / (double)root_nx3;

    id = 0;
    for (i = 0; i < root_nx1; ++i) {
        for (j = 0; j < root_nx2; ++j) {
            for (k = 0; k < root_nx3; ++k) {
                prj_block *b = &mesh->blocks[id];

                b->id = id;
                b->rank = 0;
                b->level = 0;
                b->active = 1;
                b->refine_flag = 0;
                b->can_refine = 1;
                b->xmin[0] = coord->x1min + (double)i * block_dx[0];
                b->xmax[0] = b->xmin[0] + block_dx[0];
                b->xmin[1] = coord->x2min + (double)j * block_dx[1];
                b->xmax[1] = b->xmin[1] + block_dx[1];
                b->xmin[2] = coord->x3min + (double)k * block_dx[2];
                b->xmax[2] = b->xmin[2] + block_dx[2];
                b->dx[0] = block_dx[0] / (double)PRJ_BLOCK_SIZE;
                b->dx[1] = block_dx[1] / (double)PRJ_BLOCK_SIZE;
                b->dx[2] = block_dx[2] / (double)PRJ_BLOCK_SIZE;
                prj_block_setup_geometry(b, coord);
                prj_block_update_can_refine(b, mesh);
                if (prj_block_alloc_data(b) != 0) {
                    prj_mesh_destroy(mesh);
                    return 3;
                }
                id += 1;
            }
        }
    }
    mesh->nblocks = nroot;
    prj_mesh_update_r_com(mesh);
    prj_mesh_update_max_active_level(mesh);

    for (i = 0; i < root_nx1; ++i) {
        for (j = 0; j < root_nx2; ++j) {
            for (k = 0; k < root_nx3; ++k) {
                prj_block *b = &mesh->blocks[(i * root_nx2 + j) * root_nx3 + k];
                int ox;
                int oy;
                int oz;

                for (ox = -1; ox <= 1; ++ox) {
                    for (oy = -1; oy <= 1; ++oy) {
                        for (oz = -1; oz <= 1; ++oz) {
                            int ni = i + ox;
                            int nj = j + oy;
                            int nk = k + oz;
                            int slot_index = prj_neighbor_slot_index(ox, oy, oz);

                            if (slot_index < 0) {
                                continue;
                            }
                            if (ni < 0 || ni >= root_nx1 || nj < 0 || nj >= root_nx2 || nk < 0 || nk >= root_nx3) {
                                continue;
                            }

                            id = (ni * root_nx2 + nj) * root_nx3 + nk;
                            b->slot[slot_index].id = id;
                            b->slot[slot_index].rank = 0;
                            b->slot[slot_index].xmin[0] = mesh->blocks[id].xmin[0];
                            b->slot[slot_index].xmin[1] = mesh->blocks[id].xmin[1];
                            b->slot[slot_index].xmin[2] = mesh->blocks[id].xmin[2];
                            b->slot[slot_index].xmax[0] = mesh->blocks[id].xmax[0];
                            b->slot[slot_index].xmax[1] = mesh->blocks[id].xmax[1];
                            b->slot[slot_index].xmax[2] = mesh->blocks[id].xmax[2];
                            b->slot[slot_index].dx[0] = mesh->blocks[id].dx[0];
                            b->slot[slot_index].dx[1] = mesh->blocks[id].dx[1];
                            b->slot[slot_index].dx[2] = mesh->blocks[id].dx[2];
                            prj_neighbor_compute_geometry(b, &mesh->blocks[id], &b->slot[slot_index]);
                        }
                    }
                }
            }
        }
    }


    prj_mesh_update_cell_derived_mask(mesh);
    if (prj_mesh_rebuild_morton_lookup(mesh) != 0) {
        prj_mesh_destroy(mesh);
        return 4;
    }
    return 0;
}

void prj_mesh_destroy(prj_mesh *mesh)
{
    int i;

    if (mesh == 0) {
        return;
    }

    for (i = 0; i < mesh->nblocks; ++i) {
        prj_block_free_data(&mesh->blocks[i]);
    }
    free(mesh->blocks);
    free(mesh->morton_lookup);
    mesh->blocks = 0;
    mesh->morton_lookup = 0;
    mesh->nblocks = 0;
    mesh->nblocks_max = 0;
    mesh->morton_lookup_capacity = 0;
    mesh->morton_lookup_count = 0;
    mesh->max_active_level = -1;
}
