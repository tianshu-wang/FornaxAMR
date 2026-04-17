#include <math.h>

#include "prj.h"

static prj_mesh *prj_riemann_flux_mesh = 0;

double prj_riemann_min_double(double a, double b)
{
    return a < b ? a : b;
}

double prj_riemann_max_double(double a, double b)
{
    return a > b ? a : b;
}

static double prj_abs_double(double x)
{
    return x < 0.0 ? -x : x;
}

static void prj_riemann_state(const double *W, double pressure, double gamma,
    const prj_eos *eos, double *U, double *F, double *cs)
{
    double rho;
    double v1;
    double v2;
    double v3;
    double etot;
    int v;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        U[v] = 0.0;
        F[v] = 0.0;
    }

    prj_eos_prim2cons((prj_eos *)eos, (double *)W, U);

    rho = W[PRJ_PRIM_RHO];
    v1 = W[PRJ_PRIM_V1];
    v2 = W[PRJ_PRIM_V2];
    v3 = W[PRJ_PRIM_V3];
    etot = U[PRJ_CONS_ETOT];

    *cs = sqrt(gamma * pressure / rho);

    F[PRJ_CONS_RHO] = rho * v1;
    F[PRJ_CONS_MOM1] = rho * v1 * v1 + pressure;
    F[PRJ_CONS_MOM2] = rho * v1 * v2;
    F[PRJ_CONS_MOM3] = rho * v1 * v3;
    F[PRJ_CONS_ETOT] = (etot + pressure) * v1;
    F[PRJ_CONS_YE] = rho * W[PRJ_PRIM_YE] * v1;
}

void prj_riemann_hlle(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double *flux, double v_face[3])
{
    double UL[PRJ_NVAR_CONS];
    double UR[PRJ_NVAR_CONS];
    double FL[PRJ_NVAR_CONS];
    double FR[PRJ_NVAR_CONS];
    double csL;
    double csR;
    double SL;
    double SR;
    int v;

    prj_riemann_state(WL, pL, gL, eos, UL, FL, &csL);
    prj_riemann_state(WR, pR, gR, eos, UR, FR, &csR);

    SL = prj_riemann_min_double(WL[PRJ_PRIM_V1] - csL, WR[PRJ_PRIM_V1] - csR);
    SR = prj_riemann_max_double(WL[PRJ_PRIM_V1] + csL, WR[PRJ_PRIM_V1] + csR);

    if (0.0 <= SL) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = FL[v];
        }
        if (v_face != 0) {
            v_face[0] = WL[PRJ_PRIM_V1];
            v_face[1] = WL[PRJ_PRIM_V2];
            v_face[2] = WL[PRJ_PRIM_V3];
        }
        return;
    }
    if (SR <= 0.0) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = FR[v];
        }
        if (v_face != 0) {
            v_face[0] = WR[PRJ_PRIM_V1];
            v_face[1] = WR[PRJ_PRIM_V2];
            v_face[2] = WR[PRJ_PRIM_V3];
        }
        return;
    }

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        flux[v] = (SR * FL[v] - SL * FR[v] + SL * SR * (UR[v] - UL[v])) / (SR - SL);
    }
    if (v_face != 0) {
        double rho_hll = (SR * UR[PRJ_CONS_RHO] - SL * UL[PRJ_CONS_RHO] -
            (FR[PRJ_CONS_RHO] - FL[PRJ_CONS_RHO])) / (SR - SL);
        double mom1_hll = (SR * UR[PRJ_CONS_MOM1] - SL * UL[PRJ_CONS_MOM1] -
            (FR[PRJ_CONS_MOM1] - FL[PRJ_CONS_MOM1])) / (SR - SL);
        double mom2_hll = (SR * UR[PRJ_CONS_MOM2] - SL * UL[PRJ_CONS_MOM2] -
            (FR[PRJ_CONS_MOM2] - FL[PRJ_CONS_MOM2])) / (SR - SL);
        double mom3_hll = (SR * UR[PRJ_CONS_MOM3] - SL * UL[PRJ_CONS_MOM3] -
            (FR[PRJ_CONS_MOM3] - FL[PRJ_CONS_MOM3])) / (SR - SL);
        if (rho_hll != 0.0) {
            v_face[0] = mom1_hll / rho_hll;
            v_face[1] = mom2_hll / rho_hll;
            v_face[2] = mom3_hll / rho_hll;
        } else {
            v_face[0] = 0.0;
            v_face[1] = 0.0;
            v_face[2] = 0.0;
        }
    }
}

void prj_riemann_hllc(const double *WL, const double *WR,
    double pL, double pR, double gL, double gR,
    const prj_eos *eos, double *flux, double v_face[3])
{
    double UL[PRJ_NVAR_CONS];
    double UR[PRJ_NVAR_CONS];
    double FL[PRJ_NVAR_CONS];
    double FR[PRJ_NVAR_CONS];
    double Ustar[PRJ_NVAR_CONS];
    double csL;
    double csR;
    double SL;
    double SR;
    double SM;
    double rho_star;
    double e_over_rho;
    int v;

    prj_riemann_state(WL, pL, gL, eos, UL, FL, &csL);
    prj_riemann_state(WR, pR, gR, eos, UR, FR, &csR);

    SL = prj_riemann_min_double(WL[PRJ_PRIM_V1] - csL, WR[PRJ_PRIM_V1] - csR);
    SR = prj_riemann_max_double(WL[PRJ_PRIM_V1] + csL, WR[PRJ_PRIM_V1] + csR);

    /* Toro (2009), Eq. 10.37: contact wave speed from Rankine-Hugoniot conditions. */
    SM = (pR - pL +
        UL[PRJ_CONS_MOM1] * (SL - WL[PRJ_PRIM_V1]) -
        UR[PRJ_CONS_MOM1] * (SR - WR[PRJ_PRIM_V1])) /
        (UL[PRJ_CONS_RHO] * (SL - WL[PRJ_PRIM_V1]) -
            UR[PRJ_CONS_RHO] * (SR - WR[PRJ_PRIM_V1]));

    if (0.0 <= SL) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = FL[v];
        }
        if (v_face != 0) {
            v_face[0] = WL[PRJ_PRIM_V1];
            v_face[1] = WL[PRJ_PRIM_V2];
            v_face[2] = WL[PRJ_PRIM_V3];
        }
        return;
    }
    if (SR <= 0.0) {
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = FR[v];
        }
        if (v_face != 0) {
            v_face[0] = WR[PRJ_PRIM_V1];
            v_face[1] = WR[PRJ_PRIM_V2];
            v_face[2] = WR[PRJ_PRIM_V3];
        }
        return;
    }

    if (v_face != 0) {
        /* HLLC star state: contact normal speed SM, transverse velocities preserved from upwind side. */
        v_face[0] = SM;
        if (0.0 <= SM) {
            v_face[1] = WL[PRJ_PRIM_V2];
            v_face[2] = WL[PRJ_PRIM_V3];
        } else {
            v_face[1] = WR[PRJ_PRIM_V2];
            v_face[2] = WR[PRJ_PRIM_V3];
        }
    }
    if (0.0 <= SM) {
        /* Toro (2009), Eqs. 10.33-10.36: left star state. */
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            Ustar[v] = UL[v];
        }
        rho_star = UL[PRJ_CONS_RHO] * (SL - WL[PRJ_PRIM_V1]) / (SL - SM);
        Ustar[PRJ_CONS_RHO] = rho_star;
        Ustar[PRJ_CONS_MOM1] = rho_star * SM;
        Ustar[PRJ_CONS_MOM2] = rho_star * WL[PRJ_PRIM_V2];
        Ustar[PRJ_CONS_MOM3] = rho_star * WL[PRJ_PRIM_V3];
        e_over_rho = UL[PRJ_CONS_ETOT] / UL[PRJ_CONS_RHO] +
            (SM - WL[PRJ_PRIM_V1]) * (SM + pL / (UL[PRJ_CONS_RHO] * (SL - WL[PRJ_PRIM_V1])));
        Ustar[PRJ_CONS_ETOT] = rho_star * e_over_rho;
        Ustar[PRJ_CONS_YE] = rho_star * WL[PRJ_PRIM_YE];
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = FL[v] + SL * (Ustar[v] - UL[v]);
        }
    } else {
        /* Toro (2009), Eqs. 10.33-10.36 with right star state. */
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            Ustar[v] = UR[v];
        }
        rho_star = UR[PRJ_CONS_RHO] * (SR - WR[PRJ_PRIM_V1]) / (SR - SM);
        Ustar[PRJ_CONS_RHO] = rho_star;
        Ustar[PRJ_CONS_MOM1] = rho_star * SM;
        Ustar[PRJ_CONS_MOM2] = rho_star * WR[PRJ_PRIM_V2];
        Ustar[PRJ_CONS_MOM3] = rho_star * WR[PRJ_PRIM_V3];
        e_over_rho = UR[PRJ_CONS_ETOT] / UR[PRJ_CONS_RHO] +
            (SM - WR[PRJ_PRIM_V1]) * (SM + pR / (UR[PRJ_CONS_RHO] * (SR - WR[PRJ_PRIM_V1])));
        Ustar[PRJ_CONS_ETOT] = rho_star * e_over_rho;
        Ustar[PRJ_CONS_YE] = rho_star * WR[PRJ_PRIM_YE];
        for (v = 0; v < PRJ_NVAR_CONS; ++v) {
            flux[v] = FR[v] + SR * (Ustar[v] - UR[v]);
        }
    }
}

static int prj_blocks_overlap_open(double amin, double amax, double bmin, double bmax)
{
    const double tol = 1.0e-12;
    return prj_riemann_min_double(amax, bmax) - prj_riemann_max_double(amin, bmin) > tol;
}

static double prj_overlap_length(double amin, double amax, double bmin, double bmax)
{
    double overlap = prj_riemann_min_double(amax, bmax) - prj_riemann_max_double(amin, bmin);

    return overlap > 0.0 ? overlap : 0.0;
}

void prj_riemann_set_mesh(prj_mesh *mesh)
{
    prj_riemann_flux_mesh = mesh;
}

int prj_riemann_detect_shock(const double *WL, const double *WR, double pL, double pR)
{
    double pressure_ratio;
    double compression[3];
    int best_dir;

    if (WL[PRJ_PRIM_EINT] <= 0.0 || WR[PRJ_PRIM_EINT] <= 0.0 ||
        WL[PRJ_PRIM_RHO] <= 0.0 || WR[PRJ_PRIM_RHO] <= 0.0) {
        return -1;
    }

    pressure_ratio = prj_riemann_max_double(pL, pR) /
        prj_riemann_max_double(prj_riemann_min_double(pL, pR), 1.0e-12);
    if (pressure_ratio < 1.5 || WR[PRJ_PRIM_V1] - WL[PRJ_PRIM_V1] >= 0.0) {
        return -1;
    }

    compression[0] = prj_abs_double(WR[PRJ_PRIM_V1] - WL[PRJ_PRIM_V1]);
    compression[1] = prj_abs_double(WR[PRJ_PRIM_V2] - WL[PRJ_PRIM_V2]);
    compression[2] = prj_abs_double(WR[PRJ_PRIM_V3] - WL[PRJ_PRIM_V3]);
    best_dir = X1DIR;
    if (compression[X2DIR] > compression[best_dir]) {
        best_dir = X2DIR;
    }
    if (compression[X3DIR] > compression[best_dir]) {
        best_dir = X3DIR;
    }
    if (compression[best_dir] < 2.0 * prj_riemann_max_double(
            compression[(best_dir + 1) % 3],
            compression[(best_dir + 2) % 3])) {
        return -1;
    }
    return best_dir;
}

void prj_riemann_flux_send(prj_block *block)
{
    int n;

    if (prj_riemann_flux_mesh == 0 || block == 0) {
        return;
    }

    for (n = 0; n < 56; ++n) {
        int nid = block->slot[n].id;

        if (nid >= 0 && nid < prj_riemann_flux_mesh->nblocks) {
            prj_block *neighbor = &prj_riemann_flux_mesh->blocks[nid];
            int axis = -1;
            int side = -1;
            int d;

            if (neighbor->id < 0 || neighbor->active != 1) {
                continue;
            }
            for (d = 0; d < 3; ++d) {
                if (prj_abs_double(block->xmax[d] - neighbor->xmin[d]) < 1.0e-12) {
                    if (axis >= 0) {
                        axis = -2;
                        break;
                    }
                    axis = d;
                    side = 1;
                } else if (prj_abs_double(neighbor->xmax[d] - block->xmin[d]) < 1.0e-12) {
                    if (axis >= 0) {
                        axis = -2;
                        break;
                    }
                    axis = d;
                    side = 0;
                }
            }
            if (axis < 0) {
                continue;
            }
            for (d = 0; d < 3; ++d) {
                if (d == axis) {
                    continue;
                }
                if (!prj_blocks_overlap_open(block->xmin[d], block->xmax[d], neighbor->xmin[d], neighbor->xmax[d])) {
                    axis = -1;
                    break;
                }
            }
            if (axis < 0) {
                continue;
            }
            if (neighbor->dx[axis] <= block->dx[axis]) {
                continue;
            }
            if (neighbor->rank == block->rank) {
                int ratio = (int)(neighbor->dx[axis] / block->dx[axis] + 0.5);
                int it0;
                int it1;

                if (ratio < 1) {
                    ratio = 1;
                }

                for (it0 = 0; it0 < PRJ_BLOCK_SIZE; ++it0) {
                    for (it1 = 0; it1 < PRJ_BLOCK_SIZE; ++it1) {
                        int coarse_face[3] = {0, 0, 0};
                        int tan0 = (axis + 1) % 3;
                        int tan1 = (axis + 2) % 3;
                        double coarse_min0;
                        double coarse_max0;
                        double coarse_min1;
                        double coarse_max1;
                        int fine_i0_start;
                        int fine_i0_end;
                        int fine_i1_start;
                        int fine_i1_end;
                        int r0;
                        int r1;

                        coarse_face[axis] = side == 1 ? 0 : PRJ_BLOCK_SIZE;
                        coarse_face[tan0] = it0;
                        coarse_face[tan1] = it1;
                        coarse_min0 = neighbor->xmin[tan0] + (double)it0 * neighbor->dx[tan0];
                        coarse_max0 = coarse_min0 + neighbor->dx[tan0];
                        coarse_min1 = neighbor->xmin[tan1] + (double)it1 * neighbor->dx[tan1];
                        coarse_max1 = coarse_min1 + neighbor->dx[tan1];

                        if (!prj_blocks_overlap_open(coarse_min0, coarse_max0, block->xmin[tan0], block->xmax[tan0]) ||
                            !prj_blocks_overlap_open(coarse_min1, coarse_max1, block->xmin[tan1], block->xmax[tan1])) {
                            continue;
                        }

                        fine_i0_start = (int)((coarse_min0 - block->xmin[tan0]) / block->dx[tan0]);
                        fine_i0_end = (int)((coarse_max0 - block->xmin[tan0]) / block->dx[tan0]);
                        fine_i1_start = (int)((coarse_min1 - block->xmin[tan1]) / block->dx[tan1]);
                        fine_i1_end = (int)((coarse_max1 - block->xmin[tan1]) / block->dx[tan1]);
                        if (fine_i0_start < 0) {
                            fine_i0_start = 0;
                        }
                        if (fine_i1_start < 0) {
                            fine_i1_start = 0;
                        }
                        if (fine_i0_end >= PRJ_BLOCK_SIZE) {
                            fine_i0_end = PRJ_BLOCK_SIZE - 1;
                        }
                        if (fine_i1_end >= PRJ_BLOCK_SIZE) {
                            fine_i1_end = PRJ_BLOCK_SIZE - 1;
                        }
                        {
                            int v;

                            for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                                double restricted_flux = 0.0;
                                double weight_sum = 0.0;

                                for (r0 = fine_i0_start; r0 <= fine_i0_end; ++r0) {
                                    for (r1 = fine_i1_start; r1 <= fine_i1_end; ++r1) {
                                        int fine_face[3] = {0, 0, 0};
                                        double fine_min0 = block->xmin[tan0] + (double)r0 * block->dx[tan0];
                                        double fine_max0 = fine_min0 + block->dx[tan0];
                                        double fine_min1 = block->xmin[tan1] + (double)r1 * block->dx[tan1];
                                        double fine_max1 = fine_min1 + block->dx[tan1];
                                        double overlap0 = prj_overlap_length(coarse_min0, coarse_max0, fine_min0, fine_max0);
                                        double overlap1 = prj_overlap_length(coarse_min1, coarse_max1, fine_min1, fine_max1);
                                        double weight = overlap0 * overlap1;

                                        if (weight <= 0.0) {
                                            continue;
                                        }
                                        fine_face[axis] = side == 1 ? PRJ_BLOCK_SIZE : 0;
                                        fine_face[tan0] = r0;
                                        fine_face[tan1] = r1;
                                        restricted_flux +=
                                            block->flux[axis][VIDX(v, fine_face[0], fine_face[1], fine_face[2])] * weight;
                                        weight_sum += weight;
                                    }
                                }
                                if (weight_sum > 0.0) {
                                    neighbor->flux[axis][VIDX(v, coarse_face[0], coarse_face[1], coarse_face[2])] =
                                        restricted_flux / weight_sum;
                                }
                            }
                        }
                    }
                }
            } else {
                /* MPI flux exchange buffer path not implemented yet. */
            }
        }
    }
}
