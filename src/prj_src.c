#include <math.h>

#include "prj.h"

void prj_src_geom(prj_eos *eos, double *W, double *dUdt)
{
    (void)eos;
    (void)W;
    (void)dUdt;
}

void prj_src_user(prj_eos *eos, double *W, double *dUdt)
{
    (void)eos;
    (void)W;
    (void)dUdt;
}

void prj_src_monopole_gravity(const prj_rad *rad, const prj_block *block,
    const prj_grav *grav, double *restrict W, double *restrict dUdt)
{
    int i;
    int j;
    int k;

    if (block == 0 || block->id < 0 || block->active != 1 ||
        grav == 0 || W == 0 || dUdt == 0) {
        return;
    }
    if (block->lapse == 0 || block->grav[0] == 0 || block->grav[1] == 0 ||
        block->grav[2] == 0) {
        return;
    }
    if (block->v_riemann[0] == 0 || block->v_riemann[1] == 0 || block->v_riemann[2] == 0) {
        return;
    }
#if !PRJ_USE_RADIATION_FSA
    (void)rad;
#endif

    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                int cell_idx = IDX(i, j, k);
                double g1;
                double g2;
                double g3;
                double rho;

                g1 = block->grav[0][cell_idx];
                g2 = block->grav[1][cell_idx];
                g3 = block->grav[2][cell_idx];
                rho = W[WIDX(PRJ_PRIM_RHO, i, j, k)];
                dUdt[VIDX(PRJ_CONS_MOM1, i, j, k)] += rho * g1;
                dUdt[VIDX(PRJ_CONS_MOM2, i, j, k)] += rho * g2;
                dUdt[VIDX(PRJ_CONS_MOM3, i, j, k)] += rho * g3;
                {
                    double v_avg1;
                    double v_avg2;
                    double v_avg3;

                    v_avg1 = (block->v_riemann[X1DIR][VRIDX(0, i, j, k)] +
                              block->v_riemann[X1DIR][VRIDX(0, i + 1, j, k)]) * 0.5;
                    v_avg2 = (block->v_riemann[X2DIR][VRIDX(1, i, j, k)] +
                              block->v_riemann[X2DIR][VRIDX(1, i, j + 1, k)]) * 0.5;
                    v_avg3 = (block->v_riemann[X3DIR][VRIDX(2, i, j, k)] +
                              block->v_riemann[X3DIR][VRIDX(2, i, j, k + 1)]) * 0.5;

                    dUdt[VIDX(PRJ_CONS_ETOT, i, j, k)] +=
                        rho * (v_avg1 * g1 + v_avg2 * g2 + v_avg3 * g3);
                }
#if PRJ_USE_RADIATION_M1
                {
                    double lapse = block->lapse[cell_idx];
                    int field;
                    int group;

                    double inv_c2 = 1.0 / (PRJ_CLIGHT * PRJ_CLIGHT);
                    for (field = 0; field < PRJ_NRAD; ++field) {
                        for (group = 0; group < PRJ_NEGROUP; ++group) {
                            double E = W[WIDX(PRJ_PRIM_RAD_E(field, group), i, j, k)];
                            double F1 = W[WIDX(PRJ_PRIM_RAD_F1(field, group), i, j, k)];
                            double F2 = W[WIDX(PRJ_PRIM_RAD_F2(field, group), i, j, k)];
                            double F3 = W[WIDX(PRJ_PRIM_RAD_F3(field, group), i, j, k)];

                            dUdt[VIDX(PRJ_CONS_RAD_F1(field, group), i, j, k)] += lapse * E * g1;
                            dUdt[VIDX(PRJ_CONS_RAD_F2(field, group), i, j, k)] += lapse * E * g2;
                            dUdt[VIDX(PRJ_CONS_RAD_F3(field, group), i, j, k)] += lapse * E * g3;
                            dUdt[VIDX(PRJ_CONS_RAD_E(field, group), i, j, k)] +=
                                lapse * (g1 * F1 + g2 * F2 + g3 * F3) * inv_c2;
                        }
                    }
                }
#elif PRJ_USE_RADIATION_FSA
                if (rad != 0) {
                    double lapse = block->lapse[cell_idx];
                    double inv_c = 1.0 / PRJ_CLIGHT;
                    int field;
                    int group;
                    int angle;

                    for (field = 0; field < PRJ_NRAD; ++field) {
                        for (group = 0; group < PRJ_NEGROUP; ++group) {
                            for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                                double n[3];
                                double a_dot_n;
                                double J = W[WIDX(PRJ_PRIM_RAD_I(field, group, angle), i, j, k)];

                                prj_rad_fsa_rotated_angle_dir(rad, block, angle, i, j, k, n);
                                a_dot_n = g1 * n[0] + g2 * n[1] + g3 * n[2];
                                dUdt[VIDX(PRJ_CONS_RAD_I(field, group, angle), i, j, k)] -=
                                    lapse * a_dot_n * inv_c * J;
                            }
                        }
                    }
                }
#endif
            }
        }
    }
}

/* Velocity-gradient radiation source.  The M1 branch applies the O(v/c)
 * "source-like" piece of the SR redshift (Eq. 12a/12b of the comoving-frame
 * mixed-frame moment equations), while the FSA branch applies the first
 * angular-cell-integrated intensity-equation RHS term -α (n·∇v·n) J:
 *
 *     ∂_t E_g   -= (∂_j v_i) P^{ji}_g                          (Eq. 12a piece)
 *     ∂_t F_{gj} -= (∂_j v_i) F_{gi}                            (Eq. 12b piece)
 *     ∂_t J_{g,a} -= α (n_j ∂_j v_i n_i) J_{g,a}                 (FSA piece)
 *
 * Indices: i,j ∈ {1,2,3}; sums over i (and i,j for the energy term).  ∂_j v_i is
 * built by central differencing the face-centred Riemann velocities stored in
 * block->v_riemann[face_dir][component, ...] during the most recent flux update,
 * so the velocity field used here is the same one the hydro fluxes saw. */
void prj_src_radiation_vel_grad(const prj_rad *rad, const prj_block *block,
    double *restrict W, double *restrict dUdt)
{
#if PRJ_USE_RADIATION_M1
    int i;
    int j;
    int k;

    if (block == 0 || block->id < 0 || block->active != 1 || W == 0 || dUdt == 0) {
        return;
    }
    if (block->v_riemann[0] == 0 || block->v_riemann[1] == 0 || block->v_riemann[2] == 0) {
        return;
    }

    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                double dvdx[3][3]; /* dvdx[jdir][icomp] = ∂_jdir v_icomp */
                int jdir;
                int icomp;
                int field;
                int group;
                double inv_dx[3];
                double lapse = block->lapse != 0 ? block->lapse[IDX(i, j, k)] : 1.0;

                inv_dx[0] = 1.0 / block->dx[0];
                inv_dx[1] = 1.0 / block->dx[1];
                inv_dx[2] = 1.0 / block->dx[2];

                /* ∂_jdir v_icomp at cell centre, from the two normal-direction
                 * faces bracketing the cell along jdir.  v_riemann[jdir] is
                 * laid out as [icomp * NCELLS + IDX(face)] with the face index
                 * along jdir running 0..PRJ_BLOCK_SIZE (left/right faces). */
                for (jdir = 0; jdir < 3; ++jdir) {
                    for (icomp = 0; icomp < 3; ++icomp) {
                        int il = i;
                        int jl = j;
                        int kl = k;
                        int ir = i;
                        int jr = j;
                        int kr = k;
                        double vL;
                        double vR;

                        if (jdir == X1DIR) {
                            ir = i + 1;
                        } else if (jdir == X2DIR) {
                            jr = j + 1;
                        } else {
                            kr = k + 1;
                        }
                        vL = block->v_riemann[jdir][VRIDX(icomp, il, jl, kl)];
                        vR = block->v_riemann[jdir][VRIDX(icomp, ir, jr, kr)];
                        dvdx[jdir][icomp] = (vR - vL) * inv_dx[jdir];
                    }
                }

                for (field = 0; field < PRJ_NRAD; ++field) {
                    for (group = 0; group < PRJ_NEGROUP; ++group) {
                        double E = W[WIDX(PRJ_PRIM_RAD_E(field, group), i, j, k)];
                        double F[3];
                        double P[3][3];
                        double dE_src;
                        int jj;
                        int ii;

                        F[0] = W[WIDX(PRJ_PRIM_RAD_F1(field, group), i, j, k)];
                        F[1] = W[WIDX(PRJ_PRIM_RAD_F2(field, group), i, j, k)];
                        F[2] = W[WIDX(PRJ_PRIM_RAD_F3(field, group), i, j, k)];

                        /* Closure: P^{ij} from the cell-centred (E, F). */
                        prj_rad_m1_pressure(rad, E, F[0], F[1], F[2], P);

                        /* Energy: dE/dt += sum_{ij} (∂_j v_i) P^{ji}.
                         * P is symmetric in M1 so P^{ji} = P^{ij}. */
                        dE_src = 0.0;
                        for (jj = 0; jj < 3; ++jj) {
                            for (ii = 0; ii < 3; ++ii) {
                                dE_src += dvdx[jj][ii] * P[jj][ii];
                            }
                        }
                        /* GR lapse: same α(r) factor that multiplies the
                         * spatial radiation flux and the gravity source. */
                        dUdt[VIDX(PRJ_CONS_RAD_E(field, group), i, j, k)] -= lapse * dE_src;

                        /* Flux: dF_j/dt += sum_i (∂_j v_i) F_i. */
                        {
                            int fi[3];
                            double dFj;

                            fi[0] = PRJ_CONS_RAD_F1(field, group);
                            fi[1] = PRJ_CONS_RAD_F2(field, group);
                            fi[2] = PRJ_CONS_RAD_F3(field, group);
                            for (jj = 0; jj < 3; ++jj) {
                                dFj = 0.0;
                                for (ii = 0; ii < 3; ++ii) {
                                    dFj += dvdx[jj][ii] * F[ii];
                                }
                                dUdt[VIDX(fi[jj], i, j, k)] -= lapse * dFj;
                            }
                        }
                    }
                }
            }
        }
    }
#elif PRJ_USE_RADIATION_FSA
    int i;
    int j;
    int k;

    if (rad == 0 || block == 0 || block->id < 0 || block->active != 1 ||
        W == 0 || dUdt == 0) {
        return;
    }
    if (block->v_riemann[0] == 0 || block->v_riemann[1] == 0 || block->v_riemann[2] == 0) {
        return;
    }

    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                double dvdx[3][3]; /* dvdx[jdir][icomp] = ∂_jdir v_icomp */
                double inv_dx[3];
                double lapse = block->lapse != 0 ? block->lapse[IDX(i, j, k)] : 1.0;
                int jdir;
                int icomp;
                int field;
                int group;
                int angle;

                inv_dx[0] = 1.0 / block->dx[0];
                inv_dx[1] = 1.0 / block->dx[1];
                inv_dx[2] = 1.0 / block->dx[2];

                for (jdir = 0; jdir < 3; ++jdir) {
                    for (icomp = 0; icomp < 3; ++icomp) {
                        int il = i;
                        int jl = j;
                        int kl = k;
                        int ir = i;
                        int jr = j;
                        int kr = k;
                        double vL;
                        double vR;

                        if (jdir == X1DIR) {
                            ir = i + 1;
                        } else if (jdir == X2DIR) {
                            jr = j + 1;
                        } else {
                            kr = k + 1;
                        }
                        vL = block->v_riemann[jdir][VRIDX(icomp, il, jl, kl)];
                        vR = block->v_riemann[jdir][VRIDX(icomp, ir, jr, kr)];
                        dvdx[jdir][icomp] = (vR - vL) * inv_dx[jdir];
                    }
                }

                for (field = 0; field < PRJ_NRAD; ++field) {
                    for (group = 0; group < PRJ_NEGROUP; ++group) {
                        for (angle = 0; angle < PRJ_NANGLE; ++angle) {
                            double n[3];
                            double ndvdxn = 0.0;
                            double J;
                            int jj;
                            int ii;

                            prj_rad_fsa_rotated_angle_dir(rad, block, angle, i, j, k, n);
                            for (jj = 0; jj < 3; ++jj) {
                                for (ii = 0; ii < 3; ++ii) {
                                    ndvdxn += n[jj] * dvdx[jj][ii] * n[ii];
                                }
                            }

                            J = W[WIDX(PRJ_PRIM_RAD_I(field, group, angle), i, j, k)];
                            dUdt[VIDX(PRJ_CONS_RAD_I(field, group, angle), i, j, k)] -=
                                lapse * ndvdxn * J;
                        }
                    }
                }
            }
        }
    }
#else
    (void)rad;
    (void)block;
    (void)W;
    (void)dUdt;
#endif
}

void prj_src_update(prj_eos *eos, const prj_rad *rad, const prj_grav *grav,
    const prj_block *block,
    double *restrict W, double *restrict dUdt)
{
    int v;
    int i;
    int j;
    int k;

    PRJ_SUBTIMER_START("sub_src_zero");
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    dUdt[VIDX(v, i, j, k)] = 0.0;
                }
            }
        }
    }
    PRJ_SUBTIMER_STOP("sub_src_zero");
    PRJ_SUBTIMER_START("sub_src_geom");
    prj_src_geom(eos, W, dUdt);
    PRJ_SUBTIMER_STOP("sub_src_geom");
    PRJ_SUBTIMER_START("sub_src_user");
    prj_src_user(eos, W, dUdt);
    PRJ_SUBTIMER_STOP("sub_src_user");
    PRJ_SUBTIMER_START("sub_src_gravity");
    prj_src_monopole_gravity(rad, block, grav, W, dUdt);
    PRJ_SUBTIMER_STOP("sub_src_gravity");
    PRJ_SUBTIMER_START("sub_src_rad_vel_grad");
    prj_src_radiation_vel_grad(rad, block, W, dUdt);
    PRJ_SUBTIMER_STOP("sub_src_rad_vel_grad");
}
