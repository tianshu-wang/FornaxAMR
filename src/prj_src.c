#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "prj.h"

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

void prj_src_geom(prj_eos *eos, double *W_mhd, double *W_rad,
    double *mhd_rhs, double *rad_rhs)
{
    (void)eos;
    (void)W_mhd;
    (void)W_rad;
    (void)mhd_rhs;
    (void)rad_rhs;
}

void prj_src_user(prj_eos *eos, double *W_mhd, double *W_rad,
    double *mhd_rhs, double *rad_rhs)
{
    (void)eos;
    (void)W_mhd;
    (void)W_rad;
    (void)mhd_rhs;
    (void)rad_rhs;
}

void prj_src_monopole_gravity(const prj_rad *rad, const prj_block *block,
    const prj_grav *grav, double *restrict W_mhd, double *restrict W_rad,
    double *restrict mhd_rhs, double *restrict rad_rhs)
{
    int i;
    int j;
    int k;

    if (block == 0 || block->id < 0 || block->active != 1 ||
        grav == 0 || W_mhd == 0 || mhd_rhs == 0) {
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
#if PRJ_NRAD == 0
    (void)W_rad;
    (void)rad_rhs;
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
                rho = W_mhd[WIDX(PRJ_PRIM_RHO, i, j, k)];
                mhd_rhs[MHDVIDX(PRJ_CONS_MOM1, i, j, k)] += rho * g1;
                mhd_rhs[MHDVIDX(PRJ_CONS_MOM2, i, j, k)] += rho * g2;
                mhd_rhs[MHDVIDX(PRJ_CONS_MOM3, i, j, k)] += rho * g3;
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

                    mhd_rhs[MHDVIDX(PRJ_CONS_ETOT, i, j, k)] +=
                        rho * (v_avg1 * g1 + v_avg2 * g2 + v_avg3 * g3);
                }
#if PRJ_USE_RADIATION_M1
                if (W_rad != 0 && rad_rhs != 0) {
                    double lapse = block->lapse[cell_idx];
                    int field;
                    int group;

                    double inv_c2 = 1.0 / (PRJ_CLIGHT * PRJ_CLIGHT);
                    for (field = 0; field < PRJ_NRAD; ++field) {
                        for (group = 0; group < PRJ_NEGROUP; ++group) {
                            double E = W_rad[WIDX(PRJ_RAD_PRIM_E(field, group), i, j, k)];
                            double F1 = W_rad[WIDX(PRJ_RAD_PRIM_F1(field, group), i, j, k)];
                            double F2 = W_rad[WIDX(PRJ_RAD_PRIM_F2(field, group), i, j, k)];
                            double F3 = W_rad[WIDX(PRJ_RAD_PRIM_F3(field, group), i, j, k)];

                            rad_rhs[RADVIDX(PRJ_RAD_CONS_F1(field, group), i, j, k)] += lapse * E * g1;
                            rad_rhs[RADVIDX(PRJ_RAD_CONS_F2(field, group), i, j, k)] += lapse * E * g2;
                            rad_rhs[RADVIDX(PRJ_RAD_CONS_F3(field, group), i, j, k)] += lapse * E * g3;
                            rad_rhs[RADVIDX(PRJ_RAD_CONS_E(field, group), i, j, k)] +=
                                lapse * (g1 * F1 + g2 * F2 + g3 * F3) * inv_c2;
                        }
                    }
                }
#elif PRJ_USE_RADIATION_FSA
                if (rad != 0 && W_rad != 0 && rad_rhs != 0) {
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
                                double J = W_rad[WIDX(PRJ_RAD_PRIM_I(field, group, angle), i, j, k)];

                                prj_rad_fsa_rotated_angle_dir(rad, block, angle, i, j, k, n);
                                a_dot_n = g1 * n[0] + g2 * n[1] + g3 * n[2];
                                rad_rhs[RADVIDX(PRJ_RAD_CONS_I(field, group, angle), i, j, k)] -=
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
    double *restrict W_rad, double *restrict rad_rhs)
{
#if PRJ_USE_RADIATION_M1
    int i;
    int j;
    int k;

    if (block == 0 || block->id < 0 || block->active != 1 ||
        W_rad == 0 || rad_rhs == 0) {
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
                        double E = W_rad[WIDX(PRJ_RAD_PRIM_E(field, group), i, j, k)];
                        double F[3];
                        double P[3][3];
                        double dE_src;
                        int jj;
                        int ii;

                        F[0] = W_rad[WIDX(PRJ_RAD_PRIM_F1(field, group), i, j, k)];
                        F[1] = W_rad[WIDX(PRJ_RAD_PRIM_F2(field, group), i, j, k)];
                        F[2] = W_rad[WIDX(PRJ_RAD_PRIM_F3(field, group), i, j, k)];

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
                        rad_rhs[RADVIDX(PRJ_RAD_CONS_E(field, group), i, j, k)] -= lapse * dE_src;

                        /* Flux: dF_j/dt += sum_i (∂_j v_i) F_i. */
                        {
                            int fi[3];
                            double dFj;

                            fi[0] = PRJ_RAD_CONS_F1(field, group);
                            fi[1] = PRJ_RAD_CONS_F2(field, group);
                            fi[2] = PRJ_RAD_CONS_F3(field, group);
                            for (jj = 0; jj < 3; ++jj) {
                                dFj = 0.0;
                                for (ii = 0; ii < 3; ++ii) {
                                    dFj += dvdx[jj][ii] * F[ii];
                                }
                                rad_rhs[RADVIDX(fi[jj], i, j, k)] -= lapse * dFj;
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
        W_rad == 0 || rad_rhs == 0) {
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

                            J = W_rad[WIDX(PRJ_RAD_PRIM_I(field, group, angle), i, j, k)];
                            rad_rhs[RADVIDX(PRJ_RAD_CONS_I(field, group, angle), i, j, k)] -=
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
    (void)W_rad;
    (void)rad_rhs;
#endif
}

#if PRJ_DYNAMIC_GR && !PRJ_MHD
static void prj_src_gr_fail(const char *op, int status, int i, int j, int k)
{
    fprintf(stderr, "dynamic GR source %s failed at cell (%d,%d,%d): status=%d\n",
        op, i, j, k, status);
#if defined(PRJ_ENABLE_MPI)
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
    exit(EXIT_FAILURE);
}

static void prj_src_gr_hydro_z4c(prj_eos *eos, const prj_mesh *mesh,
    const prj_block *block, int z4c_stage, double *restrict W_mhd,
    double *restrict mhd_rhs)
{
    int i;
    int j;
    int k;

    if (!prj_eos_full_dynamic_gr_enabled(mesh) || block == 0 || block->id < 0 ||
        block->active != 1 || W_mhd == 0 || mhd_rhs == 0) {
        return;
    }
    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                prj_z4c_hydro_geom geom;
                prj_eos_gr_geom egeom;
                double Wc[PRJ_NVAR_PRIM];
                double Uc[PRJ_NVAR_CONS];
                double Uloc[PRJ_NVAR_CONS];
                double beta[3];
                double beta_cov[3];
                double Scon[3];
                double Tij[3][3];
                double pressure;
                double rhoh;
                double beta2 = 0.0;
                double wlor2;
                double E;
                double energy_src = 0.0;
                int status;
                int v;
                int a;
                int b;
                int d;

                if (!prj_z4c_load_hydro_geom(mesh, block, z4c_stage, i, j, k, &geom)) {
                    prj_src_gr_fail("geometry load", -1, i, j, k);
                }
                for (a = 0; a < 3; ++a) {
                    for (b = 0; b < 3; ++b) {
                        egeom.gamma[a][b] = geom.gamma[a][b];
                    }
                }
                for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                    Wc[v] = 0.0;
                }
                for (v = 0; v < PRJ_NHYDRO; ++v) {
                    Wc[v] = W_mhd[WIDX(v, i, j, k)];
                }
                pressure = block->eosvar != 0 ?
                    block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)] : 0.0;
                status = prj_eos_gr_prim2cons(eos, &egeom, Wc, Uc, PRJ_EOS_CTX_MAIN);
                if (status != PRJ_EOS_GR_OK) {
                    prj_src_gr_fail("prim2cons", status, i, j, k);
                }
                for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                    Uloc[v] = Uc[v] / geom.sqrt_gamma;
                }
                for (a = 0; a < 3; ++a) {
                    beta[a] = Wc[PRJ_PRIM_V1 + a] / PRJ_CLIGHT;
                    beta_cov[a] = 0.0;
                    Scon[a] = 0.0;
                }
                for (a = 0; a < 3; ++a) {
                    for (b = 0; b < 3; ++b) {
                        beta_cov[a] += geom.gamma[a][b] * beta[b];
                    }
                    beta2 += beta_cov[a] * beta[a];
                }
                if (!isfinite(beta2) || beta2 < 0.0 || beta2 >= 1.0) {
                    prj_src_gr_fail("state", PRJ_EOS_GR_BAD_STATE, i, j, k);
                }
                wlor2 = 1.0 / (1.0 - beta2);
                rhoh = Wc[PRJ_PRIM_RHO] * PRJ_CLIGHT * PRJ_CLIGHT +
                    Wc[PRJ_PRIM_RHO] * Wc[PRJ_PRIM_EINT] + pressure;
                E = rhoh * wlor2 - pressure;
                for (a = 0; a < 3; ++a) {
                    for (b = 0; b < 3; ++b) {
                        Tij[a][b] = rhoh * wlor2 * beta[a] * beta[b] +
                            pressure * geom.gamma_inv[a][b];
                        Scon[a] += geom.gamma_inv[a][b] * Uloc[PRJ_CONS_MOM1 + b];
                    }
                }
                for (d = 0; d < 3; ++d) {
                    double src = -E * geom.dalpha[d];

                    for (a = 0; a < 3; ++a) {
                        src += PRJ_CLIGHT * Uloc[PRJ_CONS_MOM1 + a] *
                            geom.dbeta[d][a];
                        for (b = 0; b < 3; ++b) {
                            src += 0.5 * geom.alpha * Tij[a][b] *
                                geom.dgamma[d][a][b];
                        }
                    }
                    mhd_rhs[MHDVIDX(PRJ_CONS_MOM1 + d, i, j, k)] +=
                        geom.sqrt_gamma * src;
                    energy_src -= geom.sqrt_gamma * PRJ_CLIGHT * PRJ_CLIGHT *
                        Scon[d] * geom.dalpha[d];
                }
                for (a = 0; a < 3; ++a) {
                    for (b = 0; b < 3; ++b) {
                        energy_src += geom.sqrt_gamma * PRJ_CLIGHT *
                            geom.alpha * Tij[a][b] * geom.K_dd[a][b];
                    }
                }
                mhd_rhs[MHDVIDX(PRJ_CONS_ETOT, i, j, k)] += energy_src;
            }
        }
    }
}
#endif

void prj_src_update(prj_eos *eos, const prj_rad *rad, const prj_grav *grav,
    const prj_mesh *mesh, const prj_block *block, int z4c_stage,
    double *restrict W_mhd, double *restrict W_rad,
    double *restrict mhd_rhs, double *restrict rad_rhs)
{
    int v;
    int i;
    int j;
    int k;

    PRJ_SUBTIMER_START("sub_src_zero");
    for (v = 0; v < PRJ_NVAR_MHD_CONS; ++v) {
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    mhd_rhs[MHDVIDX(v, i, j, k)] = 0.0;
                }
            }
        }
    }
    for (v = 0; v < PRJ_NVAR_RAD_CONS; ++v) {
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    rad_rhs[RADVIDX(v, i, j, k)] = 0.0;
                }
            }
        }
    }
    PRJ_SUBTIMER_STOP("sub_src_zero");
    PRJ_SUBTIMER_START("sub_src_geom");
    prj_src_geom(eos, W_mhd, W_rad, mhd_rhs, rad_rhs);
    PRJ_SUBTIMER_STOP("sub_src_geom");
    PRJ_SUBTIMER_START("sub_src_user");
    prj_src_user(eos, W_mhd, W_rad, mhd_rhs, rad_rhs);
    PRJ_SUBTIMER_STOP("sub_src_user");
    PRJ_SUBTIMER_START("sub_src_gravity");
#if PRJ_DYNAMIC_GR && !PRJ_MHD
    if (prj_eos_full_dynamic_gr_enabled(mesh)) {
        prj_src_gr_hydro_z4c(eos, mesh, block, z4c_stage, W_mhd, mhd_rhs);
    } else
#else
    (void)mesh;
    (void)z4c_stage;
#endif
    prj_src_monopole_gravity(rad, block, grav, W_mhd, W_rad, mhd_rhs, rad_rhs);
    PRJ_SUBTIMER_STOP("sub_src_gravity");
    PRJ_SUBTIMER_START("sub_src_rad_vel_grad");
    prj_src_radiation_vel_grad(rad, block, W_rad, rad_rhs);
    PRJ_SUBTIMER_STOP("sub_src_rad_vel_grad");
}
