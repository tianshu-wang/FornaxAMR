#include <math.h>

#include "prj.h"

#if PRJ_NRAD > 0
static double prj_src_interp_lapse(const prj_grav_mono *grav_mono, double r)
{
    int idx;

    if (grav_mono == 0 || grav_mono->nbins <= 0 || grav_mono->lapse == 0) {
        return 1.0;
    }
    if (r <= grav_mono->rf[0]) {
        return grav_mono->lapse[0];
    }
    for (idx = 0; idx < grav_mono->nbins; ++idx) {
        double r0 = grav_mono->rf[idx];
        double r1 = grav_mono->rf[idx + 1];

        if (r <= r1) {
            double weight = (r - r0) / (r1 - r0);

            return (1.0 - weight) * grav_mono->lapse[idx] + weight * grav_mono->lapse[idx + 1];
        }
    }
    return grav_mono->lapse[grav_mono->nbins];
}
#endif

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

void prj_src_monopole_gravity(const prj_block *block, const prj_grav_mono *grav_mono,
    double *restrict W, double *restrict dUdt)
{
    int i;
    int j;
    int k;

    if (block == 0 || block->id < 0 || block->active != 1 ||
        grav_mono == 0 || W == 0 || dUdt == 0) {
        return;
    }
    if (block->v_riemann[0] == 0 || block->v_riemann[1] == 0 || block->v_riemann[2] == 0) {
        return;
    }

    for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
        for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
            for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                double x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                double x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                double x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                double r = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
                double accel;
                double g1;
                double g2;
                double g3;
                double rho;
                double v1;
                double v2;
                double v3;

                if (r <= 0.0) {
                    continue;
                }

                accel = prj_gravity_block_accel_at(block, i, j, k);
                g1 = accel * x1 / r;
                g2 = accel * x2 / r;
                g3 = accel * x3 / r;
                rho = W[VIDX(PRJ_PRIM_RHO, i, j, k)];
                dUdt[VIDX(PRJ_CONS_MOM1, i, j, k)] += rho * g1;
                dUdt[VIDX(PRJ_CONS_MOM2, i, j, k)] += rho * g2;
                dUdt[VIDX(PRJ_CONS_MOM3, i, j, k)] += rho * g3;
                {
                    double v_avg[3];
                    int idx000 = IDX(i, j, k);
                    int idx100 = IDX(i + 1, j, k);
                    int idx010 = IDX(i, j + 1, k);
                    int idx001 = IDX(i, j, k + 1);
                    int n;

                    for (n = 0; n < 3; ++n) {
                        v_avg[n] = (block->v_riemann[X1DIR][n * PRJ_BLOCK_NCELLS + idx000] +
                                  block->v_riemann[X1DIR][n * PRJ_BLOCK_NCELLS + idx100] +
                                  block->v_riemann[X2DIR][n * PRJ_BLOCK_NCELLS + idx000] +
                                  block->v_riemann[X2DIR][n * PRJ_BLOCK_NCELLS + idx010] +
                                  block->v_riemann[X3DIR][n * PRJ_BLOCK_NCELLS + idx000] +
                                  block->v_riemann[X3DIR][n * PRJ_BLOCK_NCELLS + idx001]) / 6.0;
                    }

                    dUdt[VIDX(PRJ_CONS_ETOT, i, j, k)] +=
                        rho * (v_avg[0] * g1 + v_avg[1] * g2 + v_avg[2] * g3);
                }
#if PRJ_NRAD > 0
                {
                    double lapse = prj_src_interp_lapse(grav_mono, r);
                    int field;
                    int group;

                    double inv_c2 = 1.0 / (PRJ_CLIGHT * PRJ_CLIGHT);
                    for (field = 0; field < PRJ_NRAD; ++field) {
                        for (group = 0; group < PRJ_NEGROUP; ++group) {
                            double E = W[VIDX(PRJ_PRIM_RAD_E(field, group), i, j, k)];
                            double F1 = W[VIDX(PRJ_PRIM_RAD_F1(field, group), i, j, k)];
                            double F2 = W[VIDX(PRJ_PRIM_RAD_F2(field, group), i, j, k)];
                            double F3 = W[VIDX(PRJ_PRIM_RAD_F3(field, group), i, j, k)];

                            dUdt[VIDX(PRJ_CONS_RAD_F1(field, group), i, j, k)] += lapse * E * g1;
                            dUdt[VIDX(PRJ_CONS_RAD_F2(field, group), i, j, k)] += lapse * E * g2;
                            dUdt[VIDX(PRJ_CONS_RAD_F3(field, group), i, j, k)] += lapse * E * g3;
                            dUdt[VIDX(PRJ_CONS_RAD_E(field, group), i, j, k)] +=
                                lapse * (g1 * F1 + g2 * F2 + g3 * F3) * inv_c2;
                        }
                    }
                }
#endif
            }
        }
    }
}

#if PRJ_NRAD > 0
/* Build the M1 radiation pressure tensor P^{ij} = E * D^{ij} for a single
 * (species, group), where D^{ij} = a δ^{ij} + b n^i n^j with the Levermore
 * closure χ(f) and n = F/|F|, f = |F|/(c E).  Falls back to the isotropic
 * limit P^{ij} = (E/3) δ^{ij} when |F| or E vanishes. */
static void prj_src_rad_m1_pressure(double E, double F1, double F2, double F3,
    double P[3][3])
{
    double F[3];
    double Fmag;
    double cE;
    double f;
    double chi;
    double a_c;
    double b_c;
    double n[3];
    int a;
    int b;

    F[0] = F1;
    F[1] = F2;
    F[2] = F3;
    Fmag = sqrt(F1 * F1 + F2 * F2 + F3 * F3);
    cE = PRJ_CLIGHT * (E > 0.0 ? E : 0.0);

    if (cE <= 0.0 || Fmag <= 0.0) {
        double third = (E > 0.0 ? E : 0.0) / 3.0;

        for (a = 0; a < 3; ++a) {
            for (b = 0; b < 3; ++b) {
                P[a][b] = (a == b) ? third : 0.0;
            }
        }
        return;
    }

    f = Fmag / cE;
    if (f > 1.0) {
        f = 1.0;
    }
    chi = (3.0 + 4.0 * f * f) / (5.0 + 2.0 * sqrt(4.0 - 3.0 * f * f));
    a_c = 0.5 * (1.0 - chi);
    b_c = 0.5 * (3.0 * chi - 1.0);
    n[0] = F[0] / Fmag;
    n[1] = F[1] / Fmag;
    n[2] = F[2] / Fmag;
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            P[a][b] = E * (a_c * (a == b ? 1.0 : 0.0) + b_c * n[a] * n[b]);
        }
    }
}
#endif

/* O(v/c) "source-like" piece of the SR redshift (Eq. 12a/12b of the
 * comoving-frame mixed-frame moment equations).  These are the bits left over
 * once the energy-space-flux divergence has been split off, and they appear as
 * pure cell-local sources:
 *
 *     ∂_t E_g   += (∂_j v_i) P^{ji}_g                          (Eq. 12a piece)
 *     ∂_t F_{gj} += (∂_j v_i) F_{gi}                            (Eq. 12b piece)
 *
 * Indices: i,j ∈ {1,2,3}; sums over i (and i,j for the energy term).  ∂_j v_i is
 * built by central differencing the face-centred Riemann velocities stored in
 * block->v_riemann[face_dir][component, ...] during the most recent flux update,
 * so the velocity field used here is the same one the hydro fluxes saw. */
void prj_src_radiation_vel_grad(const prj_block *block,
    double *restrict W, double *restrict dUdt)
{
#if PRJ_NRAD > 0
    const prj_grav_mono *grav_mono = prj_gravity_active_monopole();
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
                /* Cell-centred lapse α(r) for the GR proper-time slowdown. */
                double xc1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
                double xc2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
                double xc3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
                double r_cell = sqrt(xc1 * xc1 + xc2 * xc2 + xc3 * xc3);
                double lapse = (grav_mono != 0)
                    ? prj_src_interp_lapse(grav_mono, r_cell)
                    : 1.0;

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
                        vL = block->v_riemann[jdir][icomp * PRJ_BLOCK_NCELLS + IDX(il, jl, kl)];
                        vR = block->v_riemann[jdir][icomp * PRJ_BLOCK_NCELLS + IDX(ir, jr, kr)];
                        dvdx[jdir][icomp] = (vR - vL) * inv_dx[jdir];
                    }
                }

                for (field = 0; field < PRJ_NRAD; ++field) {
                    for (group = 0; group < PRJ_NEGROUP; ++group) {
                        double E = W[VIDX(PRJ_PRIM_RAD_E(field, group), i, j, k)];
                        double F[3];
                        double P[3][3];
                        double dE_src;
                        int jj;
                        int ii;

                        F[0] = W[VIDX(PRJ_PRIM_RAD_F1(field, group), i, j, k)];
                        F[1] = W[VIDX(PRJ_PRIM_RAD_F2(field, group), i, j, k)];
                        F[2] = W[VIDX(PRJ_PRIM_RAD_F3(field, group), i, j, k)];

                        /* Closure: P^{ij} from the cell-centred (E, F). */
                        prj_src_rad_m1_pressure(E, F[0], F[1], F[2], P);

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
                        dUdt[VIDX(PRJ_CONS_RAD_E(field, group), i, j, k)] += lapse * dE_src;

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
                                dUdt[VIDX(fi[jj], i, j, k)] += lapse * dFj;
                            }
                        }
                    }
                }
            }
        }
    }
#else
    (void)block;
    (void)W;
    (void)dUdt;
#endif
}

void prj_src_update(prj_eos *eos, const prj_block *block, double *restrict W,
    double *restrict dUdt)
{
    prj_src_geom(eos, W, dUdt);
    prj_src_user(eos, W, dUdt);
    prj_src_monopole_gravity(block, prj_gravity_active_monopole(), W, dUdt);
    prj_src_radiation_vel_grad(block, W, dUdt);
}
