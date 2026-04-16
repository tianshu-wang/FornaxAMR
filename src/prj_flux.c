#include "prj.h"

static void prj_flux_rotate_to_local(const double *Wg, int dir, double *Wl)
{
    int field;
    int group;

    Wl[PRJ_PRIM_RHO] = Wg[PRJ_PRIM_RHO];
    Wl[PRJ_PRIM_EINT] = Wg[PRJ_PRIM_EINT];
    Wl[PRJ_PRIM_YE] = Wg[PRJ_PRIM_YE];

    if (dir == X1DIR) {
        Wl[PRJ_PRIM_V1] = Wg[PRJ_PRIM_V1];
        Wl[PRJ_PRIM_V2] = Wg[PRJ_PRIM_V2];
        Wl[PRJ_PRIM_V3] = Wg[PRJ_PRIM_V3];
    } else if (dir == X2DIR) {
        Wl[PRJ_PRIM_V1] = Wg[PRJ_PRIM_V2];
        Wl[PRJ_PRIM_V2] = Wg[PRJ_PRIM_V3];
        Wl[PRJ_PRIM_V3] = Wg[PRJ_PRIM_V1];
    } else {
        Wl[PRJ_PRIM_V1] = Wg[PRJ_PRIM_V3];
        Wl[PRJ_PRIM_V2] = Wg[PRJ_PRIM_V1];
        Wl[PRJ_PRIM_V3] = Wg[PRJ_PRIM_V2];
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            Wl[PRJ_PRIM_RAD_E(field, group)] = Wg[PRJ_PRIM_RAD_E(field, group)];
            if (dir == X1DIR) {
                Wl[PRJ_PRIM_RAD_F1(field, group)] = Wg[PRJ_PRIM_RAD_F1(field, group)];
                Wl[PRJ_PRIM_RAD_F2(field, group)] = Wg[PRJ_PRIM_RAD_F2(field, group)];
                Wl[PRJ_PRIM_RAD_F3(field, group)] = Wg[PRJ_PRIM_RAD_F3(field, group)];
            } else if (dir == X2DIR) {
                Wl[PRJ_PRIM_RAD_F1(field, group)] = Wg[PRJ_PRIM_RAD_F2(field, group)];
                Wl[PRJ_PRIM_RAD_F2(field, group)] = Wg[PRJ_PRIM_RAD_F3(field, group)];
                Wl[PRJ_PRIM_RAD_F3(field, group)] = Wg[PRJ_PRIM_RAD_F1(field, group)];
            } else {
                Wl[PRJ_PRIM_RAD_F1(field, group)] = Wg[PRJ_PRIM_RAD_F3(field, group)];
                Wl[PRJ_PRIM_RAD_F2(field, group)] = Wg[PRJ_PRIM_RAD_F1(field, group)];
                Wl[PRJ_PRIM_RAD_F3(field, group)] = Wg[PRJ_PRIM_RAD_F2(field, group)];
            }
        }
    }
}

static void prj_flux_rotate_from_local(const double *Fl, int dir, double *Fg)
{
    int field;
    int group;

    Fg[PRJ_CONS_RHO] = Fl[PRJ_CONS_RHO];
    Fg[PRJ_CONS_ETOT] = Fl[PRJ_CONS_ETOT];
    Fg[PRJ_CONS_YE] = Fl[PRJ_CONS_YE];

    if (dir == X1DIR) {
        Fg[PRJ_CONS_MOM1] = Fl[PRJ_CONS_MOM1];
        Fg[PRJ_CONS_MOM2] = Fl[PRJ_CONS_MOM2];
        Fg[PRJ_CONS_MOM3] = Fl[PRJ_CONS_MOM3];
    } else if (dir == X2DIR) {
        Fg[PRJ_CONS_MOM2] = Fl[PRJ_CONS_MOM1];
        Fg[PRJ_CONS_MOM3] = Fl[PRJ_CONS_MOM2];
        Fg[PRJ_CONS_MOM1] = Fl[PRJ_CONS_MOM3];
    } else {
        Fg[PRJ_CONS_MOM3] = Fl[PRJ_CONS_MOM1];
        Fg[PRJ_CONS_MOM1] = Fl[PRJ_CONS_MOM2];
        Fg[PRJ_CONS_MOM2] = Fl[PRJ_CONS_MOM3];
    }

    for (field = 0; field < PRJ_NRAD; ++field) {
        for (group = 0; group < PRJ_NEGROUP; ++group) {
            Fg[PRJ_CONS_RAD_E(field, group)] = Fl[PRJ_CONS_RAD_E(field, group)];
            if (dir == X1DIR) {
                Fg[PRJ_CONS_RAD_F1(field, group)] = Fl[PRJ_CONS_RAD_F1(field, group)];
                Fg[PRJ_CONS_RAD_F2(field, group)] = Fl[PRJ_CONS_RAD_F2(field, group)];
                Fg[PRJ_CONS_RAD_F3(field, group)] = Fl[PRJ_CONS_RAD_F3(field, group)];
            } else if (dir == X2DIR) {
                Fg[PRJ_CONS_RAD_F2(field, group)] = Fl[PRJ_CONS_RAD_F1(field, group)];
                Fg[PRJ_CONS_RAD_F3(field, group)] = Fl[PRJ_CONS_RAD_F2(field, group)];
                Fg[PRJ_CONS_RAD_F1(field, group)] = Fl[PRJ_CONS_RAD_F3(field, group)];
            } else {
                Fg[PRJ_CONS_RAD_F3(field, group)] = Fl[PRJ_CONS_RAD_F1(field, group)];
                Fg[PRJ_CONS_RAD_F1(field, group)] = Fl[PRJ_CONS_RAD_F2(field, group)];
                Fg[PRJ_CONS_RAD_F2(field, group)] = Fl[PRJ_CONS_RAD_F3(field, group)];
            }
        }
    }
}

static void prj_flux_face_states(double *W, int dir, int i, int j, int k, double *WL, double *WR)
{
    int v;
    int il;
    int jl;
    int kl;
    int ir;
    int jr;
    int kr;

    il = i;
    jl = j;
    kl = k;
    ir = i;
    jr = j;
    kr = k;
    if (dir == X1DIR) {
        il = i - 1;
        ir = i;
    } else if (dir == X2DIR) {
        jl = j - 1;
        jr = j;
    } else {
        kl = k - 1;
        kr = k;
    }

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        double left_stencil[3];
        double right_stencil[3];
        double slope_l;
        double slope_r;
        double WLg;
        double WRg;

        if (dir == X1DIR) {
            left_stencil[0] = W[VIDX(v, il - 1, jl, kl)];
            left_stencil[1] = W[VIDX(v, il, jl, kl)];
            left_stencil[2] = W[VIDX(v, il + 1, jl, kl)];
            right_stencil[0] = W[VIDX(v, ir - 1, jr, kr)];
            right_stencil[1] = W[VIDX(v, ir, jr, kr)];
            right_stencil[2] = W[VIDX(v, ir + 1, jr, kr)];
        } else if (dir == X2DIR) {
            left_stencil[0] = W[VIDX(v, il, jl - 1, kl)];
            left_stencil[1] = W[VIDX(v, il, jl, kl)];
            left_stencil[2] = W[VIDX(v, il, jl + 1, kl)];
            right_stencil[0] = W[VIDX(v, ir, jr - 1, kr)];
            right_stencil[1] = W[VIDX(v, ir, jr, kr)];
            right_stencil[2] = W[VIDX(v, ir, jr + 1, kr)];
        } else {
            left_stencil[0] = W[VIDX(v, il, jl, kl - 1)];
            left_stencil[1] = W[VIDX(v, il, jl, kl)];
            left_stencil[2] = W[VIDX(v, il, jl, kl + 1)];
            right_stencil[0] = W[VIDX(v, ir, jr, kr - 1)];
            right_stencil[1] = W[VIDX(v, ir, jr, kr)];
            right_stencil[2] = W[VIDX(v, ir, jr, kr + 1)];
        }

        slope_l = prj_reconstruct_slope(left_stencil, 1.0);
        slope_r = prj_reconstruct_slope(right_stencil, 1.0);
        WLg = W[VIDX(v, il, jl, kl)] + 0.5 * slope_l;
        WRg = W[VIDX(v, ir, jr, kr)] - 0.5 * slope_r;
        WL[v] = WLg;
        WR[v] = WRg;
    }
}

static void prj_flux_cell_state(double *W, int i, int j, int k, double *Wc)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        Wc[v] = W[VIDX(v, i, j, k)];
    }
}

void prj_flux_update(prj_eos *eos, prj_rad *rad, prj_block *block, double *W, double *eosvar, double *flux[3])
{
    int dir;
    const prj_grav_mono *grav_mono = prj_gravity_active_monopole();

    for (dir = 0; dir < 3; ++dir) {
        int i;
        int j;
        int k;

        for (i = 0; i <= PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j <= PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k <= PRJ_BLOCK_SIZE; ++k) {
                    double WLg[PRJ_NVAR_PRIM];
                    double WRg[PRJ_NVAR_PRIM];
                    double WL[PRJ_NVAR_PRIM];
                    double WR[PRJ_NVAR_PRIM];
                    double WCLg[PRJ_NVAR_PRIM];
                    double WCRg[PRJ_NVAR_PRIM];
                    double WLLg[PRJ_NVAR_PRIM];
                    double WRRg[PRJ_NVAR_PRIM];
                    double WCL[PRJ_NVAR_PRIM];
                    double WCR[PRJ_NVAR_PRIM];
                    double WLL[PRJ_NVAR_PRIM];
                    double WRR[PRJ_NVAR_PRIM];
                    double pCLg;
                    double pCRg;
                    double pLLg;
                    double pRRg;
                    double Fl[PRJ_NVAR_CONS];
                    double Fg[PRJ_NVAR_CONS];
                    int shock_left;
                    int shock_right;
                    int v;
                    int il;
                    int jl;
                    int kl;
                    int ir;
                    int jr;
                    int kr;

                    if ((dir == X1DIR && (j >= PRJ_BLOCK_SIZE || k >= PRJ_BLOCK_SIZE)) ||
                        (dir == X2DIR && (i >= PRJ_BLOCK_SIZE || k >= PRJ_BLOCK_SIZE)) ||
                        (dir == X3DIR && (i >= PRJ_BLOCK_SIZE || j >= PRJ_BLOCK_SIZE))) {
                        continue;
                    }

                    il = i;
                    jl = j;
                    kl = k;
                    ir = i;
                    jr = j;
                    kr = k;
                    if (dir == X1DIR) {
                        il = i - 1;
                        ir = i;
                    } else if (dir == X2DIR) {
                        jl = j - 1;
                        jr = j;
                    } else {
                        kl = k - 1;
                        kr = k;
                    }

                    prj_flux_face_states(W, dir, i, j, k, WLg, WRg);
                    prj_flux_rotate_to_local(WLg, dir, WL);
                    prj_flux_rotate_to_local(WRg, dir, WR);
                    prj_flux_cell_state(W, il, jl, kl, WCLg);
                    prj_flux_cell_state(W, ir, jr, kr, WCRg);
                    if (dir == X1DIR) {
                        prj_flux_cell_state(W, il - 1, jl, kl, WLLg);
                        prj_flux_cell_state(W, ir + 1, jr, kr, WRRg);
                    } else if (dir == X2DIR) {
                        prj_flux_cell_state(W, il, jl - 1, kl, WLLg);
                        prj_flux_cell_state(W, ir, jr + 1, kr, WRRg);
                    } else {
                        prj_flux_cell_state(W, il, jl, kl - 1, WLLg);
                        prj_flux_cell_state(W, ir, jr, kr + 1, WRRg);
                    }
                    prj_flux_rotate_to_local(WCLg, dir, WCL);
                    prj_flux_rotate_to_local(WCRg, dir, WCR);
                    prj_flux_rotate_to_local(WLLg, dir, WLL);
                    prj_flux_rotate_to_local(WRRg, dir, WRR);
                    pCLg = eosvar != 0 ? eosvar[EIDX(PRJ_EOSVAR_PRESSURE, il, jl, kl)] : 0.0;
                    pCRg = eosvar != 0 ? eosvar[EIDX(PRJ_EOSVAR_PRESSURE, ir, jr, kr)] : 0.0;
                    if (dir == X1DIR) {
                        pLLg = eosvar != 0 ? eosvar[EIDX(PRJ_EOSVAR_PRESSURE, il - 1, jl, kl)] : 0.0;
                        pRRg = eosvar != 0 ? eosvar[EIDX(PRJ_EOSVAR_PRESSURE, ir + 1, jr, kr)] : 0.0;
                    } else if (dir == X2DIR) {
                        pLLg = eosvar != 0 ? eosvar[EIDX(PRJ_EOSVAR_PRESSURE, il, jl - 1, kl)] : 0.0;
                        pRRg = eosvar != 0 ? eosvar[EIDX(PRJ_EOSVAR_PRESSURE, ir, jr + 1, kr)] : 0.0;
                    } else {
                        pLLg = eosvar != 0 ? eosvar[EIDX(PRJ_EOSVAR_PRESSURE, il, jl, kl - 1)] : 0.0;
                        pRRg = eosvar != 0 ? eosvar[EIDX(PRJ_EOSVAR_PRESSURE, ir, jr, kr + 1)] : 0.0;
                    }
                    /* Detect the dominant shock direction in global coordinates.
                     * Use HLLE only on transverse fluxes; keep HLLC on the shock-aligned direction. */
                    shock_left = prj_riemann_detect_shock(WLLg, WCLg, pLLg, pCLg);
                    shock_right = prj_riemann_detect_shock(WCRg, WRRg, pCRg, pRRg);
                    double v_face_loc[3] = {0.0, 0.0, 0.0};
                    if ((shock_left >= 0 && shock_left != dir) ||
                        (shock_right >= 0 && shock_right != dir)) {
                        prj_riemann_hlle(WL, WR, eos, Fl, v_face_loc);
                    } else {
                        prj_riemann_hllc(WL, WR, eos, Fl, v_face_loc);
                    }
                    /* Unrotate (V1=normal, V2/V3=transverse) back to global (v1,v2,v3). */
                    double v_face_glob[3];
                    if (dir == X1DIR) {
                        v_face_glob[0] = v_face_loc[0];
                        v_face_glob[1] = v_face_loc[1];
                        v_face_glob[2] = v_face_loc[2];
                    } else if (dir == X2DIR) {
                        v_face_glob[1] = v_face_loc[0];
                        v_face_glob[2] = v_face_loc[1];
                        v_face_glob[0] = v_face_loc[2];
                    } else {
                        v_face_glob[2] = v_face_loc[0];
                        v_face_glob[0] = v_face_loc[1];
                        v_face_glob[1] = v_face_loc[2];
                    }
                    if (block->v_riemann[dir] != 0) {
                        block->v_riemann[dir][0 * PRJ_BLOCK_NCELLS + IDX(i, j, k)] = v_face_glob[0];
                        block->v_riemann[dir][1 * PRJ_BLOCK_NCELLS + IDX(i, j, k)] = v_face_glob[1];
                        block->v_riemann[dir][2 * PRJ_BLOCK_NCELLS + IDX(i, j, k)] = v_face_glob[2];
                    }
                    {
                        double x_face[3];
                        double dx_dir;
                        x_face[0] = block->xmin[0] + ((dir == X1DIR) ? (double)i : ((double)i + 0.5)) * block->dx[0];
                        x_face[1] = block->xmin[1] + ((dir == X2DIR) ? (double)j : ((double)j + 0.5)) * block->dx[1];
                        x_face[2] = block->xmin[2] + ((dir == X3DIR) ? (double)k : ((double)k + 0.5)) * block->dx[2];
                        dx_dir = block->dx[dir];
                        prj_rad_flux(WL, WR, eos, rad, grav_mono, x_face, dx_dir, v_face_loc[0], Fl);
                    }
                    prj_flux_rotate_from_local(Fl, dir, Fg);
                    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                        flux[dir][VIDX(v, i, j, k)] = Fg[v];
                    }
                }
            }
        }
    }
}

void prj_flux_div(double *flux[3], double area[3], double vol, double *dUdt)
{
    int v;
    int i;
    int j;
    int k;

    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    dUdt[VIDX(v, i, j, k)] =
                        -(area[X1DIR] * (flux[X1DIR][VIDX(v, i + 1, j, k)] - flux[X1DIR][VIDX(v, i, j, k)]) +
                            area[X2DIR] * (flux[X2DIR][VIDX(v, i, j + 1, k)] - flux[X2DIR][VIDX(v, i, j, k)]) +
                            area[X3DIR] * (flux[X3DIR][VIDX(v, i, j, k + 1)] - flux[X3DIR][VIDX(v, i, j, k)])) / vol;
                }
            }
        }
    }
}
