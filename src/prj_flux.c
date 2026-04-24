#include <stdio.h>

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

#if PRJ_MHD
    if (dir == X1DIR) {
        Wl[PRJ_PRIM_B1] = Wg[PRJ_PRIM_B1];
        Wl[PRJ_PRIM_B2] = Wg[PRJ_PRIM_B2];
        Wl[PRJ_PRIM_B3] = Wg[PRJ_PRIM_B3];
    } else if (dir == X2DIR) {
        Wl[PRJ_PRIM_B1] = Wg[PRJ_PRIM_B2];
        Wl[PRJ_PRIM_B2] = Wg[PRJ_PRIM_B3];
        Wl[PRJ_PRIM_B3] = Wg[PRJ_PRIM_B1];
    } else {
        Wl[PRJ_PRIM_B1] = Wg[PRJ_PRIM_B3];
        Wl[PRJ_PRIM_B2] = Wg[PRJ_PRIM_B1];
        Wl[PRJ_PRIM_B3] = Wg[PRJ_PRIM_B2];
    }
#endif

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

#if PRJ_MHD
    if (dir == X1DIR) {
        Fg[PRJ_CONS_B1] = Fl[PRJ_CONS_B1];
        Fg[PRJ_CONS_B2] = Fl[PRJ_CONS_B2];
        Fg[PRJ_CONS_B3] = Fl[PRJ_CONS_B3];
    } else if (dir == X2DIR) {
        Fg[PRJ_CONS_B2] = Fl[PRJ_CONS_B1];
        Fg[PRJ_CONS_B3] = Fl[PRJ_CONS_B2];
        Fg[PRJ_CONS_B1] = Fl[PRJ_CONS_B3];
    } else {
        Fg[PRJ_CONS_B3] = Fl[PRJ_CONS_B1];
        Fg[PRJ_CONS_B1] = Fl[PRJ_CONS_B2];
        Fg[PRJ_CONS_B2] = Fl[PRJ_CONS_B3];
    }
#endif

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

#if !PRJ_MHD
static void prj_flux_cell_state(double *W, int i, int j, int k, double *Wc)
{
    int v;

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        Wc[v] = W[VIDX(v, i, j, k)];
    }
}
#endif

static void prj_flux_face_eosvar(const double *eosvar, int var, int dir, int i, int j, int k,
    double *left_value, double *right_value)
{
    int il;
    int jl;
    int kl;
    int ir;
    int jr;
    int kr;
    double left_stencil[3];
    double right_stencil[3];
    double slope_l;
    double slope_r;

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

    if (dir == X1DIR) {
        left_stencil[0] = eosvar[EIDX(var, il - 1, jl, kl)];
        left_stencil[1] = eosvar[EIDX(var, il, jl, kl)];
        left_stencil[2] = eosvar[EIDX(var, il + 1, jl, kl)];
        right_stencil[0] = eosvar[EIDX(var, ir - 1, jr, kr)];
        right_stencil[1] = eosvar[EIDX(var, ir, jr, kr)];
        right_stencil[2] = eosvar[EIDX(var, ir + 1, jr, kr)];
    } else if (dir == X2DIR) {
        left_stencil[0] = eosvar[EIDX(var, il, jl - 1, kl)];
        left_stencil[1] = eosvar[EIDX(var, il, jl, kl)];
        left_stencil[2] = eosvar[EIDX(var, il, jl + 1, kl)];
        right_stencil[0] = eosvar[EIDX(var, ir, jr - 1, kr)];
        right_stencil[1] = eosvar[EIDX(var, ir, jr, kr)];
        right_stencil[2] = eosvar[EIDX(var, ir, jr + 1, kr)];
    } else {
        left_stencil[0] = eosvar[EIDX(var, il, jl, kl - 1)];
        left_stencil[1] = eosvar[EIDX(var, il, jl, kl)];
        left_stencil[2] = eosvar[EIDX(var, il, jl, kl + 1)];
        right_stencil[0] = eosvar[EIDX(var, ir, jr, kr - 1)];
        right_stencil[1] = eosvar[EIDX(var, ir, jr, kr)];
        right_stencil[2] = eosvar[EIDX(var, ir, jr, kr + 1)];
    }

    slope_l = prj_reconstruct_slope(left_stencil, 1.0);
    slope_r = prj_reconstruct_slope(right_stencil, 1.0);
    *left_value = eosvar[EIDX(var, il, jl, kl)] + 0.5 * slope_l;
    *right_value = eosvar[EIDX(var, ir, jr, kr)] - 0.5 * slope_r;
}

void prj_flux_update(prj_eos *eos, prj_rad *rad, prj_block *block, double *W, double *eosvar, double *flux[3])
{
    int dir;
#if PRJ_NRAD > 0
    const prj_grav_mono *grav_mono = prj_gravity_active_monopole();
#else
    (void)rad;
#endif

    for (dir = 0; dir < 3; ++dir) {
        int i_min;
        int i_max;
        int j_min;
        int j_max;
        int k_min;
        int k_max;
        int i;
        int j;
        int k;

        if (dir == X1DIR) {
            i_min = 0;
            i_max = PRJ_BLOCK_SIZE;
#if PRJ_MHD
            j_min = -1;
            j_max = PRJ_BLOCK_SIZE;
            k_min = -1;
            k_max = PRJ_BLOCK_SIZE;
#else
            j_min = 0;
            j_max = PRJ_BLOCK_SIZE - 1;
            k_min = 0;
            k_max = PRJ_BLOCK_SIZE - 1;
#endif
        } else if (dir == X2DIR) {
#if PRJ_MHD
            i_min = -1;
            i_max = PRJ_BLOCK_SIZE;
#else
            i_min = 0;
            i_max = PRJ_BLOCK_SIZE - 1;
#endif
            j_min = 0;
            j_max = PRJ_BLOCK_SIZE;
#if PRJ_MHD
            k_min = -1;
            k_max = PRJ_BLOCK_SIZE;
#else
            k_min = 0;
            k_max = PRJ_BLOCK_SIZE - 1;
#endif
        } else {
#if PRJ_MHD
            i_min = -1;
            i_max = PRJ_BLOCK_SIZE;
            j_min = -1;
            j_max = PRJ_BLOCK_SIZE;
#else
            i_min = 0;
            i_max = PRJ_BLOCK_SIZE - 1;
            j_min = 0;
            j_max = PRJ_BLOCK_SIZE - 1;
#endif
            k_min = 0;
            k_max = PRJ_BLOCK_SIZE;
        }

        for (i = i_min; i <= i_max; ++i) {
            for (j = j_min; j <= j_max; ++j) {
                for (k = k_min; k <= k_max; ++k) {
                    double WLg[PRJ_NVAR_PRIM];
                    double WRg[PRJ_NVAR_PRIM];
                    double WL[PRJ_NVAR_PRIM];
                    double WR[PRJ_NVAR_PRIM];
                    double pL;
                    double pR;
                    double gL;
                    double gR;
                    double Fl[PRJ_NVAR_CONS];
                    double Fg[PRJ_NVAR_CONS];
                    double v_face_loc[3] = {0.0, 0.0, 0.0};
                    int v;
#if !PRJ_MHD || PRJ_NRAD > 0
                    int il;
                    int jl;
                    int kl;
                    int ir;
                    int jr;
                    int kr;
#endif
#if !PRJ_MHD
                    double WCLg[PRJ_NVAR_PRIM];
                    double WCRg[PRJ_NVAR_PRIM];
                    double WLLg[PRJ_NVAR_PRIM];
                    double WRRg[PRJ_NVAR_PRIM];
                    double pCLg;
                    double pCRg;
                    double pLLg;
                    double pRRg;
                    int shock_left;
                    int shock_right;
#endif

#if !PRJ_MHD || PRJ_NRAD > 0
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
#endif

                    prj_flux_face_states(W, dir, i, j, k, WLg, WRg);
                    prj_flux_rotate_to_local(WLg, dir, WL);
                    prj_flux_rotate_to_local(WRg, dir, WR);
                    if (eosvar != 0) {
                        prj_flux_face_eosvar(eosvar, PRJ_EOSVAR_PRESSURE, dir, i, j, k, &pL, &pR);
                        prj_flux_face_eosvar(eosvar, PRJ_EOSVAR_GAMMA, dir, i, j, k, &gL, &gR);
                    } else {
                        pL = 0.0;
                        pR = 0.0;
                        gL = 0.0;
                        gR = 0.0;
                    }
#if PRJ_MHD
                    {
                        double b_face = block->Bf[dir] != 0 ? block->Bf[dir][IDX(i, j, k)] : 0.0;
                        double emf_face_loc[2] = {0.0, 0.0};

                        WL[PRJ_PRIM_B1] = b_face;
                        WR[PRJ_PRIM_B1] = b_face;
                        prj_riemann_hlld(WL, WR, pL, pR, gL, gR, eos, Fl, v_face_loc, emf_face_loc);
#if PRJ_MHD_DEBUG
                        if (dir == X1DIR &&
                            ((block->id == 1 && i == 15 && j == PRJ_BLOCK_SIZE) ||
                                (block->id == 5 && i == 15 && j == 0))) {
                            int wr_i = i;
                            int wr_j = j;
                            int wr_k = k;
                            double wr_b2_stencil[3] = {
                                W[VIDX(PRJ_PRIM_B2, wr_i - 1, wr_j, wr_k)],
                                W[VIDX(PRJ_PRIM_B2, wr_i, wr_j, wr_k)],
                                W[VIDX(PRJ_PRIM_B2, wr_i + 1, wr_j, wr_k)]
                            };
                            double wr_cell16_b2_faces[2] = {
                                block->Bf[X2DIR][IDX(wr_i + 1, wr_j, wr_k)],
                                block->Bf[X2DIR][IDX(wr_i + 1, wr_j + 1, wr_k)]
                            };
                            double wr_cell16_b3_faces[2] = {
                                block->Bf[X3DIR][IDX(wr_i + 1, wr_j, wr_k)],
                                block->Bf[X3DIR][IDX(wr_i + 1, wr_j, wr_k + 1)]
                            };

                            fprintf(stderr,
                                "[mhd-face-debug] block=%d dir=%d face=(%d,%d,%d) "
                                "b=% .17g p=(% .17g,% .17g) g=(% .17g,% .17g) "
                                "WL=(rho=% .17g eint=% .17g v=% .17g,% .17g,% .17g B=% .17g,% .17g,% .17g) "
                                "WR=(rho=% .17g eint=% .17g v=% .17g,% .17g,% .17g B=% .17g,% .17g,% .17g) "
                                "WR_B2_stencil=(% .17g,% .17g,% .17g) "
                                "WR_cell16_B2_faces=(% .17g,% .17g) "
                                "WR_cell16_B3_faces=(% .17g,% .17g) "
                                "vface=(% .17g,% .17g,% .17g) emf=(% .17g,% .17g)\n",
                                block->id, dir, i, j, k,
                                b_face, pL, pR, gL, gR,
                                WL[PRJ_PRIM_RHO], WL[PRJ_PRIM_EINT],
                                WL[PRJ_PRIM_V1], WL[PRJ_PRIM_V2], WL[PRJ_PRIM_V3],
                                WL[PRJ_PRIM_B1], WL[PRJ_PRIM_B2], WL[PRJ_PRIM_B3],
                                WR[PRJ_PRIM_RHO], WR[PRJ_PRIM_EINT],
                                WR[PRJ_PRIM_V1], WR[PRJ_PRIM_V2], WR[PRJ_PRIM_V3],
                                WR[PRJ_PRIM_B1], WR[PRJ_PRIM_B2], WR[PRJ_PRIM_B3],
                                wr_b2_stencil[0], wr_b2_stencil[1], wr_b2_stencil[2],
                                wr_cell16_b2_faces[0], wr_cell16_b2_faces[1],
                                wr_cell16_b3_faces[0], wr_cell16_b3_faces[1],
                                v_face_loc[0], v_face_loc[1], v_face_loc[2],
                                emf_face_loc[0], emf_face_loc[1]);
                        }
#endif
                        if (block->vB1[dir] != 0) {
                            block->vB1[dir][IDX(i, j, k)] = emf_face_loc[0];
                        }
                        if (block->vB2[dir] != 0) {
                            block->vB2[dir][IDX(i, j, k)] = emf_face_loc[1];
                        }
                    }
#else
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
                    shock_left = prj_riemann_detect_shock(WLLg, WCLg, pLLg, pCLg);
                    shock_right = prj_riemann_detect_shock(WCRg, WRRg, pCRg, pRRg);
                    if ((shock_left >= 0 && shock_left != dir) ||
                        (shock_right >= 0 && shock_right != dir)) {
                        prj_riemann_hlle(WL, WR, pL, pR, gL, gR, eos, Fl, v_face_loc);
                    } else {
                        prj_riemann_hllc(WL, WR, pL, pR, gL, gR, eos, Fl, v_face_loc);
                    }
#endif
                    {
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
                    }
#if PRJ_NRAD > 0
                    {
                        double x_face[3];
                        double dx_dir;
                        double chi_face[PRJ_NRAD * PRJ_NEGROUP];
                        double kappa_L[PRJ_NRAD * PRJ_NEGROUP];
                        double sigma_L[PRJ_NRAD * PRJ_NEGROUP];
                        double kappa_R[PRJ_NRAD * PRJ_NEGROUP];
                        double sigma_R[PRJ_NRAD * PRJ_NEGROUP];
                        double rho_L = W[VIDX(PRJ_PRIM_RHO, il, jl, kl)];
                        double ye_L = W[VIDX(PRJ_PRIM_YE, il, jl, kl)];
                        double rho_R = W[VIDX(PRJ_PRIM_RHO, ir, jr, kr)];
                        double ye_R = W[VIDX(PRJ_PRIM_YE, ir, jr, kr)];
                        double T_L = eosvar != 0 ? eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, il, jl, kl)] : 0.0;
                        double T_R = eosvar != 0 ? eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, ir, jr, kr)] : 0.0;
                        int idx;

                        x_face[0] = block->xmin[0] + ((dir == X1DIR) ? (double)i : ((double)i + 0.5)) * block->dx[0];
                        x_face[1] = block->xmin[1] + ((dir == X2DIR) ? (double)j : ((double)j + 0.5)) * block->dx[1];
                        x_face[2] = block->xmin[2] + ((dir == X3DIR) ? (double)k : ((double)k + 0.5)) * block->dx[2];
                        dx_dir = block->dx[dir];

                        prj_rad3_opac_lookup(rad, rho_L, T_L, ye_L, kappa_L, sigma_L, 0, 0);
                        prj_rad3_opac_lookup(rad, rho_R, T_R, ye_R, kappa_R, sigma_R, 0, 0);
                        for (idx = 0; idx < PRJ_NRAD * PRJ_NEGROUP; ++idx) {
                            double kL = kappa_L[idx];
                            double kR = kappa_R[idx];
                            double sL_o = sigma_L[idx];
                            double sR_o = sigma_R[idx];
                            double k_sum = kL + kR;
                            double s_sum = sL_o + sR_o;
                            double k_face = (k_sum > 0.0) ? (2.0 * kL * kR / k_sum) : 0.0;
                            double s_face = (s_sum > 0.0) ? (2.0 * sL_o * sR_o / s_sum) : 0.0;
                            chi_face[idx] = k_face + s_face;
                        }

                        prj_rad_flux(WL, WR, grav_mono, x_face, chi_face, dx_dir, v_face_loc[0], Fl);
                    }
#endif
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
