#ifndef PRJ_DEFS_H
#define PRJ_DEFS_H

#include <stddef.h>

/* Some strict C99 compilers do not expose M_PI from <math.h>. Define it once
   here so every translation unit that includes prj_defs.h (via prj.h) sees it. */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define PRJ_NDIM 3
#ifndef PRJ_USE_GRAVITY
#define PRJ_USE_GRAVITY 1
#endif

#ifndef PRJ_GRAV_DEBUG
#define PRJ_GRAV_DEBUG 0
#endif

#ifndef PRJ_TIMER
#define PRJ_TIMER 0
#endif

#ifndef PRJ_USE_RADIATION_M1
#define PRJ_USE_RADIATION_M1 0
#endif
#ifndef PRJ_USE_RADIATION_FSA
#define PRJ_USE_RADIATION_FSA 0
#endif
#ifndef PRJ_USE_RADIAL_FRAME_FSA
#ifdef USE_RADIAL_FRAME_FSA
#define PRJ_USE_RADIAL_FRAME_FSA USE_RADIAL_FRAME_FSA
#else
#define PRJ_USE_RADIAL_FRAME_FSA 0
#endif
#endif
#if PRJ_USE_RADIATION_M1 && PRJ_USE_RADIATION_FSA
#error "RADIATION_M1 and RADIATION_FSA are separate solvers; enable only one"
#endif

#ifndef PRJ_DUMP_SINGLE_PRECISION
#define PRJ_DUMP_SINGLE_PRECISION 1
#endif

#ifndef PRJ_MIXED_PRECISION_FLUX
#define PRJ_MIXED_PRECISION_FLUX 0
#endif

#ifndef PRJ_MIXED_PRECISION_TABLE
#define PRJ_MIXED_PRECISION_TABLE 1
#endif

#if PRJ_MIXED_PRECISION_TABLE
typedef float prj_table_real;
#else
typedef double prj_table_real;
#endif

#ifndef PRJ_MHD
#define PRJ_MHD 0
#endif

#ifndef PRJ_MHD_DEBUG
#define PRJ_MHD_DEBUG 0
#endif

#ifndef PRJ_DYNAMIC_GR
#define PRJ_DYNAMIC_GR 0
#endif
#ifndef PRJ_RAD_GR_M1_CSIGMA
#define PRJ_RAD_GR_M1_CSIGMA 0.5
#endif
#ifndef PRJ_RAD_GR_M1_EULERIAN_FBAR_BETA
#define PRJ_RAD_GR_M1_EULERIAN_FBAR_BETA 1.0e-4
#endif

#define PRJ_TIMEINT_RK2 1
#define PRJ_TIMEINT_ESSPRK 2
#define PRJ_TIMEINT_ESSPRK9_3 3
#define PRJ_TIMEINT_IMEX 4
#ifndef RK2
#define RK2 PRJ_TIMEINT_RK2
#endif
#ifndef eSSPRK9_2
#define eSSPRK9_2 PRJ_TIMEINT_ESSPRK
#endif
#ifndef eSSPRK9_3
#define eSSPRK9_3 PRJ_TIMEINT_ESSPRK9_3
#endif
#ifndef IMEX
#define IMEX PRJ_TIMEINT_IMEX
#endif
/* === Tme integration scheme (single source of truth) ===================
 * These defaults are authoritative; the build only overrides them when the
 * corresponding -D is passed on the command line (see the Makefile's
 * TIME_INTEGRATION_DEF), so editing here is enough -- no Makefile edit needed.
 *
 * Select the scheme with TIME_INTEGRATION: IMEX, RK2, eSSPRK9_2, eSSPRK9_3.
 *   - eSSPRK9_2 with a different stage count: set PRJ_TIMEINT_ESSPRK_N.
 *   - IMEX with a different Butcher tableau: set PRJ_TIMEINT_TABLEAU_NAME and
 *     PRJ_TIMEINT_TABLEAU_NSTAGES.  These two MUST agree: NSTAGES must equal the
 *     stage count of the named tableau struct in rk_tableau/, otherwise the
 *     integrator indexes a wrong-sized coefficient array (silent, no error).
 *
 * Command-line overrides win, e.g.:
 *   make TIME_INTEGRATION=RK2
 *   make TIME_INTEGRATION=IMEX PRJ_TIMEINT_TABLEAU_NAME=prj_godunov_rk2 \
 *        PRJ_TIMEINT_TABLEAU_NSTAGES=3
 */
#ifndef TIME_INTEGRATION
#define TIME_INTEGRATION RK2
#endif
#ifndef PRJ_TIMEINT_ESSPRK_N
#define PRJ_TIMEINT_ESSPRK_N 9
#endif
#ifndef PRJ_TIMEINT_TABLEAU_NAME
#define PRJ_TIMEINT_TABLEAU_NAME prj_imex_ssp2332
#endif
#ifndef PRJ_TIMEINT_TABLEAU_NSTAGES
#define PRJ_TIMEINT_TABLEAU_NSTAGES 2
#endif
#if TIME_INTEGRATION == PRJ_TIMEINT_ESSPRK
#if PRJ_TIMEINT_ESSPRK_N < 2
#error "TIME_INTEGRATION=eSSPRK<N>_2 requires N >= 2"
#endif
#elif TIME_INTEGRATION == PRJ_TIMEINT_ESSPRK9_3
/* Exact named third-order scheme. */
#elif TIME_INTEGRATION == PRJ_TIMEINT_IMEX
#if PRJ_TIMEINT_TABLEAU_NSTAGES < 1
#error "TIME_INTEGRATION=IMEX requires PRJ_TIMEINT_TABLEAU_NSTAGES >= 1"
#endif
#elif TIME_INTEGRATION != RK2
#error "Unsupported TIME_INTEGRATION value"
#endif
#if TIME_INTEGRATION == PRJ_TIMEINT_ESSPRK9_3
#define PRJ_TIMEINT_EXTRA_SAVED_STATES 1
#else
#define PRJ_TIMEINT_EXTRA_SAVED_STATES 0
#endif
#if TIME_INTEGRATION == PRJ_TIMEINT_IMEX
#define PRJ_BLOCK_NSTAGES PRJ_TIMEINT_TABLEAU_NSTAGES
#elif PRJ_TIMEINT_EXTRA_SAVED_STATES
#define PRJ_BLOCK_NSTAGES 4
#else
#define PRJ_BLOCK_NSTAGES 2
#endif
#if TIME_INTEGRATION == PRJ_TIMEINT_ESSPRK || TIME_INTEGRATION == PRJ_TIMEINT_ESSPRK9_3
#define PRJ_TIMEINT_USES_ESSPRK_STEP 1
#else
#define PRJ_TIMEINT_USES_ESSPRK_STEP 0
#endif

#if PRJ_USE_RADIATION_M1 || PRJ_USE_RADIATION_FSA
#ifndef PRJ_NRAD
#define PRJ_NRAD 3
#endif
#else
#undef PRJ_NRAD
#define PRJ_NRAD 0
#endif
#if PRJ_NRAD > 0
#ifndef NCLOSURE
#define NCLOSURE 100
#endif
#if NCLOSURE < 1
#error "NCLOSURE must be at least 1"
#endif
/* Radiation E/F can exceed the range of single-precision float; scaled float
   paths divide by RAD_SCALE before conversion and restore it afterward. */
#ifndef RAD_SCALE
#define RAD_SCALE 1e25
#endif
#endif
#ifndef PRJ_NEGROUP
#define PRJ_NEGROUP 12
#endif
#ifndef PRJ_N_ANGLE_LEV
#define PRJ_N_ANGLE_LEV 1
#endif
#if PRJ_USE_RADIATION_FSA && PRJ_N_ANGLE_LEV < 1
#error "PRJ_N_ANGLE_LEV must be at least 1 for RADIATION_FSA"
#endif
/* Fast-flavor conversion of neutrinos (default off).  Set here only. */
#ifndef DO_FFC
#define DO_FFC 0
#endif
#if DO_FFC && !PRJ_USE_RADIATION_FSA
#error "DO_FFC requires RADIATION_FSA"
#endif
#if DO_FFC && PRJ_NRAD < 3
#error "DO_FFC requires at least 3 neutrino species (nu_e, nubar_e, nu_x)"
#endif
#define PRJ_NANGLE (10 * (PRJ_N_ANGLE_LEV) * (PRJ_N_ANGLE_LEV) + 2)
#define PRJ_NARC (30 * (PRJ_N_ANGLE_LEV) * (PRJ_N_ANGLE_LEV))
#ifndef INEL_PHI_NT
#define INEL_PHI_NT 30
#endif
#ifndef PRJ_EXPE_NT
#define PRJ_EXPE_NT 1024
#endif
#if PRJ_MHD
#define PRJ_NHYDRO 9
#else
#define PRJ_NHYDRO 6
#endif
#if PRJ_USE_RADIATION_FSA
#define PRJ_NRAD_VAR (PRJ_NRAD * PRJ_NEGROUP * PRJ_NANGLE)
#else
#define PRJ_NRAD_VAR (PRJ_NRAD * PRJ_NEGROUP * (1 + PRJ_NDIM))
#endif
#define PRJ_NVAR_MHD_PRIM PRJ_NHYDRO
#define PRJ_NVAR_RAD_PRIM PRJ_NRAD_VAR
#define PRJ_NVAR_MHD_CONS PRJ_NHYDRO
#define PRJ_NVAR_RAD_CONS PRJ_NRAD_VAR
#define PRJ_NVAR_CONS (PRJ_NVAR_MHD_CONS + PRJ_NVAR_RAD_CONS)
#define PRJ_NVAR_PRIM (PRJ_NVAR_MHD_PRIM + PRJ_NVAR_RAD_PRIM)
#define PRJ_NVAR_EOSVAR 3
#ifndef PRJ_BLOCK_SIZE
#define PRJ_BLOCK_SIZE 8
#endif
#if (PRJ_BLOCK_SIZE & 1) != 0
#error "PRJ_BLOCK_SIZE must be even for AMR refinement"
#endif
#define PRJ_RECON_MC    0
#define PRJ_RECON_WENO3 1
#define PRJ_RECON_WENO7 2
#ifndef PRJ_RECON_HYDRO
#define PRJ_RECON_HYDRO PRJ_RECON_MC
#endif
#ifndef PRJ_RECON_RADIATION
#define PRJ_RECON_RADIATION PRJ_RECON_MC
#endif
/* Z4c reconstruction/finite-difference order is intentionally owned by this
 * header.  A command-line -DPRJ_RECON_Z4C_ORDER is discarded so Makefile flags
 * cannot silently change the Z4c stencil width. */
#ifdef PRJ_RECON_Z4C_ORDER
#undef PRJ_RECON_Z4C_ORDER
#endif
#define PRJ_RECON_Z4C_ORDER 4
#if PRJ_RECON_HYDRO != PRJ_RECON_MC && PRJ_RECON_HYDRO != PRJ_RECON_WENO3 && PRJ_RECON_HYDRO != PRJ_RECON_WENO7
#error "Unsupported PRJ_RECON_HYDRO value"
#endif
#if PRJ_RECON_RADIATION != PRJ_RECON_MC && PRJ_RECON_RADIATION != PRJ_RECON_WENO3 && PRJ_RECON_RADIATION != PRJ_RECON_WENO7
#error "Unsupported PRJ_RECON_RADIATION value"
#endif
#if PRJ_RECON_Z4C_ORDER != 2 && PRJ_RECON_Z4C_ORDER != 4
#error "Unsupported PRJ_RECON_Z4C_ORDER value"
#endif
#if PRJ_RECON_HYDRO == PRJ_RECON_WENO7
#define PRJ_RECON_HYDRO_NCELLS 7
#else
#define PRJ_RECON_HYDRO_NCELLS 3
#endif
#if PRJ_RECON_RADIATION == PRJ_RECON_WENO7
#define PRJ_RECON_RADIATION_NCELLS 7
#else
#define PRJ_RECON_RADIATION_NCELLS 3
#endif
#if PRJ_RECON_HYDRO_NCELLS > PRJ_RECON_RADIATION_NCELLS
#define PRJ_RECON_NCELLS PRJ_RECON_HYDRO_NCELLS
#else
#define PRJ_RECON_NCELLS PRJ_RECON_RADIATION_NCELLS
#endif
#ifndef PRJ_NGHOST
#if PRJ_RECON_HYDRO == PRJ_RECON_WENO7 || PRJ_RECON_RADIATION == PRJ_RECON_WENO7
#define PRJ_NGHOST 4
#else
#define PRJ_NGHOST 2
#endif
#endif
#if (PRJ_NGHOST & 1) != 0
#error "PRJ_NGHOST must be even for AMR prolongation alignment"
#endif
#if (PRJ_RECON_HYDRO == PRJ_RECON_WENO7 || PRJ_RECON_RADIATION == PRJ_RECON_WENO7) && PRJ_NGHOST < 4
#error "PRJ_RECON_WENO7 requires PRJ_NGHOST >= 4"
#endif
#ifndef PRJ_NGHOST_RAD
#if PRJ_RECON_RADIATION == PRJ_RECON_WENO7
#define PRJ_NGHOST_RAD 4
#else
#define PRJ_NGHOST_RAD 2
#endif
#endif
#if (PRJ_NGHOST_RAD & 1) != 0
#error "PRJ_NGHOST_RAD must be even for AMR prolongation alignment"
#endif
#if PRJ_NGHOST_RAD > PRJ_NGHOST
#error "PRJ_NGHOST_RAD must not exceed PRJ_NGHOST"
#endif
#if ((PRJ_NGHOST - PRJ_NGHOST_RAD) & 1) != 0
#error "PRJ_NGHOST - PRJ_NGHOST_RAD must be even for AMR prolongation alignment"
#endif
#define PRJ_NGHOST_Z4C PRJ_RECON_Z4C_ORDER
#if PRJ_NGHOST_Z4C < 2 || PRJ_NGHOST_Z4C > 4
#error "PRJ_NGHOST_Z4C must be in [2, 4]"
#endif
#if (PRJ_NGHOST_Z4C & 1) != 0
#error "PRJ_NGHOST_Z4C must be even for AMR prolongation alignment"
#endif
#ifndef PRJ_AMR_BUFFER_ZONE
#define PRJ_AMR_BUFFER_ZONE 2
#endif
#if PRJ_AMR_BUFFER_ZONE > PRJ_BLOCK_SIZE
#error "PRJ_AMR_BUFFER_ZONE must not exceed PRJ_BLOCK_SIZE"
#endif
#define PRJ_AMR_N 4
#define PRJ_PATH_MAX 1024
#define PRJ_BS (PRJ_BLOCK_SIZE + 2 * PRJ_NGHOST)
#define PRJ_BLOCK_NCELLS (PRJ_BS * PRJ_BS * PRJ_BS)
#define PRJ_BLOCK_NFACES ((PRJ_BS + 1) * PRJ_BS * PRJ_BS)
#define PRJ_BLOCK_NEDGES ((PRJ_BS + 1) * (PRJ_BS + 1) * PRJ_BS)
#define PRJ_BS_Z4C (PRJ_BLOCK_SIZE + 2 * PRJ_NGHOST_Z4C)
#define PRJ_BLOCK_NCELLS_Z4C (PRJ_BS_Z4C * PRJ_BS_Z4C * PRJ_BS_Z4C)

/* Untiled linear cell index within a block (k fastest).  Used for face/edge
 * fields and any per-cell SoA scratch that is swept with constant strides. */
#define LIDX(i, j, k) \
    (((i) + PRJ_NGHOST) * PRJ_BS * PRJ_BS + ((j) + PRJ_NGHOST) * PRJ_BS + ((k) + PRJ_NGHOST))
#define FACE_IDX(dir, i, j, k) \
    ((dir) == X1DIR ? \
        (((i) + PRJ_NGHOST) * PRJ_BS * PRJ_BS + ((j) + PRJ_NGHOST) * PRJ_BS + ((k) + PRJ_NGHOST)) : \
     (dir) == X2DIR ? \
        (((i) + PRJ_NGHOST) * (PRJ_BS + 1) * PRJ_BS + ((j) + PRJ_NGHOST) * PRJ_BS + ((k) + PRJ_NGHOST)) : \
        (((i) + PRJ_NGHOST) * PRJ_BS * (PRJ_BS + 1) + ((j) + PRJ_NGHOST) * (PRJ_BS + 1) + ((k) + PRJ_NGHOST)))
#define EDGE_IDX(dir, i, j, k) \
    ((dir) == X1DIR ? \
        (((i) + PRJ_NGHOST) * (PRJ_BS + 1) * (PRJ_BS + 1) + ((j) + PRJ_NGHOST) * (PRJ_BS + 1) + ((k) + PRJ_NGHOST)) : \
     (dir) == X2DIR ? \
        (((i) + PRJ_NGHOST) * PRJ_BS * (PRJ_BS + 1) + ((j) + PRJ_NGHOST) * (PRJ_BS + 1) + ((k) + PRJ_NGHOST)) : \
        (((i) + PRJ_NGHOST) * (PRJ_BS + 1) * PRJ_BS + ((j) + PRJ_NGHOST) * PRJ_BS + ((k) + PRJ_NGHOST)))
/* ---- Block memory layout (SoA at PRJ_AOSOA_W==1, AoSoA when >1) ------------
 * BIDX(nc,v,i,j,k) gives the address of component v (of an nc-component cell
 * buffer) at cell (i,j,k).  At PRJ_AOSOA_W==1 it is plain component-major SoA
 * (v*NCELLS + LIDX), reproducing the legacy layout exactly.  At W>1 it is
 * AoSoA: memory order [i][j][k_tile][var][k_lane] with the fast k axis tiled
 * into lanes of width W, so a cell's variables sit in one compact tile (better
 * locality) and the W cells of a tile vectorise per variable.
 * Per-array wrappers fix nc: VIDX for the NVAR-wide prim/cons buffers, EIDX for
 * eosvar, VRIDX for the 3-component v_riemann/v_face fields, IDX for scalar
 * (one component per cell) fields. */
/* AoSoA lane width (cells per tile along the fast k axis).  Must divide PRJ_BS
 * (12 for the default 2-ghost build, 16 for WENO7).  4 is the AVX2 double width
 * and divides both; set to 1 for the legacy plain-SoA layout. */
#ifndef PRJ_AOSOA_W
#define PRJ_AOSOA_W 4
#endif
#if (PRJ_BS % PRJ_AOSOA_W) != 0
#error "PRJ_AOSOA_W must divide PRJ_BS"
#endif
#define PRJ_AOSOA_KT (PRJ_BS / PRJ_AOSOA_W)

#if PRJ_AOSOA_W <= 1
#define BIDX(nc, v, i, j, k) ((v) * PRJ_BLOCK_NCELLS + LIDX(i, j, k))
#else
#define BIDX(nc, v, i, j, k) \
    ((((((i) + PRJ_NGHOST) * PRJ_BS + ((j) + PRJ_NGHOST)) * PRJ_AOSOA_KT \
        + ((k) + PRJ_NGHOST) / PRJ_AOSOA_W) * (nc) + (v)) * PRJ_AOSOA_W \
     + ((k) + PRJ_NGHOST) % PRJ_AOSOA_W)
#endif

#define IDX(i, j, k) BIDX(1, 0, i, j, k)
#define VIDX(v, i, j, k) BIDX(PRJ_NVAR_CONS, v, i, j, k)
#define MHDVIDX(v, i, j, k) BIDX(PRJ_NVAR_MHD_CONS, v, i, j, k)
#define RADVIDX(v, i, j, k) BIDX(PRJ_NVAR_RAD_CONS, v, i, j, k)
#define EIDX(v, i, j, k) BIDX(PRJ_NVAR_EOSVAR, v, i, j, k)
#define VRIDX(v, i, j, k) BIDX(PRJ_NDIM, v, i, j, k)
#define Z4CLIDX(i, j, k) \
    (((i) + PRJ_NGHOST_Z4C) * PRJ_BS_Z4C * PRJ_BS_Z4C + \
     ((j) + PRJ_NGHOST_Z4C) * PRJ_BS_Z4C + ((k) + PRJ_NGHOST_Z4C))
#define Z4CIDX(v, i, j, k) ((v) * PRJ_BLOCK_NCELLS_Z4C + Z4CLIDX(i, j, k))
#if PRJ_USE_RADIATION_FSA && PRJ_USE_RADIAL_FRAME_FSA
#define PRJ_FSA_ROT_IDX(row, col, i, j, k) BIDX(9, 3 * (row) + (col), i, j, k)
#define PRJ_FSA_ANG_GEOM_IDX(arc, d, i, j, k) \
    BIDX(3 * PRJ_NARC, 3 * (arc) + (d), i, j, k)
#endif
#define WIDX(v, i, j, k) \
    (((size_t)(v) * (size_t)PRJ_BLOCK_NSTAGES * (size_t)PRJ_BLOCK_NCELLS) + \
     (size_t)IDX(i, j, k))
#define PRJ_BLOCK_STAGE_W(W, stage) \
    (&(W)[(size_t)(stage) * (size_t)PRJ_BLOCK_NCELLS])
#define PRJ_BLOCK_STAGE_BF(BF, stage) \
    (&(BF)[(size_t)(stage) * (size_t)PRJ_BLOCK_NFACES])
#define PRJ_BLOCK_STAGE_Z4C(Z, stage) \
    (&(Z)[(size_t)(stage) * (size_t)PRJ_NZ4C * (size_t)PRJ_BLOCK_NCELLS_Z4C])

/* Legacy SoA per-variable plane base (valid only at PRJ_AOSOA_W==1).  Still used
 * by the pencil-reconstruction scratch; removed from block-array access in the
 * reconstruction rewrite. */
#define VPLANE(v) ((v) * PRJ_BLOCK_NCELLS)

enum prj_cons_var {
    PRJ_CONS_RHO = 0,
    PRJ_CONS_MOM1 = 1,
    PRJ_CONS_MOM2 = 2,
    PRJ_CONS_MOM3 = 3,
    PRJ_CONS_ETOT = 4,
    PRJ_CONS_YE = 5
#if PRJ_MHD
    ,
    PRJ_CONS_B1 = 6,
    PRJ_CONS_B2 = 7,
    PRJ_CONS_B3 = 8
#endif
};

enum prj_prim_var {
    PRJ_PRIM_RHO = 0,
    PRJ_PRIM_V1 = 1,
    PRJ_PRIM_V2 = 2,
    PRJ_PRIM_V3 = 3,
    PRJ_PRIM_EINT = 4,
    PRJ_PRIM_YE = 5
#if PRJ_MHD
    ,
    PRJ_PRIM_B1 = 6,
    PRJ_PRIM_B2 = 7,
    PRJ_PRIM_B3 = 8
#endif
};

enum prj_z4c_var {
    PRJ_Z4C_CHI = 0,
    PRJ_Z4C_GXX = 1,
    PRJ_Z4C_GXY = 2,
    PRJ_Z4C_GXZ = 3,
    PRJ_Z4C_GYY = 4,
    PRJ_Z4C_GYZ = 5,
    PRJ_Z4C_GZZ = 6,
    PRJ_Z4C_KHAT = 7,
    PRJ_Z4C_AXX = 8,
    PRJ_Z4C_AXY = 9,
    PRJ_Z4C_AXZ = 10,
    PRJ_Z4C_AYY = 11,
    PRJ_Z4C_AYZ = 12,
    PRJ_Z4C_AZZ = 13,
    PRJ_Z4C_GAMX = 14,
    PRJ_Z4C_GAMY = 15,
    PRJ_Z4C_GAMZ = 16,
    PRJ_Z4C_THETA = 17,
    PRJ_Z4C_ALPHA = 18,
    PRJ_Z4C_BETAX = 19,
    PRJ_Z4C_BETAY = 20,
    PRJ_Z4C_BETAZ = 21,
    PRJ_NZ4C = 22
};

enum prj_z4c_tmunu_var {
    PRJ_TMUNU_SXX = 0,
    PRJ_TMUNU_SXY = 1,
    PRJ_TMUNU_SXZ = 2,
    PRJ_TMUNU_SYY = 3,
    PRJ_TMUNU_SYZ = 4,
    PRJ_TMUNU_SZZ = 5,
    PRJ_TMUNU_E = 6,
    PRJ_TMUNU_SX = 7,
    PRJ_TMUNU_SY = 8,
    PRJ_TMUNU_SZ = 9,
    PRJ_NTMUNU = 10
};

#if PRJ_USE_RADIATION_FSA
#define PRJ_RAD_GROUP_STRIDE PRJ_NANGLE
#else
#define PRJ_RAD_GROUP_STRIDE (1 + PRJ_NDIM)
#endif
#define PRJ_CONS_RAD_E(field, group) (PRJ_NHYDRO + (((field) * PRJ_NEGROUP + (group)) * PRJ_RAD_GROUP_STRIDE))
#define PRJ_CONS_RAD_F1(field, group) (PRJ_CONS_RAD_E(field, group) + 1)
#define PRJ_CONS_RAD_F2(field, group) (PRJ_CONS_RAD_E(field, group) + 2)
#define PRJ_CONS_RAD_F3(field, group) (PRJ_CONS_RAD_E(field, group) + 3)
#define PRJ_PRIM_RAD_E(field, group) PRJ_CONS_RAD_E(field, group)
#define PRJ_PRIM_RAD_F1(field, group) PRJ_CONS_RAD_F1(field, group)
#define PRJ_PRIM_RAD_F2(field, group) PRJ_CONS_RAD_F2(field, group)
#define PRJ_PRIM_RAD_F3(field, group) PRJ_CONS_RAD_F3(field, group)
#define PRJ_CONS_RAD_I(field, group, angle) (PRJ_NHYDRO + (((field) * PRJ_NEGROUP + (group)) * PRJ_NANGLE + (angle)))
#define PRJ_PRIM_RAD_I(field, group, angle) PRJ_CONS_RAD_I(field, group, angle)
#define PRJ_RAD_PRIM_E(field, group) ((((field) * PRJ_NEGROUP + (group)) * PRJ_RAD_GROUP_STRIDE))
#define PRJ_RAD_PRIM_F1(field, group) (PRJ_RAD_PRIM_E(field, group) + 1)
#define PRJ_RAD_PRIM_F2(field, group) (PRJ_RAD_PRIM_E(field, group) + 2)
#define PRJ_RAD_PRIM_F3(field, group) (PRJ_RAD_PRIM_E(field, group) + 3)
#define PRJ_RAD_PRIM_I(field, group, angle) ((((field) * PRJ_NEGROUP + (group)) * PRJ_NANGLE + (angle)))
#define PRJ_RAD_CONS_E(field, group) PRJ_RAD_PRIM_E(field, group)
#define PRJ_RAD_CONS_F1(field, group) PRJ_RAD_PRIM_F1(field, group)
#define PRJ_RAD_CONS_F2(field, group) PRJ_RAD_PRIM_F2(field, group)
#define PRJ_RAD_CONS_F3(field, group) PRJ_RAD_PRIM_F3(field, group)
#define PRJ_RAD_CONS_I(field, group, angle) PRJ_RAD_PRIM_I(field, group, angle)

enum prj_dir {
    X1DIR = 0,
    X2DIR = 1,
    X3DIR = 2
};

enum prj_face_pos {
    FACE_L = 0,
    FACE_R = 1
};

enum prj_bc_type {
    PRJ_BC_OUTFLOW = 0,
    PRJ_BC_REFLECT = 1,
    PRJ_BC_USER = 2
};

enum prj_amr_estimator {
    PRJ_AMR_ESTIMATOR_LOEHNER = 0,
    PRJ_AMR_ESTIMATOR_VELOCITY = 1,
    PRJ_AMR_ESTIMATOR_PRESSURE_SCALE_HEIGHT = 2,
    PRJ_AMR_ESTIMATOR_FRACTIONAL_JUMP = 3
};

#if PRJ_MHD
enum prj_mhd_init_type {
    PRJ_MHD_INIT_UNIFORM = 0,
    PRJ_MHD_INIT_DIPOLE_CORE = 1
};

enum prj_mhd_fidelity {
    PRJ_MHD_FIDELITY_NONE = 0,
    PRJ_MHD_FIDELITY_COARSER = 1,
    PRJ_MHD_FIDELITY_SAME = 2,
    PRJ_MHD_FIDELITY_FINER = 3
};
#endif

enum prj_lohner_var {
    PRJ_LOHNER_VAR_DENSITY = 0,
    PRJ_LOHNER_VAR_PRESSURE = 1,
    PRJ_LOHNER_VAR_TEMPERATURE = 2
};

enum prj_fractional_jump_var {
    PRJ_FRACTIONAL_JUMP_VAR_DENSITY = 0,
    PRJ_FRACTIONAL_JUMP_VAR_PRESSURE = 1
};

enum prj_eos_kind {
    PRJ_EOS_KIND_IDEAL = 0,
    PRJ_EOS_KIND_TABLE = 1
};

enum {
    PRJ_BOUNDARY_FILL_NONRECON = 0,
    PRJ_BOUNDARY_FILL_RECON = 1,
    PRJ_BOUNDARY_FILL_ALL = 2,
    PRJ_BOUNDARY_FILL_SAME_LEVEL = 3,
    PRJ_BOUNDARY_FILL_RESTRICTION = 4,
    PRJ_BOUNDARY_FILL_PROLONGATION = 5
};

#define PRJ_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define PRJ_MIN(a,b) (((a) < (b)) ? (a) : (b))

#endif
