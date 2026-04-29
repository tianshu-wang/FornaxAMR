#ifndef PRJ_DEFS_H
#define PRJ_DEFS_H

#define PRJ_NDIM 3
#ifndef PRJ_USE_GRAVITY
#define PRJ_USE_GRAVITY 1
#endif

#ifndef PRJ_USE_RADIATION
#define PRJ_USE_RADIATION 0
#endif

#ifndef PRJ_DUMP_SINGLE_PRECISION
#define PRJ_DUMP_SINGLE_PRECISION 1
#endif

#ifndef PRJ_MHD
#define PRJ_MHD 0
#endif

#if PRJ_USE_RADIATION
#ifndef PRJ_NRAD
#define PRJ_NRAD 3
#endif
#else
#undef PRJ_NRAD
#define PRJ_NRAD 0
#endif
#ifndef PRJ_NEGROUP
#define PRJ_NEGROUP 3
#endif
#if PRJ_MHD
#define PRJ_NHYDRO 9
#else
#define PRJ_NHYDRO 6
#endif
#define PRJ_NRAD_VAR (PRJ_NRAD * PRJ_NEGROUP * (1 + PRJ_NDIM))
#define PRJ_NVAR_CONS (PRJ_NHYDRO + PRJ_NRAD_VAR)
#define PRJ_NVAR_PRIM (PRJ_NHYDRO + PRJ_NRAD_VAR)
#define PRJ_NVAR_EOSVAR 3
#ifndef PRJ_BLOCK_SIZE
#define PRJ_BLOCK_SIZE 16
#endif
#define PRJ_NGHOST 2
#define PRJ_AMR_N 4
#define PRJ_PATH_MAX 1024
#define PRJ_BS (PRJ_BLOCK_SIZE + 2 * PRJ_NGHOST)
#define PRJ_BLOCK_NCELLS (PRJ_BS * PRJ_BS * PRJ_BS)

#define IDX(i, j, k) \
    (((i) + PRJ_NGHOST) * PRJ_BS * PRJ_BS + ((j) + PRJ_NGHOST) * PRJ_BS + ((k) + PRJ_NGHOST))
#define VIDX(v, i, j, k) ((v) * PRJ_BLOCK_NCELLS + IDX(i, j, k))
#define EIDX(v, i, j, k) ((v) * PRJ_BLOCK_NCELLS + IDX(i, j, k))

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

#define PRJ_RAD_GROUP_STRIDE (1 + PRJ_NDIM)
#define PRJ_CONS_RAD_E(field, group) (PRJ_NHYDRO + (((field) * PRJ_NEGROUP + (group)) * PRJ_RAD_GROUP_STRIDE))
#define PRJ_CONS_RAD_F1(field, group) (PRJ_CONS_RAD_E(field, group) + 1)
#define PRJ_CONS_RAD_F2(field, group) (PRJ_CONS_RAD_E(field, group) + 2)
#define PRJ_CONS_RAD_F3(field, group) (PRJ_CONS_RAD_E(field, group) + 3)
#define PRJ_PRIM_RAD_E(field, group) PRJ_CONS_RAD_E(field, group)
#define PRJ_PRIM_RAD_F1(field, group) PRJ_CONS_RAD_F1(field, group)
#define PRJ_PRIM_RAD_F2(field, group) PRJ_CONS_RAD_F2(field, group)
#define PRJ_PRIM_RAD_F3(field, group) PRJ_CONS_RAD_F3(field, group)

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
    PRJ_AMR_ESTIMATOR_DENSITY_JUMP = 3
};

#if PRJ_MHD
enum prj_mhd_init_type {
    PRJ_MHD_INIT_UNIFORM = 0,
    PRJ_MHD_INIT_DIPOLE_CORE = 1
};
#endif

enum prj_lohner_var {
    PRJ_LOHNER_VAR_DENSITY = 0,
    PRJ_LOHNER_VAR_PRESSURE = 1,
    PRJ_LOHNER_VAR_TEMPERATURE = 2,
    PRJ_LOHNER_VAR_LOG_DENSITY = 3
};

enum prj_eos_kind {
    PRJ_EOS_KIND_IDEAL = 0,
    PRJ_EOS_KIND_TABLE = 1
};

#endif
