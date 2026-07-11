#include "prj.h"

/*
 * Classical explicit RK4 encoded as an additive IMEX tableau.  The implicit
 * half is zero, so this is useful for decoupled dynamics such as Z4c runs with
 * radiation disabled.
 */

const double prj_rk4_a_ex[4 * 4] = {
    0.0,       0.0,       0.0, 0.0,
    1.0 / 2.0, 0.0,       0.0, 0.0,
    0.0,       1.0 / 2.0, 0.0, 0.0,
    0.0,       0.0,       1.0, 0.0
};

const double prj_rk4_a_im[4 * 4] = {
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0
};

const double prj_rk4_b_ex[4] = {
    1.0 / 6.0,
    1.0 / 3.0,
    1.0 / 3.0,
    1.0 / 6.0
};

const double prj_rk4_b_im[4] = {
    0.0,
    0.0,
    0.0,
    0.0
};

const prj_timeint_imex_tableau prj_rk4 = {
    4,
    prj_rk4_a_ex,
    prj_rk4_a_im,
    prj_rk4_b_ex,
    prj_rk4_b_im,
    1.0
};
