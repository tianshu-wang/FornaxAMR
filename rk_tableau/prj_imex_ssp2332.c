#include "prj.h"

/*
 * prj_imex_ssp2332 -- the SSP2(3,3,2) additive IMEX Runge-Kutta method.
 *
 * Explicit A:                Implicit A (tilde):
 *   [ 0    0    0 ]            [ 1/5   0    0   ]
 *   [ 1/2  0    0 ]            [ 1/10  1/5  0   ]
 *   [ 1/2  1/2  0 ]            [ 1/3   1/3  1/3 ]
 * b_ex = [1/3, 1/3, 1/3]     b_im = [1/3, 1/3, 1/3]
 *
 * SSP coefficient r = 2.
 */

const double prj_imex_ssp2332_a_ex[3 * 3] = {
    0.0,       0.0,       0.0,
    1.0 / 2.0, 0.0,       0.0,
    1.0 / 2.0, 1.0 / 2.0, 0.0
};

const double prj_imex_ssp2332_a_im[3 * 3] = {
    1.0 / 5.0,  0.0,       0.0,
    1.0 / 10.0, 1.0 / 5.0, 0.0,
    1.0 / 3.0,  1.0 / 3.0, 1.0 / 3.0
};

const double prj_imex_ssp2332_b_ex[3] = {
    1.0 / 3.0,
    1.0 / 3.0,
    1.0 / 3.0
};

const double prj_imex_ssp2332_b_im[3] = {
    1.0 / 3.0,
    1.0 / 3.0,
    1.0 / 3.0
};

const prj_timeint_imex_tableau prj_imex_ssp2332 = {
    3,
    prj_imex_ssp2332_a_ex,
    prj_imex_ssp2332_a_im,
    prj_imex_ssp2332_b_ex,
    prj_imex_ssp2332_b_im,
    2.0
};
