#include "prj.h"

const double prj_2s2pe2pi2plinKInf_a_ex[2 * 2] = {
    0.0, 0.0,
    1.0, 0.0
};

const double prj_2s2pe2pi2plinKInf_a_im[2 * 2] = {
    0.7192175840741879, 0.0,
    0.11776435415794313, 0.163018061767869
};

const double prj_2s2pe2pi2plinKInf_b_ex[2] = {
    0.5,
    0.5
};

const double prj_2s2pe2pi2plinKInf_b_im[2] = {
    0.49999999999999994,
    0.5
};

const prj_timeint_imex_tableau prj_2s2pe2pi2plinKInf = {
    2,
    prj_2s2pe2pi2plinKInf_a_ex,
    prj_2s2pe2pi2plinKInf_a_im,
    prj_2s2pe2pi2plinKInf_b_ex,
    prj_2s2pe2pi2plinKInf_b_im,
    1.0
};
