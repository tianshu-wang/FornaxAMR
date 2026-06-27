#include <string.h>

#include "prj_rk_tableau.h"

void prj_2s2pe2pi2plinKInf(prj_timeint_imex_tableau *tableau)
{
    if (tableau == 0) {
        return;
    }

    memset(tableau, 0, sizeof(*tableau));
    tableau->nstages = 2;
    tableau->a_ex[0][0] = 0.0;
    tableau->a_ex[0][1] = 0.0;
    tableau->a_ex[1][0] = 1.0;
    tableau->a_ex[1][1] = 0.0;

    tableau->a_im[0][0] = 0.7192175840741879;
    tableau->a_im[0][1] = 0.0;
    tableau->a_im[1][0] = 0.11776435415794313;
    tableau->a_im[1][1] = 0.163018061767869;

    tableau->b_ex[0] = 0.5;
    tableau->b_ex[1] = 0.5;
    tableau->b_im[0] = 0.49999999999999994;
    tableau->b_im[1] = 0.5;
    tableau->r = 1.0;
}
