#include <math.h>
#include <string.h>

#include "prj_types.h"
#include "prj_eos.h"

#define PRJ_TABULATED_EOS_PATH "../eos_tmp/SFHoEOS__ye__0.035_0.56_50__logT_-4.793_2.176_500__logrho_-8.699_15.5_500_extend.dat"

static int prj_almost_equal(double a, double b)
{
    return fabs(a - b) < 1.0e-12;
}

int main(void)
{
    prj_eos eos;
    double from_rty[PRJ_EOS_NQUANT];
    double from_rey[PRJ_EOS_NQUANT];
    double rho = 2.0;
    double T = 3.0;
    double ye = 0.1;

    eos.filename[0] = '\0';
    prj_eos_init(&eos);

    prj_eos_rty(&eos, rho, T, ye, from_rty);
    prj_eos_rey(&eos, rho, from_rty[PRJ_EOS_EINT], ye, from_rey);

    if (!prj_almost_equal(from_rty[PRJ_EOS_TEMPERATURE], from_rey[PRJ_EOS_TEMPERATURE])) {
        return 1;
    }
    if (!prj_almost_equal(from_rty[PRJ_EOS_EINT], from_rey[PRJ_EOS_EINT])) {
        return 1;
    }
    if (!prj_almost_equal(from_rty[PRJ_EOS_PRESSURE], from_rey[PRJ_EOS_PRESSURE])) {
        return 1;
    }
    if (!prj_almost_equal(from_rty[PRJ_EOS_GAMMA], from_rey[PRJ_EOS_GAMMA])) {
        return 1;
    }

    memset(&eos, 0, sizeof(eos));
    strncpy(eos.filename, PRJ_TABULATED_EOS_PATH, sizeof(eos.filename) - 1);
    eos.filename[sizeof(eos.filename) - 1] = '\0';
    prj_eos_init(&eos);
    if (eos.table_loaded != 1 || eos.table == 0) {
        return 2;
    }

    rho = 1.0e10;
    T = 1.0;
    ye = 0.2;
    prj_eos_rty(&eos, rho, T, ye, from_rty);
    prj_eos_rey(&eos, rho, from_rty[PRJ_EOS_EINT], ye, from_rey);

    if (!(from_rty[PRJ_EOS_PRESSURE] > 0.0) ||
        !(from_rty[PRJ_EOS_EINT] > 0.0) ||
        !(from_rty[PRJ_EOS_GAMMA] > 0.0) ||
        !(from_rty[PRJ_EOS_TEMPERATURE] > 0.0)) {
        return 3;
    }
    if (fabs((from_rty[PRJ_EOS_TEMPERATURE] - from_rey[PRJ_EOS_TEMPERATURE]) / from_rty[PRJ_EOS_TEMPERATURE]) > 5.0e-6) {
        return 4;
    }
    if (fabs((from_rty[PRJ_EOS_PRESSURE] - from_rey[PRJ_EOS_PRESSURE]) / from_rty[PRJ_EOS_PRESSURE]) > 5.0e-6) {
        return 5;
    }
    if (fabs((from_rty[PRJ_EOS_GAMMA] - from_rey[PRJ_EOS_GAMMA]) / from_rty[PRJ_EOS_GAMMA]) > 5.0e-6) {
        return 6;
    }

    return 0;
}
