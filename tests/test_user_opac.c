#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "prj.h"

#if PRJ_NRAD > 0 && PRJ_OPAC_PROVIDER == PRJ_PROVIDER_USER
static int close_rel(double a, double b)
{
    return fabs(a - b) <= 2.0e-13 * fmax(1.0, fmax(fabs(a), fabs(b)));
}
#endif

int main(void)
{
#if PRJ_NRAD > 0 && PRJ_OPAC_PROVIDER == PRJ_PROVIDER_USER
    prj_rad rad;
    double k[PRJ_RAD3_OPAC_NGROUPS], s[PRJ_RAD3_OPAC_NGROUPS];
    double d[PRJ_RAD3_OPAC_NGROUPS], e[PRJ_RAD3_OPAC_NGROUPS];
    double dk0[PRJ_RAD3_OPAC_NGROUPS], dk1[PRJ_RAD3_OPAC_NGROUPS], dk2[PRJ_RAD3_OPAC_NGROUPS];
    double ds0[PRJ_RAD3_OPAC_NGROUPS], ds1[PRJ_RAD3_OPAC_NGROUPS], ds2[PRJ_RAD3_OPAC_NGROUPS];
    double dd0[PRJ_RAD3_OPAC_NGROUPS], dd1[PRJ_RAD3_OPAC_NGROUPS], dd2[PRJ_RAD3_OPAC_NGROUPS];
    double de0[PRJ_RAD3_OPAC_NGROUPS], de1[PRJ_RAD3_OPAC_NGROUPS], de2[PRJ_RAD3_OPAC_NGROUPS];
    double lkt[PRJ_RAD3_OPAC_NGROUPS], lky[PRJ_RAD3_OPAC_NGROUPS];
    double let[PRJ_RAD3_OPAC_NGROUPS], ley[PRJ_RAD3_OPAC_NGROUPS];
    const double rho = 2.0, temp = 3.0, ye = 0.25;
    int q;

    memset(&rad, 0, sizeof(rad));
    for (q = 0; q < PRJ_NRAD; ++q) {
        rad.emin[q] = 1.0;
        rad.emax[q] = 9.0;
    }
    prj_rad3_opac_init(&rad);
    assert(rad.freq_grid_valid[0] == 1);
    assert(rad.absopac[0] == 0 && rad.emis[0] == 0 && rad.scaopac[0] == 0);
    prj_rad3_opac_lookup_derivs(&rad, rho, temp, ye, k, s, d, e,
        dk0, dk1, dk2, ds0, ds1, ds2, dd0, dd1, dd2, de0, de1, de2);
    for (q = 0; q < PRJ_RAD3_OPAC_NGROUPS; ++q) {
        int species = q / PRJ_NEGROUP;
        int group = q % PRJ_NEGROUP;
        double scale = (double)(species + 1) * (double)(group + 1);
        assert(close_rel(k[q], scale * rho * temp * (1.0 + ye)));
        assert(close_rel(dk0[q], scale * temp * (1.0 + ye)));
        assert(close_rel(dk1[q], scale * rho * (1.0 + ye)));
        assert(close_rel(dk2[q], scale * rho * temp));
        assert(close_rel(ds0[q], 2.0) && close_rel(ds1[q], 3.0) && close_rel(ds2[q], 4.0));
        assert(close_rel(dd0[q], 1.0) && close_rel(dd1[q], -1.0) && close_rel(dd2[q], 1.0));
        assert(close_rel(de0[q], 2.0 * rho));
        assert(close_rel(de1[q], 2.0 * temp));
        assert(close_rel(de2[q], 1.0));
    }
    prj_rad3_opac_lookup_ke(&rad, rho, temp, ye, k, e, lkt, lky, let, ley);
    for (q = 0; q < PRJ_RAD3_OPAC_NGROUPS; ++q) {
        assert(close_rel(lkt[q], 1.0));
        assert(close_rel(lky[q], 1.0 / (1.0 + ye)));
        assert(close_rel(let[q], temp * (2.0 * temp) / e[q]));
        assert(close_rel(ley[q], 1.0 / e[q]));
    }
    prj_rad3_opac_free(&rad);
#endif
    printf("user opacity tests passed\n");
    return 0;
}
