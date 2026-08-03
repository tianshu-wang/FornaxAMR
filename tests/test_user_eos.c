#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "prj.h"

static int close_rel(double a, double b, double tol)
{
    return fabs(a - b) <= tol * fmax(1.0, fmax(fabs(a), fabs(b)));
}

int main(void)
{
#if PRJ_EOS_PROVIDER == PRJ_PROVIDER_USER
    prj_eos eos;
    double q[PRJ_EOS_NQUANT];
    double inverse[PRJ_EOS_NQUANT];
    double eint, pressure;
    double de[3], dp[3];
    const double rho = 2.3;
    const double temperature = 4.7;
    const double ye = 0.21;
    const double h[3] = {1.0e-6 * rho, 1.0e-6 * temperature, 1.0e-6};
    int d;

    memset(&eos, 0, sizeof(eos));
    prj_eos_init(&eos, 0);
    prj_eos_rty(&eos, rho, temperature, ye, q, PRJ_EOS_CTX_MAIN);
    prj_eos_rey(&eos, rho, q[PRJ_EOS_EINT], ye, inverse, PRJ_EOS_CTX_MAIN);
    assert(close_rel(inverse[PRJ_EOS_TEMPERATURE], temperature, 2.0e-15));
    assert(close_rel(inverse[PRJ_EOS_PRESSURE], q[PRJ_EOS_PRESSURE], 2.0e-15));
    assert(prj_eos_rty_derivs(&eos, rho, temperature, ye, &eint, &pressure,
        &de[0], &de[1], &de[2], &dp[0], &dp[1], &dp[2], PRJ_EOS_CTX_MAIN));
    for (d = 0; d < 3; ++d) {
        double xp[3] = {rho, temperature, ye};
        double xm[3] = {rho, temperature, ye};
        double qp[PRJ_EOS_NQUANT];
        double qm[PRJ_EOS_NQUANT];
        double de_fd;
        double dp_fd;

        xp[d] += h[d];
        xm[d] -= h[d];
        prj_eos_rty(&eos, xp[0], xp[1], xp[2], qp, PRJ_EOS_CTX_MAIN);
        prj_eos_rty(&eos, xm[0], xm[1], xm[2], qm, PRJ_EOS_CTX_MAIN);
        de_fd = (qp[PRJ_EOS_EINT] - qm[PRJ_EOS_EINT]) / (2.0 * h[d]);
        dp_fd = (qp[PRJ_EOS_PRESSURE] - qm[PRJ_EOS_PRESSURE]) / (2.0 * h[d]);
        assert(close_rel(de[d], de_fd, 2.0e-9));
        assert(close_rel(dp[d], dp_fd, 2.0e-9));
    }
#endif
    printf("user EOS tests passed\n");
    return 0;
}
