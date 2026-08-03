#include "prj.h"
#include "../../problems/prj_problem_ideal_eos.h"

PRJ_DEFINE_IDEAL_EOS(5.0 / 3.0)

void prj_problem_user_provider(prj_sim *sim, prj_mpi *mpi)
{
    (void)sim;
    (void)mpi;
}

static void fixture_derivatives(double drho, double dtemp, double dye,
    double derivatives[3])
{
    derivatives[0] = drho;
    derivatives[1] = dtemp;
    derivatives[2] = dye;
}

void opac(const prj_rad *rad, int species, int group, double rho,
    double temperature, double ye, double *value, double derivatives[3])
{
    double scale = (double)(species + 1) * (double)(group + 1);
    (void)rad;
    *value = scale * rho * temperature * (1.0 + ye);
    fixture_derivatives(scale * temperature * (1.0 + ye),
        scale * rho * (1.0 + ye), scale * rho * temperature, derivatives);
}

void emis(const prj_rad *rad, int species, int group, double rho,
    double temperature, double ye, double *value, double derivatives[3])
{
    (void)rad;
    *value = rho * rho + temperature * temperature + ye + species + group;
    fixture_derivatives(2.0 * rho, 2.0 * temperature, 1.0, derivatives);
}

void scat(const prj_rad *rad, int species, int group, double rho,
    double temperature, double ye, double *value, double derivatives[3])
{
    (void)rad;
    *value = 2.0 * rho + 3.0 * temperature + 4.0 * ye + species + group;
    fixture_derivatives(2.0, 3.0, 4.0, derivatives);
}

void delta(const prj_rad *rad, int species, int group, double rho,
    double temperature, double ye, double *value, double derivatives[3])
{
    (void)rad;
    (void)species;
    (void)group;
    *value = rho - temperature + ye;
    fixture_derivatives(1.0, -1.0, 1.0, derivatives);
}
