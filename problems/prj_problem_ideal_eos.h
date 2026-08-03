#ifndef PRJ_PROBLEM_IDEAL_EOS_H
#define PRJ_PROBLEM_IDEAL_EOS_H

#define PRJ_IDEAL_EOS_ENERGY_SCALE 0.95655684e18

#define PRJ_DEFINE_IDEAL_EOS(ADIABATIC_INDEX) \
void eos_rty(double rho, double temperature, double ye, \
    double *eint, double *pressure, double *gamma, double *eta, \
    double deint[3], double dpressure[3]) \
{ \
    const double g = (ADIABATIC_INDEX); \
    (void)ye; \
    *eint = PRJ_IDEAL_EOS_ENERGY_SCALE * temperature / (g - 1.0); \
    *pressure = rho * PRJ_IDEAL_EOS_ENERGY_SCALE * temperature; \
    *gamma = g; \
    *eta = 0.0; \
    deint[0] = 0.0; \
    deint[1] = PRJ_IDEAL_EOS_ENERGY_SCALE / (g - 1.0); \
    deint[2] = 0.0; \
    dpressure[0] = PRJ_IDEAL_EOS_ENERGY_SCALE * temperature; \
    dpressure[1] = rho * PRJ_IDEAL_EOS_ENERGY_SCALE; \
    dpressure[2] = 0.0; \
} \
void eos_rey(double rho, double eint, double ye, \
    double *temperature, double *pressure, double *gamma) \
{ \
    const double g = (ADIABATIC_INDEX); \
    (void)ye; \
    *temperature = (g - 1.0) * eint / PRJ_IDEAL_EOS_ENERGY_SCALE; \
    *pressure = (g - 1.0) * rho * eint; \
    *gamma = g; \
}

#endif
