#include <math.h>
#include <stdio.h>

#include "prj.h"

int main(void)
{
    prj_z4c_params params;
    int status;

    prj_z4c_init_params(&params);
    if (params.wave_radii_cm != 1.0e9 || params.wave_lmax != 4) {
        fprintf(stderr, "unexpected Z4c waveform defaults\n");
        return 1;
    }
    status = prj_z4c_wave_self_test();
    if (status != 0) {
        fprintf(stderr, "Z4c waveform self-test failed: %d\n", status);
        return 1;
    }
    return 0;
}
