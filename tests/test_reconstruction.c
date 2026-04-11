#include <math.h>

#include "prj_reconstruct.h"

int main(void)
{
    double stencil[3] = {1.0, 3.0, 5.0};
    double slope = prj_reconstruct_slope(stencil, 2.0);

    return fabs(slope - 1.0) < 1.0e-12 ? 0 : 1;
}
