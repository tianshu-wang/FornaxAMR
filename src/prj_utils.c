#include "prj.h"

void prj_fill(double *data, size_t n, double value)
{
    size_t i;

    if (data == 0) {
        return;
    }

    for (i = 0; i < n; ++i) {
        data[i] = value;
    }
}
