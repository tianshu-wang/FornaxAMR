#include "prj.h"

int main(void)
{
    double data[4];

    prj_fill(data, 4U, 3.0);
    return data[0] == 3.0 && data[3] == 3.0 ? 0 : 1;
}
