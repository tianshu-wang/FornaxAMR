#include "prj.h"

int main(int argc, char **argv)
{
    prj_mpi mpi;

    prj_mpi_init(&argc, &argv, &mpi);
    prj_mpi_finalize();
    return 0;
}
