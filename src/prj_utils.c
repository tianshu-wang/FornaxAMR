#include "prj.h"

#include <stdio.h>
#include <stdlib.h>
#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

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

/* Report a failed allocation and tear the job down.  Kept out of line so the
 * fast path in the allocation wrappers stays a single null check. */
static void prj_alloc_abort(const char *what, size_t bytes, const char *file, int line)
{
    int rank = -1;

#if defined(PRJ_ENABLE_MPI)
    int mpi_ready = 0;

    MPI_Initialized(&mpi_ready);
    if (mpi_ready) {
        int finalized = 0;

        MPI_Finalized(&finalized);
        if (!finalized) {
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        }
    }
#endif

    fprintf(stderr, "[ALLOC-FAIL] rank %d: %s of %zu bytes failed at %s:%d\n",
        rank, what, bytes, file, line);
    fflush(stderr);

#if defined(PRJ_ENABLE_MPI)
    {
        int mpi_ready = 0;

        MPI_Initialized(&mpi_ready);
        if (mpi_ready) {
            int finalized = 0;

            MPI_Finalized(&finalized);
            if (!finalized) {
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }
#endif
    abort();
}

void *prj_malloc_impl(size_t size, const char *file, int line)
{
    void *p = malloc(size);

    if (p == 0 && size != 0) {
        prj_alloc_abort("malloc", size, file, line);
    }
    return p;
}

void *prj_calloc_impl(size_t nmemb, size_t size, const char *file, int line)
{
    void *p = calloc(nmemb, size);

    if (p == 0 && nmemb != 0 && size != 0) {
        prj_alloc_abort("calloc", nmemb * size, file, line);
    }
    return p;
}

void *prj_realloc_impl(void *ptr, size_t size, const char *file, int line)
{
    void *p = realloc(ptr, size);

    if (p == 0 && size != 0) {
        prj_alloc_abort("realloc", size, file, line);
    }
    return p;
}
