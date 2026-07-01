#ifndef PRJ_UTILS_H
#define PRJ_UTILS_H

#include <stddef.h>

void prj_fill(double *data, size_t n, double value);

/* Checked allocation helpers.  These wrap the standard allocators one-to-one
 * (same arguments, same semantics) but abort the whole job with a diagnostic
 * message -- rank, requested byte count, and call site -- if the allocation
 * fails.  Call sites use the prj_malloc/prj_calloc/prj_realloc macros so the
 * originating file and line are captured automatically; the *_impl functions
 * are the actual implementations and are not meant to be called directly. */
void *prj_malloc_impl(size_t size, const char *file, int line);
void *prj_calloc_impl(size_t nmemb, size_t size, const char *file, int line);
void *prj_realloc_impl(void *ptr, size_t size, const char *file, int line);

#define prj_malloc(size) prj_malloc_impl((size), __FILE__, __LINE__)
#define prj_calloc(nmemb, size) prj_calloc_impl((nmemb), (size), __FILE__, __LINE__)
#define prj_realloc(ptr, size) prj_realloc_impl((ptr), (size), __FILE__, __LINE__)

/* Trilinear interpolation that also returns the three partial derivatives with
 * respect to the normalized cell coordinates (d0, d1, d2), each in [0, 1].
 * Corner values are indexed by bit pattern: v[i0 + 2*i1 + 4*i2] is the corner at
 * (d0=i0, d1=i1, d2=i2).  Callers multiply each returned partial by the inverse
 * grid spacing of the corresponding axis to get a coordinate derivative. */
static inline double prj_trilinear_with_deriv(const double v[8],
    double d0, double d1, double d2,
    double *dfdd0, double *dfdd1, double *dfdd2)
{
    double o0 = 1.0 - d0;
    double o1 = 1.0 - d1;
    double o2 = 1.0 - d2;

    *dfdd0 = o1 * o2 * (v[1] - v[0]) + d1 * o2 * (v[3] - v[2]) +
             o1 * d2 * (v[5] - v[4]) + d1 * d2 * (v[7] - v[6]);
    *dfdd1 = o0 * o2 * (v[2] - v[0]) + d0 * o2 * (v[3] - v[1]) +
             o0 * d2 * (v[6] - v[4]) + d0 * d2 * (v[7] - v[5]);
    *dfdd2 = o0 * o1 * (v[4] - v[0]) + d0 * o1 * (v[5] - v[1]) +
             o0 * d1 * (v[6] - v[2]) + d0 * d1 * (v[7] - v[3]);

    return o0 * o1 * o2 * v[0] + d0 * o1 * o2 * v[1] +
           o0 * d1 * o2 * v[2] + d0 * d1 * o2 * v[3] +
           o0 * o1 * d2 * v[4] + d0 * o1 * d2 * v[5] +
           o0 * d1 * d2 * v[6] + d0 * d1 * d2 * v[7];
}

#endif
