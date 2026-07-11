#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <hdf5.h>

enum { NZ4C = 22 };

static const char *const z4c_names[NZ4C] = {
    "chi",
    "gxx", "gxy", "gxz", "gyy", "gyz", "gzz",
    "Khat",
    "Axx", "Axy", "Axz", "Ayy", "Ayz", "Azz",
    "Gamx", "Gamy", "Gamz",
    "Theta",
    "alpha",
    "betax", "betay", "betaz"
};

static int read_attr_int(hid_t file, const char *name, int fallback)
{
    hid_t attr;
    int value = fallback;

    if (H5Aexists(file, name) <= 0) {
        return fallback;
    }
    attr = H5Aopen(file, name, H5P_DEFAULT);
    if (attr < 0) {
        return fallback;
    }
    (void)H5Aread(attr, H5T_NATIVE_INT, &value);
    H5Aclose(attr);
    return value;
}

static double read_attr_double(hid_t file, const char *name, double fallback)
{
    hid_t attr;
    double value = fallback;

    if (H5Aexists(file, name) <= 0) {
        return fallback;
    }
    attr = H5Aopen(file, name, H5P_DEFAULT);
    if (attr < 0) {
        return fallback;
    }
    (void)H5Aread(attr, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attr);
    return value;
}

int main(int argc, char **argv)
{
    const char *filename;
    hid_t file;
    hid_t dset;
    hid_t space;
    hsize_t dims[3] = {0, 0, 0};
    int ndims;
    double *z4c;
    hsize_t nblocks;
    hsize_t ncells;
    int var;

    if (argc != 2) {
        fprintf(stderr, "usage: %s output/dump_00000.h5\n", argv[0]);
        return 2;
    }
    filename = argv[1];
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "%s: failed to open file\n", filename);
        return 1;
    }
    if (H5Lexists(file, "Z4c", H5P_DEFAULT) <= 0) {
        fprintf(stderr, "%s: no Z4c dataset\n", filename);
        H5Fclose(file);
        return 1;
    }
    dset = H5Dopen2(file, "Z4c", H5P_DEFAULT);
    space = H5Dget_space(dset);
    ndims = H5Sget_simple_extent_dims(space, dims, 0);
    if (ndims != 3 || dims[1] != NZ4C) {
        fprintf(stderr, "%s: unexpected Z4c shape\n", filename);
        H5Sclose(space);
        H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }
    nblocks = dims[0];
    ncells = dims[2];
    z4c = (double *)calloc((size_t)nblocks * (size_t)NZ4C * (size_t)ncells, sizeof(*z4c));
    if (z4c == 0) {
        fprintf(stderr, "%s: allocation failed\n", filename);
        H5Sclose(space);
        H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }
    if (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, z4c) < 0) {
        fprintf(stderr, "%s: failed to read Z4c\n", filename);
        free(z4c);
        H5Sclose(space);
        H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }

    printf("# file=%s step=%d time=%.17e shape=[%llu,%d,%llu]\n",
        filename, read_attr_int(file, "step", -1),
        read_attr_double(file, "time", NAN),
        (unsigned long long)nblocks, NZ4C, (unsigned long long)ncells);
    for (var = 0; var < NZ4C; ++var) {
        double minv = DBL_MAX;
        double maxv = -DBL_MAX;
        double maxabs = 0.0;
        unsigned long long finite = 0;
        hsize_t b;

        for (b = 0; b < nblocks; ++b) {
            hsize_t c;

            for (c = 0; c < ncells; ++c) {
                size_t idx = ((size_t)b * (size_t)NZ4C + (size_t)var) * (size_t)ncells +
                    (size_t)c;
                double value = z4c[idx];

                if (!isfinite(value)) {
                    continue;
                }
                if (value < minv) minv = value;
                if (value > maxv) maxv = value;
                if (fabs(value) > maxabs) maxabs = fabs(value);
                finite += 1ULL;
            }
        }
        if (finite == 0ULL) {
            printf("%7s finite=0\n", z4c_names[var]);
        } else {
            printf("%7s finite=%llu min=%.9e max=%.9e maxabs=%.9e\n",
                z4c_names[var], finite, minv, maxv, maxabs);
        }
    }

    free(z4c);
    H5Sclose(space);
    H5Dclose(dset);
    H5Fclose(file);
    return 0;
}
