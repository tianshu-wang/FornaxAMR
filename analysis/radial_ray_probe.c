#include <hdf5.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define META_COLS 13
#define HYDRO_VARS 6
#define EOS_VARS 3

typedef struct {
    double r;
    double perp;
    double rho;
    double vr;
    double speed;
    double eint;
    double ye;
    double pressure;
    double temperature;
    double gamma;
    double dx;
    int block;
    int cell;
    int level;
} sample;

static void fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}

static int compare_sample(const void *a, const void *b)
{
    const sample *sa = (const sample *)a;
    const sample *sb = (const sample *)b;

    if (sa->r < sb->r) return -1;
    if (sa->r > sb->r) return 1;
    if (sa->perp < sb->perp) return -1;
    if (sa->perp > sb->perp) return 1;
    return 0;
}

static void read_dataset_batch(hid_t dataset, int block0, int block_count,
    int nvars, int ncells, float *buffer)
{
    hid_t file_space = H5Dget_space(dataset);
    hsize_t start[3] = {(hsize_t)block0, 0, 0};
    hsize_t count[3] = {(hsize_t)block_count, (hsize_t)nvars, (hsize_t)ncells};
    hid_t mem_space = H5Screate_simple(3, count, 0);

    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, 0, count, 0);
    if (H5Dread(dataset, H5T_NATIVE_FLOAT, mem_space, file_space,
            H5P_DEFAULT, buffer) < 0) {
        fail("failed to read dataset batch");
    }
    H5Sclose(mem_space);
    H5Sclose(file_space);
}

int main(int argc, char **argv)
{
    const int batch_blocks = 128;
    const char *filename;
    double ray[3];
    double ray_norm;
    double rmin;
    double rmax;
    double max_perp_dx;
    hid_t file;
    hid_t data;
    hid_t eos;
    hid_t metadata_dset;
    hid_t attr;
    int nblocks;
    int block_size;
    int ncells;
    double x_com[3];
    double *metadata;
    float *hydro_buffer;
    float *eos_buffer;
    sample *samples = 0;
    size_t nsamples = 0;
    size_t capacity = 0;
    int block0;
    int d;

    if (argc != 8) {
        fprintf(stderr,
            "usage: %s dump.h5 ray_x ray_y ray_z rmin_cm rmax_cm max_perp_dx\n",
            argv[0]);
        return 2;
    }
    filename = argv[1];
    ray[0] = strtod(argv[2], 0);
    ray[1] = strtod(argv[3], 0);
    ray[2] = strtod(argv[4], 0);
    rmin = strtod(argv[5], 0);
    rmax = strtod(argv[6], 0);
    max_perp_dx = strtod(argv[7], 0);

    ray_norm = sqrt(ray[0] * ray[0] + ray[1] * ray[1] + ray[2] * ray[2]);
    if (!(ray_norm > 0.0) || !(rmax > rmin) || !(max_perp_dx > 0.0)) {
        fail("invalid ray or radial range");
    }
    for (d = 0; d < 3; ++d) ray[d] /= ray_norm;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) fail("failed to open dump");
    attr = H5Aopen(file, "nblocks", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &nblocks);
    H5Aclose(attr);
    attr = H5Aopen(file, "block_size", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &block_size);
    H5Aclose(attr);
    attr = H5Aopen(file, "x_com", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, x_com);
    H5Aclose(attr);
    ncells = block_size * block_size * block_size;

    metadata = (double *)malloc((size_t)nblocks * META_COLS * sizeof(*metadata));
    hydro_buffer = (float *)malloc(
        (size_t)batch_blocks * HYDRO_VARS * ncells * sizeof(*hydro_buffer));
    eos_buffer = (float *)malloc(
        (size_t)batch_blocks * EOS_VARS * ncells * sizeof(*eos_buffer));
    if (metadata == 0 || hydro_buffer == 0 || eos_buffer == 0) fail("allocation failed");

    metadata_dset = H5Dopen2(file, "MetaData", H5P_DEFAULT);
    {
        hid_t file_space = H5Dget_space(metadata_dset);
        hsize_t start[2] = {0, 0};
        hsize_t count[2] = {(hsize_t)nblocks, META_COLS};
        hid_t mem_space = H5Screate_simple(2, count, 0);

        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, 0, count, 0);
        H5Dread(metadata_dset, H5T_NATIVE_DOUBLE, mem_space, file_space,
            H5P_DEFAULT, metadata);
        H5Sclose(mem_space);
        H5Sclose(file_space);
    }
    H5Dclose(metadata_dset);

    data = H5Dopen2(file, "Data", H5P_DEFAULT);
    eos = H5Dopen2(file, "Eos", H5P_DEFAULT);
    for (block0 = 0; block0 < nblocks; block0 += batch_blocks) {
        int block_count = nblocks - block0;
        int local_block;

        if (block_count > batch_blocks) block_count = batch_blocks;
        read_dataset_batch(data, block0, block_count, HYDRO_VARS, ncells, hydro_buffer);
        read_dataset_batch(eos, block0, block_count, EOS_VARS, ncells, eos_buffer);

        for (local_block = 0; local_block < block_count; ++local_block) {
            int block = block0 + local_block;
            const double *meta = &metadata[(size_t)block * META_COLS];
            size_t hydro_base = (size_t)local_block * HYDRO_VARS * ncells;
            size_t eos_base = (size_t)local_block * EOS_VARS * ncells;
            int cell;

            if ((int)meta[3] != 1) continue;
            for (cell = 0; cell < ncells; ++cell) {
                int i = cell / (block_size * block_size);
                int j = (cell / block_size) % block_size;
                int k = cell % block_size;
                int ijk[3] = {i, j, k};
                double xyz[3];
                double radius2 = 0.0;
                double along = 0.0;
                double perp2;
                double velocity[3];
                double speed2 = 0.0;
                sample s;

                for (d = 0; d < 3; ++d) {
                    xyz[d] = meta[4 + d] + ((double)ijk[d] + 0.5) * meta[10 + d] -
                        x_com[d];
                    radius2 += xyz[d] * xyz[d];
                    along += xyz[d] * ray[d];
                }
                if (along <= 0.0 || radius2 < rmin * rmin || radius2 > rmax * rmax) {
                    continue;
                }
                perp2 = radius2 - along * along;
                if (perp2 < 0.0) perp2 = 0.0;
                if (sqrt(perp2) > max_perp_dx * meta[10]) continue;

                for (d = 0; d < 3; ++d) {
                    velocity[d] = hydro_buffer[
                        hydro_base + (size_t)(d + 1) * ncells + cell];
                    speed2 += velocity[d] * velocity[d];
                }
                s.r = sqrt(radius2);
                s.perp = sqrt(perp2);
                s.rho = hydro_buffer[hydro_base + cell];
                s.vr = (velocity[0] * xyz[0] + velocity[1] * xyz[1] +
                    velocity[2] * xyz[2]) / s.r;
                s.speed = sqrt(speed2);
                s.eint = hydro_buffer[hydro_base + (size_t)4 * ncells + cell];
                s.ye = hydro_buffer[hydro_base + (size_t)5 * ncells + cell];
                s.pressure = eos_buffer[eos_base + cell];
                s.temperature = eos_buffer[eos_base + (size_t)ncells + cell];
                s.gamma = eos_buffer[eos_base + (size_t)2 * ncells + cell];
                s.dx = meta[10];
                s.block = block;
                s.cell = cell;
                s.level = (int)meta[2];

                if (nsamples == capacity) {
                    size_t new_capacity = capacity == 0 ? 256 : 2 * capacity;
                    sample *new_samples = (sample *)realloc(
                        samples, new_capacity * sizeof(*samples));

                    if (new_samples == 0) fail("sample allocation failed");
                    samples = new_samples;
                    capacity = new_capacity;
                }
                samples[nsamples++] = s;
            }
        }
    }
    H5Dclose(eos);
    H5Dclose(data);
    H5Fclose(file);

    qsort(samples, nsamples, sizeof(*samples), compare_sample);
    printf("r_cm,perp_cm,block,cell,level,dx_cm,rho,vr,speed,eint,Ye,"
           "pressure,temperature,gamma,cs_proxy,mach_proxy\n");
    for (size_t idx = 0; idx < nsamples; ++idx) {
        const sample *s = &samples[idx];
        double cs = s->rho > 0.0 && s->pressure > 0.0 && s->gamma > 0.0 ?
            sqrt(s->gamma * s->pressure / s->rho) : 0.0;

        printf("%.9e,%.9e,%d,%d,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,"
               "%.9e,%.9e,%.9e,%.9e,%.9e\n",
            s->r, s->perp, s->block, s->cell, s->level, s->dx,
            s->rho, s->vr, s->speed, s->eint, s->ye, s->pressure,
            s->temperature, s->gamma, cs, cs > 0.0 ? fabs(s->vr) / cs : 0.0);
    }

    free(samples);
    free(eos_buffer);
    free(hydro_buffer);
    free(metadata);
    return 0;
}
