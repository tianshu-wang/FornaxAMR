#include <hdf5.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define META_COLS 13
#define NVAR_MATCH 6
#define KEEP 64

typedef struct {
    double score;
    double radiation_score;
    double flux_direction_score;
    int block;
    int cell;
    double rho;
    double v[3];
    double ye;
} match;

static void fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}

static void read_attr_int(hid_t file, const char *name, int *value)
{
    hid_t attr = H5Aopen(file, name, H5P_DEFAULT);

    if (attr < 0 || H5Aread(attr, H5T_NATIVE_INT, value) < 0) {
        fail("failed to read integer attribute");
    }
    H5Aclose(attr);
}

static void read_attr_double3(hid_t file, const char *name, double value[3])
{
    hid_t attr = H5Aopen(file, name, H5P_DEFAULT);

    if (attr < 0 || H5Aread(attr, H5T_NATIVE_DOUBLE, value) < 0) {
        fail("failed to read vector attribute");
    }
    H5Aclose(attr);
}

static void keep_match(match best[KEEP], int *nkeep, const match *candidate)
{
    int worst = 0;
    int n;

    if (*nkeep < KEEP) {
        best[*nkeep] = *candidate;
        *nkeep += 1;
        return;
    }
    for (n = 1; n < KEEP; ++n) {
        if (best[n].score > best[worst].score) {
            worst = n;
        }
    }
    if (candidate->score < best[worst].score) {
        best[worst] = *candidate;
    }
}

static int compare_match(const void *a, const void *b)
{
    const match *ma = (const match *)a;
    const match *mb = (const match *)b;

    if (ma->score < mb->score) return -1;
    if (ma->score > mb->score) return 1;
    return 0;
}

static int read_raw_u(const char *filename, double *u, int nvar)
{
    FILE *fp = fopen(filename, "r");
    char line[512];
    int *found;
    int nfound = 0;

    if (fp == 0) return 1;
    found = (int *)calloc((size_t)nvar, sizeof(*found));
    if (found == 0) fail("raw-input allocation failed");
    while (fgets(line, sizeof(line), fp) != 0) {
        int idx;
        double value;

        if (sscanf(line, " [%d] = %lf,", &idx, &value) == 2 &&
            idx >= 0 && idx < nvar) {
            u[idx] = value;
            if (!found[idx]) {
                found[idx] = 1;
                nfound += 1;
            }
        }
    }
    free(found);
    fclose(fp);
    return nfound == nvar ? 0 : 1;
}

int main(int argc, char **argv)
{
    const int batch_blocks = 128;
    const char *filename;
    double target_rho;
    double target_v[3];
    double target_ye;
    double target_u[153] = {0};
    int have_target_u = 0;
    hid_t file;
    hid_t data;
    hid_t metadata_dset;
    int nblocks;
    int block_size;
    int ncells;
    double x_com[3];
    double *metadata;
    float *buffer;
    match best[KEEP];
    int nkeep = 0;
    int active_count = 0;
    int block0;
    int n;

    if (argc != 7 && argc != 8) {
        fprintf(stderr, "usage: %s dump.h5 rho v1 v2 v3 Ye [crash.err]\n", argv[0]);
        return 2;
    }
    filename = argv[1];
    target_rho = strtod(argv[2], 0);
    target_v[0] = strtod(argv[3], 0);
    target_v[1] = strtod(argv[4], 0);
    target_v[2] = strtod(argv[5], 0);
    target_ye = strtod(argv[6], 0);
    if (argc == 8) {
        if (read_raw_u(argv[7], target_u, 153) != 0) {
            fail("failed to parse raw u from .err file");
        }
        have_target_u = 1;
    }

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) fail("failed to open dump");
    read_attr_int(file, "nblocks", &nblocks);
    read_attr_int(file, "block_size", &block_size);
    read_attr_double3(file, "x_com", x_com);
    ncells = block_size * block_size * block_size;

    metadata = (double *)malloc((size_t)nblocks * META_COLS * sizeof(*metadata));
    buffer = (float *)malloc((size_t)batch_blocks * NVAR_MATCH * ncells * sizeof(*buffer));
    if (metadata == 0 || buffer == 0) fail("allocation failed");

    metadata_dset = H5Dopen2(file, "MetaData", H5P_DEFAULT);
    if (metadata_dset < 0) fail("failed to open MetaData");
    {
        hid_t file_space = H5Dget_space(metadata_dset);
        hsize_t start[2] = {0, 0};
        hsize_t count[2] = {(hsize_t)nblocks, META_COLS};
        hid_t mem_space = H5Screate_simple(2, count, 0);

        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, 0, count, 0);
        if (H5Dread(metadata_dset, H5T_NATIVE_DOUBLE, mem_space, file_space,
                H5P_DEFAULT, metadata) < 0) {
            fail("failed to read metadata");
        }
        H5Sclose(mem_space);
        H5Sclose(file_space);
    }
    H5Dclose(metadata_dset);
    for (n = 0; n < nblocks; ++n) {
        if ((int)metadata[(size_t)n * META_COLS + 3] == 1) {
            active_count += 1;
        }
    }
    fprintf(stderr, "nblocks=%d active=%d block_size=%d\n",
        nblocks, active_count, block_size);

    data = H5Dopen2(file, "Data", H5P_DEFAULT);
    if (data < 0) fail("failed to open Data");
    for (block0 = 0; block0 < nblocks; block0 += batch_blocks) {
        int block_count = nblocks - block0;
        hid_t file_space;
        hid_t mem_space;
        hsize_t start[3] = {(hsize_t)block0, 0, 0};
        hsize_t count[3];
        int local_block;

        if (block_count > batch_blocks) block_count = batch_blocks;
        count[0] = (hsize_t)block_count;
        count[1] = NVAR_MATCH;
        count[2] = (hsize_t)ncells;
        file_space = H5Dget_space(data);
        mem_space = H5Screate_simple(3, count, 0);
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, 0, count, 0);
        if (H5Dread(data, H5T_NATIVE_FLOAT, mem_space, file_space,
                H5P_DEFAULT, buffer) < 0) {
            fail("failed to read Data batch");
        }
        H5Sclose(mem_space);
        H5Sclose(file_space);

        for (local_block = 0; local_block < block_count; ++local_block) {
            int block = block0 + local_block;
            const double *meta = &metadata[(size_t)block * META_COLS];
            size_t block_base = (size_t)local_block * NVAR_MATCH * ncells;
            int cell;

            if ((int)meta[3] != 1) continue;
            for (cell = 0; cell < ncells; ++cell) {
                double rho = buffer[block_base + cell];
                double v[3];
                double ye;
                double logrho;
                double score;
                match candidate;
                int d;

                if (!(rho > 0.0)) continue;
                for (d = 0; d < 3; ++d) {
                    v[d] = buffer[block_base + (size_t)(d + 1) * ncells + cell];
                }
                ye = buffer[block_base + (size_t)5 * ncells + cell];
                logrho = log10(rho / target_rho);
                score = (logrho / 0.08) * (logrho / 0.08) +
                    ((v[0] - target_v[0]) / 4.0e8) * ((v[0] - target_v[0]) / 4.0e8) +
                    ((v[1] - target_v[1]) / 4.0e8) * ((v[1] - target_v[1]) / 4.0e8) +
                    ((v[2] - target_v[2]) / 4.0e8) * ((v[2] - target_v[2]) / 4.0e8) +
                    ((ye - target_ye) / 0.025) * ((ye - target_ye) / 0.025);
                candidate.score = score;
                candidate.radiation_score = 0.0;
                candidate.flux_direction_score = 0.0;
                candidate.block = block;
                candidate.cell = cell;
                candidate.rho = rho;
                candidate.v[0] = v[0];
                candidate.v[1] = v[1];
                candidate.v[2] = v[2];
                candidate.ye = ye;
                keep_match(best, &nkeep, &candidate);
            }
        }
    }

    if (have_target_u) {
        float row[153];
        int candidate_idx;

        for (candidate_idx = 0; candidate_idx < nkeep; ++candidate_idx) {
            match *m = &best[candidate_idx];
            hid_t file_space = H5Dget_space(data);
            hsize_t start[3] = {
                (hsize_t)m->block, 0, (hsize_t)m->cell
            };
            hsize_t count[3] = {1, 153, 1};
            hid_t mem_space = H5Screate_simple(3, count, 0);
            double radiation_score = 0.0;
            double flux_score = 0.0;
            double flux_weight = 0.0;
            int field;
            int group;

            H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, 0, count, 0);
            if (H5Dread(data, H5T_NATIVE_FLOAT, mem_space, file_space,
                    H5P_DEFAULT, row) < 0) {
                fail("failed to read candidate radiation state");
            }
            H5Sclose(mem_space);
            H5Sclose(file_space);

            for (field = 0; field < 3; ++field) {
                for (group = 0; group < 12; ++group) {
                    int eidx = 9 + 4 * (field * 12 + group);
                    double candidate_e = (double)row[eidx] / 1.0e25;
                    double target_e = target_u[eidx];
                    double candidate_log = copysign(log1p(fabs(candidate_e)), candidate_e);
                    double target_log = copysign(log1p(fabs(target_e)), target_e);
                    double target_flux[3];
                    double candidate_flux[3];
                    double target_flux2 = 0.0;
                    double candidate_flux2 = 0.0;
                    double dot = 0.0;
                    int component;

                    radiation_score += (candidate_log - target_log) *
                        (candidate_log - target_log);
                    for (component = 0; component < 3; ++component) {
                        target_flux[component] = target_u[eidx + 1 + component];
                        candidate_flux[component] =
                            (double)row[eidx + 1 + component] / 1.0e25;
                        target_flux2 += target_flux[component] * target_flux[component];
                        candidate_flux2 += candidate_flux[component] *
                            candidate_flux[component];
                        dot += target_flux[component] * candidate_flux[component];
                    }
                    if (target_e > 0.0 && candidate_e > 0.0 &&
                        isfinite(target_flux2) && isfinite(candidate_flux2) &&
                        isfinite(dot) &&
                        target_flux2 > 0.0 && candidate_flux2 > 0.0) {
                        double weight = sqrt(target_e);
                        double cosine = dot / sqrt(target_flux2 * candidate_flux2);

                        if (cosine > 1.0) cosine = 1.0;
                        if (cosine < -1.0) cosine = -1.0;
                        flux_score += weight * (1.0 - cosine);
                        flux_weight += weight;
                    }
                }
            }
            m->radiation_score = radiation_score / 36.0;
            m->flux_direction_score = flux_weight > 0.0 ?
                flux_score / flux_weight : 0.0;
            m->score += m->radiation_score + 20.0 * m->flux_direction_score;
        }
    }
    H5Dclose(data);
    H5Fclose(file);

    qsort(best, (size_t)nkeep, sizeof(*best), compare_match);
    printf("rank,score,radiation_score,flux_direction_score,"
           "block,i,j,k,level,dx,x,y,z,radius,rho,v1,v2,v3,speed,vr,Ye\n");
    for (n = 0; n < nkeep; ++n) {
        const match *m = &best[n];
        const double *meta = &metadata[(size_t)m->block * META_COLS];
        int i = m->cell / (block_size * block_size);
        int j = (m->cell / block_size) % block_size;
        int k = m->cell % block_size;
        double xyz[3];
        double radius2 = 0.0;
        double speed2 = 0.0;
        double radial_momentum = 0.0;
        double radius;
        int d;

        for (d = 0; d < 3; ++d) {
            int ijk = d == 0 ? i : (d == 1 ? j : k);
            xyz[d] = meta[4 + d] + ((double)ijk + 0.5) * meta[10 + d] - x_com[d];
            radius2 += xyz[d] * xyz[d];
            speed2 += m->v[d] * m->v[d];
            radial_momentum += m->v[d] * xyz[d];
        }
        radius = sqrt(radius2);
        printf("%d,%.9e,%.9e,%.9e,%d,%d,%d,%d,%d,%.9e,%.9e,%.9e,%.9e,%.9e,"
               "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
            n, m->score, m->radiation_score, m->flux_direction_score,
            m->block, i, j, k, (int)meta[2], meta[10],
            xyz[0], xyz[1], xyz[2], radius, m->rho,
            m->v[0], m->v[1], m->v[2], sqrt(speed2),
            radius > 0.0 ? radial_momentum / radius : 0.0, m->ye);
    }

    free(buffer);
    free(metadata);
    return 0;
}
