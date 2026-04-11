#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#define PRJ_IO_METADATA_SIZE 638

static void prj_io_dataset_name(int var, char *name, size_t size)
{
    if (var == PRJ_PRIM_RHO) {
        snprintf(name, size, "density");
    } else if (var == PRJ_PRIM_V1) {
        snprintf(name, size, "v1");
    } else if (var == PRJ_PRIM_V2) {
        snprintf(name, size, "v2");
    } else if (var == PRJ_PRIM_V3) {
        snprintf(name, size, "v3");
    } else if (var == PRJ_PRIM_EINT) {
        snprintf(name, size, "eint");
    } else if (var == PRJ_PRIM_YE) {
        snprintf(name, size, "ye");
    } else {
        int rad_idx = var - PRJ_NHYDRO;
        int field = rad_idx / (PRJ_NEGROUP * PRJ_RAD_GROUP_STRIDE);
        int local = rad_idx % (PRJ_NEGROUP * PRJ_RAD_GROUP_STRIDE);
        int group = local / PRJ_RAD_GROUP_STRIDE;
        int component = local % PRJ_RAD_GROUP_STRIDE;

        if (component == 0) {
            snprintf(name, size, "rad%d_g%02d_E", field, group);
        } else {
            snprintf(name, size, "rad%d_g%02d_F%d", field, group, component);
        }
    }
}

static void prj_io_fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1);
}

static int prj_io_is_local_owner(const prj_block *block)
{
    prj_mpi *mpi = prj_mpi_current();

    return block != 0 && block->id >= 0 && (mpi == 0 || block->rank == mpi->rank);
}

static hid_t prj_io_create_file(const char *filename)
{
#if defined(PRJ_ENABLE_MPI)
    prj_mpi *mpi = prj_mpi_current();

    if (mpi != 0 && mpi->totrank > 1) {
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        hid_t file;

        H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
        file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        H5Pclose(fapl);
        return file;
    }
#endif
    return H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

static hid_t prj_io_open_file_readonly(const char *filename)
{
#if defined(PRJ_ENABLE_MPI)
    prj_mpi *mpi = prj_mpi_current();

    if (mpi != 0 && mpi->totrank > 1) {
        hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
        hid_t file;

        H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
        file = H5Fopen(filename, H5F_ACC_RDONLY, fapl);
        H5Pclose(fapl);
        return file;
    }
#endif
    return H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
}

static hid_t prj_io_data_xfer_plist(void)
{
#if defined(PRJ_ENABLE_MPI)
    prj_mpi *mpi = prj_mpi_current();

    if (mpi != 0 && mpi->totrank > 1) {
        hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);

        H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
        return dxpl;
    }
#endif
    return H5P_DEFAULT;
}

static void prj_io_close_dxpl(hid_t dxpl)
{
    if (dxpl != H5P_DEFAULT) {
        H5Pclose(dxpl);
    }
}

static void prj_io_write_attr_double(hid_t obj, const char *name, double value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(obj, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attr);
    H5Sclose(space);
}

static void prj_io_write_attr_int(hid_t obj, const char *name, int value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(obj, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr, H5T_NATIVE_INT, &value);
    H5Aclose(attr);
    H5Sclose(space);
}

static void prj_io_write_attr_int3(hid_t obj, const char *name, const int values[3])
{
    hsize_t dims[1] = {3};
    hid_t space = H5Screate_simple(1, dims, dims);
    hid_t attr = H5Acreate2(obj, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr, H5T_NATIVE_INT, values);
    H5Aclose(attr);
    H5Sclose(space);
}

static void prj_io_write_attr_double6(hid_t obj, const char *name, const prj_coord *coord)
{
    double values[6];
    hsize_t dims[1] = {6};
    hid_t space;
    hid_t attr;

    values[0] = coord->x1min;
    values[1] = coord->x1max;
    values[2] = coord->x2min;
    values[3] = coord->x2max;
    values[4] = coord->x3min;
    values[5] = coord->x3max;
    space = H5Screate_simple(1, dims, dims);
    attr = H5Acreate2(obj, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, values);
    H5Aclose(attr);
    H5Sclose(space);
}

static void prj_io_read_attr_double(hid_t obj, const char *name, double *value)
{
    hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);

    H5Aread(attr, H5T_NATIVE_DOUBLE, value);
    H5Aclose(attr);
}

static void prj_io_read_attr_int(hid_t obj, const char *name, int *value)
{
    hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);

    H5Aread(attr, H5T_NATIVE_INT, value);
    H5Aclose(attr);
}

static void prj_io_read_attr_int3(hid_t obj, const char *name, int values[3])
{
    hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);

    H5Aread(attr, H5T_NATIVE_INT, values);
    H5Aclose(attr);
}

static void prj_io_read_attr_double6(hid_t obj, const char *name, prj_coord *coord)
{
    double values[6];
    hid_t attr = H5Aopen(obj, name, H5P_DEFAULT);

    H5Aread(attr, H5T_NATIVE_DOUBLE, values);
    H5Aclose(attr);
    coord->x1min = values[0];
    coord->x1max = values[1];
    coord->x2min = values[2];
    coord->x2max = values[3];
    coord->x3min = values[4];
    coord->x3max = values[5];
}

static void prj_io_fill_metadata(const prj_block *block, double *metadata_row)
{
    int n;
    int idx;

    for (n = 0; n < PRJ_IO_METADATA_SIZE; ++n) {
        metadata_row[n] = 0.0;
    }
    metadata_row[0] = (double)block->id;
    metadata_row[1] = (double)block->rank;
    metadata_row[2] = (double)block->level;
    metadata_row[3] = (double)block->active;
    metadata_row[4] = block->xmin[0];
    metadata_row[5] = block->xmin[1];
    metadata_row[6] = block->xmin[2];
    metadata_row[7] = block->xmax[0];
    metadata_row[8] = block->xmax[1];
    metadata_row[9] = block->xmax[2];
    metadata_row[10] = block->dx[0];
    metadata_row[11] = block->dx[1];
    metadata_row[12] = block->dx[2];
    metadata_row[13] = (double)block->parent;
    for (n = 0; n < 8; ++n) {
        metadata_row[14 + n] = (double)block->children[n];
    }
    idx = 22;
    for (n = 0; n < 56; ++n) {
        int d;

        metadata_row[idx++] = (double)block->slot[n].id;
        metadata_row[idx++] = (double)block->slot[n].rank;
        for (d = 0; d < 3; ++d) {
            metadata_row[idx++] = block->slot[n].xmin[d];
        }
        for (d = 0; d < 3; ++d) {
            metadata_row[idx++] = block->slot[n].xmax[d];
        }
        for (d = 0; d < 3; ++d) {
            metadata_row[idx++] = block->slot[n].dx[d];
        }
    }
}

static void prj_io_unpack_metadata(prj_block *block, const double *metadata_row)
{
    int n;
    int idx;

    block->id = (int)metadata_row[0];
    block->rank = (int)metadata_row[1];
    block->level = (int)metadata_row[2];
    block->active = (int)metadata_row[3];
    block->xmin[0] = metadata_row[4];
    block->xmin[1] = metadata_row[5];
    block->xmin[2] = metadata_row[6];
    block->xmax[0] = metadata_row[7];
    block->xmax[1] = metadata_row[8];
    block->xmax[2] = metadata_row[9];
    block->dx[0] = metadata_row[10];
    block->dx[1] = metadata_row[11];
    block->dx[2] = metadata_row[12];
    block->parent = (int)metadata_row[13];
    for (n = 0; n < 8; ++n) {
        block->children[n] = (int)metadata_row[14 + n];
    }
    idx = 22;
    for (n = 0; n < 56; ++n) {
        int d;

        block->slot[n].id = (int)metadata_row[idx++];
        block->slot[n].rank = (int)metadata_row[idx++];
        for (d = 0; d < 3; ++d) {
            block->slot[n].xmin[d] = metadata_row[idx++];
        }
        for (d = 0; d < 3; ++d) {
            block->slot[n].xmax[d] = metadata_row[idx++];
        }
        for (d = 0; d < 3; ++d) {
            block->slot[n].dx[d] = metadata_row[idx++];
        }
    }
}

void prj_io_write_restart(const prj_mesh *mesh, double time, int step)
{
    char filename[64];
    hid_t file;
    hid_t space_data;
    hid_t space_meta;
    hid_t dset_data;
    hid_t dset_meta;
    hsize_t dims_data[3];
    hsize_t dims_meta[2];
    int bidx;

    snprintf(filename, sizeof(filename), "output/restart_%08d.h5", step);
    dims_data[0] = (hsize_t)mesh->nblocks;
    dims_data[1] = (hsize_t)PRJ_NVAR_PRIM;
    dims_data[2] = (hsize_t)(PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE);
    dims_meta[0] = (hsize_t)mesh->nblocks;
    dims_meta[1] = (hsize_t)PRJ_IO_METADATA_SIZE;
    file = prj_io_create_file(filename);
    prj_io_write_attr_double(file, "time", time);
    prj_io_write_attr_int(file, "step", step);
    prj_io_write_attr_int(file, "nblocks", mesh->nblocks);
    prj_io_write_attr_int(file, "nvar_prim", PRJ_NVAR_PRIM);
    prj_io_write_attr_int(file, "block_size", PRJ_BLOCK_SIZE);
    prj_io_write_attr_int(file, "max_level", mesh->max_level);
    prj_io_write_attr_int3(file, "root_nx", mesh->root_nx);
    prj_io_write_attr_double6(file, "coord", &mesh->coord);

    space_data = H5Screate_simple(3, dims_data, dims_data);
    dset_data = H5Dcreate2(file, "Data", H5T_NATIVE_DOUBLE, space_data, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    space_meta = H5Screate_simple(2, dims_meta, dims_meta);
    dset_meta = H5Dcreate2(file, "MetaData", H5T_NATIVE_DOUBLE, space_meta, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];

        if (block->id < 0 || !prj_io_is_local_owner(block)) {
            continue;
        }
        {
            hsize_t start_meta[2] = {(hsize_t)bidx, 0};
            hsize_t count_meta[2] = {1, (hsize_t)PRJ_IO_METADATA_SIZE};
            hid_t mem_meta = H5Screate_simple(2, count_meta, count_meta);
            hid_t file_meta = H5Dget_space(dset_meta);
            double metadata_row[PRJ_IO_METADATA_SIZE];
            hid_t dxpl = prj_io_data_xfer_plist();

            prj_io_fill_metadata(block, metadata_row);
            H5Sselect_hyperslab(file_meta, H5S_SELECT_SET, start_meta, 0, count_meta, 0);
            H5Dwrite(dset_meta, H5T_NATIVE_DOUBLE, mem_meta, file_meta, dxpl, metadata_row);
            prj_io_close_dxpl(dxpl);
            H5Sclose(file_meta);
            H5Sclose(mem_meta);
        }
        if (block->active == 1 && block->W != 0) {
            hsize_t start_data[3] = {(hsize_t)bidx, 0, 0};
            hsize_t count_data[3] = {1, (hsize_t)PRJ_NVAR_PRIM, (hsize_t)(PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE)};
            hid_t mem_data = H5Screate_simple(3, count_data, count_data);
            hid_t file_data = H5Dget_space(dset_data);
            double *buffer = (double *)calloc((size_t)PRJ_NVAR_PRIM * (size_t)(PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE), sizeof(*buffer));
            hid_t dxpl = prj_io_data_xfer_plist();
            int v;
            int i;
            int j;
            int k;

            if (buffer == 0) {
                prj_io_fail("prj_io_write_restart: allocation failed");
            }
            for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                    for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                        for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                            size_t cell = (size_t)i * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE + (size_t)j * PRJ_BLOCK_SIZE + (size_t)k;
                            buffer[(size_t)v * count_data[2] + cell] = block->W[VIDX(v, i, j, k)];
                        }
                    }
                }
            }
            H5Sselect_hyperslab(file_data, H5S_SELECT_SET, start_data, 0, count_data, 0);
            H5Dwrite(dset_data, H5T_NATIVE_DOUBLE, mem_data, file_data, dxpl, buffer);
            prj_io_close_dxpl(dxpl);
            free(buffer);
            H5Sclose(file_data);
            H5Sclose(mem_data);
        }
    }
    H5Dclose(dset_data);
    H5Sclose(space_data);
    H5Dclose(dset_meta);
    H5Sclose(space_meta);
    H5Fclose(file);
}

void prj_io_read_restart(prj_mesh *mesh, const prj_eos *eos, const char *filename, double *time, int *step)
{
    hid_t file;
    hid_t dset_data;
    hid_t dset_meta;
    int nblocks;
    int block_size;
    int nvar_prim;
    int root_nx[3];
    int max_level;
    prj_coord coord;
    double *metadata;
    int bidx;
    prj_bc bc = {
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW,
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW,
        PRJ_BC_OUTFLOW, PRJ_BC_OUTFLOW
    };

    file = prj_io_open_file_readonly(filename);
    prj_io_read_attr_double(file, "time", time);
    prj_io_read_attr_int(file, "step", step);
    prj_io_read_attr_int(file, "nblocks", &nblocks);
    prj_io_read_attr_int(file, "nvar_prim", &nvar_prim);
    prj_io_read_attr_int(file, "block_size", &block_size);
    prj_io_read_attr_int(file, "max_level", &max_level);
    prj_io_read_attr_int3(file, "root_nx", root_nx);
    prj_io_read_attr_double6(file, "coord", &coord);
    if (nvar_prim != PRJ_NVAR_PRIM || block_size != PRJ_BLOCK_SIZE) {
        prj_io_fail("prj_io_read_restart: incompatible restart metadata");
    }

    if (prj_mesh_init(mesh, root_nx[0], root_nx[1], root_nx[2], max_level, &coord) != 0) {
        prj_io_fail("prj_io_read_restart: mesh init failed");
    }
    mesh->nblocks = nblocks;

    metadata = (double *)calloc((size_t)nblocks * PRJ_IO_METADATA_SIZE, sizeof(*metadata));
    if (metadata == 0) {
        prj_io_fail("prj_io_read_restart: allocation failed");
    }

    dset_meta = H5Dopen2(file, "MetaData", H5P_DEFAULT);
    H5Dread(dset_meta, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata);
    H5Dclose(dset_meta);

    for (bidx = 0; bidx < nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        const double *meta_row = &metadata[(size_t)bidx * PRJ_IO_METADATA_SIZE];

        prj_io_unpack_metadata(block, meta_row);
        if (block->id >= 0 && prj_io_is_local_owner(block) && block->W == 0) {
            prj_block_alloc_data(block);
        }
        prj_block_setup_geometry(block, &coord);
        if (block->id >= 0 && block->active == 1 && prj_io_is_local_owner(block) && block->W != 0) {
            hsize_t start_data[3] = {(hsize_t)bidx, 0, 0};
            hsize_t count_data[3] = {1, (hsize_t)PRJ_NVAR_PRIM, (hsize_t)(PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE)};
            hid_t mem_data = H5Screate_simple(3, count_data, count_data);
            hid_t file_data;
            double *buffer = (double *)calloc((size_t)PRJ_NVAR_PRIM * (size_t)(PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE), sizeof(*buffer));
            hid_t dxpl = prj_io_data_xfer_plist();
            int v;
            int i;
            int j;
            int k;

            if (buffer == 0) {
                prj_io_fail("prj_io_read_restart: allocation failed");
            }
            dset_data = H5Dopen2(file, "Data", H5P_DEFAULT);
            file_data = H5Dget_space(dset_data);
            H5Sselect_hyperslab(file_data, H5S_SELECT_SET, start_data, 0, count_data, 0);
            H5Dread(dset_data, H5T_NATIVE_DOUBLE, mem_data, file_data, dxpl, buffer);
            prj_io_close_dxpl(dxpl);
            H5Sclose(file_data);
            H5Dclose(dset_data);
            for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                    for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                        for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                            size_t cell = (size_t)i * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE + (size_t)j * PRJ_BLOCK_SIZE + (size_t)k;
                            block->W[VIDX(v, i, j, k)] = buffer[(size_t)v * count_data[2] + cell];
                        }
                    }
                }
            }
            free(buffer);
            H5Sclose(mem_data);
        }
    }
    H5Fclose(file);
    prj_boundary_fill_ghosts(mesh, &bc, 1);
    (void)eos;
    free(metadata);
}

void prj_io_write_dump(const prj_mesh *mesh, const char *basename, int step)
{
    char filename[256];
    hid_t file;
    hid_t dset_level;
    hid_t dset_coord;
    hid_t dset_var[PRJ_NVAR_PRIM];
    hid_t space_level;
    hid_t space_coord;
    hid_t space_var;
    hsize_t dims_level[1];
    hsize_t dims_coord[2];
    hsize_t dims_var[4];
    int local_active = 0;
    int nactive = 0;
    int bidx;
    int active_idx = 0;

    snprintf(filename, sizeof(filename), "%s_%05d.h5", basename, step);
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        if (mesh->blocks[bidx].id >= 0 && mesh->blocks[bidx].active == 1 && prj_io_is_local_owner(&mesh->blocks[bidx])) {
            local_active += 1;
        }
    }
    nactive = (int)prj_mpi_global_sum((double)local_active);
    dims_level[0] = (hsize_t)nactive;
    dims_coord[0] = (hsize_t)nactive;
    dims_coord[1] = 3;
    dims_var[0] = (hsize_t)nactive;
    dims_var[1] = PRJ_BLOCK_SIZE;
    dims_var[2] = PRJ_BLOCK_SIZE;
    dims_var[3] = PRJ_BLOCK_SIZE;
    file = prj_io_create_file(filename);
    prj_io_write_attr_int3(file, "root_nx", mesh->root_nx);
    prj_io_write_attr_double6(file, "coord", &mesh->coord);
    space_level = H5Screate_simple(1, dims_level, dims_level);
    dset_level = H5Dcreate2(file, "level", H5T_NATIVE_INT, space_level, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    space_coord = H5Screate_simple(2, dims_coord, dims_coord);
    dset_coord = H5Dcreate2(file, "coordinate", H5T_NATIVE_DOUBLE, space_coord, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    space_var = H5Screate_simple(4, dims_var, dims_var);
    for (bidx = 0; bidx < PRJ_NVAR_PRIM; ++bidx) {
        char name[32];

        prj_io_dataset_name(bidx, name, sizeof(name));
        dset_var[bidx] = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, space_var, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    {
        int offset = 0;
#if defined(PRJ_ENABLE_MPI)
        prj_mpi *mpi = prj_mpi_current();

        if (mpi != 0 && mpi->totrank > 1) {
            MPI_Exscan(&local_active, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (mpi->rank == 0) {
                offset = 0;
            }
        }
#endif
        active_idx = offset;
    }
    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        const prj_block *block = &mesh->blocks[bidx];
        int var;

        if (block->id < 0 || block->active != 1 || !prj_io_is_local_owner(block) || block->W == 0) {
            continue;
        }
        {
            hsize_t start_level[1] = {(hsize_t)active_idx};
            hsize_t count_level[1] = {1};
            hid_t mem_level = H5Screate_simple(1, count_level, count_level);
            hid_t file_level = H5Dget_space(dset_level);
            int level_value = block->level;
            hid_t dxpl = prj_io_data_xfer_plist();

            H5Sselect_hyperslab(file_level, H5S_SELECT_SET, start_level, 0, count_level, 0);
            H5Dwrite(dset_level, H5T_NATIVE_INT, mem_level, file_level, dxpl, &level_value);
            prj_io_close_dxpl(dxpl);
            H5Sclose(file_level);
            H5Sclose(mem_level);
        }
        {
            hsize_t start_coord[2] = {(hsize_t)active_idx, 0};
            hsize_t count_coord[2] = {1, 3};
            hid_t mem_coord = H5Screate_simple(2, count_coord, count_coord);
            hid_t file_coord = H5Dget_space(dset_coord);
            double coord_value[3];
            hid_t dxpl = prj_io_data_xfer_plist();

            coord_value[0] = block->xmin[0];
            coord_value[1] = block->xmin[1];
            coord_value[2] = block->xmin[2];
            H5Sselect_hyperslab(file_coord, H5S_SELECT_SET, start_coord, 0, count_coord, 0);
            H5Dwrite(dset_coord, H5T_NATIVE_DOUBLE, mem_coord, file_coord, dxpl, coord_value);
            prj_io_close_dxpl(dxpl);
            H5Sclose(file_coord);
            H5Sclose(mem_coord);
        }
        for (var = 0; var < PRJ_NVAR_PRIM; ++var) {
            hsize_t start_var[4] = {(hsize_t)active_idx, 0, 0, 0};
            hsize_t count_var[4] = {1, PRJ_BLOCK_SIZE, PRJ_BLOCK_SIZE, PRJ_BLOCK_SIZE};
            hid_t mem_var = H5Screate_simple(4, count_var, count_var);
            hid_t file_var = H5Dget_space(dset_var[var]);
            double *buffer = (double *)calloc((size_t)(PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE), sizeof(*buffer));
            hid_t dxpl = prj_io_data_xfer_plist();
            int i;
            int j;
            int k;

            if (buffer == 0) {
                prj_io_fail("prj_io_write_dump: allocation failed");
            }
            for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
                for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                    for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                        size_t cell = (size_t)i * PRJ_BLOCK_SIZE * PRJ_BLOCK_SIZE + (size_t)j * PRJ_BLOCK_SIZE + (size_t)k;
                        buffer[cell] = block->W[VIDX(var, i, j, k)];
                    }
                }
            }
            H5Sselect_hyperslab(file_var, H5S_SELECT_SET, start_var, 0, count_var, 0);
            H5Dwrite(dset_var[var], H5T_NATIVE_DOUBLE, mem_var, file_var, dxpl, buffer);
            prj_io_close_dxpl(dxpl);
            free(buffer);
            H5Sclose(file_var);
            H5Sclose(mem_var);
        }
        active_idx += 1;
    }
    H5Dclose(dset_level);
    H5Sclose(space_level);
    H5Dclose(dset_coord);
    H5Sclose(space_coord);
    for (bidx = 0; bidx < PRJ_NVAR_PRIM; ++bidx) {
        H5Dclose(dset_var[bidx]);
    }
    H5Sclose(space_var);
    H5Fclose(file);
}
