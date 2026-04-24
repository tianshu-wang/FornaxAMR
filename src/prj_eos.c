#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif

#include "prj.h"

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#define PRJ_EOS_NUMEL 16
#define PRJ_EOS_EXT_NT 500
#define PRJ_EOS_EXT_NR 500
#define PRJ_EOS_EXT_NY 50
#define PRJ_EOS_EXT_R1 (-8.699)
#define PRJ_EOS_EXT_R2 (15.5)
#define PRJ_EOS_EXT_T1 (-4.793)
#define PRJ_EOS_EXT_T2 (2.176)
#define PRJ_EOS_EXT_Y1 (0.035)
#define PRJ_EOS_EXT_Y2 (0.56)
/* 1 MeV per baryon in cgs specific-energy units (erg g^-1). */
#define PRJ_EOS_ENERGY_SCALE 0.95655684e18
#define PRJ_EOS_PRESSURE_SCALE 1.60217733e33

enum {
    PRJ_EOS_REC_EINT = 1,
    /* Stored in-memory as log(pressure) after table initialization. */
    PRJ_EOS_REC_PRESSURE = 2,
    PRJ_EOS_REC_GAMMA = 12
};

static double prj_eos_gamma_value(void)
{
    return 5.0 / 3.0;
}

static double prj_eos_exp10(double x)
{
    return exp(M_LN10 * x);
}

static size_t prj_eos_tab_index(const prj_eos *eos, int rec, int jy, int jr, int jt)
{
    return (size_t)(rec - 1) * (size_t)eos->ny * (size_t)eos->nr * (size_t)eos->nt +
        (size_t)(jy - 1) * (size_t)eos->nr * (size_t)eos->nt +
        (size_t)(jr - 1) * (size_t)eos->nt +
        (size_t)(jt - 1);
}

static double prj_eos_tab_elem(const prj_eos *eos, int rec, int jy, int jr, int jt)
{
    return eos->table[prj_eos_tab_index(eos, rec, jy, jr, jt)];
}

static double prj_eos_clamp_double(double x, double lo, double hi)
{
    if (x < lo) {
        return lo;
    }
    if (x > hi) {
        return hi;
    }
    return x;
}

static void prj_eos_table_range_fail(const char *name, double value, double lo, double hi)
{
    prj_mpi *mpi = prj_mpi_current();

    fprintf(stderr,
        "tabulated EOS lookup out of range for %s: value=%.17e allowed=[%.17e, %.17e]\n",
        name, value, lo, hi);
#if defined(PRJ_ENABLE_MPI)
    if (mpi != 0 && mpi->totrank > 1) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
#endif
    exit(1);
}

static void prj_eos_convert_pressure_slab_to_log(prj_eos *eos, size_t slab_size)
{
    size_t offset;
    size_t idx;

    if (eos == 0 || eos->table == 0 || slab_size == 0) {
        return;
    }
    offset = (size_t)(PRJ_EOS_REC_PRESSURE - 1) * slab_size;
    for (idx = 0; idx < slab_size; ++idx) {
        double p = eos->table[offset + idx];

        eos->table[offset + idx] = log(p > 0.0 ? p : 1.0e-300);
    }
}

static void prj_eos_print_fill_neighbors(const prj_block *block, double x1, double x2, double x3)
{
    int n;
    const double tol = 1.0e-12;
    int found = 0;

    if (block == 0) {
        return;
    }

    fprintf(stderr, "candidate fill neighbors for x=(%.17e, %.17e, %.17e):\n", x1, x2, x3);
    for (n = 0; n < 56; ++n) {
        const prj_neighbor *slot = &block->slot[n];

        if (slot->id < 0) {
            continue;
        }
        if (x1 >= slot->xmin[0] - tol && x1 <= slot->xmax[0] + tol &&
            x2 >= slot->xmin[1] - tol && x2 <= slot->xmax[1] + tol &&
            x3 >= slot->xmin[2] - tol && x3 <= slot->xmax[2] + tol) {
            fprintf(stderr,
                "  slot=%d neighbor_id=%d current_rank=%d neighbor_rank=%d "
                "xmin=(%.17e, %.17e, %.17e) xmax=(%.17e, %.17e, %.17e) "
                "dx=(%.17e, %.17e, %.17e)\n",
                n, slot->id, block->rank, slot->rank,
                slot->xmin[0], slot->xmin[1], slot->xmin[2],
                slot->xmax[0], slot->xmax[1], slot->xmax[2],
                slot->dx[0], slot->dx[1], slot->dx[2]);
            found = 1;
        }
    }
    if (found == 0) {
        fprintf(stderr, "  no neighbor slot contains this cell center\n");
    }
}

static void prj_eos_fail_zero_rho_in_block(const char *caller, const prj_block *block,
    int i, int j, int k, double *W)
{
    double x1;
    double x2;
    double x3;

    if (block == 0 || W == 0) {
        fprintf(stderr, "%s: rho is zero and block/state context is unavailable\n", caller);
        exit(1);
    }

    x1 = block->xmin[0] + ((double)i + 0.5) * block->dx[0];
    x2 = block->xmin[1] + ((double)j + 0.5) * block->dx[1];
    x3 = block->xmin[2] + ((double)k + 0.5) * block->dx[2];
    fprintf(stderr,
        "%s: rho=0 before EOS call for current_block id=%d current_rank=%d level=%d cell=(%d,%d,%d) "
        "x=(%.17e, %.17e, %.17e)\n",
        caller, block->id, block->rank, block->level, i, j, k, x1, x2, x3);
    fprintf(stderr,
        "  block xmin=(%.17e, %.17e, %.17e) xmax=(%.17e, %.17e, %.17e) "
        "dx=(%.17e, %.17e, %.17e)\n",
        block->xmin[0], block->xmin[1], block->xmin[2],
        block->xmax[0], block->xmax[1], block->xmax[2],
        block->dx[0], block->dx[1], block->dx[2]);
    fprintf(stderr,
        "  primitive state: rho=%.17e v=(%.17e, %.17e, %.17e) eint=%.17e ye=%.17e\n",
        W[VIDX(PRJ_PRIM_RHO, i, j, k)],
        W[VIDX(PRJ_PRIM_V1, i, j, k)],
        W[VIDX(PRJ_PRIM_V2, i, j, k)],
        W[VIDX(PRJ_PRIM_V3, i, j, k)],
        W[VIDX(PRJ_PRIM_EINT, i, j, k)],
        W[VIDX(PRJ_PRIM_YE, i, j, k)]);
    prj_eos_print_fill_neighbors(block, x1, x2, x3);
    exit(1);
}

static void prj_eos_fill_cell(prj_eos *eos, prj_block *block, double *W, int i, int j, int k)
{
    double eos_quantities[PRJ_EOS_NQUANT];

    if (block == 0 || W == 0 || block->eosvar == 0) {
        return;
    }
    if (W[VIDX(PRJ_PRIM_RHO, i, j, k)] == 0.0) {
        prj_eos_fail_zero_rho_in_block("prj_eos_fill_block", block, i, j, k, W);
    }

    prj_eos_rey(eos,
        W[VIDX(PRJ_PRIM_RHO, i, j, k)],
        W[VIDX(PRJ_PRIM_EINT, i, j, k)],
        W[VIDX(PRJ_PRIM_YE, i, j, k)],
        eos_quantities);
    block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)] = eos_quantities[PRJ_EOS_PRESSURE];
    block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)] = eos_quantities[PRJ_EOS_TEMPERATURE];
    block->eosvar[EIDX(PRJ_EOSVAR_GAMMA, i, j, k)] = eos_quantities[PRJ_EOS_GAMMA];
    if (block->eos_done != 0) {
        block->eos_done[IDX(i, j, k)] = 1;
    }
}

static void prj_eos_table_check_rty_inputs(const prj_eos *eos, double rho, double T, double ye)
{
    if (rho < eos->rho_min || rho > eos->rho_max) {
        prj_eos_table_range_fail("rho", rho, eos->rho_min, eos->rho_max);
    }
    if (T < eos->temp_min || T > eos->temp_max) {
        prj_eos_table_range_fail("T", T, eos->temp_min, eos->temp_max);
    }
    if (ye < eos->y1c || ye > eos->y2c) {
        prj_eos_table_range_fail("ye", ye, eos->y1c, eos->y2c);
    }
}

static int prj_eos_prepare_table(prj_eos *eos)
{
    prj_mpi *mpi;
    FILE *f;
    size_t slab_size;
    size_t expected_bytes;
    int rec;
    int status = 0;

    if (eos == 0 || eos->kind != PRJ_EOS_KIND_TABLE || eos->filename[0] == '\0') {
        return 0;
    }
    if (eos->table_loaded == 1 && eos->table != 0) {
        return 0;
    }

    eos->nt = PRJ_EOS_EXT_NT;
    eos->nr = PRJ_EOS_EXT_NR;
    eos->ny = PRJ_EOS_EXT_NY;
    eos->r1 = PRJ_EOS_EXT_R1;
    eos->r2 = PRJ_EOS_EXT_R2;
    eos->t1 = PRJ_EOS_EXT_T1;
    eos->t2 = PRJ_EOS_EXT_T2;
    eos->y1c = PRJ_EOS_EXT_Y1;
    eos->y2c = PRJ_EOS_EXT_Y2;
    eos->dlogrho = (eos->r2 - eos->r1) / (double)(eos->nr - 1);
    eos->dlogT = (eos->t2 - eos->t1) / (double)(eos->nt - 1);
    eos->dYe = (eos->y2c - eos->y1c) / (double)(eos->ny - 1);
    eos->rho_min = prj_eos_exp10(eos->r1);
    eos->rho_max = prj_eos_exp10(eos->r2);
    eos->temp_min = prj_eos_exp10(eos->t1);
    eos->temp_max = prj_eos_exp10(eos->t2);
    eos->inv_dlogrho = eos->dlogrho > 0.0 ? 1.0 / eos->dlogrho : 0.0;
    eos->inv_dlogT = eos->dlogT > 0.0 ? 1.0 / eos->dlogT : 0.0;
    eos->inv_dYe = eos->dYe > 0.0 ? 1.0 / eos->dYe : 0.0;
    eos->ln10_t1 = M_LN10 * eos->t1;
    eos->ln10_dlogT = M_LN10 * eos->dlogT;
    slab_size = (size_t)eos->nt * (size_t)eos->nr * (size_t)eos->ny;
    expected_bytes = (size_t)PRJ_EOS_NUMEL * slab_size * sizeof(*eos->table);
    mpi = prj_mpi_current();

#if defined(PRJ_ENABLE_MPI)
    if (mpi != 0 && mpi->totrank > 1) {
        if (mpi->rank == 0) {
            eos->table = (double *)malloc(expected_bytes);
            if (eos->table == 0) {
                status = 1;
            } else {
                f = fopen(eos->filename, "rb");
                if (f == 0) {
                    status = 1;
                } else {
                    for (rec = 0; rec < PRJ_EOS_NUMEL; ++rec) {
                        if (fread(&eos->table[(size_t)rec * slab_size], sizeof(double), slab_size, f) != slab_size) {
                            status = 1;
                            break;
                        }
                    }
                    if (status == 0) {
                        prj_eos_convert_pressure_slab_to_log(eos, slab_size);
                    }
                    fclose(f);
                }
            }
        }

        MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (status != 0) {
            if (mpi->rank == 0) {
                free(eos->table);
                eos->table = 0;
            }
            return 1;
        }

        if (mpi->rank != 0) {
            eos->table = (double *)malloc(expected_bytes);
            if (eos->table == 0) {
                status = 1;
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &status, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (status != 0) {
            free(eos->table);
            eos->table = 0;
            return 1;
        }

        MPI_Bcast(eos->table, (int)(expected_bytes / sizeof(*eos->table)), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        eos->table_loaded = 1;
        eos->table_is_mmap = 0;
        eos->table_bytes = expected_bytes;
        return 0;
    }
#endif

    eos->table = (double *)malloc(expected_bytes);
    if (eos->table == 0) {
        return 1;
    }
    f = fopen(eos->filename, "rb");
    if (f == 0) {
        free(eos->table);
        eos->table = 0;
        return 1;
    }
    for (rec = 0; rec < PRJ_EOS_NUMEL; ++rec) {
        if (fread(&eos->table[(size_t)rec * slab_size], sizeof(double), slab_size, f) != slab_size) {
            fclose(f);
            free(eos->table);
            eos->table = 0;
            return 1;
        }
    }
    fclose(f);
    prj_eos_convert_pressure_slab_to_log(eos, slab_size);

    eos->table_loaded = 1;
    eos->table_is_mmap = 0;
    eos->table_bytes = expected_bytes;
    return 0;
}

void prj_eos_init(prj_eos *eos)
{
    prj_mpi *mpi;

    if (eos == 0) {
        return;
    }
    if (eos->kind != PRJ_EOS_KIND_TABLE || eos->filename[0] == '\0') {
        eos->table_loaded = 0;
        eos->table_is_mmap = 0;
        eos->nt = 0;
        eos->nr = 0;
        eos->ny = 0;
        eos->r1 = 0.0;
        eos->r2 = 0.0;
        eos->t1 = 0.0;
        eos->t2 = 0.0;
        eos->y1c = 0.0;
        eos->y2c = 0.0;
        eos->dlogrho = 0.0;
        eos->dlogT = 0.0;
        eos->dYe = 0.0;
        eos->rho_min = 0.0;
        eos->rho_max = 0.0;
        eos->temp_min = 0.0;
        eos->temp_max = 0.0;
        eos->inv_dlogrho = 0.0;
        eos->inv_dlogT = 0.0;
        eos->inv_dYe = 0.0;
        eos->ln10_t1 = 0.0;
        eos->ln10_dlogT = 0.0;
        eos->table_bytes = 0;
        eos->table = 0;
        return;
    }
    (void)prj_eos_prepare_table(eos);
    mpi = prj_mpi_current();
    if (eos->table_loaded == 1 && (mpi == 0 || mpi->rank == 0)) {
        fprintf(stderr, "EOS table initialization finished (%zu bytes loaded)\n", eos->table_bytes);
    }
}

static void prj_eos_table_interp_base(const prj_eos *eos, double rho, double T, double ye,
    int *restrict jy, int *restrict jyp, int *restrict jr, int *restrict jrp, int *restrict jt, int *restrict jtp,
    double *restrict dye, double *restrict drho, double *restrict dtemp)
{
    double rl;
    double tl;
    
    prj_eos_table_check_rty_inputs(eos, rho, T, ye);
    rl = log10(rho);
    tl = log10(T);

    *jr = 1 + (int)((rl - eos->r1) * eos->inv_dlogrho);
    *jt = 1 + (int)((tl - eos->t1) * eos->inv_dlogT);
    *jy = 1 + (int)((ye - eos->y1c) * eos->inv_dYe);
    *jr = *jr < 1 ? 1 : (*jr >= eos->nr ? eos->nr - 1 : *jr);
    *jt = *jt < 1 ? 1 : (*jt >= eos->nt ? eos->nt - 1 : *jt);
    *jy = *jy < 1 ? 1 : (*jy >= eos->ny ? eos->ny - 1 : *jy);
    *jrp = *jr + 1;
    *jtp = *jt + 1;
    *jyp = *jy + 1;
    *drho = (rl - (eos->r1 + (double)(*jr - 1) * eos->dlogrho)) * eos->inv_dlogrho;
    *dtemp = (tl - (eos->t1 + (double)(*jt - 1) * eos->dlogT)) * eos->inv_dlogT;
    *dye = (ye - (eos->y1c + (double)(*jy - 1) * eos->dYe)) * eos->inv_dYe;
    *drho = prj_eos_clamp_double(*drho, 0.0, 1.0);
    *dtemp = prj_eos_clamp_double(*dtemp, 0.0, 1.0);
    *dye = prj_eos_clamp_double(*dye, 0.0, 1.0);
}

static double prj_eos_table_interp_trilinear(const prj_eos *eos, int rec,
    int jy, int jyp, int jr, int jrp, int jt, int jtp, double dye, double drho, double dtemp)
{
    double c000 = prj_eos_tab_elem(eos, rec, jy, jr, jt);
    double c100 = prj_eos_tab_elem(eos, rec, jyp, jr, jt);
    double c010 = prj_eos_tab_elem(eos, rec, jy, jrp, jt);
    double c110 = prj_eos_tab_elem(eos, rec, jyp, jrp, jt);
    double c001 = prj_eos_tab_elem(eos, rec, jy, jr, jtp);
    double c101 = prj_eos_tab_elem(eos, rec, jyp, jr, jtp);
    double c011 = prj_eos_tab_elem(eos, rec, jy, jrp, jtp);
    double c111 = prj_eos_tab_elem(eos, rec, jyp, jrp, jtp);
    double c00 = (1.0 - dye) * c000 + dye * c100;
    double c10 = (1.0 - dye) * c010 + dye * c110;
    double c01 = (1.0 - dye) * c001 + dye * c101;
    double c11 = (1.0 - dye) * c011 + dye * c111;
    double c0 = (1.0 - drho) * c00 + drho * c10;
    double c1 = (1.0 - drho) * c01 + drho * c11;

    return (1.0 - dtemp) * c0 + dtemp * c1;
}

static double prj_eos_pressure_log_interp(const prj_eos *eos,
    int jy, int jyp, int jr, int jrp, int jt, int jtp, double dye, double drho, double dtemp)
{
    return prj_eos_table_interp_trilinear(eos, PRJ_EOS_REC_PRESSURE,
        jy, jyp, jr, jrp, jt, jtp, dye, drho, dtemp);
}

static double prj_eos_rey_slice_eint(const prj_eos *eos, int jy, int jyp, int jr, int jrp, int jt,
    double coeff0, double coeff1, double coeff2, double coeff3)
{
    return coeff0 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy, jr, jt) +
        coeff1 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jr, jt) +
        coeff2 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy, jrp, jt) +
        coeff3 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jrp, jt);
}

void prj_eos_rty(prj_eos *eos, double rho, double T, double ye, double *eos_quantities)
{
    double gamma;
    double eint;

    if (eos_quantities == 0) {
        return;
    }
    if (eos != 0 && eos->kind == PRJ_EOS_KIND_TABLE &&
        eos->filename[0] != '\0' && prj_eos_prepare_table(eos) == 0 && eos->table_loaded == 1) {
        int jy;
        int jyp;
        int jr;
        int jrp;
        int jt;
        int jtp;
        double dye;
        double drho;
        double dtemp;
        double pressure_log;

        prj_eos_table_interp_base(eos, rho, T, ye, &jy, &jyp, &jr, &jrp, &jt, &jtp, &dye, &drho, &dtemp);
        eint = prj_eos_table_interp_trilinear(eos, PRJ_EOS_REC_EINT,
            jy, jyp, jr, jrp, jt, jtp, dye, drho, dtemp) * PRJ_EOS_ENERGY_SCALE;
        pressure_log = prj_eos_pressure_log_interp(eos,
            jy, jyp, jr, jrp, jt, jtp, dye, drho, dtemp);
        gamma = prj_eos_table_interp_trilinear(eos, PRJ_EOS_REC_GAMMA,
            jy, jyp, jr, jrp, jt, jtp, dye, drho, dtemp);

        eos_quantities[PRJ_EOS_EINT] = eint;
        eos_quantities[PRJ_EOS_PRESSURE] = exp(pressure_log) * PRJ_EOS_PRESSURE_SCALE;
        eos_quantities[PRJ_EOS_GAMMA] = gamma;
        eos_quantities[PRJ_EOS_TEMPERATURE] = T;
        return;
    }

    gamma = prj_eos_gamma_value();
    eint = PRJ_EOS_ENERGY_SCALE * T / (gamma - 1.0);
    eos_quantities[PRJ_EOS_EINT] = eint;
    eos_quantities[PRJ_EOS_PRESSURE] = rho * PRJ_EOS_ENERGY_SCALE * T;
    eos_quantities[PRJ_EOS_GAMMA] = gamma;
    eos_quantities[PRJ_EOS_TEMPERATURE] = T;
}

void prj_eos_rey(prj_eos *eos, double rho, double eint, double ye, double *eos_quantities)
{
    double gamma;
    double T;

    if (eos_quantities == 0) {
        return;
    }
    if (eos != 0 && eos->kind == PRJ_EOS_KIND_TABLE &&
        eos->filename[0] != '\0' && prj_eos_prepare_table(eos) == 0 && eos->table_loaded == 1) {
        double e_table;
        double rl;
        double jyf;
        double jrf;
        int jy;
        int jyp;
        int jr;
        int jrp;
        double dye;
        double drho;
        double coeff0;
        double coeff1;
        double coeff2;
        double coeff3;
        int jt;
        int jtp;
        double elo;
        double ehi;
        double dtemp;
        double pressure_log;
        int jt_lo;
        int jt_hi;
        int jt_mid;

        e_table = eint / PRJ_EOS_ENERGY_SCALE;
        if (rho < eos->rho_min || rho > eos->rho_max) {
            prj_eos_table_range_fail("rho", rho, eos->rho_min, eos->rho_max);
        }
        if (ye < eos->y1c || ye > eos->y2c) {
            prj_eos_table_range_fail("ye", ye, eos->y1c, eos->y2c);
        }
        rl = log10(rho);
        jrf = (rl - eos->r1) * eos->inv_dlogrho;
        jyf = (ye - eos->y1c) * eos->inv_dYe;
        jr = 1 + (int)jrf;
        jy = 1 + (int)jyf;
        jr = jr < 1 ? 1 : (jr >= eos->nr ? eos->nr - 1 : jr);
        jy = jy < 1 ? 1 : (jy >= eos->ny ? eos->ny - 1 : jy);
        jrp = jr + 1;
        jyp = jy + 1;
        drho = prj_eos_clamp_double((rl - (eos->r1 + (double)(jr - 1) * eos->dlogrho)) * eos->inv_dlogrho, 0.0, 1.0);
        dye = prj_eos_clamp_double((ye - (eos->y1c + (double)(jy - 1) * eos->dYe)) * eos->inv_dYe, 0.0, 1.0);
        coeff0 = (1.0 - drho) * (1.0 - dye);
        coeff1 = (1.0 - drho) * dye;
        coeff2 = drho * (1.0 - dye);
        coeff3 = drho * dye;

        elo = prj_eos_rey_slice_eint(eos, jy, jyp, jr, jrp, 1, coeff0, coeff1, coeff2, coeff3);
        ehi = prj_eos_rey_slice_eint(eos, jy, jyp, jr, jrp, eos->nt, coeff0, coeff1, coeff2, coeff3);
        if (e_table < elo || e_table > ehi) {
            prj_eos_table_range_fail("eint/PRJ_EOS_ENERGY_SCALE", e_table, elo, ehi);
        }
        if (e_table == elo) {
            jt = 1;
            ehi = prj_eos_rey_slice_eint(eos, jy, jyp, jr, jrp, 2, coeff0, coeff1, coeff2, coeff3);
        } else if (e_table == ehi) {
            jt = eos->nt - 1;
            elo = prj_eos_rey_slice_eint(eos, jy, jyp, jr, jrp, jt, coeff0, coeff1, coeff2, coeff3);
            ehi = prj_eos_rey_slice_eint(eos, jy, jyp, jr, jrp, eos->nt, coeff0, coeff1, coeff2, coeff3);
        } else {
            jt_lo = 1;
            jt_hi = eos->nt;
            while (jt_hi - jt_lo > 1) {
                jt_mid = jt_lo + (jt_hi - jt_lo) / 2;
                if (prj_eos_rey_slice_eint(eos, jy, jyp, jr, jrp, jt_mid, coeff0, coeff1, coeff2, coeff3) <= e_table) {
                    jt_lo = jt_mid;
                } else {
                    jt_hi = jt_mid;
                }
            }
            jt = jt_lo;
            elo = prj_eos_rey_slice_eint(eos, jy, jyp, jr, jrp, jt, coeff0, coeff1, coeff2, coeff3);
            ehi = prj_eos_rey_slice_eint(eos, jy, jyp, jr, jrp, jt + 1, coeff0, coeff1, coeff2, coeff3);
        }
        jtp = jt + 1;
        if (fabs(ehi - elo) > 0.0) {
            dtemp = prj_eos_clamp_double((e_table - elo) / (ehi - elo), 0.0, 1.0);
        } else {
            dtemp = 0.0;
        }
        T = exp(eos->ln10_t1 + ((double)(jt - 1) + dtemp) * eos->ln10_dlogT);
        pressure_log = prj_eos_pressure_log_interp(eos,
            jy, jyp, jr, jrp, jt, jtp, dye, drho, dtemp);
        gamma = prj_eos_table_interp_trilinear(eos, PRJ_EOS_REC_GAMMA,
            jy, jyp, jr, jrp, jt, jtp, dye, drho, dtemp);
        eos_quantities[PRJ_EOS_EINT] = eint;
        eos_quantities[PRJ_EOS_PRESSURE] = exp(pressure_log) * PRJ_EOS_PRESSURE_SCALE;
        eos_quantities[PRJ_EOS_GAMMA] = gamma;
        eos_quantities[PRJ_EOS_TEMPERATURE] = T;
        return;
    }

    gamma = prj_eos_gamma_value();
    T = (gamma - 1.0) * eint / PRJ_EOS_ENERGY_SCALE;
    eos_quantities[PRJ_EOS_EINT] = eint;
    eos_quantities[PRJ_EOS_PRESSURE] = rho * PRJ_EOS_ENERGY_SCALE * T;
    eos_quantities[PRJ_EOS_GAMMA] = gamma;
    eos_quantities[PRJ_EOS_TEMPERATURE] = T;
}

void prj_eos_fill_block(prj_eos *eos, prj_block *block, double *W)
{
    int i;
    int j;
    int k;

    if (block == 0 || W == 0 || block->eosvar == 0) {
        return;
    }

    for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
        for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
            for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                if (block->eos_done != 0 && block->eos_done[IDX(i, j, k)] != 0) {
                    continue;
                }
                prj_eos_fill_cell(eos, block, W, i, j, k);
            }
        }
    }
}

void prj_eos_fill_active_cells(prj_mesh *mesh, prj_eos *eos, int stage)
{
    int bidx;

    if (mesh == 0) {
        return;
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        prj_mpi *mpi = prj_mpi_current();
        double *W = stage == 2 ? block->W1 : block->W;
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1 || W == 0 || block->eosvar == 0) {
            continue;
        }
        if (mpi != 0 && block->rank != mpi->rank) {
            continue;
        }
        if (block->eos_done != 0) {
            memset(block->eos_done, 0, (size_t)PRJ_BLOCK_NCELLS * sizeof(*block->eos_done));
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    prj_eos_fill_cell(eos, block, W, i, j, k);
                }
            }
        }
    }
}

void prj_eos_fill_mesh(prj_mesh *mesh, prj_eos *eos, int stage)
{
    int bidx;

    if (mesh == 0) {
        return;
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        prj_mpi *mpi = prj_mpi_current();
        double *W = stage == 2 ? block->W1 : block->W;

        if (block->id < 0 || block->active != 1 || W == 0 || block->eosvar == 0) {
            continue;
        }
        if (mpi != 0 && block->rank != mpi->rank) {
            continue;
        }
        prj_eos_fill_block(eos, block, W);
    }
}

void prj_eos_fill_ghost_cons(prj_mesh *mesh, prj_eos *eos, int stage)
{
    int bidx;

    if (mesh == 0) {
        return;
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        prj_mpi *mpi = prj_mpi_current();
        double *W = stage == 2 ? block->W1 : block->W;
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1 || W == 0 || block->U == 0) {
            continue;
        }
        if (mpi != 0 && block->rank != mpi->rank) {
            continue;
        }
        for (i = -PRJ_NGHOST; i < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++i) {
            for (j = -PRJ_NGHOST; j < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++j) {
                for (k = -PRJ_NGHOST; k < PRJ_BLOCK_SIZE + PRJ_NGHOST; ++k) {
                    double Wc[PRJ_NVAR_PRIM];
                    double Uc[PRJ_NVAR_CONS];
                    int v;

                    if (i >= 0 && i < PRJ_BLOCK_SIZE &&
                        j >= 0 && j < PRJ_BLOCK_SIZE &&
                        k >= 0 && k < PRJ_BLOCK_SIZE) {
                        continue;
                    }
                    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
                        Wc[v] = W[VIDX(v, i, j, k)];
                    }
                    prj_eos_prim2cons(eos, Wc, Uc);
                    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
                        block->U[VIDX(v, i, j, k)] = Uc[v];
                    }
                }
            }
        }
    }
}

void prj_eos_prim2cons(prj_eos *eos, double *W, double *U)
{
    double rho;
    double eint;

    (void)eos;

    if (W == 0 || U == 0) {
        return;
    }

    rho = W[PRJ_PRIM_RHO];
    eint = W[PRJ_PRIM_EINT];

    U[PRJ_CONS_RHO] = rho;
    U[PRJ_CONS_MOM1] = rho * W[PRJ_PRIM_V1];
    U[PRJ_CONS_MOM2] = rho * W[PRJ_PRIM_V2];
    U[PRJ_CONS_MOM3] = rho * W[PRJ_PRIM_V3];
    U[PRJ_CONS_ETOT] = rho * eint +
        prj_eos_kinetic_energy_density_prim(W) +
        prj_eos_magnetic_energy_density_prim(W);
    U[PRJ_CONS_YE] = rho * W[PRJ_PRIM_YE];
#if PRJ_MHD
    U[PRJ_CONS_B1] = W[PRJ_PRIM_B1];
    U[PRJ_CONS_B2] = W[PRJ_PRIM_B2];
    U[PRJ_CONS_B3] = W[PRJ_PRIM_B3];
#endif
    prj_rad_prim2cons(W, U);
}

void prj_eos_cons2prim(prj_eos *eos, double *U, double *W)
{
    double rho;
    double v1;
    double v2;
    double v3;
    double kinetic;
    double magnetic;

    (void)eos;

    if (U == 0 || W == 0) {
        return;
    }

    rho = U[PRJ_CONS_RHO];
    if (rho == 0.0) {
        W[PRJ_PRIM_RHO] = 0.0;
        W[PRJ_PRIM_V1] = 0.0;
        W[PRJ_PRIM_V2] = 0.0;
        W[PRJ_PRIM_V3] = 0.0;
        W[PRJ_PRIM_EINT] = 0.0;
        W[PRJ_PRIM_YE] = 0.0;
#if PRJ_MHD
        W[PRJ_PRIM_B1] = U[PRJ_CONS_B1];
        W[PRJ_PRIM_B2] = U[PRJ_CONS_B2];
        W[PRJ_PRIM_B3] = U[PRJ_CONS_B3];
#endif
        prj_rad_cons2prim(U, W);
        return;
    }

    v1 = U[PRJ_CONS_MOM1] / rho;
    v2 = U[PRJ_CONS_MOM2] / rho;
    v3 = U[PRJ_CONS_MOM3] / rho;
    kinetic = prj_eos_kinetic_energy_density_cons(U) / rho;
    magnetic = prj_eos_magnetic_energy_density_cons(U) / rho;

    W[PRJ_PRIM_RHO] = rho;
    W[PRJ_PRIM_V1] = v1;
    W[PRJ_PRIM_V2] = v2;
    W[PRJ_PRIM_V3] = v3;
    W[PRJ_PRIM_EINT] = U[PRJ_CONS_ETOT] / rho - kinetic - magnetic;
    W[PRJ_PRIM_YE] = U[PRJ_CONS_YE] / rho;
#if PRJ_MHD
    W[PRJ_PRIM_B1] = U[PRJ_CONS_B1];
    W[PRJ_PRIM_B2] = U[PRJ_CONS_B2];
    W[PRJ_PRIM_B3] = U[PRJ_CONS_B3];
#endif
    prj_rad_cons2prim(U, W);
}

double prj_eos_kinetic_energy_density_prim(const double *W)
{
    if (W == 0) {
        return 0.0;
    }

    return 0.5 * W[PRJ_PRIM_RHO] * (
        W[PRJ_PRIM_V1] * W[PRJ_PRIM_V1] +
        W[PRJ_PRIM_V2] * W[PRJ_PRIM_V2] +
        W[PRJ_PRIM_V3] * W[PRJ_PRIM_V3]);
}

double prj_eos_kinetic_energy_density_cons(const double *U)
{
    double rho;

    if (U == 0) {
        return 0.0;
    }

    rho = U[PRJ_CONS_RHO];
    if (rho <= 0.0) {
        return 0.0;
    }

    return 0.5 * (
        U[PRJ_CONS_MOM1] * U[PRJ_CONS_MOM1] +
        U[PRJ_CONS_MOM2] * U[PRJ_CONS_MOM2] +
        U[PRJ_CONS_MOM3] * U[PRJ_CONS_MOM3]) / rho;
}

double prj_eos_magnetic_energy_density_prim(const double *W)
{
#if PRJ_MHD
    if (W == 0) {
        return 0.0;
    }

    return 0.5 * (
        W[PRJ_PRIM_B1] * W[PRJ_PRIM_B1] +
        W[PRJ_PRIM_B2] * W[PRJ_PRIM_B2] +
        W[PRJ_PRIM_B3] * W[PRJ_PRIM_B3]);
#else
    (void)W;
    return 0.0;
#endif
}

double prj_eos_magnetic_energy_density_cons(const double *U)
{
#if PRJ_MHD
    if (U == 0) {
        return 0.0;
    }

    return 0.5 * (
        U[PRJ_CONS_B1] * U[PRJ_CONS_B1] +
        U[PRJ_CONS_B2] * U[PRJ_CONS_B2] +
        U[PRJ_CONS_B3] * U[PRJ_CONS_B3]);
#else
    (void)U;
    return 0.0;
#endif
}
