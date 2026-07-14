#include <float.h>
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
#if PRJ_MIXED_PRECISION_TABLE
#define PRJ_TABLE_MPI_TYPE MPI_FLOAT
#else
#define PRJ_TABLE_MPI_TYPE MPI_DOUBLE
#endif
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
    PRJ_EOS_REC_PRESSURE = 2,
    PRJ_EOS_REC_GAMMA = 12
};

#define PRJ_EOS_COMPACT_NUMEL 4

static int prj_eos_rec_to_compact(int rec)
{
    switch (rec) {
        case 1: return 0;
        case 2: return 1;
        case 12: return 2;
        case 15: return 3;
        default: return 0;
    }
}

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
    return (size_t)prj_eos_rec_to_compact(rec) * (size_t)eos->ny * (size_t)eos->nr * (size_t)eos->nt +
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

static const char *prj_eos_call_ctx_label(enum prj_eos_call_ctx ctx)
{
    switch (ctx) {
    case PRJ_EOS_CTX_AMR:
        return "amr";
    case PRJ_EOS_CTX_MAIN:
    default:
        return "main loop";
    }
}

static void prj_eos_table_range_fail(const char *name, double value, double lo, double hi,
    double rho, double ye, enum prj_eos_call_ctx ctx)
{
    fprintf(stderr,
        "tabulated EOS lookup out of range for %s: value=%.17e allowed=[%.17e, %.17e] "
        "(cell rho=%.17e ye=%.17e, called from %s)\n",
        name, value, lo, hi, rho, ye, prj_eos_call_ctx_label(ctx));
#if defined(PRJ_ENABLE_MPI)
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    exit(1);
}

static void prj_eos_missing_table_filename_fail(void)
{
    fprintf(stderr, "tabulated EOS requires eos_file in the parameter file\n");
#if defined(PRJ_ENABLE_MPI)
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    exit(1);
}

static int prj_eos_compact_src_rec(int idx)
{
    switch (idx) {
        case 0: return 1;
        case 1: return 2;
        case 2: return 12;
        case 3: return 15;
        default: return 1;
    }
}

static int prj_eos_read_compact_table(FILE *f, prj_eos *eos, size_t slab_size)
{
    double *slab;
    int i;

    if (f == 0 || eos == 0 || eos->table == 0 || slab_size == 0) {
        return 1;
    }
    slab = (double *)prj_malloc(slab_size * sizeof(*slab));
    if (slab == 0) {
        return 1;
    }

    for (i = 0; i < PRJ_EOS_COMPACT_NUMEL; ++i) {
        int rec = prj_eos_compact_src_rec(i);
        size_t dst_off = (size_t)i * slab_size;
        long file_off = (long)((size_t)(rec - 1) * slab_size * sizeof(double));
        size_t idx;

        if (fseek(f, file_off, SEEK_SET) != 0 ||
            fread(slab, sizeof(*slab), slab_size, f) != slab_size) {
            free(slab);
            return 1;
        }
        for (idx = 0; idx < slab_size; ++idx) {
            double value = slab[idx];

            if (rec == PRJ_EOS_REC_PRESSURE) {
                value = log(value > 0.0 ? value : 1.0e-300);
            }
            eos->table[dst_off + idx] = (prj_table_real)value;
        }
    }

    free(slab);
    return 0;
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
        W[WIDX(PRJ_PRIM_RHO, i, j, k)],
        W[WIDX(PRJ_PRIM_V1, i, j, k)],
        W[WIDX(PRJ_PRIM_V2, i, j, k)],
        W[WIDX(PRJ_PRIM_V3, i, j, k)],
        W[WIDX(PRJ_PRIM_EINT, i, j, k)],
        W[WIDX(PRJ_PRIM_YE, i, j, k)]);
    prj_eos_print_fill_neighbors(block, x1, x2, x3);
    exit(1);
}

static void prj_eos_fill_cell(prj_eos *eos, prj_block *block, double *W, int i, int j, int k,
    enum prj_eos_call_ctx ctx)
{
    double eos_quantities[PRJ_EOS_NQUANT];

    if (block == 0 || W == 0 || block->eosvar == 0) {
        return;
    }
    if (W[WIDX(PRJ_PRIM_RHO, i, j, k)] == 0.0) {
        prj_eos_fail_zero_rho_in_block("prj_eos_fill_block", block, i, j, k, W);
    }

    prj_eos_rey(eos,
        W[WIDX(PRJ_PRIM_RHO, i, j, k)],
        W[WIDX(PRJ_PRIM_EINT, i, j, k)],
        W[WIDX(PRJ_PRIM_YE, i, j, k)],
        eos_quantities,
        ctx);
    block->eosvar[EIDX(PRJ_EOSVAR_PRESSURE, i, j, k)] = eos_quantities[PRJ_EOS_PRESSURE];
    block->eosvar[EIDX(PRJ_EOSVAR_TEMPERATURE, i, j, k)] = eos_quantities[PRJ_EOS_TEMPERATURE];
    block->eosvar[EIDX(PRJ_EOSVAR_GAMMA, i, j, k)] = eos_quantities[PRJ_EOS_GAMMA];
}

static void prj_eos_table_check_rty_inputs(const prj_eos *eos, double rho, double T, double ye,
    enum prj_eos_call_ctx ctx)
{
    if (rho < eos->rho_min || rho > eos->rho_max) {
        prj_eos_table_range_fail("rho", rho, eos->rho_min, eos->rho_max, rho, ye, ctx);
    }
    if (T < eos->temp_min || T > eos->temp_max) {
        prj_eos_table_range_fail("T", T, eos->temp_min, eos->temp_max, rho, ye, ctx);
    }
    if (ye < eos->y1c || ye > eos->y2c) {
        prj_eos_table_range_fail("ye", ye, eos->y1c, eos->y2c, rho, ye, ctx);
    }
}

static int prj_eos_prepare_table(prj_eos *eos, const prj_mpi *mpi)
{
    FILE *f;
    size_t slab_size;
    size_t table_count;
    size_t table_bytes;
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
    table_count = (size_t)PRJ_EOS_COMPACT_NUMEL * slab_size;
    table_bytes = table_count * sizeof(*eos->table);
#if defined(PRJ_ENABLE_MPI)
    if (mpi != 0 && mpi->totrank > 1) {
        if (mpi->rank == 0) {
            eos->table = (prj_table_real *)prj_malloc(table_bytes);
            if (eos->table == 0) {
                status = 1;
            } else {
                f = fopen(eos->filename, "rb");
                if (f == 0) {
                    status = 1;
                } else {
                    status = prj_eos_read_compact_table(f, eos, slab_size);
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
            eos->table = (prj_table_real *)prj_malloc(table_bytes);
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

        MPI_Bcast(eos->table, (int)table_count, PRJ_TABLE_MPI_TYPE, 0, MPI_COMM_WORLD);
        eos->table_bytes = table_bytes;
        eos->table_loaded = 1;
        eos->table_is_mmap = 0;
        return 0;
    }
#endif

    eos->table = (prj_table_real *)prj_malloc(table_bytes);
    if (eos->table == 0) {
        return 1;
    }
    f = fopen(eos->filename, "rb");
    if (f == 0) {
        free(eos->table);
        eos->table = 0;
        return 1;
    }
    status = prj_eos_read_compact_table(f, eos, slab_size);
    fclose(f);
    if (status != 0) {
        free(eos->table);
        eos->table = 0;
        return 1;
    }

    eos->table_bytes = table_bytes;
    eos->table_loaded = 1;
    eos->table_is_mmap = 0;
    return 0;
}

void prj_eos_init(prj_eos *eos, const prj_mpi *mpi)
{
    if (eos == 0) {
        return;
    }
    if (eos->kind == PRJ_EOS_KIND_TABLE && eos->filename[0] == '\0') {
        prj_eos_missing_table_filename_fail();
    }
    if (eos->kind != PRJ_EOS_KIND_TABLE) {
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
    (void)prj_eos_prepare_table(eos, mpi);
    if (eos->table_loaded == 1 && (mpi == 0 || mpi->rank == 0)) {
        fprintf(stderr, "EOS table initialization finished\n");
    }
}

static void prj_eos_table_interp_base(const prj_eos *eos, double rho, double T, double ye,
    int *restrict jy, int *restrict jyp, int *restrict jr, int *restrict jrp, int *restrict jt, int *restrict jtp,
    double *restrict dye, double *restrict drho, double *restrict dtemp, enum prj_eos_call_ctx ctx)
{
    double rl;
    double tl;

    prj_eos_table_check_rty_inputs(eos, rho, T, ye, ctx);
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

static double prj_eos_rey_slice_eint(const prj_eos *eos,
    int base_eint, int off_y_r, int off_yp_r, int off_y_rp, int off_yp_rp,
    int jt, double coeff0, double coeff1, double coeff2, double coeff3)
{
    int jt_off = jt - 1;
    const prj_table_real *t = eos->table;
    return coeff0 * t[base_eint + off_y_r   + jt_off] +
           coeff1 * t[base_eint + off_yp_r  + jt_off] +
           coeff2 * t[base_eint + off_y_rp  + jt_off] +
           coeff3 * t[base_eint + off_yp_rp + jt_off];
}

void prj_eos_rty(prj_eos *eos, double rho, double T, double ye, double *eos_quantities,
    enum prj_eos_call_ctx ctx)
{
    double gamma;
    double eint;

    if (eos_quantities == 0) {
        return;
    }
    if (eos != 0 && eos->kind == PRJ_EOS_KIND_TABLE &&
        eos->filename[0] != '\0' && prj_eos_prepare_table(eos, 0) == 0 && eos->table_loaded == 1) {
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

        prj_eos_table_interp_base(eos, rho, T, ye, &jy, &jyp, &jr, &jrp, &jt, &jtp, &dye, &drho, &dtemp, ctx);
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

/* Like prj_eos_rty() but returns only the internal energy (skipping the pressure
 * and gamma interpolations) plus its derivatives d(eint)/d(lnT) and
 * d(eint)/d(Ye).  Used in the radiation implicit-solver residual.  The table's
 * temperature axis is log10(T); the chain rule d(lgT)/d(lnT) = 1/M_LN10
 * converts the log10-T slope to a natural-log-T slope. */
double prj_eos_rty_eint(prj_eos *eos, double rho, double T, double ye,
    double *deint_dlnT, double *deint_dYe, enum prj_eos_call_ctx ctx)
{
    if (eos != 0 && eos->kind == PRJ_EOS_KIND_TABLE &&
        eos->filename[0] != '\0' && prj_eos_prepare_table(eos, 0) == 0 && eos->table_loaded == 1) {
        int jy;
        int jyp;
        int jr;
        int jrp;
        int jt;
        int jtp;
        double dye;
        double drho;
        double dtemp;
        double v[8];
        double dfd_dye;
        double dfd_drho;
        double dfd_dtemp;
        double e_raw;

        prj_eos_table_interp_base(eos, rho, T, ye, &jy, &jyp, &jr, &jrp, &jt, &jtp, &dye, &drho, &dtemp, ctx);
        /* Corner bit order for prj_trilinear_with_deriv: d0=dye, d1=drho, d2=dtemp. */
        v[0] = prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy,  jr,  jt);
        v[1] = prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jr,  jt);
        v[2] = prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy,  jrp, jt);
        v[3] = prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jrp, jt);
        v[4] = prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy,  jr,  jtp);
        v[5] = prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jr,  jtp);
        v[6] = prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy,  jrp, jtp);
        v[7] = prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jrp, jtp);
        e_raw = prj_trilinear_with_deriv(v, dye, drho, dtemp, &dfd_dye, &dfd_drho, &dfd_dtemp);

        *deint_dYe = PRJ_EOS_ENERGY_SCALE * dfd_dye * eos->inv_dYe;
        *deint_dlnT = PRJ_EOS_ENERGY_SCALE * dfd_dtemp * eos->inv_dlogT / M_LN10;
        return e_raw * PRJ_EOS_ENERGY_SCALE;
    }

    {
        double eint = PRJ_EOS_ENERGY_SCALE * T / (prj_eos_gamma_value() - 1.0);

        /* eint is linear in T, so d(eint)/d(lnT) = T d(eint)/dT = eint. */
        *deint_dlnT = eint;
        *deint_dYe = 0.0;
        return eint;
    }
}

double prj_eos_rty_geteta(prj_eos *eos, double rho, double T, double ye,
    enum prj_eos_call_ctx ctx)
{
    if (eos != 0 && eos->kind == PRJ_EOS_KIND_TABLE &&
        eos->filename[0] != '\0' && prj_eos_prepare_table(eos, 0) == 0 && eos->table_loaded == 1) {
        int jy, jyp, jr, jrp, jt, jtp;
        double dye, drho, dtemp;
        double eta_raw;

        prj_eos_table_interp_base(eos, rho, T, ye,
            &jy, &jyp, &jr, &jrp, &jt, &jtp, &dye, &drho, &dtemp, ctx);
        eta_raw = prj_eos_table_interp_trilinear(eos, 15,
            jy, jyp, jr, jrp, jt, jtp, dye, drho, dtemp);
        if (T > 0.0) {
            return eta_raw / T;
        }
        return 0.0;
    }

    return 0.0;
}

void prj_eos_rey(prj_eos *eos, double rho, double eint, double ye, double *eos_quantities,
    enum prj_eos_call_ctx ctx)
{
    double gamma;
    double T;

    if (eos_quantities == 0) {
        return;
    }
    if (eos != 0 && eos->kind == PRJ_EOS_KIND_TABLE &&
        eos->filename[0] != '\0' && prj_eos_prepare_table(eos, 0) == 0 && eos->table_loaded == 1) {
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
            prj_eos_table_range_fail("rho", rho, eos->rho_min, eos->rho_max, rho, ye, ctx);
        }
        if (ye < eos->y1c || ye > eos->y2c) {
            prj_eos_table_range_fail("ye", ye, eos->y1c, eos->y2c, rho, ye, ctx);
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

        int base_eint = prj_eos_rec_to_compact(PRJ_EOS_REC_EINT) * eos->ny * eos->nr * eos->nt;
        int off_y_r   = (jy  - 1) * eos->nr * eos->nt + (jr  - 1) * eos->nt;
        int off_yp_r  = (jyp - 1) * eos->nr * eos->nt + (jr  - 1) * eos->nt;
        int off_y_rp  = (jy  - 1) * eos->nr * eos->nt + (jrp - 1) * eos->nt;
        int off_yp_rp = (jyp - 1) * eos->nr * eos->nt + (jrp - 1) * eos->nt;

        elo = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r, off_y_rp, off_yp_rp, 1, coeff0, coeff1, coeff2, coeff3);
        ehi = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r, off_y_rp, off_yp_rp, eos->nt, coeff0, coeff1, coeff2, coeff3);
        if (e_table < elo || e_table > ehi) {
            prj_eos_table_range_fail("eint/PRJ_EOS_ENERGY_SCALE", e_table, elo, ehi, rho, ye, ctx);
        }
        if (e_table == elo) {
            jt = 1;
            ehi = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r, off_y_rp, off_yp_rp, 2, coeff0, coeff1, coeff2, coeff3);
        } else if (e_table == ehi) {
            jt = eos->nt - 1;
            elo = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r, off_y_rp, off_yp_rp, jt, coeff0, coeff1, coeff2, coeff3);
            ehi = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r, off_y_rp, off_yp_rp, eos->nt, coeff0, coeff1, coeff2, coeff3);
        } else {
            jt_lo = 1;
            jt_hi = eos->nt;
            while (jt_hi - jt_lo > 1) {
                jt_mid = jt_lo + (jt_hi - jt_lo) / 2;
                if (prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r, off_y_rp, off_yp_rp, jt_mid, coeff0, coeff1, coeff2, coeff3) <= e_table) {
                    jt_lo = jt_mid;
                } else {
                    jt_hi = jt_mid;
                }
            }
            jt = jt_lo;
            elo = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r, off_y_rp, off_yp_rp, jt, coeff0, coeff1, coeff2, coeff3);
            ehi = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r, off_y_rp, off_yp_rp, jt + 1, coeff0, coeff1, coeff2, coeff3);
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

static int prj_eos_pressure_try(prj_eos *eos, double rho, double eint, double ye,
    double *pressure)
{
    double gamma;

    if (pressure == 0 || !isfinite(rho) || !isfinite(eint) || !isfinite(ye) ||
        rho <= 0.0 || eint < 0.0) {
        return 0;
    }
    if (eos != 0 && eos->kind == PRJ_EOS_KIND_TABLE && eos->filename[0] != '\0') {
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
        int base_eint;
        int off_y_r;
        int off_yp_r;
        int off_y_rp;
        int off_yp_rp;

        if (prj_eos_prepare_table(eos, 0) != 0 || eos->table_loaded != 1) {
            return 0;
        }
        if (rho < eos->rho_min || rho > eos->rho_max ||
            ye < eos->y1c || ye > eos->y2c) {
            return 0;
        }
        e_table = eint / PRJ_EOS_ENERGY_SCALE;
        rl = log10(rho);
        jrf = (rl - eos->r1) * eos->inv_dlogrho;
        jyf = (ye - eos->y1c) * eos->inv_dYe;
        jr = 1 + (int)jrf;
        jy = 1 + (int)jyf;
        jr = jr < 1 ? 1 : (jr >= eos->nr ? eos->nr - 1 : jr);
        jy = jy < 1 ? 1 : (jy >= eos->ny ? eos->ny - 1 : jy);
        jrp = jr + 1;
        jyp = jy + 1;
        drho = prj_eos_clamp_double(
            (rl - (eos->r1 + (double)(jr - 1) * eos->dlogrho)) * eos->inv_dlogrho,
            0.0, 1.0);
        dye = prj_eos_clamp_double(
            (ye - (eos->y1c + (double)(jy - 1) * eos->dYe)) * eos->inv_dYe,
            0.0, 1.0);
        coeff0 = (1.0 - drho) * (1.0 - dye);
        coeff1 = (1.0 - drho) * dye;
        coeff2 = drho * (1.0 - dye);
        coeff3 = drho * dye;
        base_eint = prj_eos_rec_to_compact(PRJ_EOS_REC_EINT) * eos->ny * eos->nr * eos->nt;
        off_y_r = (jy - 1) * eos->nr * eos->nt + (jr - 1) * eos->nt;
        off_yp_r = (jyp - 1) * eos->nr * eos->nt + (jr - 1) * eos->nt;
        off_y_rp = (jy - 1) * eos->nr * eos->nt + (jrp - 1) * eos->nt;
        off_yp_rp = (jyp - 1) * eos->nr * eos->nt + (jrp - 1) * eos->nt;

        elo = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r,
            off_y_rp, off_yp_rp, 1, coeff0, coeff1, coeff2, coeff3);
        ehi = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r,
            off_y_rp, off_yp_rp, eos->nt, coeff0, coeff1, coeff2, coeff3);
        if (e_table < elo || e_table > ehi) {
            return 0;
        }
        if (e_table == elo) {
            jt = 1;
            ehi = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r,
                off_y_rp, off_yp_rp, 2, coeff0, coeff1, coeff2, coeff3);
        } else if (e_table == ehi) {
            jt = eos->nt - 1;
            elo = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r,
                off_y_rp, off_yp_rp, jt, coeff0, coeff1, coeff2, coeff3);
            ehi = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r,
                off_y_rp, off_yp_rp, eos->nt, coeff0, coeff1, coeff2, coeff3);
        } else {
            jt_lo = 1;
            jt_hi = eos->nt;
            while (jt_hi - jt_lo > 1) {
                jt_mid = jt_lo + (jt_hi - jt_lo) / 2;
                if (prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r,
                        off_y_rp, off_yp_rp, jt_mid, coeff0, coeff1, coeff2,
                        coeff3) <= e_table) {
                    jt_lo = jt_mid;
                } else {
                    jt_hi = jt_mid;
                }
            }
            jt = jt_lo;
            elo = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r,
                off_y_rp, off_yp_rp, jt, coeff0, coeff1, coeff2, coeff3);
            ehi = prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r,
                off_y_rp, off_yp_rp, jt + 1, coeff0, coeff1, coeff2, coeff3);
        }
        jtp = jt + 1;
        if (fabs(ehi - elo) > 0.0) {
            dtemp = prj_eos_clamp_double((e_table - elo) / (ehi - elo), 0.0, 1.0);
        } else {
            dtemp = 0.0;
        }
        pressure_log = prj_eos_pressure_log_interp(eos,
            jy, jyp, jr, jrp, jt, jtp, dye, drho, dtemp);
        *pressure = exp(pressure_log) * PRJ_EOS_PRESSURE_SCALE;
        return isfinite(*pressure) && *pressure >= 0.0;
    }

    gamma = prj_eos_gamma_value();
    *pressure = (gamma - 1.0) * rho * eint;
    return isfinite(*pressure) && *pressure >= 0.0;
}

double prj_eos_low_temp_eint(prj_eos *eos, double rho, double ye, enum prj_eos_call_ctx ctx)
{
    double rl;
    double drho;
    double dye;
    int jr;
    int jrp;
    int jy;
    int jyp;
    double coeff0;
    double coeff1;
    double coeff2;
    double coeff3;
    int base_eint;
    int off_y_r;
    int off_yp_r;
    int off_y_rp;
    int off_yp_rp;

    /* The internal-energy offset only exists for the tabulated EOS; for the
     * ideal-gas EOS the floor is applied directly so there is no boundary. */
    if (eos == 0 || eos->kind != PRJ_EOS_KIND_TABLE || eos->filename[0] == '\0' ||
        prj_eos_prepare_table(eos, 0) != 0 || eos->table_loaded != 1) {
        return 0.0;
    }

    if (rho < eos->rho_min || rho > eos->rho_max) {
        prj_eos_table_range_fail("rho", rho, eos->rho_min, eos->rho_max, rho, ye, ctx);
    }
    if (ye < eos->y1c || ye > eos->y2c) {
        prj_eos_table_range_fail("ye", ye, eos->y1c, eos->y2c, rho, ye, ctx);
    }

    rl = log10(rho);
    jr = 1 + (int)((rl - eos->r1) * eos->inv_dlogrho);
    jy = 1 + (int)((ye - eos->y1c) * eos->inv_dYe);
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

    base_eint = prj_eos_rec_to_compact(PRJ_EOS_REC_EINT) * eos->ny * eos->nr * eos->nt;
    off_y_r   = (jy  - 1) * eos->nr * eos->nt + (jr  - 1) * eos->nt;
    off_yp_r  = (jyp - 1) * eos->nr * eos->nt + (jr  - 1) * eos->nt;
    off_y_rp  = (jy  - 1) * eos->nr * eos->nt + (jrp - 1) * eos->nt;
    off_yp_rp = (jyp - 1) * eos->nr * eos->nt + (jrp - 1) * eos->nt;

    /* jt = 1 is the lowest-temperature slice of the table. */
    return prj_eos_rey_slice_eint(eos, base_eint, off_y_r, off_yp_r, off_y_rp,
        off_yp_rp, 1, coeff0, coeff1, coeff2, coeff3) * PRJ_EOS_ENERGY_SCALE;
}

void prj_eos_fill_block(prj_eos *eos, prj_block *block, double *W, enum prj_eos_call_ctx ctx)
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
                if (block->cell_derived_done != 0 && block->cell_derived_done[IDX(i, j, k)] != 0) {
                    continue;
                }
                prj_eos_fill_cell(eos, block, W, i, j, k, ctx);
            }
        }
    }
}

void prj_eos_fill_active_cells(prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi, int stage,
    enum prj_eos_call_ctx ctx)
{
    int bidx;

    if (mesh == 0) {
        return;
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *W = prj_block_prim_stage(block, prj_stage_slot_from_stage_arg(stage));
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1 || W == 0 || block->eosvar == 0) {
            continue;
        }
        if (mpi != 0 && block->rank != mpi->rank) {
            continue;
        }
        for (i = 0; i < PRJ_BLOCK_SIZE; ++i) {
            for (j = 0; j < PRJ_BLOCK_SIZE; ++j) {
                for (k = 0; k < PRJ_BLOCK_SIZE; ++k) {
                    prj_eos_fill_cell(eos, block, W, i, j, k, ctx);
                }
            }
        }
    }
}

void prj_eos_fill_mesh(prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi, int stage,
    enum prj_eos_call_ctx ctx)
{
    int bidx;

    if (mesh == 0) {
        return;
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *W = prj_block_prim_stage(block, prj_stage_slot_from_stage_arg(stage));

        if (block->id < 0 || block->active != 1 || W == 0 || block->eosvar == 0) {
            continue;
        }
        if (mpi != 0 && block->rank != mpi->rank) {
            continue;
        }
        prj_eos_fill_block(eos, block, W, ctx);
    }
}

void prj_eos_fill_ghost_cons(prj_mesh *mesh, prj_eos *eos, const prj_mpi *mpi, int stage,
    enum prj_eos_call_ctx ctx)
{
    int bidx;

    /* This routine rebuilds conserved ghost values without an EOS lookup on the
     * non-GR path; ctx is still passed through for dynamic-GR recovery errors. */

    if (mesh == 0) {
        return;
    }

    for (bidx = 0; bidx < mesh->nblocks; ++bidx) {
        prj_block *block = &mesh->blocks[bidx];
        double *W = prj_block_prim_stage(block, prj_stage_slot_from_stage_arg(stage));
        int i;
        int j;
        int k;

        if (block->id < 0 || block->active != 1 || W == 0 ||
            !prj_block_has_cons_storage(block)) {
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
#if PRJ_NRAD > 0
                    int in_rad_zone =
                        i >= -PRJ_NGHOST_RAD && i < PRJ_BLOCK_SIZE + PRJ_NGHOST_RAD &&
                        j >= -PRJ_NGHOST_RAD && j < PRJ_BLOCK_SIZE + PRJ_NGHOST_RAD &&
                        k >= -PRJ_NGHOST_RAD && k < PRJ_BLOCK_SIZE + PRJ_NGHOST_RAD;
#endif

                    if (i >= 0 && i < PRJ_BLOCK_SIZE &&
                        j >= 0 && j < PRJ_BLOCK_SIZE &&
                        k >= 0 && k < PRJ_BLOCK_SIZE) {
                        continue;
                    }
                    for (v = 0; v < PRJ_NHYDRO; ++v) {
                        Wc[v] = W[WIDX(v, i, j, k)];
                    }
#if PRJ_NRAD > 0
                    if (in_rad_zone) {
                        for (v = PRJ_NHYDRO; v < PRJ_NVAR_PRIM; ++v) {
                            Wc[v] = W[WIDX(v, i, j, k)];
                        }
                    } else {
                        for (v = PRJ_NHYDRO; v < PRJ_NVAR_PRIM; ++v) {
                            Wc[v] = 0.0;
                        }
                    }
#endif
                    prj_eos_cell_prim2cons(eos, mesh, block,
                        prj_stage_slot_from_stage_arg(stage), i, j, k, Wc, Uc, ctx);
                    for (v = 0; v < PRJ_NHYDRO; ++v) {
                        prj_block_set_cons_value(block, v, i, j, k, Uc[v]);
                    }
#if PRJ_NRAD > 0
                    if (in_rad_zone) {
                        for (v = PRJ_NHYDRO; v < PRJ_NVAR_CONS; ++v) {
                            prj_block_set_cons_value(block, v, i, j, k, Uc[v]);
                        }
                    }
#endif
                }
            }
        }
    }
}

typedef struct prj_eos_gr_recovery {
    prj_eos *eos;
    double gu[3][3];
    double Bcon[3];
    double Bcov[3];
    double Scov[3];
    double D;
    double tau;
    double Bsq;
    double Ssq;
    double SB;
    double c2;
    double ye;
} prj_eos_gr_recovery;

static double prj_eos_gr_det3(const double a[3][3])
{
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

static void prj_eos_gr_inv3(const double a[3][3], double inv[3][3], double det)
{
    double idet = 1.0 / det;

    inv[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * idet;
    inv[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * idet;
    inv[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * idet;
    inv[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * idet;
    inv[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * idet;
    inv[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * idet;
    inv[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * idet;
    inv[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * idet;
    inv[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * idet;
}

static int prj_eos_gr_metric_prepare(const prj_eos_gr_geom *geom, double g[3][3],
    double gu[3][3], double *det, double *sqrt_det)
{
    double scale = 0.0;
    double sym_tol;
    double minor2;
    int a;
    int b;

    if (geom == 0 || det == 0 || sqrt_det == 0) {
        return 0;
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            g[a][b] = geom->gamma[a][b];
            if (!isfinite(g[a][b])) {
                return 0;
            }
            if (fabs(g[a][b]) > scale) {
                scale = fabs(g[a][b]);
            }
        }
    }
    if (scale <= 0.0) {
        return 0;
    }
    sym_tol = 1.0e-12 * scale;
    for (a = 0; a < 3; ++a) {
        for (b = a + 1; b < 3; ++b) {
            if (fabs(g[a][b] - g[b][a]) > sym_tol) {
                return 0;
            }
        }
    }
    minor2 = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    *det = prj_eos_gr_det3(g);
    if (g[0][0] <= 0.0 || minor2 <= 0.0 || *det <= 0.0 || !isfinite(*det)) {
        return 0;
    }
    prj_eos_gr_inv3(g, gu, *det);
    *sqrt_det = sqrt(*det);
    return isfinite(*sqrt_det) && *sqrt_det > 0.0;
}

static void prj_eos_gr_lower(const double g[3][3], const double vcon[3], double vcov[3])
{
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        vcov[a] = 0.0;
        for (b = 0; b < 3; ++b) {
            vcov[a] += g[a][b] * vcon[b];
        }
    }
}

static void prj_eos_gr_raise(const double gu[3][3], const double vcov[3], double vcon[3])
{
    int a;
    int b;

    for (a = 0; a < 3; ++a) {
        vcon[a] = 0.0;
        for (b = 0; b < 3; ++b) {
            vcon[a] += gu[a][b] * vcov[b];
        }
    }
}

static double prj_eos_gr_dot_cov_con(const double vcov[3], const double wcon[3])
{
    return vcov[0] * wcon[0] + vcov[1] * wcon[1] + vcov[2] * wcon[2];
}

static int prj_eos_gr_valid_array(const double *a, int n)
{
    int i;

    if (a == 0) {
        return 0;
    }
    for (i = 0; i < n; ++i) {
        if (!isfinite(a[i])) {
            return 0;
        }
    }
    return 1;
}

static double prj_eos_gr_lorentz_minus_one(double beta2, double sqrt_one_minus_beta2)
{
    return beta2 / (sqrt_one_minus_beta2 * (1.0 + sqrt_one_minus_beta2));
}

static int prj_eos_gr_recovery_residual(const prj_eos_gr_recovery *rec, double x,
    double *residual, double *rho_out, double *eint_out, double beta_con_out[3])
{
    double z;
    double denom;
    double sb_over_z;
    double beta2;
    double sqrt_one_minus_beta2;
    double wlor;
    double wlor2;
    double wlor_m1;
    double rho;
    double p;
    double eps;
    double eint;
    double pressure;

    if (rec == 0 || !isfinite(x) || x < 0.0 || rec->D <= 0.0) {
        return 0;
    }
    z = rec->D + x;
    denom = z + rec->Bsq;
    if (!isfinite(z) || !isfinite(denom) || z <= 0.0 || denom <= 0.0) {
        return 0;
    }
    sb_over_z = rec->SB / z;
    beta2 = (rec->Ssq + rec->SB * rec->SB * (2.0 * z + rec->Bsq) / (z * z)) /
        (denom * denom);
    if (!isfinite(beta2) || beta2 < 0.0 || beta2 >= 1.0) {
        return 0;
    }
    sqrt_one_minus_beta2 = sqrt(1.0 - beta2);
    wlor = 1.0 / sqrt_one_minus_beta2;
    wlor2 = wlor * wlor;
    wlor_m1 = prj_eos_gr_lorentz_minus_one(beta2, sqrt_one_minus_beta2);
    rho = rec->D / wlor;
    p = x + rec->Bsq - rec->tau -
        0.5 * (rec->Bsq / wlor2 + sb_over_z * sb_over_z);
    if (!isfinite(rho) || !isfinite(p) || rho <= 0.0 || p <= 0.0) {
        return 0;
    }
    eps = x / (rec->D * wlor) - wlor_m1 / wlor - p / rho;
    if (!isfinite(eps) || eps < 0.0) {
        return 0;
    }
    eint = eps * rec->c2;
    if (!isfinite(eint) ||
        !prj_eos_pressure_try(rec->eos, rho, eint, rec->ye, &pressure)) {
        return 0;
    }
    if (residual != 0) {
        *residual = pressure / rec->c2 - p;
    }
    if (rho_out != 0) {
        *rho_out = rho;
    }
    if (eint_out != 0) {
        *eint_out = eint;
    }
    if (beta_con_out != 0) {
        double beta_cov[3];

        beta_cov[0] = (rec->Scov[0] + sb_over_z * rec->Bcov[0]) / denom;
        beta_cov[1] = (rec->Scov[1] + sb_over_z * rec->Bcov[1]) / denom;
        beta_cov[2] = (rec->Scov[2] + sb_over_z * rec->Bcov[2]) / denom;
        prj_eos_gr_raise(rec->gu, beta_cov, beta_con_out);
    }
    return residual == 0 || isfinite(*residual);
}

static int prj_eos_gr_recovery_residual_ye(const prj_eos_gr_recovery *rec, double ye,
    double x, double *residual, double *rho_out, double *eint_out, double beta_con_out[3])
{
    prj_eos_gr_recovery trial = *rec;

    if (!isfinite(ye)) {
        return 0;
    }
    trial.ye = ye;
    if (!prj_eos_gr_recovery_residual(&trial, x, residual, rho_out, eint_out,
            beta_con_out)) {
        return 0;
    }
    return 1;
}

/* Bisect for the recovery root on [xlo, xhi], assuming the residual is
   monotonically decreasing (positive below the root, negative above it) on the
   valid interval. xhi must be a valid sample with residual <= 0. xlo may lie in
   the invalid region below the validity boundary (e.g. where p <= 0): such a
   point sits below the root, so an invalid midpoint simply raises the lower
   bound. This brackets roots that fall in the first valid geometric bin, where
   the sign change straddles the invalid->valid boundary. */
static int prj_eos_gr_bisect_recovery(const prj_eos_gr_recovery *rec, double ye,
    double xlo, double xhi, double *xroot)
{
    double scale;
    int iter;

    if (xroot == 0) {
        return 0;
    }
    scale = fmax(fmax(fabs(rec->tau), rec->Bsq), fmax(fabs(xlo), fabs(xhi)));
    if (scale <= 0.0) {
        scale = 1.0;
    }
    for (iter = 0; iter < 200; ++iter) {
        double xm = 0.5 * (xlo + xhi);
        double fm;

        if (!prj_eos_gr_recovery_residual_ye(rec, ye, xm, &fm, 0, 0, 0)) {
            xlo = xm;
            continue;
        }
        if (fabs(fm) <= 1.0e-14 * scale ||
            fabs(xhi - xlo) <= 1.0e-14 * (1.0 + fabs(xm))) {
            *xroot = xm;
            return 1;
        }
        if (fm > 0.0) {
            xlo = xm;
        } else {
            xhi = xm;
        }
    }
    *xroot = 0.5 * (xlo + xhi);
    return 1;
}

static int prj_eos_gr_solve_recovery(const prj_eos_gr_recovery *rec, double ye,
    double *xroot)
{
    double x_scale;
    double x;
    double prev_x = 0.0;
    double prev_f = 0.0;
    double last_x = 0.0;
    int have_prev = 0;
    int saw_valid = 0;
    int scan;

    if (rec == 0 || xroot == 0 || rec->D <= 0.0 || !isfinite(ye)) {
        return PRJ_EOS_GR_BAD_STATE;
    }
    x_scale = fabs(rec->tau) + sqrt(fmax(0.0, rec->Ssq)) + rec->Bsq;
    if (rec->Bsq > 0.0) {
        x_scale += fabs(rec->SB) / sqrt(rec->Bsq);
    }
    if (!isfinite(x_scale) || x_scale <= 0.0) {
        x_scale = rec->D * DBL_EPSILON;
    }
    if (x_scale <= 0.0) {
        return PRJ_EOS_GR_BAD_STATE;
    }
    x = x_scale * 1.0e-14;
    if (x < DBL_MIN) {
        x = DBL_MIN;
    }
    for (scan = 0; scan < 420; ++scan) {
        double f;
        double xtest = scan == 0 ? 0.0 : x;

        if (prj_eos_gr_recovery_residual_ye(rec, ye, xtest, &f, 0, 0, 0)) {
            double tol = 1.0e-14 * fmax(1.0, x_scale);

            saw_valid = 1;
            if (fabs(f) <= tol) {
                *xroot = xtest;
                return PRJ_EOS_GR_OK;
            }
            if (have_prev &&
                ((prev_f <= 0.0 && f >= 0.0) || (prev_f >= 0.0 && f <= 0.0))) {
                /* Sign change between two valid samples. */
                return prj_eos_gr_bisect_recovery(rec, ye, prev_x, xtest,
                    xroot) ? PRJ_EOS_GR_OK : PRJ_EOS_GR_NO_CONVERGE;
            }
            if (!have_prev && f < 0.0) {
                /* First valid sample already overshot the (decreasing) root:
                   the root sits in the first valid geometric bin, so bracket
                   against the previous, still-invalid lower point. */
                return prj_eos_gr_bisect_recovery(rec, ye, last_x, xtest,
                    xroot) ? PRJ_EOS_GR_OK : PRJ_EOS_GR_NO_CONVERGE;
            }
            prev_x = xtest;
            prev_f = f;
            have_prev = 1;
        }
        last_x = xtest;
        if (scan > 0) {
            if (x > DBL_MAX / 1.35) {
                break;
            }
            x *= 1.35;
        }
    }
    return saw_valid ? PRJ_EOS_GR_NO_CONVERGE : PRJ_EOS_GR_BAD_STATE;
}

int prj_eos_gr_prim2cons(prj_eos *eos, const prj_eos_gr_geom *geom,
    const double *W, double *U, enum prj_eos_call_ctx ctx)
{
    double g[3][3];
    double gu[3][3];
    double det;
    double sqrt_det;
    double Utmp[PRJ_NVAR_CONS];
    double rho;
    double vcon[3];
    double beta_con[3];
    double beta_cov[3];
    double beta2;
    double sqrt_one_minus_beta2;
    double wlor;
    double wlor2;
    double wlor_m1;
    double eint;
    double ye;
    double pressure;
    double c = PRJ_CLIGHT;
    double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    double Bcon[3] = {0.0, 0.0, 0.0};
    double Bcov[3] = {0.0, 0.0, 0.0};
    double Bsq = 0.0;
    double Bbeta = 0.0;
    double w;
    int v;
    int d;

    (void)ctx;
    (void)gu;
    if (W == 0 || U == 0 || geom == 0) {
        return PRJ_EOS_GR_NULL_ARG;
    }
    if (!prj_eos_gr_metric_prepare(geom, g, gu, &det, &sqrt_det)) {
        return PRJ_EOS_GR_BAD_METRIC;
    }
    if (!prj_eos_gr_valid_array(W, PRJ_NVAR_PRIM)) {
        return PRJ_EOS_GR_BAD_STATE;
    }
    rho = W[PRJ_PRIM_RHO];
    vcon[0] = W[PRJ_PRIM_V1];
    vcon[1] = W[PRJ_PRIM_V2];
    vcon[2] = W[PRJ_PRIM_V3];
    eint = W[PRJ_PRIM_EINT];
    ye = W[PRJ_PRIM_YE];
    if (rho <= 0.0 || eint < 0.0 ||
        !prj_eos_pressure_try(eos, rho, eint, ye, &pressure)) {
        return PRJ_EOS_GR_BAD_STATE;
    }
    for (d = 0; d < 3; ++d) {
        beta_con[d] = vcon[d] / c;
    }
    prj_eos_gr_lower(g, beta_con, beta_cov);
    beta2 = prj_eos_gr_dot_cov_con(beta_cov, beta_con);
    if (!isfinite(beta2) || beta2 < 0.0 || beta2 >= 1.0) {
        return PRJ_EOS_GR_BAD_STATE;
    }
    sqrt_one_minus_beta2 = sqrt(1.0 - beta2);
    wlor = 1.0 / sqrt_one_minus_beta2;
    wlor2 = wlor * wlor;
    wlor_m1 = prj_eos_gr_lorentz_minus_one(beta2, sqrt_one_minus_beta2);
#if PRJ_MHD
    Bcon[0] = W[PRJ_PRIM_B1];
    Bcon[1] = W[PRJ_PRIM_B2];
    Bcon[2] = W[PRJ_PRIM_B3];
#endif
    prj_eos_gr_lower(g, Bcon, Bcov);
    Bsq = prj_eos_gr_dot_cov_con(Bcov, Bcon);
    Bbeta = prj_eos_gr_dot_cov_con(Bcov, beta_con);
    if (!isfinite(Bsq) || !isfinite(Bbeta) || Bsq < 0.0) {
        return PRJ_EOS_GR_BAD_STATE;
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        Utmp[v] = 0.0;
    }
    w = rho * c2 + rho * eint + pressure;
    Utmp[PRJ_CONS_RHO] = rho * wlor;
    for (d = 0; d < 3; ++d) {
        Utmp[PRJ_CONS_MOM1 + d] =
            ((w + Bsq) * wlor2 * beta_cov[d] - Bbeta * Bcov[d]) / c;
    }
    Utmp[PRJ_CONS_ETOT] = (rho * eint + pressure) * wlor2 +
        rho * c2 * wlor * wlor_m1 - pressure + Bsq -
        0.5 * (Bbeta * Bbeta + Bsq / wlor2);
    Utmp[PRJ_CONS_YE] = Utmp[PRJ_CONS_RHO] * ye;
#if PRJ_MHD
    Utmp[PRJ_CONS_B1] = Bcon[0];
    Utmp[PRJ_CONS_B2] = Bcon[1];
    Utmp[PRJ_CONS_B3] = Bcon[2];
#endif
    prj_rad_prim2cons(W, Utmp);
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        if (!isfinite(Utmp[v])) {
            return PRJ_EOS_GR_BAD_STATE;
        }
        U[v] = sqrt_det * Utmp[v];
    }
    return PRJ_EOS_GR_OK;
}

int prj_eos_gr_cons2prim(prj_eos *eos, const prj_eos_gr_geom *geom,
    const double *U, double *W, enum prj_eos_call_ctx ctx)
{
    double g[3][3];
    double gu[3][3];
    double det;
    double sqrt_det;
    double Uloc[PRJ_NVAR_CONS];
    double Wtmp[PRJ_NVAR_PRIM];
    double c = PRJ_CLIGHT;
    double c2 = PRJ_CLIGHT * PRJ_CLIGHT;
    double ye;
    double xroot;
    double rho;
    double eint;
    double beta_con[3];
    prj_eos_gr_recovery rec;
    double Scon[3];
    int status;
    int v;
    int d;

    (void)ctx;
    if (U == 0 || W == 0 || geom == 0) {
        return PRJ_EOS_GR_NULL_ARG;
    }
    if (!prj_eos_gr_metric_prepare(geom, g, gu, &det, &sqrt_det)) {
        return PRJ_EOS_GR_BAD_METRIC;
    }
    if (!prj_eos_gr_valid_array(U, PRJ_NVAR_CONS)) {
        return PRJ_EOS_GR_BAD_STATE;
    }
    for (v = 0; v < PRJ_NVAR_CONS; ++v) {
        Uloc[v] = U[v] / sqrt_det;
    }
    memset(&rec, 0, sizeof(rec));
    rec.eos = eos;
    rec.c2 = c2;
    rec.D = Uloc[PRJ_CONS_RHO];
    rec.tau = Uloc[PRJ_CONS_ETOT] / c2;
    if (!isfinite(rec.D) || !isfinite(rec.tau) || rec.D <= 0.0) {
        return PRJ_EOS_GR_BAD_STATE;
    }
    ye = Uloc[PRJ_CONS_YE] / rec.D;
    if (!isfinite(ye)) {
        return PRJ_EOS_GR_BAD_STATE;
    }
    rec.ye = ye;
    for (d = 0; d < 3; ++d) {
        int mom = PRJ_CONS_MOM1 + d;

        rec.Scov[d] = Uloc[mom] / c;
        rec.Bcon[d] = 0.0;
        rec.gu[d][0] = gu[d][0];
        rec.gu[d][1] = gu[d][1];
        rec.gu[d][2] = gu[d][2];
    }
#if PRJ_MHD
    rec.Bcon[0] = Uloc[PRJ_CONS_B1] / c;
    rec.Bcon[1] = Uloc[PRJ_CONS_B2] / c;
    rec.Bcon[2] = Uloc[PRJ_CONS_B3] / c;
#endif
    prj_eos_gr_lower(g, rec.Bcon, rec.Bcov);
    prj_eos_gr_raise(gu, rec.Scov, Scon);
    rec.Bsq = prj_eos_gr_dot_cov_con(rec.Bcov, rec.Bcon);
    rec.Ssq = prj_eos_gr_dot_cov_con(rec.Scov, Scon);
    rec.SB = prj_eos_gr_dot_cov_con(rec.Scov, rec.Bcon);
    if (!isfinite(rec.Bsq) || !isfinite(rec.Ssq) || !isfinite(rec.SB) ||
        rec.Bsq < 0.0 || rec.Ssq < 0.0) {
        return PRJ_EOS_GR_BAD_STATE;
    }
    status = prj_eos_gr_solve_recovery(&rec, ye, &xroot);
    if (status != PRJ_EOS_GR_OK) {
        return status;
    }
    if (!prj_eos_gr_recovery_residual_ye(&rec, ye, xroot, 0, &rho, &eint, beta_con)) {
        return PRJ_EOS_GR_NO_CONVERGE;
    }
    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        Wtmp[v] = 0.0;
    }
    Wtmp[PRJ_PRIM_RHO] = rho;
    Wtmp[PRJ_PRIM_V1] = c * beta_con[0];
    Wtmp[PRJ_PRIM_V2] = c * beta_con[1];
    Wtmp[PRJ_PRIM_V3] = c * beta_con[2];
    Wtmp[PRJ_PRIM_EINT] = eint;
    Wtmp[PRJ_PRIM_YE] = ye;
#if PRJ_MHD
    Wtmp[PRJ_PRIM_B1] = Uloc[PRJ_CONS_B1];
    Wtmp[PRJ_PRIM_B2] = Uloc[PRJ_CONS_B2];
    Wtmp[PRJ_PRIM_B3] = Uloc[PRJ_CONS_B3];
#endif
    prj_rad_cons2prim(Uloc, Wtmp);
    if (!prj_eos_gr_valid_array(Wtmp, PRJ_NVAR_PRIM)) {
        return PRJ_EOS_GR_BAD_STATE;
    }
    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        W[v] = Wtmp[v];
    }
    return PRJ_EOS_GR_OK;
}

int prj_eos_full_dynamic_gr_enabled(const prj_mesh *mesh)
{
#if PRJ_DYNAMIC_GR && !PRJ_MHD
    return mesh != 0 && mesh->use_full_dynamic_gr != 0;
#else
    (void)mesh;
    return 0;
#endif
}

static const char *prj_eos_gr_status_name(int status)
{
    switch (status) {
    case PRJ_EOS_GR_OK:
        return "ok";
    case PRJ_EOS_GR_NULL_ARG:
        return "null_arg";
    case PRJ_EOS_GR_BAD_METRIC:
        return "bad_metric";
    case PRJ_EOS_GR_BAD_STATE:
        return "bad_state";
    case PRJ_EOS_GR_NO_CONVERGE:
        return "no_converge";
    default:
        return "missing_geometry";
    }
}

static void prj_eos_gr_cell_fail(const char *op, int status,
    int i, int j, int k, enum prj_eos_call_ctx ctx)
{
    fprintf(stderr, "dynamic GR %s failed at cell (%d,%d,%d): %s (ctx=%d)\n",
        op, i, j, k, prj_eos_gr_status_name(status), (int)ctx);
#if defined(PRJ_ENABLE_MPI)
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
    exit(EXIT_FAILURE);
}

static void prj_eos_gr_load_cell_geom(const prj_mesh *mesh, const prj_block *block,
    int z4c_stage, int i, int j, int k, prj_eos_gr_geom *geom,
    enum prj_eos_call_ctx ctx)
{
    prj_z4c_hydro_geom zgeom;
    int a, b;

    if (!prj_z4c_load_hydro_geom(mesh, block, z4c_stage, i, j, k, &zgeom)) {
        prj_eos_gr_cell_fail("geometry load", -1, i, j, k, ctx);
    }
    for (a = 0; a < 3; ++a) {
        for (b = 0; b < 3; ++b) {
            geom->gamma[a][b] = zgeom.gamma[a][b];
        }
    }
}

void prj_eos_cell_prim2cons(prj_eos *eos, const prj_mesh *mesh,
    const prj_block *block, int z4c_stage, int i, int j, int k,
    const double *W, double *U, enum prj_eos_call_ctx ctx)
{
    prj_eos_gr_geom geom;
    double Ugr[PRJ_NVAR_CONS];
    int status;
    int v;

    if (!prj_eos_full_dynamic_gr_enabled(mesh)) {
        prj_eos_prim2cons(eos, (double *)W, U);
        return;
    }

    prj_eos_prim2cons(eos, (double *)W, U);
    prj_eos_gr_load_cell_geom(mesh, block, z4c_stage, i, j, k, &geom, ctx);
    status = prj_eos_gr_prim2cons(eos, &geom, W, Ugr, ctx);
    if (status != PRJ_EOS_GR_OK) {
        prj_eos_gr_cell_fail("prim2cons", status, i, j, k, ctx);
    }
    for (v = 0; v < PRJ_NHYDRO; ++v) {
        U[v] = Ugr[v];
    }
}

void prj_eos_cell_cons2prim(prj_eos *eos, const prj_mesh *mesh,
    const prj_block *block, int z4c_stage, int i, int j, int k,
    const double *U, double *W, enum prj_eos_call_ctx ctx)
{
    prj_eos_gr_geom geom;
    double Wgr[PRJ_NVAR_PRIM];
    int status;
    int v;

    if (!prj_eos_full_dynamic_gr_enabled(mesh)) {
        prj_eos_cons2prim(eos, (double *)U, W);
        return;
    }

    for (v = 0; v < PRJ_NVAR_PRIM; ++v) {
        W[v] = 0.0;
    }
    prj_rad_cons2prim(U, W);
    prj_eos_gr_load_cell_geom(mesh, block, z4c_stage, i, j, k, &geom, ctx);
    status = prj_eos_gr_cons2prim(eos, &geom, U, Wgr, ctx);
    if (status != PRJ_EOS_GR_OK) {
        prj_eos_gr_cell_fail("cons2prim", status, i, j, k, ctx);
    }
    for (v = 0; v < PRJ_NHYDRO; ++v) {
        W[v] = Wgr[v];
    }
}

void prj_eos_prim2cons(prj_eos *eos, double *W, double *U)
{
    double rho;
    double v1;
    double v2;
    double v3;
    double eint;
#if PRJ_MHD
    double b1;
    double b2;
    double b3;
#endif

    (void)eos;

    if (W == 0 || U == 0) {
        return;
    }

    rho = W[PRJ_PRIM_RHO];
    v1 = W[PRJ_PRIM_V1];
    v2 = W[PRJ_PRIM_V2];
    v3 = W[PRJ_PRIM_V3];
    eint = W[PRJ_PRIM_EINT];
#if PRJ_MHD
    b1 = W[PRJ_PRIM_B1];
    b2 = W[PRJ_PRIM_B2];
    b3 = W[PRJ_PRIM_B3];
#endif

    U[PRJ_CONS_RHO] = rho;
    U[PRJ_CONS_MOM1] = rho * v1;
    U[PRJ_CONS_MOM2] = rho * v2;
    U[PRJ_CONS_MOM3] = rho * v3;
    U[PRJ_CONS_ETOT] = rho * eint + 0.5 * rho * (v1 * v1 + v2 * v2 + v3 * v3)
#if PRJ_MHD
        + 0.5 * (b1 * b1 + b2 * b2 + b3 * b3)
#endif
        ;
    U[PRJ_CONS_YE] = rho * W[PRJ_PRIM_YE];
#if PRJ_MHD
    U[PRJ_CONS_B1] = b1;
    U[PRJ_CONS_B2] = b2;
    U[PRJ_CONS_B3] = b3;
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
#if PRJ_MHD
    double magnetic;
#endif

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
        W[PRJ_PRIM_B1] = 0.0;
        W[PRJ_PRIM_B2] = 0.0;
        W[PRJ_PRIM_B3] = 0.0;
#endif
        prj_rad_cons2prim(U, W);
        return;
    }

    v1 = U[PRJ_CONS_MOM1] / rho;
    v2 = U[PRJ_CONS_MOM2] / rho;
    v3 = U[PRJ_CONS_MOM3] / rho;
    kinetic = 0.5 * (v1 * v1 + v2 * v2 + v3 * v3);
#if PRJ_MHD
    magnetic = 0.5 * (U[PRJ_CONS_B1] * U[PRJ_CONS_B1] +
        U[PRJ_CONS_B2] * U[PRJ_CONS_B2] +
        U[PRJ_CONS_B3] * U[PRJ_CONS_B3]) / rho;
#endif

    W[PRJ_PRIM_RHO] = rho;
    W[PRJ_PRIM_V1] = v1;
    W[PRJ_PRIM_V2] = v2;
    W[PRJ_PRIM_V3] = v3;
    W[PRJ_PRIM_EINT] = U[PRJ_CONS_ETOT] / rho - kinetic
#if PRJ_MHD
        - magnetic
#endif
        ;
    W[PRJ_PRIM_YE] = U[PRJ_CONS_YE] / rho;
#if PRJ_MHD
    W[PRJ_PRIM_B1] = U[PRJ_CONS_B1];
    W[PRJ_PRIM_B2] = U[PRJ_CONS_B2];
    W[PRJ_PRIM_B3] = U[PRJ_CONS_B3];
#endif
    prj_rad_cons2prim(U, W);
}
