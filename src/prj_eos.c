#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prj.h"

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

static double prj_eos_gamma_value(void)
{
    return 5.0 / 3.0;
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

static double prj_eos_tab_pressure_log(const prj_eos *eos, int jy, int jr, int jt)
{
    double p = prj_eos_tab_elem(eos, PRJ_EOS_REC_PRESSURE, jy, jr, jt);

    return log(p > 0.0 ? p : 1.0e-300);
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

static int prj_eos_prepare_table(prj_eos *eos)
{
    FILE *f;
    size_t slab_size;
    size_t expected_bytes;
    int rec;

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
    slab_size = (size_t)eos->nt * (size_t)eos->nr * (size_t)eos->ny;
    expected_bytes = (size_t)PRJ_EOS_NUMEL * slab_size * sizeof(*eos->table);

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

    eos->table_loaded = 1;
    eos->table_is_mmap = 0;
    eos->table_bytes = expected_bytes;
    return 0;
}

void prj_eos_init(prj_eos *eos)
{
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
        eos->table_bytes = 0;
        eos->table = 0;
        return;
    }
    (void)prj_eos_prepare_table(eos);
}

static void prj_eos_table_interp_base(const prj_eos *eos, double rho, double T, double ye,
    int *restrict jy, int *restrict jyp, int *restrict jr, int *restrict jrp, int *restrict jt, int *restrict jtp,
    double *restrict dye, double *restrict drho, double *restrict dtemp)
{
    double rl;
    double tl;
    double dy;
    double dr;
    double dt;
    double yval;
    double tval;

    yval = prj_eos_clamp_double(ye, eos->y1c, eos->y2c);
    tval = prj_eos_clamp_double(T, pow(10.0, eos->t1), pow(10.0, eos->t2));
    rl = log10(prj_eos_clamp_double(rho, pow(10.0, eos->r1), pow(10.0, eos->r2)));
    tl = log10(tval);
    dr = (eos->r2 - eos->r1) / (double)(eos->nr - 1);
    dt = (eos->t2 - eos->t1) / (double)(eos->nt - 1);
    dy = (eos->y2c - eos->y1c) / (double)(eos->ny - 1);

    *jr = 1 + (int)((rl - eos->r1) / dr);
    *jt = 1 + (int)((tl - eos->t1) / dt);
    *jy = 1 + (int)((yval - eos->y1c) / dy);
    *jr = *jr < 1 ? 1 : (*jr >= eos->nr ? eos->nr - 1 : *jr);
    *jt = *jt < 1 ? 1 : (*jt >= eos->nt ? eos->nt - 1 : *jt);
    *jy = *jy < 1 ? 1 : (*jy >= eos->ny ? eos->ny - 1 : *jy);
    *jrp = *jr + 1;
    *jtp = *jt + 1;
    *jyp = *jy + 1;
    *drho = (rl - (eos->r1 + (double)(*jr - 1) * dr)) / dr;
    *dtemp = (tl - (eos->t1 + (double)(*jt - 1) * dt)) / dt;
    *dye = (yval - (eos->y1c + (double)(*jy - 1) * dy)) / dy;
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
    double c000 = prj_eos_tab_pressure_log(eos, jy, jr, jt);
    double c100 = prj_eos_tab_pressure_log(eos, jyp, jr, jt);
    double c010 = prj_eos_tab_pressure_log(eos, jy, jrp, jt);
    double c110 = prj_eos_tab_pressure_log(eos, jyp, jrp, jt);
    double c001 = prj_eos_tab_pressure_log(eos, jy, jr, jtp);
    double c101 = prj_eos_tab_pressure_log(eos, jyp, jr, jtp);
    double c011 = prj_eos_tab_pressure_log(eos, jy, jrp, jtp);
    double c111 = prj_eos_tab_pressure_log(eos, jyp, jrp, jtp);
    double c00 = (1.0 - dye) * c000 + dye * c100;
    double c10 = (1.0 - dye) * c010 + dye * c110;
    double c01 = (1.0 - dye) * c001 + dye * c101;
    double c11 = (1.0 - dye) * c011 + dye * c111;
    double c0 = (1.0 - drho) * c00 + drho * c10;
    double c1 = (1.0 - drho) * c01 + drho * c11;

    return (1.0 - dtemp) * c0 + dtemp * c1;
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
        double dr;
        double dy;
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
        double tlo;
        double thi;
        double dtemp;
        double pressure_log;

        e_table = eint / PRJ_EOS_ENERGY_SCALE;
        ye = prj_eos_clamp_double(ye, eos->y1c, eos->y2c);
        rl = log10(prj_eos_clamp_double(rho, pow(10.0, eos->r1), pow(10.0, eos->r2)));
        dr = (eos->r2 - eos->r1) / (double)(eos->nr - 1);
        dy = (eos->y2c - eos->y1c) / (double)(eos->ny - 1);
        jrf = (rl - eos->r1) / dr;
        jyf = (ye - eos->y1c) / dy;
        jr = 1 + (int)jrf;
        jy = 1 + (int)jyf;
        jr = jr < 1 ? 1 : (jr >= eos->nr ? eos->nr - 1 : jr);
        jy = jy < 1 ? 1 : (jy >= eos->ny ? eos->ny - 1 : jy);
        jrp = jr + 1;
        jyp = jy + 1;
        drho = prj_eos_clamp_double((rl - (eos->r1 + (double)(jr - 1) * dr)) / dr, 0.0, 1.0);
        dye = prj_eos_clamp_double((ye - (eos->y1c + (double)(jy - 1) * dy)) / dy, 0.0, 1.0);
        coeff0 = (1.0 - drho) * (1.0 - dye);
        coeff1 = (1.0 - drho) * dye;
        coeff2 = drho * (1.0 - dye);
        coeff3 = drho * dye;

        {
            double temp_guess = prj_eos_clamp_double((prj_eos_gamma_value() - 1.0) * eint,
                pow(10.0, eos->t1), pow(10.0, eos->t2));

            jt = 1 + (int)((log10(temp_guess) - eos->t1) /
                ((eos->t2 - eos->t1) / (double)(eos->nt - 1)));
        }
        jt = jt < 1 ? 1 : (jt >= eos->nt ? eos->nt - 1 : jt);

        elo = coeff0 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy, jr, jt) +
            coeff1 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jr, jt) +
            coeff2 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy, jrp, jt) +
            coeff3 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jrp, jt);
        while (elo > e_table && jt > 1) {
            jt -= 1;
            elo = coeff0 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy, jr, jt) +
                coeff1 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jr, jt) +
                coeff2 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy, jrp, jt) +
                coeff3 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jrp, jt);
        }
        ehi = elo;
        while (ehi < e_table && jt < eos->nt - 1) {
            jt += 1;
            ehi = coeff0 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy, jr, jt) +
                coeff1 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jr, jt) +
                coeff2 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy, jrp, jt) +
                coeff3 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jrp, jt);
        }
        if (jt > 1) {
            jt -= 1;
            elo = coeff0 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy, jr, jt) +
                coeff1 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jr, jt) +
                coeff2 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy, jrp, jt) +
                coeff3 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jrp, jt);
            ehi = coeff0 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy, jr, jt + 1) +
                coeff1 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jr, jt + 1) +
                coeff2 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jy, jrp, jt + 1) +
                coeff3 * prj_eos_tab_elem(eos, PRJ_EOS_REC_EINT, jyp, jrp, jt + 1);
        }
        jtp = jt + 1;
        if (fabs(ehi - elo) > 0.0) {
            dtemp = prj_eos_clamp_double((e_table - elo) / (ehi - elo), 0.0, 1.0);
        } else {
            dtemp = 0.0;
        }
        tlo = eos->t1 + (double)(jt - 1) * (eos->t2 - eos->t1) / (double)(eos->nt - 1);
        thi = eos->t1 + (double)(jtp - 1) * (eos->t2 - eos->t1) / (double)(eos->nt - 1);
        T = pow(10.0, tlo + (thi - tlo) * dtemp);
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

void prj_eos_prim2cons(prj_eos *eos, double *W, double *U)
{
    double rho;
    double v1;
    double v2;
    double v3;
    double eint;

    (void)eos;

    if (W == 0 || U == 0) {
        return;
    }

    rho = W[PRJ_PRIM_RHO];
    v1 = W[PRJ_PRIM_V1];
    v2 = W[PRJ_PRIM_V2];
    v3 = W[PRJ_PRIM_V3];
    eint = W[PRJ_PRIM_EINT];

    U[PRJ_CONS_RHO] = rho;
    U[PRJ_CONS_MOM1] = rho * v1;
    U[PRJ_CONS_MOM2] = rho * v2;
    U[PRJ_CONS_MOM3] = rho * v3;
    U[PRJ_CONS_ETOT] = rho * eint + 0.5 * rho * (v1 * v1 + v2 * v2 + v3 * v3);
    U[PRJ_CONS_YE] = rho * W[PRJ_PRIM_YE];
    prj_rad_prim2cons(W, U);
}

void prj_eos_cons2prim(prj_eos *eos, double *U, double *W)
{
    double rho;
    double v1;
    double v2;
    double v3;
    double kinetic;

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
        prj_rad_cons2prim(U, W);
        return;
    }

    v1 = U[PRJ_CONS_MOM1] / rho;
    v2 = U[PRJ_CONS_MOM2] / rho;
    v3 = U[PRJ_CONS_MOM3] / rho;
    kinetic = 0.5 * (v1 * v1 + v2 * v2 + v3 * v3);

    W[PRJ_PRIM_RHO] = rho;
    W[PRJ_PRIM_V1] = v1;
    W[PRJ_PRIM_V2] = v2;
    W[PRJ_PRIM_V3] = v3;
    W[PRJ_PRIM_EINT] = U[PRJ_CONS_ETOT] / rho - kinetic;
    W[PRJ_PRIM_YE] = U[PRJ_CONS_YE] / rho;
    prj_rad_cons2prim(U, W);
}
