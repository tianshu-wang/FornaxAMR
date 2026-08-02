#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "prj.h"

int main(void)
{
#if PRJ_RAD_TEST_PROBLEM != 0
    assert(PRJ_USE_INELASTIC_SCATTERING == 0);
#endif
#if PRJ_NRAD > 0 && PRJ_RAD_MICROPHYSICS != PRJ_RAD_MICROPHYSICS_TABLE
    prj_rad rad;
    double k[PRJ_RAD3_OPAC_NGROUPS], s[PRJ_RAD3_OPAC_NGROUPS];
    double d[PRJ_RAD3_OPAC_NGROUPS], e[PRJ_RAD3_OPAC_NGROUPS];
    double dkt[PRJ_RAD3_OPAC_NGROUPS], dky[PRJ_RAD3_OPAC_NGROUPS];
    double det[PRJ_RAD3_OPAC_NGROUPS], dey[PRJ_RAD3_OPAC_NGROUPS];
    double lte_sum = 0.0;
    double test_temp;
    prj_rad3_opac_interp_result interp;
    int q;
    memset(&rad, 0, sizeof(rad));
#if PRJ_RAD_MICROPHYSICS == PRJ_RAD_MICROPHYSICS_CONSTANT_ABSORPTION
    rad.emin[0] = 1.0e-8;
    rad.emax[0] = 1.0e-2;
#else
    rad.emin[0] = 1.0e-14;
    rad.emax[0] = 1.0e-5;
#endif
    prj_rad3_opac_init(&rad);
#if PRJ_EOS_T3
    test_temp = 1000.0;
#elif PRJ_RAD_MICROPHYSICS == PRJ_RAD_MICROPHYSICS_CONSTANT_ABSORPTION
    test_temp = 8.617333262145e-5;
#else
    test_temp = 8.617333262145e-8;
#endif
    prj_rad3_opac_lookup(&rad, 1.0, test_temp, 0.5, k, s, d, e);
    prj_rad3_opac_lookup_ke(&rad, 1.0, test_temp, 0.5,
        k, e, dkt, dky, det, dey);
    prj_rad3_opac_lookup_interp(&rad, 1.0, test_temp, 0.5, &interp);
    for (q = 0; q < PRJ_RAD3_OPAC_NGROUPS; ++q) {
        assert(isfinite(k[q]) && isfinite(s[q]) && isfinite(d[q]) && isfinite(e[q]));
        assert(isfinite(dkt[q]) && isfinite(dky[q]) && isfinite(det[q]) && isfinite(dey[q]));
        assert(k[q] >= 0.0 && s[q] >= 0.0 && e[q] >= 0.0);
        assert(dkt[q] == 0.0 && dky[q] == 0.0 && dey[q] == 0.0);
        if (k[q] > 0.0) lte_sum += e[q] / (2.99792458e10 * k[q]);
    }
#if PRJ_RAD_MICROPHYSICS == PRJ_RAD_MICROPHYSICS_CONSTANT_ABSORPTION
    assert(k[0] == 1.0 && s[0] == 0.0);
#elif PRJ_RAD_MICROPHYSICS == PRJ_RAD_MICROPHYSICS_CONSTANT_SCATTERING
    assert(k[0] == 0.0 && s[0] == 2.5e-6 && e[0] == 0.0);
#endif
    {
        double dk[3], ds[3], dd[3], de[3];
        assert(prj_rad3_opac_interp_group_derivs(&interp, 0, dk, ds, dd, de));
        for (q = 0; q < 3; ++q)
            assert(isfinite(dk[q]) && isfinite(ds[q]) && isfinite(dd[q]) && isfinite(de[q]));
#if PRJ_RAD_MICROPHYSICS == PRJ_RAD_MICROPHYSICS_CONSTANT_ABSORPTION
        assert(fabs(dk[0] - 1.0) < 1.0e-14);
#else
        assert(fabs(dk[0]) < 1.0e-14);
#endif
    }
#if PRJ_RAD_MICROPHYSICS == PRJ_RAD_MICROPHYSICS_PICKET_FENCE
    assert(k[0] == 2.0 && k[1] == 20.0);
#endif
#if PRJ_RAD_MICROPHYSICS != PRJ_RAD_MICROPHYSICS_CONSTANT_SCATTERING
    {
        double temp_k = PRJ_EOS_T3 ? test_temp : test_temp / 8.617333262145e-11;
        double expected = 7.5657e-15 * pow(temp_k, 4.0) / RAD_SCALE;
        assert(fabs(lte_sum / expected - 1.0) < 5.0e-3);
    }
#else
    (void)lte_sum;
#endif
    prj_rad3_opac_free(&rad);

    {
        prj_mesh mesh;
        prj_block block;
        prj_bc bc;
        size_t count = (size_t)PRJ_NVAR_PRIM * PRJ_BLOCK_NSTAGES * PRJ_BLOCK_NCELLS;
        double *storage = calloc(count, sizeof(*storage));
        int fv = PRJ_PRIM_RAD_F1(0, 0);
        assert(storage != 0);
        memset(&mesh, 0, sizeof(mesh)); memset(&block, 0, sizeof(block)); memset(&bc, 0, sizeof(bc));
        mesh.coord.x1min = 0.0; mesh.coord.x1max = 1.0;
        mesh.coord.x2min = 0.0; mesh.coord.x2max = 1.0;
        mesh.coord.x3min = 0.0; mesh.coord.x3max = 1.0;
        block.xmin[0] = block.xmin[1] = block.xmin[2] = 0.0;
        block.xmax[0] = block.xmax[1] = block.xmax[2] = 1.0;
        block.W_mhd = storage;
        block.W_rad = storage + (size_t)PRJ_NVAR_MHD_PRIM * PRJ_BLOCK_NSTAGES * PRJ_BLOCK_NCELLS;
        bc.bc_x1_inner = PRJ_BC_REFLECT;
        bc.bc_x1_outer = bc.bc_x2_inner = bc.bc_x2_outer =
            bc.bc_x3_inner = bc.bc_x3_outer = PRJ_BC_PERIODIC;
        storage[WIDX(fv, 0, 0, 0)] = 7.0;
        prj_boundary_physical(&mesh, &bc, &block, 0, PRJ_BOUNDARY_PHYS_FACE_ONLY);
        assert(storage[WIDX(fv, -1, 0, 0)] == -7.0);
        free(storage);
    }
#endif
#if PRJ_EOS_T3
    {
        prj_eos eos;
        double q[PRJ_EOS_NQUANT];
        memset(&eos, 0, sizeof(eos));
        prj_eos_rty(&eos, 1.0, 1000.0, 0.5, q, PRJ_EOS_CTX_MAIN);
        assert(fabs(q[PRJ_EOS_EINT] - 7.5657e-3) < 1.0e-12);
        prj_eos_rey(&eos, 1.0, q[PRJ_EOS_EINT], 0.5, q, PRJ_EOS_CTX_MAIN);
        assert(fabs(q[PRJ_EOS_TEMPERATURE] - 1000.0) < 1.0e-9);
    }
#endif
#if PRJ_RAD_TEST_PROBLEM == 3
    {
        prj_mesh mesh;
        prj_block block;
        size_t nrad = (size_t)PRJ_NVAR_RAD_CONS * PRJ_BLOCK_NCELLS;
        double *rhs = calloc(nrad, sizeof(*rhs));
        size_t n;
        double total = 0.0;
        assert(rhs != 0);
        memset(&mesh, 0, sizeof(mesh)); memset(&block, 0, sizeof(block));
        block.dx[0] = (3.0 / 11.0) / (19.0 * PRJ_BLOCK_SIZE);
        prj_problem_user_source(&mesh, &block, 0, 0, 0, rhs);
        for (n = 0; n < nrad; ++n) total += rhs[n];
        assert(total > 0.0);
        memset(rhs, 0, nrad * sizeof(*rhs));
        mesh.time_seconds = 10.0001 / (2.99792458e10 * 11.0);
        prj_problem_user_source(&mesh, &block, 0, 0, 0, rhs);
        for (n = 0; n < nrad; ++n) assert(rhs[n] == 0.0);
        free(rhs);
    }
#endif
    return 0;
}
