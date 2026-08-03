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
    double test_rho = 1.0;
    prj_rad3_opac_interp_result interp;
    int q;
    memset(&rad, 0, sizeof(rad));
#if PRJ_RAD_MICROPHYSICS == PRJ_RAD_MICROPHYSICS_CONSTANT_ABSORPTION
    rad.emin[0] = 1.0e-8;
    rad.emax[0] = 1.0e-2;
#elif PRJ_RAD_MICROPHYSICS == PRJ_RAD_MICROPHYSICS_BENCHMARK_LTE && PRJ_RAD_FERMI_DIRAC
    rad.emin[0] = 1.0;
    rad.emax[0] = 50.0;
#else
    rad.emin[0] = 1.0e-14;
    rad.emax[0] = 1.0e-5;
#endif
    prj_rad3_opac_init(&rad);
#if PRJ_EOS_T3
    test_temp = 1000.0;
#elif PRJ_RAD_MICROPHYSICS == PRJ_RAD_MICROPHYSICS_CONSTANT_ABSORPTION
    test_temp = 8.617333262145e-5;
#elif PRJ_RAD_MICROPHYSICS == PRJ_RAD_MICROPHYSICS_BENCHMARK_LTE && PRJ_RAD_FERMI_DIRAC
    test_temp = 5.0;
#else
    test_temp = 8.617333262145e-8;
#endif
#if PRJ_RAD_MICROPHYSICS == PRJ_RAD_MICROPHYSICS_BENCHMARK_LTE
    test_rho = PRJ_RAD_TEST_RHO_CUTOFF + 9.0e14;
#endif
    prj_rad3_opac_lookup(&rad, test_rho, test_temp, 0.5, k, s, d, e);
    prj_rad3_opac_lookup_ke(&rad, test_rho, test_temp, 0.5,
        k, e, dkt, dky, det, dey);
    prj_rad3_opac_lookup_interp(&rad, test_rho, test_temp, 0.5, &interp);
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
#elif PRJ_RAD_MICROPHYSICS == PRJ_RAD_MICROPHYSICS_BENCHMARK_LTE
    assert(fabs(k[0] - PRJ_RAD_TEST_ABSORPTION) < 1.0e-14);
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
        double lte_tol = 5.0e-3;
#if PRJ_RAD_FERMI_DIRAC
        expected *= 7.0 / 8.0;
        lte_tol = 1.5e-2;
#endif
        assert(fabs(lte_sum / expected - 1.0) < lte_tol);
    }
#else
    (void)lte_sum;
#endif
#if PRJ_SPHERICAL_1D && PRJ_NRAD > 0
    {
        prj_block block;
        double div[PRJ_NVAR_CONS];
        int d;
        memset(&block, 0, sizeof(block));
        block.dx[0] = block.dx[1] = block.dx[2] = 1.0;
        block.area[0] = block.area[1] = block.area[2] = 1.0;
        block.vol = 1.0;
        for (d = 0; d < 3; ++d)
            block.flux[d] = calloc((size_t)PRJ_NVAR_CONS * PRJ_BLOCK_NFACES,
                sizeof(*block.flux[d]));
        block.flux[0][VIDX(PRJ_CONS_RAD_E(0, 0), 1, 0, 0)] = 1.0;
        prj_flux_div(&block, 0, 0, 0, div);
        assert(fabs(div[PRJ_CONS_RAD_E(0, 0)] + 3.0) < 1.0e-14);
        for (d = 0; d < 3; ++d) free(block.flux[d]);
        assert(fabs(prj_src_spherical_average_inv_r(0.0, 1.0) - 1.5) < 1.0e-14);
    }
#endif
#if PRJ_USER_STATIC_METRIC
    {
        double pos[3] = {1.0e6, 0.0, 0.0};
        double alpha, grr, grad;
        prj_problem_static_metric(pos, &alpha, &grr, &grad);
        assert(alpha > 0.0 && alpha < 1.0 && grr > 1.0 && grad > 0.0);
    }
#endif
    {
        prj_eos eos;
        double q[PRJ_EOS_NQUANT];
        memset(&eos, 0, sizeof(eos));
        prj_eos_rty(&eos, 2.0, 1.0, 0.5, q, PRJ_EOS_CTX_MAIN);
        assert(fabs(q[PRJ_EOS_GAMMA] - PRJ_IDEAL_GAMMA) < 1.0e-14);
    }
#if PRJ_PRESCRIBED_MATTER
    {
        prj_eos eos;
        double u[PRJ_NVAR_CONS] = {0.0};
        double hydro[PRJ_NHYDRO];
        double final_temp;
        int v;
        memset(&eos, 0, sizeof(eos));
        u[PRJ_CONS_RHO] = test_rho;
        u[PRJ_CONS_ETOT] = test_rho * 1.0e18;
        u[PRJ_CONS_YE] = 0.5 * test_rho;
        for (v = 0; v < PRJ_NHYDRO; ++v) hydro[v] = u[v];
        prj_rad_energy_update(&rad, &eos, u, 1.0e-12, 1.0, &final_temp, 0);
        for (v = 0; v < PRJ_NHYDRO; ++v) assert(u[v] == hydro[v]);
        assert(isfinite(final_temp) && u[PRJ_CONS_RAD_E(0, 0)] >= 0.0);
    }
#endif
    {
        double ux = 0.69;
        double lor = sqrt(1.0 + ux * ux);
        double beta = ux / lor;
        double f = 4.0 * beta / (3.0 + beta * beta);
        double recovered = (2.0 - sqrt(4.0 - 3.0 * f * f)) / f;
        assert(fabs(recovered - beta) < 1.0e-14);
    }
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
