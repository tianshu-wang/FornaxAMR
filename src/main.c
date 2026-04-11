#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <errno.h>

#include "prj.h"

static prj_problem_init_fn prj_select_problem(const char *name)
{
    if (strcmp(name, "sedov") == 0) {
        return prj_problem_sedov;
    }
    if (strcmp(name, "cc") == 0) {
        return prj_problem_cc;
    }
    if (strcmp(name, "sedov_offcenter") == 0) {
        return prj_problem_sedov_offcenter;
    }
    if (strcmp(name, "shock1d") == 0) {
        return prj_problem_shock1d;
    }
    return prj_problem_general;
}

static void prj_copy_file(const char *src, const char *dst)
{
    FILE *fin = fopen(src, "rb");
    FILE *fout;
    char buffer[4096];
    size_t nread;

    if (fin == 0) {
        return;
    }
    fout = fopen(dst, "wb");
    if (fout == 0) {
        fclose(fin);
        return;
    }
    while ((nread = fread(buffer, 1, sizeof(buffer), fin)) > 0) {
        fwrite(buffer, 1, nread, fout);
    }
    fclose(fin);
    fclose(fout);
}

int main(int argc, char *argv[])
{
    prj_sim sim;
    prj_problem_init_fn init_fn = prj_problem_general;
    const char *restart_file = 0;
    int resolution = -1;
    int max_level_override = -1;
    int i;

    memset(&sim, 0, sizeof(sim));
    sim.mesh.max_level = -1;
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--problem") == 0 && i + 1 < argc) {
            init_fn = prj_select_problem(argv[++i]);
        } else if (strcmp(argv[i], "--resolution") == 0 && i + 1 < argc) {
            resolution = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--max_level") == 0 && i + 1 < argc) {
            max_level_override = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--restart") == 0 && i + 1 < argc) {
            restart_file = argv[++i];
        }
    }
    if (resolution > 0) {
        sim.mesh.root_nx[0] = resolution;
        sim.mesh.root_nx[1] = resolution;
        sim.mesh.root_nx[2] = resolution;
    }
    if (max_level_override >= 0) {
        sim.mesh.max_level = max_level_override;
    }

    init_fn(&sim);
    mkdir("output", 0777);
    if (restart_file != 0) {
        prj_mesh_destroy(&sim.mesh);
        prj_io_read_restart(&sim.mesh, &sim.eos, restart_file, &sim.time, &sim.step);
    }
    prj_rad_init(&sim.rad);
    prj_gravity_init(&sim);
    prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
    prj_io_write_dump(&sim.mesh, sim.output_dir, sim.step);

    while (sim.time < sim.t_end && sim.step < sim.max_steps) {
        sim.dt = prj_timeint_calc_dt(&sim.mesh, &sim.eos, sim.cfl);
        if (sim.time + sim.dt > sim.t_end) {
            sim.dt = sim.t_end - sim.time;
        }
        prj_timeint_step(&sim.mesh, &sim.coord, &sim.bc, &sim.eos, &sim.rad, sim.dt);
        sim.time += sim.dt;
        sim.step += 1;
        if (sim.amr_interval > 0 && sim.step % sim.amr_interval == 0) {
            prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
            prj_amr_adapt(&sim.mesh, &sim.eos);
            prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
        }
        if (sim.restart_interval > 0 && sim.step % sim.restart_interval == 0) {
            prj_io_write_restart(&sim.mesh, sim.time, sim.step);
        }
        if (sim.output_interval > 0 && sim.step % sim.output_interval == 0) {
            prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
            prj_io_write_dump(&sim.mesh, sim.output_dir, sim.step);
        }
        printf("step=%d  t=%.6e  dt=%.6e  blocks=%d\n",
            sim.step, sim.time, sim.dt, prj_mesh_count_active(&sim.mesh));
    }

    prj_io_write_restart(&sim.mesh, sim.time, sim.step);
    {
        char final_restart[64];

        snprintf(final_restart, sizeof(final_restart), "output/restart_%08d.h5", sim.step);
        prj_copy_file(final_restart, "output/final.h5");
    }
    prj_gravity_free(&sim.monopole_grav);
    prj_mesh_destroy(&sim.mesh);
    return 0;
}
