#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>

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

static const char *prj_eos_label(const prj_sim *sim)
{
    if (sim == 0) {
        return "unknown";
    }
    if (sim->eos.kind == PRJ_EOS_KIND_TABLE && sim->eos.filename[0] != '\0') {
        return sim->eos.filename;
    }
    if (sim->eos.kind == PRJ_EOS_KIND_TABLE) {
        return "table";
    }
    return "ideal";
}

static const char *prj_amr_label(const prj_sim *sim)
{
    if (sim == 0) {
        return "off";
    }
    if (sim->mesh.max_level > 0 && sim->amr_interval > 0) {
        return "on";
    }
    return "off";
}

static const char *prj_amr_estimator_label(const prj_sim *sim)
{
    if (sim == 0) {
        return "unknown";
    }
    {
        static char label[128];
        int offset = 0;
        int i;
        int any = 0;

        label[0] = '\0';
        for (i = 0; i < PRJ_AMR_N; ++i) {
            const char *name;

            if (sim->mesh.amr_criterion_set[i] == 0) {
                continue;
            }
            if (sim->mesh.amr_estimator[i] == PRJ_AMR_ESTIMATOR_PRESSURE_SCALE_HEIGHT) {
                name = "pressure_scale_height";
            } else if (sim->mesh.amr_estimator[i] == PRJ_AMR_ESTIMATOR_PRESSURE_JUMP) {
                name = "pressure_jump";
            } else if (sim->mesh.amr_estimator[i] == PRJ_AMR_ESTIMATOR_VELOCITY) {
                name = "velocity";
            } else {
                name = "lohner";
            }
            offset += snprintf(label + offset, sizeof(label) - (size_t)offset,
                "%s%s", any ? "," : "", name);
            any = 1;
            if (offset >= (int)sizeof(label) - 1) {
                break;
            }
        }
        if (any != 0) {
            return label;
        }
    }
    return "none";
}

static void prj_print_config(const prj_sim *sim, int rank)
{
    if (rank != 0 || sim == 0) {
        return;
    }

    fprintf(stderr, "config:\n");
    fprintf(stderr, "mpi: %s\n",
#if defined(PRJ_ENABLE_MPI)
        "on"
#else
        "off"
#endif
    );
    fprintf(stderr, "radiation: %s\n",
#if PRJ_USE_RADIATION
        "on"
#else
        "off"
#endif
    );
    fprintf(stderr, "gravity: %s\n",
#if PRJ_USE_GRAVITY
        "on"
#else
        "off"
#endif
    );
    fprintf(stderr, "eos: %s\n",
        prj_eos_label(sim)
    );
    fprintf(stderr, "amr: %s\n",
        prj_amr_label(sim)
    );
    fprintf(stderr, "amr estimator: %s\n",
        prj_amr_estimator_label(sim)
    );
}

int main(int argc, char *argv[])
{
    prj_sim sim;
    prj_mpi mpi;
    prj_problem_init_fn init_fn;
    int init_with_mpi = 0;
    char *param_file = 0;
    double saved_amr_refine_thresh[PRJ_AMR_N];
    double saved_amr_derefine_thresh[PRJ_AMR_N];
    double saved_amr_eps;
    int saved_amr_estimator[PRJ_AMR_N];
    int saved_amr_criterion_set[PRJ_AMR_N];
    int saved_use_amr_angle_resolution;
    double saved_amr_angle_resolution_limit;
    int resolution = -1;
    int max_level_override = -1;
    double next_output_time = -1.0;
    double next_restart_time = -1.0;
    int i;

    memset(&sim, 0, sizeof(sim));
    prj_io_parser(&sim, 0);
    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--param") == 0 && i + 1 < argc) {
            param_file = argv[++i];
        } else if (strcmp(argv[i], "--resolution") == 0 && i + 1 < argc) {
            resolution = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--max_level") == 0 && i + 1 < argc) {
            max_level_override = atoi(argv[++i]);
        }
    }
    if (param_file != 0) {
        prj_io_parser(&sim, param_file);
    } else {
        fprintf(stderr, "param_file is missing\n");
        return 1;
    }
    init_fn = prj_select_problem(sim.problem_name);
    if (sim.restart_from_file != 0 && sim.restart_file_name[0] == '\0') {
        fprintf(stderr, "restart_from_file is enabled but restart_file_name is empty\n");
        return 1;
    }
    if (resolution > 0) {
        sim.mesh.root_nx[0] = resolution;
        sim.mesh.root_nx[1] = resolution;
        sim.mesh.root_nx[2] = resolution;
    }
    if (max_level_override >= 0) {
        sim.mesh.max_level = max_level_override;
    }

    init_with_mpi = (init_fn == prj_problem_cc);
    prj_mpi_init(&argc, &argv, &mpi);
    prj_print_config(&sim, mpi.rank);
    init_fn(&sim);
    if (!init_with_mpi) {
        prj_mpi_decompose(&sim.mesh);
        prj_mpi_prepare(&sim.mesh, &mpi);
    }
    if (mpi.rank == 0) {
        mkdir("output", 0777);
    }
    if (sim.restart_from_file != 0) {
        for (i = 0; i < PRJ_AMR_N; ++i) {
            saved_amr_refine_thresh[i] = sim.mesh.amr_refine_thresh[i];
            saved_amr_derefine_thresh[i] = sim.mesh.amr_derefine_thresh[i];
            saved_amr_estimator[i] = sim.mesh.amr_estimator[i];
            saved_amr_criterion_set[i] = sim.mesh.amr_criterion_set[i];
        }
        saved_amr_eps = sim.mesh.amr_eps;
        saved_use_amr_angle_resolution = sim.mesh.use_amr_angle_resolution;
        saved_amr_angle_resolution_limit = sim.mesh.amr_angle_resolution_limit;
        prj_mesh_destroy(&sim.mesh);
        prj_io_read_restart(&sim.mesh, &sim.eos, sim.restart_file_name, &sim.time, &sim.step);
        for (i = 0; i < PRJ_AMR_N; ++i) {
            sim.mesh.amr_refine_thresh[i] = saved_amr_refine_thresh[i];
            sim.mesh.amr_derefine_thresh[i] = saved_amr_derefine_thresh[i];
            sim.mesh.amr_estimator[i] = saved_amr_estimator[i];
            sim.mesh.amr_criterion_set[i] = saved_amr_criterion_set[i];
        }
        sim.mesh.amr_eps = saved_amr_eps;
        sim.mesh.use_amr_angle_resolution = saved_use_amr_angle_resolution;
        sim.mesh.amr_angle_resolution_limit = saved_amr_angle_resolution_limit;
        prj_print_config(&sim, mpi.rank);
    }
    prj_mesh_mark_base_blocks(&sim.mesh);
    prj_rad_init(&sim.rad);
 #if PRJ_USE_GRAVITY
    prj_gravity_init(&sim);
 #endif
    prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
    prj_eos_fill_mesh(&sim.mesh, &sim.eos, 1);
    prj_io_write_dump(&sim.mesh, sim.output_dir, sim.step);
    if (sim.output_dt >= 0.0) {
        next_output_time = sim.time + sim.output_dt;
    }
    if (sim.restart_dt >= 0.0) {
        next_restart_time = sim.time + sim.restart_dt;
    }

    struct timeval wall_start;
    gettimeofday(&wall_start, 0);

    while (sim.time < sim.t_end && sim.step < sim.max_steps) {
        int write_output = 0;
        int write_restart = 0;
        double next_event_time = -1.0;

        sim.dt = prj_timeint_calc_dt(&sim.mesh, &sim.eos, sim.cfl);
        if (sim.output_dt >= 0.0 && next_output_time >= 0.0) {
            next_event_time = next_output_time;
        }
        if (sim.restart_dt >= 0.0 && next_restart_time >= 0.0 &&
            (next_event_time < 0.0 || next_restart_time < next_event_time)) {
            next_event_time = next_restart_time;
        }
        if (next_event_time >= 0.0 && sim.time + sim.dt > next_event_time) {
            sim.dt = 1.000001 * (next_event_time - sim.time);
        }
        if (sim.time + sim.dt > sim.t_end) {
            sim.dt = sim.t_end - sim.time;
        }
        prj_timeint_step(&sim.mesh, &sim.coord, &sim.bc, &sim.eos, &sim.rad, sim.dt);
        sim.time += sim.dt;
        sim.step += 1;
        if (sim.amr_interval > 0 && sim.step % sim.amr_interval == 0) {
            prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
            prj_eos_fill_mesh(&sim.mesh, &sim.eos, 1);
            prj_amr_adapt(&sim.mesh, &sim.eos);
            prj_mpi_rebalance(&sim.mesh);
            prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
            prj_eos_fill_mesh(&sim.mesh, &sim.eos, 1);
#if PRJ_USE_GRAVITY
            prj_gravity_rebuild_grid(&sim);
#endif
        }
        if (sim.output_interval > 0 && sim.step % sim.output_interval == 0) {
            write_output = 1;
        }
        if (sim.restart_interval > 0 && sim.step % sim.restart_interval == 0) {
            write_restart = 1;
        }
        if (sim.output_dt >= 0.0 && next_output_time >= 0.0 && sim.time >= next_output_time) {
            write_output = 1;
            if (sim.output_dt > 0.0) {
                do {
                    next_output_time += sim.output_dt;
                } while (sim.time >= next_output_time);
            } else {
                next_output_time = -1.0;
            }
        }
        if (sim.restart_dt >= 0.0 && next_restart_time >= 0.0 && sim.time >= next_restart_time) {
            write_restart = 1;
            if (sim.restart_dt > 0.0) {
                do {
                    next_restart_time += sim.restart_dt;
                } while (sim.time >= next_restart_time);
            } else {
                next_restart_time = -1.0;
            }
        }
        if (write_output) {
            prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
            prj_eos_fill_mesh(&sim.mesh, &sim.eos, 1);
            prj_io_write_dump(&sim.mesh, sim.output_dir, sim.step);
        }
        if (write_restart) {
            prj_io_write_restart(&sim.mesh, sim.time, sim.step);
        }
        if (mpi.rank == 0) {
            struct timeval wall_now;
            double wall_elapsed;

            gettimeofday(&wall_now, 0);
            wall_elapsed = (double)(wall_now.tv_sec - wall_start.tv_sec) +
                1.0e-6 * (double)(wall_now.tv_usec - wall_start.tv_usec);
            fprintf(stderr, "step=%d  t=%.6e  dt=%.6e  blocks=%d  wall=%.3fs\n",
                sim.step, sim.time, sim.dt, prj_mesh_count_active(&sim.mesh), wall_elapsed);
        }
    }

    prj_io_write_restart(&sim.mesh, sim.time, sim.step);
    if (mpi.rank == 0) {
        char final_restart[64];

        snprintf(final_restart, sizeof(final_restart), "output/restart_%08d.h5", sim.step);
        prj_copy_file(final_restart, "output/final.h5");
    }
#if PRJ_USE_GRAVITY
    prj_gravity_free(&sim.monopole_grav);
#endif
    prj_mesh_destroy(&sim.mesh);
    prj_mpi_finalize();
    return 0;
}
