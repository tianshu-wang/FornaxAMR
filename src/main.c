#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>

#include <errno.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#if PRJ_TIMER
static void prj_write_timer_report(const prj_timer *timer, int rank)
{
    char filename[64];
    FILE *fp;

    snprintf(filename, sizeof(filename), "timer-rank%d.txt", rank);
    fp = fopen(filename, "w");
    if (fp == 0) {
        fprintf(stderr, "failed to open %s for timer report: %s\n", filename, strerror(errno));
        return;
    }
    prj_timer_report(timer, fp, rank);
    fclose(fp);
}
#endif

#if PRJ_TIMER
#define PRJ_TIMER_START(timer, name) prj_timer_start((timer), (name))
#define PRJ_TIMER_STOP(timer, name) prj_timer_stop((timer), (name))
#else
#define PRJ_TIMER_START(timer, name) ((void)(timer), (void)(name))
#define PRJ_TIMER_STOP(timer, name) ((void)(timer), (void)(name))
#endif

static prj_problem_init_fn prj_select_problem(const char *name)
{
    if (strcmp(name, "sedov") == 0) {
        return prj_problem_sedov;
    }
    if (strcmp(name, "magnetized_sedov") == 0) {
        return prj_problem_sedov;
    }
    if (strcmp(name, "cc") == 0) {
        return prj_problem_cc;
    }
    if (strcmp(name, "magnetized_cc") == 0) {
        return prj_problem_cc;
    }
    if (strcmp(name, "ccsn") == 0) {
        return prj_problem_ccsn;
    }
    if (strcmp(name, "magnetized_ccsn") == 0) {
        return prj_problem_ccsn;
    }
    if (strcmp(name, "sedov_offcenter") == 0) {
        return prj_problem_sedov_offcenter;
    }
    if (strcmp(name, "shock1d") == 0) {
        return prj_problem_shock1d;
    }
    return prj_problem_general;
}

static void prj_prepare_restart_problem(prj_sim *sim, prj_problem_init_fn init_fn)
{
    static const char *cc_eos_path =
        "../eos_tmp/SFHoEOS__ye__0.035_0.56_50__logT_-4.793_2.176_500__logrho_-8.699_15.5_500_extend.dat";

    if (sim == 0) {
        return;
    }
    if ((init_fn == prj_problem_cc || init_fn == prj_problem_ccsn) &&
        sim->eos.kind == PRJ_EOS_KIND_TABLE &&
        sim->eos.filename[0] == '\0') {
        strncpy(sim->eos.filename, cc_eos_path, sizeof(sim->eos.filename) - 1);
        sim->eos.filename[sizeof(sim->eos.filename) - 1] = '\0';
    }
    prj_eos_init(&sim->eos);
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

static double prj_last_event_time(double next_event_time, double event_dt)
{
    if (event_dt < 0.0 || next_event_time < 0.0) {
        return -1.0;
    }
    return next_event_time - event_dt;
}

static double prj_next_event_time(double last_event_time, double event_dt, double current_time)
{
    double next_event_time;

    if (event_dt < 0.0) {
        return -1.0;
    }
    if (last_event_time < 0.0) {
        return current_time + event_dt;
    }
    if (event_dt == 0.0) {
        return -1.0;
    }

    next_event_time = last_event_time + event_dt;
    while (current_time >= next_event_time) {
        next_event_time += event_dt;
    }
    return next_event_time;
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
    if (sim->mesh.max_level != 0 && sim->amr_interval > 0) {
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
            } else if (sim->mesh.amr_estimator[i] == PRJ_AMR_ESTIMATOR_DENSITY_JUMP) {
                name = "density_jump";
            } else if (sim->mesh.amr_estimator[i] == PRJ_AMR_ESTIMATOR_VELOCITY) {
                name = "velocity";
            } else {
                if (sim->mesh.amr_lohner_var[i] == PRJ_LOHNER_VAR_LOG_DENSITY) {
                    name = "lohner_log_density";
                } else if (sim->mesh.amr_lohner_var[i] == PRJ_LOHNER_VAR_DENSITY) {
                    name = "lohner_density";
                } else if (sim->mesh.amr_lohner_var[i] == PRJ_LOHNER_VAR_TEMPERATURE) {
                    name = "lohner_temperature";
                } else {
                    name = "lohner_pressure";
                }
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
    fprintf(stderr, "max_level: %d\n", sim->mesh.max_level);
    fprintf(stderr, "min_dx: %.6e\n", sim->mesh.min_dx);
}

int main(int argc, char *argv[])
{
    prj_sim sim;
    prj_mpi mpi;
    prj_timer timer;
    prj_problem_init_fn init_fn;
    int init_with_mpi = 0;
    char *param_file = 0;
    double saved_amr_refine_thresh[PRJ_AMR_N];
    double saved_amr_derefine_thresh[PRJ_AMR_N];
    double saved_amr_lohner_eps[PRJ_AMR_N];
    int saved_amr_estimator[PRJ_AMR_N];
    int saved_amr_lohner_var[PRJ_AMR_N];
    int saved_amr_criterion_set[PRJ_AMR_N];
    int saved_use_amr_angle_resolution;
    double saved_amr_angle_resolution_limit;
    double saved_min_dx;
    int resolution = -1;
    int max_level_override = -1;
    double next_output_time = -1.0;
    double next_restart_time = -1.0;
    double last_output_time = -1.0;
    double last_restart_time = -1.0;
    int i;

    memset(&sim, 0, sizeof(sim));
    prj_timer_init(&timer);
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

    init_with_mpi = (init_fn == prj_problem_cc || init_fn == prj_problem_ccsn);
    prj_mpi_init(&argc, &argv, &mpi);
    prj_print_config(&sim, mpi.rank);
    if (sim.restart_from_file == 0) {
        init_fn(&sim);
        if (!init_with_mpi) {
            prj_mpi_decompose(&sim.mesh);
            prj_mpi_prepare(&sim.mesh, &mpi);
        }
    } else {
        prj_prepare_restart_problem(&sim, init_fn);
    }
    if (mpi.rank == 0) {
        mkdir("output", 0777);
    }
    sim.dump_count = 0;
    if (sim.restart_from_file != 0) {
        for (i = 0; i < PRJ_AMR_N; ++i) {
            saved_amr_refine_thresh[i] = sim.mesh.amr_refine_thresh[i];
            saved_amr_derefine_thresh[i] = sim.mesh.amr_derefine_thresh[i];
            saved_amr_estimator[i] = sim.mesh.amr_estimator[i];
            saved_amr_lohner_var[i] = sim.mesh.amr_lohner_var[i];
            saved_amr_lohner_eps[i] = sim.mesh.amr_lohner_eps[i];
            saved_amr_criterion_set[i] = sim.mesh.amr_criterion_set[i];
        }
        saved_use_amr_angle_resolution = sim.mesh.use_amr_angle_resolution;
        saved_amr_angle_resolution_limit = sim.mesh.amr_angle_resolution_limit;
        saved_min_dx = sim.mesh.min_dx;
        prj_io_read_restart(&sim.mesh, &sim.eos, sim.restart_file_name, &sim.time, &sim.step, &sim.dump_count,
            &last_output_time, &last_restart_time, &sim.dt);
        for (i = 0; i < PRJ_AMR_N; ++i) {
            sim.mesh.amr_refine_thresh[i] = saved_amr_refine_thresh[i];
            sim.mesh.amr_derefine_thresh[i] = saved_amr_derefine_thresh[i];
            sim.mesh.amr_estimator[i] = saved_amr_estimator[i];
            sim.mesh.amr_lohner_var[i] = saved_amr_lohner_var[i];
            sim.mesh.amr_lohner_eps[i] = saved_amr_lohner_eps[i];
            sim.mesh.amr_criterion_set[i] = saved_amr_criterion_set[i];
        }
        sim.mesh.use_amr_angle_resolution = saved_use_amr_angle_resolution;
        sim.mesh.amr_angle_resolution_limit = saved_amr_angle_resolution_limit;
        sim.mesh.min_dx = saved_min_dx;
        next_output_time = prj_next_event_time(last_output_time, sim.output_dt, sim.time);
        next_restart_time = prj_next_event_time(last_restart_time, sim.restart_dt, sim.time);
        prj_print_config(&sim, mpi.rank);
    }
    if (sim.restart_from_file == 0) {
        prj_mesh_mark_base_blocks(&sim.mesh);
    }
    prj_rad_init(&sim.rad);
 #if PRJ_USE_GRAVITY
    prj_gravity_init(&sim);
 #endif

    prj_eos_fill_active_cells(&sim.mesh, &sim.eos, 1);
    prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
    prj_eos_fill_mesh(&sim.mesh, &sim.eos, 1);
#if PRJ_MHD
    prj_boundary_fill_bf(&sim.mesh, &sim.bc, 0);
#endif
#if PRJ_USE_GRAVITY
    prj_gravity_monopole_reduce(&sim.mesh, 1);
    prj_gravity_monopole_integrate(&sim.mesh);
#endif

    if (sim.restart_from_file == 0) {
        prj_io_write_dump(&sim.mesh, sim.output_dir, sim.dump_count, sim.step, sim.time);
        sim.dump_count += 1;
        if (sim.output_dt >= 0.0) {
            next_output_time = sim.time + sim.output_dt;
        }
        if (sim.restart_dt >= 0.0) {
            next_restart_time = sim.time + sim.restart_dt;
        }
    } else {
        if (sim.output_dt >= 0.0 && next_output_time < 0.0) {
            next_output_time = sim.time + sim.output_dt;
        }
        if (sim.restart_dt >= 0.0 && next_restart_time < 0.0) {
            next_restart_time = sim.time + sim.restart_dt;
        }
    }

    struct timeval wall_start;
    gettimeofday(&wall_start, 0);

    while (sim.time < sim.t_end && sim.step < sim.max_steps) {
        int write_output = 0;
        int write_restart = 0;
        double next_event_time = -1.0;
        double dt_step;
        double dt_src = 1.0e100;

        PRJ_TIMER_START(&timer, "main_loop");
        {
            double dt_new = prj_timeint_calc_dt(&sim.mesh, &sim.eos, sim.cfl);

            if (sim.dt > 0.0) {
                double dt_limit = sim.dt_factor * sim.dt;

                sim.dt = dt_new <= dt_limit ? dt_new : dt_limit;
            } else {
                sim.dt = dt_new;
            }
        }
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
        dt_step = sim.dt;
        prj_timeint_step(&sim.mesh, &sim.coord, &sim.bc, &sim.eos, &sim.rad, dt_step, &dt_src,
#if PRJ_TIMER
            &timer
#else
            0
#endif
        );
        dt_src = prj_mpi_min_dt(dt_src);
        if (sim.cfl * dt_src < sim.dt) {
            sim.dt = sim.cfl * dt_src;
        }
        sim.time += dt_step;
        sim.step += 1;
        if (sim.amr_interval > 0 && sim.step % sim.amr_interval == 0) {
            PRJ_TIMER_START(&timer, "amr");
#if PRJ_USE_GRAVITY
            if (prj_amr_criteria_need_gravity(&sim.mesh)) {
                prj_gravity_monopole_reduce(&sim.mesh, 1);
                prj_gravity_monopole_integrate(&sim.mesh);
            }
#endif
            prj_eos_fill_ghost_cons(&sim.mesh, &sim.eos, 1);
            int block_changed = prj_amr_adapt(&sim.mesh, &sim.eos);
            prj_mpi_rebalance(&sim.mesh);
            if (block_changed) {
#if PRJ_USE_GRAVITY
                prj_gravity_rebuild_grid(&sim);
#endif
            }
            PRJ_TIMER_STOP(&timer, "amr");
            if (block_changed) {
                PRJ_TIMER_START(&timer, "ghost_fill_post_amr");
                prj_eos_fill_active_cells(&sim.mesh, &sim.eos, 1);
                prj_boundary_fill_ghosts(&sim.mesh, &sim.bc, 1);
                prj_eos_fill_mesh(&sim.mesh, &sim.eos, 1);
            #if PRJ_MHD
                prj_boundary_fill_bf(&sim.mesh, &sim.bc, 0);
            #endif
                PRJ_TIMER_STOP(&timer, "ghost_fill_post_amr");
            #if PRJ_USE_GRAVITY
                prj_gravity_monopole_reduce(&sim.mesh, 1);
                prj_gravity_monopole_integrate(&sim.mesh);
            #endif
            }
        }
#if PRJ_MHD && PRJ_MHD_DEBUG
        prj_mhd_debug_check_divb(&sim.mesh, 0);
#endif
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
            prj_io_write_dump(&sim.mesh, sim.output_dir, sim.dump_count, sim.step, sim.time);
            sim.dump_count += 1;
        }
        if (write_restart) {
            prj_io_write_restart(&sim.mesh, sim.time, sim.step, sim.dump_count,
                prj_last_event_time(next_output_time, sim.output_dt),
                prj_last_event_time(next_restart_time, sim.restart_dt), sim.dt);
        }
        if (mpi.rank == 0) {
            struct timeval wall_now;
            double wall_elapsed;
            double min_cell_size;

            long wall_days;
            long wall_hours;
            long wall_minutes;
            double wall_seconds;
            long total_sec;

            gettimeofday(&wall_now, 0);
            wall_elapsed = (double)(wall_now.tv_sec - wall_start.tv_sec) +
                1.0e-6 * (double)(wall_now.tv_usec - wall_start.tv_usec);
            min_cell_size = prj_mesh_min_cell_size(&sim.mesh);
            total_sec = (long)wall_elapsed;
            wall_days = total_sec / 86400;
            wall_hours = (total_sec % 86400) / 3600;
            wall_minutes = (total_sec % 3600) / 60;
            wall_seconds = wall_elapsed - (double)(total_sec - total_sec % 60);
            fprintf(stderr,
                "step=%d  t=%.6e  dt=%.6e  blocks=%d  max_active_level=%d  min_cell_size=%.6e  wall=%ldd %ldh %ldm %.3fs\n",
                sim.step, sim.time, dt_step, prj_mesh_count_active(&sim.mesh), sim.mesh.max_active_level,
                min_cell_size,
                wall_days, wall_hours, wall_minutes, wall_seconds);
        }
        PRJ_TIMER_STOP(&timer, "main_loop");
    }

#if PRJ_TIMER
    prj_write_timer_report(&timer, mpi.rank);
#endif

    prj_io_write_restart(&sim.mesh, sim.time, sim.step, sim.dump_count,
        prj_last_event_time(next_output_time, sim.output_dt),
        prj_last_event_time(next_restart_time, sim.restart_dt), sim.dt);
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
