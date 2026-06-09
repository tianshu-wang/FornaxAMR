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

#ifndef PRJ_BRANCH
#define PRJ_BRANCH "unknown"
#endif
#ifndef PRJ_COMMIT
#define PRJ_COMMIT "unknown"
#endif

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

static void prj_prepare_restart_problem(prj_sim *sim, prj_mpi *mpi)
{
    if (sim == 0) {
        return;
    }
    prj_eos_init(&sim->eos, mpi);
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

static void prj_amr_estimator_label(const prj_sim *sim, char *label, size_t label_size)
{
    int written = 0;
    int i;
    int any = 0;

    if (label == 0 || label_size == 0) {
        return;
    }
    label[0] = '\0';
    if (sim == 0) {
        snprintf(label, label_size, "unknown");
        return;
    }
    for (i = 0; i < PRJ_AMR_N; ++i) {
        const char *name;

        if (sim->mesh.amr_criterion_set[i] == 0) {
            continue;
        }
        if (sim->mesh.amr_estimator[i] == PRJ_AMR_ESTIMATOR_PRESSURE_SCALE_HEIGHT) {
            name = "pressure_scale_height";
        } else if (sim->mesh.amr_estimator[i] == PRJ_AMR_ESTIMATOR_FRACTIONAL_JUMP) {
            if (sim->mesh.amr_fractional_jump_var[i] == PRJ_FRACTIONAL_JUMP_VAR_PRESSURE) {
                name = "fractional_jump_pressure";
            } else {
                name = "fractional_jump_density";
            }
        } else if (sim->mesh.amr_estimator[i] == PRJ_AMR_ESTIMATOR_VELOCITY) {
            name = "velocity";
        } else {
            if (sim->mesh.amr_lohner_var[i] == PRJ_LOHNER_VAR_DENSITY) {
                name = "lohner_density";
            } else if (sim->mesh.amr_lohner_var[i] == PRJ_LOHNER_VAR_TEMPERATURE) {
                name = "lohner_temperature";
            } else {
                name = "lohner_pressure";
            }
        }
        written += snprintf(label + written, label_size - (size_t)written,
            "%s%s", any ? "," : "", name);
        any = 1;
        if (written >= (int)label_size - 1) {
            break;
        }
    }
    if (any == 0) {
        snprintf(label, label_size, "none");
    }
}

static void prj_print_config(const prj_sim *sim, int rank)
{
    char amr_estimator_label[128];

    if (rank != 0 || sim == 0) {
        return;
    }
    prj_amr_estimator_label(sim, amr_estimator_label, sizeof(amr_estimator_label));

    fprintf(stderr, "config:\n");
    fprintf(stderr, "git branch: %s\n", PRJ_BRANCH);
    fprintf(stderr, "git commit: %s\n", PRJ_COMMIT);
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
        amr_estimator_label
    );
    fprintf(stderr, "use_BJ: %s\n", sim->mesh.use_BJ != 0 ? "on" : "off");
    fprintf(stderr, "max_level: %d\n", sim->mesh.max_level);
    fprintf(stderr, "max_blocks: %d\n", sim->mesh.max_blocks);
    fprintf(stderr, "min_dx: %.6e\n", sim->mesh.min_dx);
    fprintf(stderr, "x_com_err_tol: %.6e\n", sim->x_com_err_tol);
    fprintf(stderr, "amr_init_scale_factor: %.6e\n", sim->mesh.amr_init_scale_factor);
    fprintf(stderr, "timer_interval: %d\n", sim->timer_interval);
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
    int saved_amr_fractional_jump_var[PRJ_AMR_N];
    int saved_amr_criterion_set[PRJ_AMR_N];
    int saved_use_amr_angular_resolution_limit;
    int saved_use_BJ;
    double saved_min_dx;
    double restart_x_com[3] = {0.0, 0.0, 0.0};
    int resolution = -1;
    int max_level_override = -1;
    int restart_latest_id = -1;
    double next_output_time = -1.0;
    double next_restart_time = -1.0;
    double last_output_time = -1.0;
    double last_restart_time = -1.0;
    int i;

    memset(&sim, 0, sizeof(sim));
    prj_timer_init(&timer);
#if PRJ_TIMER
    prj_timer_set_current(&timer);
#endif
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
    if (sim.restart_from_latest != 0) {
        if (prj_io_find_latest_restart("output", sim.restart_file_name,
                sizeof(sim.restart_file_name), &restart_latest_id) != 0) {
            fprintf(stderr, "restart_from_latest is enabled but no restart_*.h5 found in output/\n");
            return 1;
        }
        sim.restart_from_file = 1;
    }
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

    init_with_mpi = (init_fn == prj_problem_cc || init_fn == prj_problem_ccsn || init_fn == prj_problem_sedov);
    prj_mpi_init(&argc, &argv, &mpi);
    if (sim.restart_from_latest != 0 && mpi.rank == 0) {
        fprintf(stderr, "restart_from_latest selected %s (id %d)\n", sim.restart_file_name, restart_latest_id);
    }
    if (sim.restart_from_file == 0) {
        prj_print_config(&sim, mpi.rank);
        prj_eos_init(&sim.eos, &mpi);
        init_fn(&sim, &mpi);
        if (!init_with_mpi) {
            prj_mpi_decompose(&sim.mesh, &mpi);
            prj_mpi_prepare(&sim.mesh, &mpi);
        }
    } else {
        prj_prepare_restart_problem(&sim, &mpi);
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
            saved_amr_fractional_jump_var[i] = sim.mesh.amr_fractional_jump_var[i];
            saved_amr_lohner_eps[i] = sim.mesh.amr_lohner_eps[i];
            saved_amr_criterion_set[i] = sim.mesh.amr_criterion_set[i];
        }
        saved_use_amr_angular_resolution_limit = sim.mesh.use_amr_angular_resolution_limit;
        saved_use_BJ = sim.mesh.use_BJ;
        saved_min_dx = sim.mesh.min_dx;
        prj_io_read_restart(&sim.mesh, &sim.eos, &mpi, sim.restart_file_name, &sim.time, &sim.step, &sim.dump_count,
            &last_output_time, &last_restart_time, &sim.dt, restart_x_com);
        for (i = 0; i < PRJ_AMR_N; ++i) {
            sim.mesh.amr_refine_thresh[i] = saved_amr_refine_thresh[i];
            sim.mesh.amr_derefine_thresh[i] = saved_amr_derefine_thresh[i];
            sim.mesh.amr_estimator[i] = saved_amr_estimator[i];
            sim.mesh.amr_lohner_var[i] = saved_amr_lohner_var[i];
            sim.mesh.amr_fractional_jump_var[i] = saved_amr_fractional_jump_var[i];
            sim.mesh.amr_lohner_eps[i] = saved_amr_lohner_eps[i];
            sim.mesh.amr_criterion_set[i] = saved_amr_criterion_set[i];
        }
        sim.mesh.use_amr_angular_resolution_limit = saved_use_amr_angular_resolution_limit;
        sim.mesh.use_BJ = saved_use_BJ;
        sim.mesh.min_dx = saved_min_dx;
        prj_mesh_update_min_allowable_cell_size(&sim.mesh);
        next_output_time = prj_next_event_time(last_output_time, sim.output_dt, sim.time);
        next_restart_time = prj_next_event_time(last_restart_time, sim.restart_dt, sim.time);
        prj_print_config(&sim, mpi.rank);
    }
    prj_rad_init(&sim.rad);
 #if PRJ_USE_GRAVITY
    prj_gravity_init(&sim, &mpi);
    if (sim.restart_from_file != 0) {
        sim.grav.x_com[0] = restart_x_com[0];
        sim.grav.x_com[1] = restart_x_com[1];
        sim.grav.x_com[2] = restart_x_com[2];
        sim.grav.x_com_new[0] = restart_x_com[0];
        sim.grav.x_com_new[1] = restart_x_com[1];
        sim.grav.x_com_new[2] = restart_x_com[2];
    }
 #endif

    if (sim.restart_from_file == 0 &&
        (init_fn == prj_problem_cc || init_fn == prj_problem_ccsn) &&
        sim.perturbation_gaussian_norm != 0.0) {
        prj_set_perturbation(&sim.mesh, &sim.eos, &mpi,
            sim.perturbation_gaussian_norm, sim.perturbation_seed);
    }
    prj_eos_fill_active_cells(&sim.mesh, &sim.eos, &mpi, 1);
    prj_boundary_fill_ghosts_and_bf(&sim.mesh, &mpi, &sim.bc, 1, 0, &sim.eos, 0,
        &sim.rad, PRJ_BOUNDARY_TIMER_SCOPE_NONE);
    prj_eos_fill_mesh(&sim.mesh, &sim.eos, &mpi, 1);
    prj_flux_fill_transport_opacity_halo(&sim.mesh, &sim.rad, &mpi, 1);
#if PRJ_USE_GRAVITY
    prj_gravity_monopole_reduce(&sim.mesh, &sim.grav, &mpi, 1);
    prj_gravity_monopole_integrate(&sim.mesh, &sim.grav, &mpi);
#endif

    if (sim.restart_from_file == 0) {
        prj_io_write_dump(&sim.mesh, &sim.grav, &mpi, sim.dump_count, sim.step, sim.time);
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

        PRJ_TIMER_BARRIER_START(&timer, &mpi, "main_loop");

#if PRJ_USE_GRAVITY
        PRJ_TIMER_BARRIER_START(&timer, &mpi, "gravity_update_x_com");
        prj_gravity_update_center_of_mass(&sim.mesh, &sim.grav, &mpi, sim.x_com_err_tol);
        PRJ_TIMER_BARRIER_STOP(&timer, &mpi, "gravity_update_x_com");
#endif
        {
            double dt_new;

            PRJ_TIMER_BARRIER_START(&timer, &mpi, "calc_dt");
            dt_new = prj_timeint_calc_dt(&sim.mesh, &sim.eos, &mpi, sim.cfl);
            PRJ_TIMER_BARRIER_STOP(&timer, &mpi, "calc_dt");

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
        dt_step = sim.dt;
        if (next_event_time >= 0.0 && sim.time + dt_step > next_event_time) {
            dt_step = 1.000001 * (next_event_time - sim.time);
        }
        if (sim.time + dt_step > sim.t_end) {
            dt_step = sim.t_end - sim.time;
        }
        prj_timeint_step(&sim.mesh, &sim.coord, &sim.bc, &sim.eos, &sim.rad, &sim.grav, &mpi,
            dt_step, &dt_src,
#if PRJ_TIMER
            &timer
#else
            0
#endif
        );

        PRJ_TIMER_BARRIER_START(&timer, &mpi, "min_dt_src");
        dt_src = prj_mpi_min_dt(&mpi, dt_src);
        PRJ_TIMER_BARRIER_STOP(&timer, &mpi, "min_dt_src");

        if (sim.cfl * dt_src < sim.dt) {
            sim.dt = sim.cfl * dt_src;
        }
        sim.time += dt_step;
        sim.step += 1;
        if (sim.amr_interval > 0 && sim.step % sim.amr_interval == 0) {
            PRJ_TIMER_BARRIER_START(&timer, &mpi, "amr");
            /* The pressure_scale_height criterion reads per-cell b->grav, which
             * is already current from the just-completed step's stage-2 ghost
             * fill (prj_boundary_fill_ghosts_and_bf runs reduce+integrate). The
             * mesh/density/x_com are unchanged since, so no pre-tag rebuild is
             * needed here; the post-adapt rebuild below handles grid changes. */
            prj_eos_fill_ghost_cons(&sim.mesh, &sim.eos, &mpi, 1);
            double E_injected_before = sim.eos.E_injected;
            int block_changed = prj_amr_adapt(&sim.mesh, &sim.eos, &mpi, &sim.grav);
            if (mpi.rank == 0 && sim.eos.E_injected != E_injected_before) {
                fprintf(stderr, "E_injected changed after amr_adapt: %.6e -> %.6e (delta=%.6e)\n",
                    E_injected_before, sim.eos.E_injected, sim.eos.E_injected - E_injected_before);
            }
            if (block_changed) {
                prj_mpi_rebalance(&sim.mesh, &mpi);
#if PRJ_USE_GRAVITY
                prj_gravity_rebuild_grid(&sim, &mpi);
#endif
            }
            if (block_changed) {
                prj_eos_fill_active_cells(&sim.mesh, &sim.eos, &mpi, 1);
                prj_boundary_fill_ghosts_and_bf(&sim.mesh, &mpi, &sim.bc, 1, 0,
                    &sim.eos, 0, &sim.rad, PRJ_BOUNDARY_TIMER_SCOPE_NONE);
                prj_eos_fill_mesh(&sim.mesh, &sim.eos, &mpi, 1);
                prj_flux_fill_transport_opacity_halo(&sim.mesh, &sim.rad, &mpi, 1);
            #if PRJ_USE_GRAVITY
                prj_gravity_monopole_reduce(&sim.mesh, &sim.grav, &mpi, 1);
                prj_gravity_monopole_integrate(&sim.mesh, &sim.grav, &mpi);
            #endif
            }
            PRJ_TIMER_BARRIER_STOP(&timer, &mpi, "amr");
        }
#if PRJ_MHD && PRJ_MHD_DEBUG
        PRJ_TIMER_BARRIER_START(&timer, &mpi, "mhd_debug");
        prj_mhd_debug_check_divb(&sim.mesh, &mpi, 0);
        PRJ_TIMER_BARRIER_STOP(&timer, &mpi, "mhd_debug");
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
            PRJ_TIMER_BARRIER_START(&timer, &mpi, "write_output");
            prj_io_write_dump(&sim.mesh, &sim.grav, &mpi, sim.dump_count, sim.step, sim.time);
            sim.dump_count += 1;
            PRJ_TIMER_BARRIER_STOP(&timer, &mpi, "write_output");
        }
        if (write_restart) {
            PRJ_TIMER_BARRIER_START(&timer, &mpi, "write_restart");
            prj_io_write_restart(&sim.mesh, &mpi, sim.time, sim.step, sim.dump_count,
                prj_last_event_time(next_output_time, sim.output_dt),
                prj_last_event_time(next_restart_time, sim.restart_dt), sim.dt,
                sim.grav.x_com);
            PRJ_TIMER_BARRIER_STOP(&timer, &mpi, "write_restart");
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
#if PRJ_TIMER
        if (sim.timer_interval > 0 && sim.step % sim.timer_interval == 0) {
            PRJ_TIMER_BARRIER_START(&timer, &mpi, "timer_write");
            prj_write_timer_report(&timer, mpi.rank);
            PRJ_TIMER_BARRIER_STOP(&timer, &mpi, "timer_write");
        }
#endif

        PRJ_TIMER_BARRIER_STOP(&timer, &mpi, "main_loop");
    }

    prj_io_write_restart(&sim.mesh, &mpi, sim.time, sim.step, sim.dump_count,
        prj_last_event_time(next_output_time, sim.output_dt),
        prj_last_event_time(next_restart_time, sim.restart_dt), sim.dt,
        sim.grav.x_com);
    if (mpi.rank == 0) {
        char final_restart[64];

        snprintf(final_restart, sizeof(final_restart), "output/restart_%08d.h5", sim.step);
        prj_copy_file(final_restart, "output/final.h5");
    }
#if PRJ_TIMER
    prj_write_timer_report(&timer, mpi.rank);
#endif
#if PRJ_USE_GRAVITY
    prj_gravity_free(&sim.grav);
#endif
    prj_mesh_destroy(&sim.mesh);
    prj_mpi_finalize(&mpi);
    return 0;
}
