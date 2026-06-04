#ifndef PRJ_TIMER_H
#define PRJ_TIMER_H

#include <stdio.h>

#include "prj_defs.h"

#ifndef PRJ_TIMER_MAX_ENTRIES
#define PRJ_TIMER_MAX_ENTRIES 256
#endif

#ifndef PRJ_TIMER_MAX_DEPTH
#define PRJ_TIMER_MAX_DEPTH 128
#endif

typedef struct prj_timer_entry {
    const char *name;
    double total;
    double inclusive_total;
    double start;
    double inclusive_start;
    int count;
    int active;
} prj_timer_entry;

typedef struct prj_timer {
    prj_timer_entry entry[PRJ_TIMER_MAX_ENTRIES];
    int stack[PRJ_TIMER_MAX_DEPTH];
    int stack_depth;
    int nentry;
} prj_timer;

double prj_timer_now(void);
void prj_timer_init(prj_timer *timer);
void prj_timer_reset(prj_timer *timer);
void prj_timer_set_current(prj_timer *timer);
prj_timer *prj_timer_current(void);
int prj_timer_start(prj_timer *timer, const char *name);
int prj_timer_stop(prj_timer *timer, const char *name);
int prj_timer_start_cached(prj_timer *timer, const char *name, int *cache_idx);
int prj_timer_stop_cached(prj_timer *timer, const char *name, int *cache_idx);
double prj_timer_total(const prj_timer *timer, const char *name);
double prj_timer_inclusive_total(const prj_timer *timer, const char *name);
int prj_timer_count(const prj_timer *timer, const char *name);
void prj_timer_report(const prj_timer *timer, FILE *stream, int rank);

struct prj_mpi;
void prj_mpi_barrier(const struct prj_mpi *mpi);

#if PRJ_TIMER
#define PRJ_TIMER_START(timer, name) do { \
        static int prj_timer_cache_idx = -1; \
        prj_timer_start_cached((timer), (name), &prj_timer_cache_idx); \
    } while (0)
#define PRJ_TIMER_STOP(timer, name) do { \
        static int prj_timer_cache_idx = -1; \
        prj_timer_stop_cached((timer), (name), &prj_timer_cache_idx); \
    } while (0)
#define PRJ_TIMER_CURRENT_START(name) PRJ_TIMER_START(prj_timer_current(), (name))
#define PRJ_TIMER_CURRENT_STOP(name) PRJ_TIMER_STOP(prj_timer_current(), (name))
#define PRJ_TIMER_BARRIER_START(timer, mpi, name) do { \
        PRJ_TIMER_START((timer), (name)); \
        prj_mpi_barrier((mpi)); \
    } while (0)
#define PRJ_TIMER_BARRIER_STOP(timer, mpi, name) do { \
        (void)(mpi); \
        PRJ_TIMER_STOP((timer), (name)); \
    } while (0)
#else
#define PRJ_TIMER_START(timer, name) ((void)(timer), (void)(name))
#define PRJ_TIMER_STOP(timer, name) ((void)(timer), (void)(name))
#define PRJ_TIMER_CURRENT_START(name) ((void)(name))
#define PRJ_TIMER_CURRENT_STOP(name) ((void)(name))
#define PRJ_TIMER_BARRIER_START(timer, mpi, name) \
        ((void)(timer), (void)(mpi), (void)(name))
#define PRJ_TIMER_BARRIER_STOP(timer, mpi, name) \
        ((void)(timer), (void)(mpi), (void)(name))
#endif

#endif
