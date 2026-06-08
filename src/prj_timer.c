#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "prj.h"

#if PRJ_TIMER
static prj_timer *prj_timer_active = 0;
#endif

double prj_timer_now(void)
{
#if defined(CLOCK_MONOTONIC)
    struct timespec ts;

    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
        return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
    }
#endif
    struct timeval tv;

    gettimeofday(&tv, 0);
    return (double)tv.tv_sec + 1.0e-6 * (double)tv.tv_usec;
}

void prj_timer_init(prj_timer *timer)
{
    if (timer == 0) {
        return;
    }
    timer->nentry = 0;
    timer->stack_depth = 0;
}

void prj_timer_reset(prj_timer *timer)
{
    int i;

    if (timer == 0) {
        return;
    }
    for (i = 0; i < timer->nentry; ++i) {
        timer->entry[i].total = 0.0;
        timer->entry[i].inclusive_total = 0.0;
        timer->entry[i].start = 0.0;
        timer->entry[i].inclusive_start = 0.0;
        timer->entry[i].count = 0;
        timer->entry[i].active = 0;
    }
    timer->stack_depth = 0;
}

void prj_timer_set_current(prj_timer *timer)
{
#if PRJ_TIMER
    prj_timer_active = timer;
#else
    (void)timer;
#endif
}

prj_timer *prj_timer_current(void)
{
#if PRJ_TIMER
    return prj_timer_active;
#else
    return 0;
#endif
}

static int prj_timer_find(const prj_timer *timer, const char *name)
{
    int i;

    if (timer == 0 || name == 0) {
        return -1;
    }
    for (i = 0; i < timer->nentry; ++i) {
        if (timer->entry[i].name != 0 && strcmp(timer->entry[i].name, name) == 0) {
            return i;
        }
    }
    return -1;
}

static int prj_timer_find_or_add(prj_timer *timer, const char *name)
{
    int idx;

    idx = prj_timer_find(timer, name);
    if (idx >= 0) {
        return idx;
    }
    if (timer == 0 || name == 0 || timer->nentry >= PRJ_TIMER_MAX_ENTRIES) {
        return -1;
    }
    idx = timer->nentry;
    timer->entry[idx].name = name;
    timer->entry[idx].total = 0.0;
    timer->entry[idx].inclusive_total = 0.0;
    timer->entry[idx].start = 0.0;
    timer->entry[idx].inclusive_start = 0.0;
    timer->entry[idx].count = 0;
    timer->entry[idx].active = 0;
    timer->nentry += 1;
    return idx;
}

static int prj_timer_start_idx(prj_timer *timer, int idx)
{
    double now;

    if (timer == 0 || idx < 0 || idx >= timer->nentry ||
        timer->entry[idx].active != 0 ||
        timer->stack_depth >= PRJ_TIMER_MAX_DEPTH) {
        return -1;
    }

    now = prj_timer_now();
    if (timer->stack_depth > 0) {
        int parent_idx = timer->stack[timer->stack_depth - 1];
        prj_timer_entry *parent = &timer->entry[parent_idx];

        if (parent->active != 0 && parent->start > 0.0) {
            parent->total += now - parent->start;
            parent->start = 0.0;
        }
    }
    timer->entry[idx].start = now;
    timer->entry[idx].inclusive_start = now;
    timer->entry[idx].active = 1;
    timer->stack[timer->stack_depth] = idx;
    timer->stack_depth += 1;
    return 0;
}

int prj_timer_start(prj_timer *timer, const char *name)
{
    int idx = prj_timer_find_or_add(timer, name);

    if (idx < 0) {
        return -1;
    }
    return prj_timer_start_idx(timer, idx);
}

static int prj_timer_cache_lookup(prj_timer *timer, const char *name, int *cache_idx, int add)
{
    int idx;

    if (timer == 0 || name == 0) {
        return -1;
    }
    if (cache_idx != 0 && *cache_idx >= 0 && *cache_idx < timer->nentry) {
        const char *cached_name = timer->entry[*cache_idx].name;

        if (cached_name == name || (cached_name != 0 && strcmp(cached_name, name) == 0)) {
            return *cache_idx;
        }
    }
    idx = add != 0 ? prj_timer_find_or_add(timer, name) : prj_timer_find(timer, name);
    if (cache_idx != 0) {
        *cache_idx = idx;
    }
    return idx;
}

int prj_timer_start_cached(prj_timer *timer, const char *name, int *cache_idx)
{
    int idx = prj_timer_cache_lookup(timer, name, cache_idx, 1);

    if (idx < 0) {
        return -1;
    }
    return prj_timer_start_idx(timer, idx);
}

static int prj_timer_stop_idx(prj_timer *timer, int idx)
{
    double now;

    if (timer == 0 || idx < 0 || idx >= timer->nentry ||
        timer->entry[idx].active == 0 ||
        timer->stack_depth <= 0 ||
        timer->stack[timer->stack_depth - 1] != idx) {
        return -1;
    }
    now = prj_timer_now();
    timer->entry[idx].total += now - timer->entry[idx].start;
    timer->entry[idx].inclusive_total += now - timer->entry[idx].inclusive_start;
    timer->entry[idx].start = 0.0;
    timer->entry[idx].inclusive_start = 0.0;
    timer->entry[idx].count += 1;
    timer->entry[idx].active = 0;
    timer->stack_depth -= 1;
    if (timer->stack_depth > 0) {
        int parent_idx = timer->stack[timer->stack_depth - 1];
        prj_timer_entry *parent = &timer->entry[parent_idx];

        if (parent->active != 0) {
            parent->start = now;
        }
    }
    return 0;
}

int prj_timer_stop(prj_timer *timer, const char *name)
{
    int idx = prj_timer_find(timer, name);

    return prj_timer_stop_idx(timer, idx);
}

int prj_timer_stop_cached(prj_timer *timer, const char *name, int *cache_idx)
{
    int idx = prj_timer_cache_lookup(timer, name, cache_idx, 0);

    return prj_timer_stop_idx(timer, idx);
}

int prj_timer_add(prj_timer *timer, const char *name, double elapsed)
{
    return prj_timer_add_cached(timer, name, 0, elapsed);
}

int prj_timer_add_cached(prj_timer *timer, const char *name, int *cache_idx, double elapsed)
{
    int idx;

    if (timer == 0 || name == 0 || elapsed < 0.0) {
        return -1;
    }
    idx = prj_timer_cache_lookup(timer, name, cache_idx, 1);
    if (idx < 0) {
        return -1;
    }
    timer->entry[idx].total += elapsed;
    timer->entry[idx].inclusive_total += elapsed;
    timer->entry[idx].count += 1;
    if (timer->stack_depth > 0) {
        int parent_idx = timer->stack[timer->stack_depth - 1];
        prj_timer_entry *parent = &timer->entry[parent_idx];

        if (parent_idx != idx && parent->active != 0 && parent->start > 0.0) {
            parent->start += elapsed;
        }
    }
    return 0;
}

double prj_timer_total(const prj_timer *timer, const char *name)
{
    int idx = prj_timer_find(timer, name);
    double total;

    if (idx < 0) {
        return 0.0;
    }
    total = timer->entry[idx].total;
    if (timer->entry[idx].active != 0 &&
        timer->stack_depth > 0 &&
        timer->stack[timer->stack_depth - 1] == idx &&
        timer->entry[idx].start > 0.0) {
        total += prj_timer_now() - timer->entry[idx].start;
    }
    return total;
}

double prj_timer_inclusive_total(const prj_timer *timer, const char *name)
{
    int idx = prj_timer_find(timer, name);
    double total;

    if (idx < 0) {
        return 0.0;
    }
    total = timer->entry[idx].inclusive_total;
    if (timer->entry[idx].active != 0 && timer->entry[idx].inclusive_start > 0.0) {
        total += prj_timer_now() - timer->entry[idx].inclusive_start;
    }
    return total;
}

int prj_timer_count(const prj_timer *timer, const char *name)
{
    int idx = prj_timer_find(timer, name);

    if (idx < 0) {
        return 0;
    }
    return timer->entry[idx].count;
}

void prj_timer_report(const prj_timer *timer, FILE *stream, int rank)
{
    int i;

    if (timer == 0) {
        return;
    }
    if (stream == 0) {
        stream = stderr;
    }
    fprintf(stream, "timer report rank=%d entries=%d\n", rank, timer->nentry);
    for (i = 0; i < timer->nentry; ++i) {
        const prj_timer_entry *entry = &timer->entry[i];
        double total = prj_timer_total(timer, entry->name);
        double inclusive = prj_timer_inclusive_total(timer, entry->name);

        fprintf(stream, "timer rank=%d name=%s count=%d active=%d total=%.6e avg=%.6e inclusive=%.6e inclusive_avg=%.6e\n",
            rank, entry->name != 0 ? entry->name : "(null)", entry->count, entry->active,
            total, entry->count > 0 ? total / (double)entry->count : 0.0,
            inclusive, entry->count > 0 ? inclusive / (double)entry->count : 0.0);
    }
}
