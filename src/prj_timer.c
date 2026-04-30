#include <string.h>
#include <sys/time.h>

#include "prj.h"

double prj_timer_now(void)
{
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
}

void prj_timer_reset(prj_timer *timer)
{
    int i;

    if (timer == 0) {
        return;
    }
    for (i = 0; i < timer->nentry; ++i) {
        timer->entry[i].total = 0.0;
        timer->entry[i].start = 0.0;
        timer->entry[i].count = 0;
        timer->entry[i].active = 0;
    }
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
    timer->entry[idx].start = 0.0;
    timer->entry[idx].count = 0;
    timer->entry[idx].active = 0;
    timer->nentry += 1;
    return idx;
}

int prj_timer_start(prj_timer *timer, const char *name)
{
    int idx = prj_timer_find_or_add(timer, name);

    if (idx < 0 || timer->entry[idx].active != 0) {
        return -1;
    }
    timer->entry[idx].start = prj_timer_now();
    timer->entry[idx].active = 1;
    return 0;
}

int prj_timer_stop(prj_timer *timer, const char *name)
{
    int idx = prj_timer_find(timer, name);
    double now;

    if (idx < 0 || timer->entry[idx].active == 0) {
        return -1;
    }
    now = prj_timer_now();
    timer->entry[idx].total += now - timer->entry[idx].start;
    timer->entry[idx].start = 0.0;
    timer->entry[idx].count += 1;
    timer->entry[idx].active = 0;
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
    if (timer->entry[idx].active != 0) {
        total += prj_timer_now() - timer->entry[idx].start;
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

        fprintf(stream, "timer rank=%d name=%s count=%d active=%d total=%.6e avg=%.6e\n",
            rank, entry->name != 0 ? entry->name : "(null)", entry->count, entry->active,
            total, entry->count > 0 ? total / (double)entry->count : 0.0);
    }
}
