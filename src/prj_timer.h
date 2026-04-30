#ifndef PRJ_TIMER_H
#define PRJ_TIMER_H

#include <stdio.h>

#ifndef PRJ_TIMER_MAX_ENTRIES
#define PRJ_TIMER_MAX_ENTRIES 256
#endif

typedef struct prj_timer_entry {
    const char *name;
    double total;
    double start;
    int count;
    int active;
} prj_timer_entry;

typedef struct prj_timer {
    prj_timer_entry entry[PRJ_TIMER_MAX_ENTRIES];
    int nentry;
} prj_timer;

double prj_timer_now(void);
void prj_timer_init(prj_timer *timer);
void prj_timer_reset(prj_timer *timer);
int prj_timer_start(prj_timer *timer, const char *name);
int prj_timer_stop(prj_timer *timer, const char *name);
double prj_timer_total(const prj_timer *timer, const char *name);
int prj_timer_count(const prj_timer *timer, const char *name);
void prj_timer_report(const prj_timer *timer, FILE *stream, int rank);

#endif
