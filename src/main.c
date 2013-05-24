/* main.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

/* Global variables */

Particle *part_dpd;
Particle *part_mon;
System sys;


/* Global functions */

double  calc_energy(int i);
double  calc_energy_list(int i);
double  energy_ij(int i, int j);
void    initialize(void);
void    init_param(void);
void    init_stats(void);
void    input(void);
int     mod(int, int);
void    monte_carlo(void);
void    new_list(void);
void    output(void);
double  ran3(void);
void    random_move(int i);
void    sample(void);
void    setup_coords(void);
double  total_energy(void);
void    write_log(void);
void    write_mon(void);
Vector  vdist(Vector, Vector);
double  vmag(Vector);

main() {
    int i;
    double run_time;
    clock_t begin, end;

    begin = clock();
    initialize();

    srand(time(NULL));
    if (sys.calc_list == 1) new_list();
    for (i = 0; i <= sys.nsteps; i++) {
        monte_carlo();
        if (i % sys.freq_sample == 0) sample();
    }

    end = clock();
    final_stats();
    output();
    write_mon();

    run_time = (double) (end - begin) / CLOCKS_PER_SEC;
    printf("\nRun time of %lf seconds.\n\n", run_time);
}

