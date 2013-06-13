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

int     accept_move(void);
void    calc_bond_length(void);
void    calc_energy_brute(void);
double  calc_energy_dpd(int i);
double  calc_energy_mon(int i);
void    calc_cm(void);
void    calc_pressure(void);
void    calc_re(void);
void    calc_rg(void);
void    check_bond(int i);
int     check_cell(Vector, Vector);
void    check_wall(Vector);
double  energy_c(Vector);
double  energy_fene(int i, int j);
void    initialize(void);
void    init_param(void);
void    init_part(void);
void    init_stats(void);
void    init_wall(void);
void    input(void);
int     mod(int, int);
void    monitor_mem(void);
void    monte_carlo(void);
void    new_list(void);
void    periodic_bc(Vector *);
void    print_stats(void);
double  ran3(void);
void    random_move_dpd(int i);
void    random_move_mon(int i);
void    sample(void);
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

  for (i = 0; i <= sys.n_steps; i++) {
    monte_carlo();
    if (i % sys.freq_sample == 0) sample();
    if (i % sys.freq_monitor == 0) monitor();
  }

  final_stats();
  print_stats();
  write_mon();

  end = clock();
  run_time = (double) (end - begin) / CLOCKS_PER_SEC;
  printf("\nRun time of %lf seconds.\n\n", run_time);
}
