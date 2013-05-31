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
double  calc_energy_dpd_debug(int i, int j);
double  calc_energy_mon(int i);
double  calc_energy_mon_debug(int i);
void    calc_pressure(void);
void    calc_re(void);
void    check_bond(int i);
int     check_cell(Vector, Vector);
double  energy_c(Vector);
double  energy_fene(int i, int j);
void    initialize(void);
void    init_param(void);
void    init_stats(void);
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
    //printf("beginning step %d\n", i);
    monte_carlo();
    if (i % sys.freq_sample == 0) sample();
    if (i % sys.freq_monitor == 0) monitor();
  }

  end = clock();
  final_stats();
  print_stats();
  write_mon();

  run_time = (double) (end - begin) / CLOCKS_PER_SEC;
  printf("\nRun time of %lf seconds.\n\n", run_time);
}
