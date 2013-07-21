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
void    calc_energy_brute(void);
double  calc_energy_dpd(int i);
double  calc_energy_mon(int i);
void    calc_cm(void);
void    calc_nseg(void);
void    calc_pressure(void);
void    calc_q(void);
void    calc_re(void);
void    calc_rg(void);
void    check_bond(int i);
void    check_cell(Vector, Vector);
int     check_pore(Vector);
int     check_side(Vector);
void    check_wall(Vector);
double  energy_c(Vector);
double  energy_fene(Vector);
void    initialize(void);
void    init_energy(void);
void    init_monitor(void);
void    init_param(void);
void    init_polymer(void);
void    init_pore(void);
void    init_solvent(void);
void    init_stats(void);
void    init_wall(void);
void    input(void);
int     mod(int, int);
void    monte_carlo(void);
void    move_monomer(int i);
void    move_solvent(int i);
void    new_list(void);
void    periodic_bc(Vector *);
void    print_stats(void);
double  ran3(void);
void    sample(void);
void    update_monitor(void);
void    update_stats(void);
void    write_log(void);
void    write_mon(void);
Vector  vdist(Vector, Vector);
double  vmag(Vector);

main() {
  double run_time;
  clock_t begin, end;

  begin = clock();
  initialize();
  srand(time(NULL));

  for (sys.step = 0; sys.step <= sys.n_steps; sys.step++) {
    monte_carlo();
    if (sys.step % sys.freq_sample == 0) sample();
  }

  final_stats();
  print_stats();
  write_mon();

  end = clock();
  run_time = (double) (end - begin) / CLOCKS_PER_SEC;
  printf("\nRun time of %lf seconds.\n\n", run_time);
}
