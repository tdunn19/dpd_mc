/* init.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void initialize(void) {
  input();
  init_param();
  init_stats();
  monitor_mem();
  write_log();
  setup_coords();
}

void init_param(void) {
  int fact;
  int i, j;

  sys.n_dpd = sys.density * sys.volume;
  sys.length = pow(sys.volume, 1.0/3.0);
  sys.monitor_step = 0;
  sys.n_accept_dpd = 0;
  sys.n_accept_mon = 0;
  sys.n_attempt_dpd = 0;
  sys.n_attempt_mon = 0;

  // Monomer-monomer bonds
  // sys.r_0 = r_max - r_eq, where r_max is the maximum
  // bond length and r_eq is the equilibrium bond length
  sys.r_max = 2.0;
  sys.r_eq = 0.7;
  sys.r_0 = sys.r_max - sys.r_eq;
  sys.k_fene = 40;

  if (sys.calc_list == 1) {
    // Determine cell size and the number of cells
    fact = (int) (sys.length / sys.r_c + 1e-9);
    sys.r_cell = (double) sys.length / fact;
    sys.n_cell = (long) (sys.length / sys.r_cell + 1e-9);

    // Allocate memory for the head of chain array
    sys.hoc = (int ***) malloc((sys.n_cell+1)*sizeof(int **));
    sys.hoc_copy = (int ***) malloc((sys.n_cell+1)*sizeof(int **));

    for (i = 0; i < sys.n_cell; i++) {
      sys.hoc[i] = (int **) malloc((sys.n_cell+1)*sizeof(int *));
      sys.hoc_copy[i] = (int **) malloc((sys.n_cell+1)*sizeof(int *));

      for (j = 0; j < sys.n_cell; j++) {
        sys.hoc[i][j] = (int *) malloc((sys.n_cell+1)*sizeof(int));
        sys.hoc_copy[i][j] = (int *) malloc((sys.n_cell+1)*sizeof(int));
      }
    }
  }
}

void init_stats(void) {
  sys.n_stats = 11;
  sys.stats = (Stats *) calloc(sys.n_stats, sizeof(Stats));

  sys.stats[0].name = "Pressure               ";
  sys.stats[1].name = "Energy                 ";
  sys.stats[2].name = "Re2                    ";
  sys.stats[3].name = "Re2x                   ";
  sys.stats[4].name = "Re2y                   ";
  sys.stats[5].name = "Re2z                   ";
  sys.stats[6].name = "Rg2                    ";
  sys.stats[7].name = "Rg2x                   ";
  sys.stats[8].name = "Rg2y                   ";
  sys.stats[9].name = "Rg2z                   ";
  sys.stats[10].name = "Bond_length            ";
}

void setup_coords(void) {
  int i;

  part_dpd = (Particle *) calloc(sys.n_dpd, sizeof(Particle));
  part_mon = (Particle *) calloc(sys.n_mon, sizeof(Particle));

  for (i = 0; i < sys.n_dpd; i++) {
    part_dpd[i].r.x = sys.length*ran3();
    part_dpd[i].r.y = sys.length*ran3();
    part_dpd[i].r.z = sys.length*ran3();
    part_dpd[i].ro.x = part_dpd[i].r.x;
    part_dpd[i].ro.y = part_dpd[i].r.y;
    part_dpd[i].ro.z = part_dpd[i].r.z;
  }

  for (i = 0; i < sys.n_mon; i++) {
    part_mon[i].r.x = sys.length/2;
    part_mon[i].r.y = sys.length/2;
    part_mon[i].r.z = 0.5 * i;
    part_mon[i].ro.x = part_mon[i].r.x;
    part_mon[i].ro.y = part_mon[i].r.y;
    part_mon[i].ro.z = part_mon[i].r.z;
  }

  if (sys.calc_list) {
    new_list();

    for (i = 0; i < sys.n_dpd; i++) {
      part_dpd[i].E = calc_energy_dpd(i);
      part_dpd[i].Eo = part_dpd[i].E;
    }

    for (i = 0; i < sys.n_mon; i++) {
      part_mon[i].E = calc_energy_mon(i);
      part_mon[i].Eo = part_mon[i].E;
    }
  } else {
    calc_energy_brute();

    for (i = 0; i < sys.n_dpd; i++) {
      part_dpd[i].Eo = part_dpd[i].E;
    }

    for (i = 0; i < sys.n_mon; i++) {
      part_mon[i].Eo = part_mon[i].E;
    }
  }

  sys.energy = total_energy();
}

void monitor_mem(void) {
  int i;
  int nsize;

  nsize = (sys.nsteps / sys.freq_monitor) + 1;

  sys.mon.energy = (double *) calloc(nsize, sizeof(double));
  sys.mon.re2 = (double *) calloc(nsize, sizeof(double));
  sys.mon.rex = (double *) calloc(nsize, sizeof(double));
  sys.mon.rey = (double *) calloc(nsize, sizeof(double));
  sys.mon.rez = (double *) calloc(nsize, sizeof(double));
  sys.mon.rg2 = (double *) calloc(nsize, sizeof(double));
  sys.mon.rgx = (double *) calloc(nsize, sizeof(double));
  sys.mon.rgy = (double *) calloc(nsize, sizeof(double));
  sys.mon.rgz = (double *) calloc(nsize, sizeof(double));
  sys.mon.cmx = (double *) calloc(nsize, sizeof(double));
  sys.mon.cmy = (double *) calloc(nsize, sizeof(double));
  sys.mon.cmz = (double *) calloc(nsize, sizeof(double));
  sys.mon.bond_length = (double *) calloc(nsize, sizeof(double));
}

