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
  init_wall();
  init_part();
  monitor_mem();
  write_log();
}

void init_param(void) {
  int fact;
  int i, j;

  sys.monitor_step = 0;
  sys.n_accept_dpd = 0;
  sys.n_accept_mon = 0;
  sys.n_attempt_dpd = 0;
  sys.n_attempt_mon = 0;

  // Wall particle spacing
  sys.r_wall = pow(sys.density_w, -1.0/3.0);

  sys.volume = sys.length.x * sys.length.y * sys.length.z;
  // Account for the volume occupied by the wall
  sys.volume -= sys.length.x * sys.length.y * sys.n_layers * sys.r_wall;

  sys.n_dpd = sys.density_s * sys.volume;

  // Monomer-monomer bonds
  sys.r_max = 2.0;
  sys.r_eq = 0.7;
  sys.r_0 = sys.r_max - sys.r_eq;
  sys.k_fene = 40;

  if (sys.calc_list == 1) {
    // Determine cell size and the number of cells
    fact = (int) (sys.length.x / sys.r_c + 1e-9);
    sys.r_cell.x = (double) sys.length.x / fact;
    sys.n_cell.x = (int) (sys.length.x / sys.r_cell.x + 1e-9);

    fact = (int) (sys.length.y / sys.r_c + 1e-9);
    sys.r_cell.y = (double) sys.length.y / fact;
    sys.n_cell.y = (int) (sys.length.y / sys.r_cell.y + 1e-9);

    fact = (int) (sys.length.z / sys.r_c + 1e-9);
    sys.r_cell.z = (double) sys.length.z / fact;
    sys.n_cell.z = (int) (sys.length.z / sys.r_cell.z + 1e-9);

    // Allocate memory for the head of chain array
    sys.hoc = (int ***) malloc((sys.n_cell.x+1)*sizeof(int **));
    sys.hoc_copy = (int ***) malloc((sys.n_cell.x+1)*sizeof(int **));

    for (i = 0; i < sys.n_cell.x; i++) {
      sys.hoc[i] = (int **) malloc((sys.n_cell.y+1)*sizeof(int *));
      sys.hoc_copy[i] = (int **) malloc((sys.n_cell.y+1)*sizeof(int *));

      for (j = 0; j < sys.n_cell.y; j++) {
        sys.hoc[i][j] = (int *) malloc((sys.n_cell.z+1)*sizeof(int));
        sys.hoc_copy[i][j] = (int *) malloc((sys.n_cell.z+1)*sizeof(int));
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

void init_wall(void) {
  int i, j, k, n, n_col, n_row;

  n_col = (int) sys.length.y / sys.r_wall + 1;
  n_row = (int) sys.length.x / sys.r_wall + 1;
  n = sys.n_dpd;
  sys.n_wall = sys.n_layers * n_col * n_row;

  part_dpd = (Particle *) calloc(sys.n_wall+sys.n_dpd, sizeof(Particle));

  for (k = 0; k < sys.n_layers; k++) {
    for (j = 0; j < n_col; j++) {
      for (i = 0; i < n_row; i++) {
        part_dpd[n].r.x = i * sys.r_wall;
        part_dpd[n].r.y = j * sys.r_wall;
        part_dpd[n].r.z = sys.length.z / 2 - k * sys.r_wall;

        n++;
      }
    }
  }

  sys.wall_max_z = sys.length.z / 2;
  sys.wall_min_z = sys.length.z / 2 - k * sys.r_wall;
}

void init_part(void) {
  int i;
  double x, y, z, d;

  for (i = 0; i < sys.n_dpd; i++) {
    // Repeat until the particle is outside the wall
    do {
      part_dpd[i].r.x = sys.length.x*ran3();
      part_dpd[i].r.y = sys.length.y*ran3();
      part_dpd[i].r.z = sys.length.z*ran3();
      part_dpd[i].ro.x = part_dpd[i].r.x;
      part_dpd[i].ro.y = part_dpd[i].r.y;
      part_dpd[i].ro.z = part_dpd[i].r.z;

      // Check for wall overlap
      check_wall(part_dpd[i].r);
    } while (!sys.wall_overlap);
  }

  part_mon = (Particle *) calloc(sys.n_mon, sizeof(Particle));

  for (i = 0; i < sys.n_mon; i++) {
    part_mon[i].r.x = sys.length.x/2;
    part_mon[i].r.y = sys.length.y/2;
    part_mon[i].r.z = sys.pol_init_z + i*sys.pol_init_bl;

    periodic_bc_r(&part_mon[i].r);

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
  int nsize;

  nsize = (sys.n_steps / sys.freq_monitor) + 1;

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

