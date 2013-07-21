/* init.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void initialize(void) {
  input();
  init_param();
  write_log();
  init_polymer();
  init_solvent();
  init_pore();
  init_wall();
  init_energy();
  init_monitor();
  init_stats();
}

void init_param(void) {
  int i, j, fact;

  sys.monitor_step = 0;
  sys.n_accept_mon = 0;
  sys.n_accept_solvent = 0;
  sys.n_attempt_mon = 0;
  sys.n_attempt_solvent = 0;

  // System flags
  sys.bond_break = 0;
  sys.wall_overlap = 0;
  sys.new_window = 0;

  // Wall particle spacing
  sys.r_wall = pow(sys.density_w, -1.0/3.0);

  // Pore particle spacing
  sys.n_pore_1d.x = (int) (2*sys.pore_radius / sys.r_wall + 0.5) + 1;
  sys.r_pore = 2*sys.pore_radius / (sys.n_pore_1d.x-1);
  sys.n_pore_1d.y = sys.n_pore_1d.x;

  // Pore bounds
  sys.pore_max.x = sys.length.x/2 + sys.pore_radius;
  sys.pore_min.x = sys.length.x/2 - sys.pore_radius;
  sys.pore_max.y = sys.length.y/2 + sys.pore_radius;
  sys.pore_min.y = sys.length.y/2 - sys.pore_radius;
  sys.pore_max.z = sys.length.z/2;
  sys.n_pore_1d.z = (int) (sys.pore_length / sys.r_wall) + 1;
  sys.pore_min.z = sys.pore_max.z - (sys.n_pore_1d.z-1)*sys.r_wall;
  sys.pore_length = sys.pore_max.z - sys.pore_min.z;

  // Wall bounds
  sys.n_wall_1d.z = sys.n_layers + 1;
  sys.wall_max.z = sys.length.z/2;
  sys.wall_min.z = sys.length.z/2 - sys.n_layers*sys.r_wall;
  sys.n_wall_1d.x = (int) (sys.pore_min.x / sys.r_wall);
  sys.wall_min.x = sys.pore_min.x - (sys.n_wall_1d.x)*sys.r_wall;
  sys.wall_max.x = sys.pore_max.x + (sys.n_wall_1d.x)*sys.r_wall;
  sys.n_wall_1d.x = 2*sys.n_wall_1d.x +  sys.n_pore_1d.x;
  sys.n_wall_1d.y = (int) (sys.pore_min.y / sys.r_wall);
  sys.wall_min.y = sys.pore_min.y - (sys.n_wall_1d.y)*sys.r_wall;
  sys.wall_max.y = sys.pore_max.y + (sys.n_wall_1d.y)*sys.r_wall;
  sys.n_wall_1d.y = 2*sys.n_wall_1d.y +  sys.n_pore_1d.y;

  // Adjust the system length to keep a consistent wall density
  sys.length.x += sys.r_wall - (sys.wall_min.x+sys.length.x) + sys.wall_max.x;
  sys.length.y += sys.r_wall - (sys.wall_min.y+sys.length.y) + sys.wall_max.y;

  // Number of pore particles
  sys.n_pore = 2 * (sys.n_pore_1d.x + sys.n_pore_1d.y - 2);
  sys.n_pore *= sys.n_pore_1d.z;

  // Number of wall particles
  sys.n_wall = sys.n_wall_1d.x * sys.n_wall_1d.y;
  sys.n_wall -= sys.n_pore_1d.x * sys.n_pore_1d.y;
  sys.n_wall *= sys.n_wall_1d.z;
  // Wall particles surrounding the pore (as part of layers)
  sys.n_wall += (sys.n_pore_1d.z-sys.n_layers-1) * (2*sys.n_layers)
    * (sys.n_pore_1d.x + sys.n_pore_1d.y + 2*sys.n_layers);

  // Wall volume
  sys.wall_volume = sys.length.x * sys.length.y;
  sys.wall_volume -= (sys.pore_max.x - sys.pore_min.x)*(sys.pore_max.y - sys.pore_min.y);
  sys.wall_volume *= sys.n_layers * sys.r_wall;

  // Pore volume
  sys.pore_volume = sys.pore_max.x - sys.pore_min.x + 2*sys.n_layers*sys.r_wall;
  sys.pore_volume *= sys.pore_max.y - sys.pore_min.y + 2*sys.n_layers*sys.r_wall;
  sys.pore_volume *= sys.pore_max.z - sys.pore_min.z;

  // Solvent volume and number of solvent particles
  sys.volume = sys.length.x * sys.length.y * sys.length.z;
  sys.volume -= sys.wall_volume + sys.pore_volume;
  sys.volume_cis = sys.length.x * sys.length.y * sys.length.z / 2;
  sys.volume_trans = sys.volume_cis - sys.wall_volume - sys.pore_volume;
  // sys.n_solvent = (int) (sys.density_s * sys.volume);
  sys.n_solvent_cis = (int) (sys.density_cis * sys.volume_cis);
  sys.n_solvent_trans = (int) (sys.density_trans * sys.volume_trans);
  sys.n_solvent = sys.n_solvent_cis + sys.n_solvent_trans;
  // Debugging
  printf("sys.volume = %lf\n", sys.volume);
  printf("sys.volume_cis = %lf\n", sys.volume_cis);
  printf("sys.volume_trans = %lf\n", sys.volume_trans);

  // Total dpd particles
  sys.n_dpd = sys.n_solvent + sys.n_pore + sys.n_wall;
  part_dpd = (Particle *) calloc(sys.n_dpd, sizeof(Particle));

  // Monomer-monomer bonds
  sys.r_max = 2.0;
  sys.r_max2 = sys.r_max * sys.r_max;
  sys.r_eq = 0.7;
  sys.r_0 = sys.r_max - sys.r_eq;
  sys.k_fene = 40;
  sys.r_c2 = sys.r_c * sys.r_c;

  // Windows and bins in the nanopore
  sys.window_width = 2.0 / (sys.n_wins+1);
  sys.bin_width = sys.window_width / sys.n_bins;
  sys.Q_init = sys.iQ_init * sys.window_width / 2;
  sys.Q_min = (sys.iQ_init-1) * sys.window_width / 2;
  sys.Q_max = (sys.iQ_init+1) * sys.window_width / 2;
  if (sys.Q_min < 0) {
    sys.Q_min = 0;
  }
  if (sys.Q_max > 1) {
    sys.Q_max = 1.0;
  }

  sys.bin_count = (int *) malloc((sys.n_bins)*sizeof(int));
  for (i = 0; i < sys.n_bins; i++) {
    sys.bin_count[i] = 0;
  }

  if (sys.calc_list) {
    // Determine cell size and the number of cells
    fact = (int) (sys.length.x / sys.r_c + 1e-9);
    sys.r_cell.x = (double) (sys.length.x / fact);
    sys.n_cell_1d.x = (int) (sys.length.x / sys.r_cell.x + 1e-9);

    fact = (int) (sys.length.y / sys.r_c + 1e-9);
    sys.r_cell.y = (double) (sys.length.y / fact);
    sys.n_cell_1d.y = (int) (sys.length.y / sys.r_cell.y + 1e-9);

    fact = (int) (sys.length.z / sys.r_c + 1e-9);
    sys.r_cell.z = (double) (sys.length.z / fact);
    sys.n_cell_1d.z = (int) (sys.length.z / sys.r_cell.z + 1e-9);

    // Allocate memory for the head of chain array
    sys.hoc = (int ***) malloc((sys.n_cell_1d.x+1)*sizeof(int **));
    sys.hoc_copy = (int ***) malloc((sys.n_cell_1d.x+1)*sizeof(int **));

    for (i = 0; i < sys.n_cell_1d.x; i++) {
      sys.hoc[i] = (int **) malloc((sys.n_cell_1d.y+1)*sizeof(int *));
      sys.hoc_copy[i] = (int **) malloc((sys.n_cell_1d.y+1)*sizeof(int *));

      for (j = 0; j < sys.n_cell_1d.y; j++) {
        sys.hoc[i][j] = (int *) malloc((sys.n_cell_1d.z+1)*sizeof(int));
        sys.hoc_copy[i][j] = (int *) malloc((sys.n_cell_1d.z+1)*sizeof(int));
      }
    }
  }
}

void init_polymer(void) {
  int i, count;
  double z_min, z_max, z_init;

  part_mon = (Particle *) calloc(sys.n_mon, sizeof(Particle));

  for (i = 0; i < sys.n_mon; i++) {
    part_mon[i].r.x = sys.pore_max.x - sys.pore_radius;
    part_mon[i].r.y = sys.pore_max.y - sys.pore_radius;
  }

  z_max = sys.pore_max.z + sys.bl_init*(sys.n_mon-1);
  z_min = sys.pore_min.z;

  count = 0;
  do {
    z_init = (z_max + z_min) / 2;

    for (i = 0; i < sys.n_mon; i++) {
      part_mon[i].r.z = z_init - i * sys.bl_init;
    }

    calc_q();

    if (sys.Q <= sys.Q_min) {
      z_min = z_init;
    } else if (sys.Q >= sys.Q_max) {
      z_max = z_init;
    }

    count++;
    if (count > 1000) {
      printf("Problem positioning the polymer\n");
      exit(0);
    }
  } while (sys.Q <= sys.Q_min || sys.Q >= sys.Q_max);

  for (i = 0; i < sys.n_mon; i++) {

    part_mon[i].ro.x = part_mon[i].r.x;
    part_mon[i].ro.y = part_mon[i].r.y;
    part_mon[i].ro.z = part_mon[i].r.z;

    periodic_bc_r(&part_mon[i].r);
    if (part_mon[i].r.z != part_mon[i].ro.z) {
      printf("System size (L = %lf, %lf, %lf) too small.\n",
        sys.length.x, sys.length.y, sys.length.z);
      exit(0);
    }
  }
}

void init_solvent(void) {
  int i, j;
  FILE *fp;
  printf("Beginning cis solvent\n\n");

  for (i = 0; i < sys.n_solvent_cis; i++) {
    // Repeat until the particle is outside the wall/pore on the cis side
    do {
      part_dpd[i].r.x = sys.length.x*ran3();
      part_dpd[i].r.y = sys.length.y*ran3();
      part_dpd[i].r.z = sys.length.z*ran3();

      // Check for wall overlap
      check_wall(part_dpd[i].r);
    } while (sys.wall_overlap || check_side(part_dpd[i].r) || check_pore(part_dpd[i].r));
    printf("dpd[%d].r=%lf,%lf,%lf\n",
      i,part_dpd[i].r.x,part_dpd[i].r.y,part_dpd[i].r.z);

    part_dpd[i].ro.x = part_dpd[i].r.x;
    part_dpd[i].ro.y = part_dpd[i].r.y;
    part_dpd[i].ro.z = part_dpd[i].r.z;
  }

  printf("Beginning trans solvent\n\n");
  for (j = i; j < sys.n_solvent; j++) {
    // Repeat until the particle is outside the wall/pore on the trans side
    do {
      part_dpd[j].r.x = sys.length.x*ran3();
      part_dpd[j].r.y = sys.length.y*ran3();
      part_dpd[j].r.z = sys.length.z*ran3();

      // Check for wall overlap
      check_wall(part_dpd[j].r);
    } while (sys.wall_overlap || !check_side(part_dpd[j].r) || check_pore(part_dpd[j].r));
    printf("dpd[%d].r=%lf,%lf,%lf\n",
      j,part_dpd[j].r.x,part_dpd[j].r.y,part_dpd[j].r.z);

    part_dpd[j].ro.x = part_dpd[j].r.x;
    part_dpd[j].ro.y = part_dpd[j].r.y;
    part_dpd[j].ro.z = part_dpd[j].r.z;
  }
}

void init_pore(void) {
  int i, j, k, n;
  Vector r;

  n = sys.n_solvent + sys.n_wall;
  r.x = sys.pore_min.x;
  r.y = sys.pore_min.y;
  r.z = sys.pore_min.z;

  for (k = 0; k < sys.n_pore_1d.z; k++) {
    for (j = 0; j < sys.n_pore_1d.y; j++) {
      for (i = 0; i < sys.n_pore_1d.x; i++) {
        if (i == 0 || i == sys.n_pore_1d.x-1) {
          part_dpd[n].r = r;
          n++;
        } else if (j == 0 || j == sys.n_pore_1d.y-1) {
          part_dpd[n].r = r;
          n++;
        }
        r.x += sys.r_pore;
      }
      r.x = sys.pore_min.x;
      r.y += sys.r_pore;
    }
    r.y = sys.pore_min.y;
    r.z += sys.r_wall;
  }
}

void init_wall(void) {
  int i, j, k, n;
  Vector r;

  r.x = sys.wall_min.x;
  r.y = sys.wall_min.y;
  r.z = sys.wall_min.z;
  n = sys.n_solvent;

  // Generate the wall in the xy plane
  for (k = 0; k < sys.n_wall_1d.z; k++) {
    for (j = 0; j < sys.n_wall_1d.y; j++) {
      for (i = 0; i < sys.n_wall_1d.x; i++) {

        if (r.x < sys.pore_min.x || r.x > sys.pore_max.x+1e-9) {
          part_dpd[n].r = r;
          n++;
        } else if (r.y < sys.pore_min.y || r.y > sys.pore_max.y+1e-9) {
          part_dpd[n].r = r;
          n++;
        }

        if (r.x >= sys.pore_min.x && r.x < sys.pore_max.x) {
          r.x += sys.r_pore;
        } else {
          r.x += sys.r_wall;
        }
      }

      r.x = sys.wall_min.x;

      if (r.y >= sys.pore_min.y && r.y < sys.pore_max.y) {
        r.y += sys.r_pore;
      } else {
        r.y += sys.r_wall;
      }
    }

    r.y = sys.wall_min.y;

    r.z += sys.r_wall;
  }

  // Add layers to the nanopore
  if (sys.n_layers > 0) {
    r.x = sys.pore_min.x - sys.n_layers*sys.r_wall;
    r.y = sys.pore_min.y - sys.n_layers*sys.r_wall;
    r.z = sys.pore_min.z;

    for (k = sys.n_wall_1d.z; k < sys.n_pore_1d.z; k++) {
      for (j = -sys.n_layers; j < sys.n_pore_1d.y+sys.n_layers; j++) {
        for (i = -sys.n_layers; i < sys.n_pore_1d.x+sys.n_layers; i++) {
          if (i < 0 || i > sys.n_pore_1d.x-1) {
            part_dpd[n].r = r;
            n++;
          } else if (j < 0 || j > sys.n_pore_1d.y-1) {
            part_dpd[n].r = r;
            n++;
          }

          if (i < 0 || i >= sys.n_pore_1d.x-1) {
            r.x += sys.r_wall;
          } else {
            r.x += sys.r_pore;
          }
        }

        r.x = sys.pore_min.x - sys.n_layers*sys.r_wall;

        if (j < 0 || j >= sys.n_pore_1d.y-1) {
          r.y += sys.r_wall;
        } else {
          r.y += sys.r_pore;
        }
      }

      r.y = sys.pore_min.y - sys.n_layers*sys.r_wall;
      r.z += sys.r_wall;
    }
  }
}

void init_energy(void) {
  int i;

  sys.energy = 0;

  if (sys.calc_list) {
    new_list();

    for (i = 0; i < sys.n_dpd; i++) {
      part_dpd[i].E = calc_energy_dpd(i);
      part_dpd[i].Eo = part_dpd[i].E;
      sys.energy += part_dpd[i].E;
    }

    for (i = 0; i < sys.n_mon; i++) {
      part_mon[i].E = calc_energy_mon(i);
      part_mon[i].Eo = part_mon[i].E;
      sys.energy += part_dpd[i].E;

    }
  } else {
    calc_energy_brute();

    for (i = 0; i < sys.n_dpd; i++) {
      part_dpd[i].Eo = part_dpd[i].E;
      sys.energy += part_dpd[i].E;
    }

    for (i = 0; i < sys.n_mon; i++) {
      part_mon[i].Eo = part_mon[i].E;
      sys.energy += part_dpd[i].E;
    }
  }
  sys.energy_old = sys.energy;
}

void init_monitor(void) {
  int n_size;

  n_size = (sys.n_steps / sys.freq_sample) + 1;

  sys.mon.energy = (double *) calloc(n_size, sizeof(double));
  sys.mon.Q = (double *) calloc(n_size, sizeof(double));
  sys.mon.re2 = (double *) calloc(n_size, sizeof(double));
  sys.mon.rex = (double *) calloc(n_size, sizeof(double));
  sys.mon.rey = (double *) calloc(n_size, sizeof(double));
  sys.mon.rez = (double *) calloc(n_size, sizeof(double));
  sys.mon.re2_cis = (double *) calloc(n_size, sizeof(double));
  sys.mon.rex_cis = (double *) calloc(n_size, sizeof(double));
  sys.mon.rey_cis = (double *) calloc(n_size, sizeof(double));
  sys.mon.rez_cis = (double *) calloc(n_size, sizeof(double));
  sys.mon.re2_trans = (double *) calloc(n_size, sizeof(double));
  sys.mon.rex_trans = (double *) calloc(n_size, sizeof(double));
  sys.mon.rey_trans = (double *) calloc(n_size, sizeof(double));
  sys.mon.rez_trans = (double *) calloc(n_size, sizeof(double));
  sys.mon.rg2 = (double *) calloc(n_size, sizeof(double));
  sys.mon.rg2x = (double *) calloc(n_size, sizeof(double));
  sys.mon.rg2y = (double *) calloc(n_size, sizeof(double));
  sys.mon.rg2z = (double *) calloc(n_size, sizeof(double));
  sys.mon.rg2_cis = (double *) calloc(n_size, sizeof(double));
  sys.mon.rg2x_cis = (double *) calloc(n_size, sizeof(double));
  sys.mon.rg2y_cis = (double *) calloc(n_size, sizeof(double));
  sys.mon.rg2z_cis = (double *) calloc(n_size, sizeof(double));
  sys.mon.rg2_trans = (double *) calloc(n_size, sizeof(double));
  sys.mon.rg2x_trans = (double *) calloc(n_size, sizeof(double));
  sys.mon.rg2y_trans = (double *) calloc(n_size, sizeof(double));
  sys.mon.rg2z_trans = (double *) calloc(n_size, sizeof(double));
  sys.mon.cmx = (double *) calloc(n_size, sizeof(double));
  sys.mon.cmy = (double *) calloc(n_size, sizeof(double));
  sys.mon.cmz = (double *) calloc(n_size, sizeof(double));
  sys.mon.Q = (double *) calloc(n_size, sizeof(double));
}

void init_stats(void) {
  sys.n_stats = 10;
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
}

