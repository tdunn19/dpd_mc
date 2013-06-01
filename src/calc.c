#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void sample(void) {
  int i;

  P.now = 0;
  BL.now = 0;
  Etot.now = sys.energy;

  calc_pressure();
  calc_re();
  calc_bond_length();

  for (i = 0; i < sys.n_stats; i++) {
    sys.stats[i].sum += sys.stats[i].now;
    sys.stats[i].sumsq += sys.stats[i].now*sys.stats[i].now;
    sys.stats[i].num += 1;
  }
}

void monitor(void) {
  sys.mon.energy[sys.monitor_step] = Etot.now;
  sys.mon.re2[sys.monitor_step] = RE2.now;
  sys.mon.rex[sys.monitor_step] = sqrt(RE2x.now);
  sys.mon.rey[sys.monitor_step] = sqrt(RE2y.now);
  sys.mon.rez[sys.monitor_step] = sqrt(RE2z.now);
  sys.mon.bond_length[sys.monitor_step] = BL.now;

  sys.monitor_step += 1;
}

void calc_pressure(void) {
  // Calculate the pressure using the truncated virial expansion
  int i, j;
  double r_ij, fact;
  Vector dr;

  // Contribution from solvent-solvent forces
  for (i = 0; i < sys.n_dpd-1; i++) {
    for (j = i+1; j < sys.n_dpd; j++) {
      dr = vdist(part_dpd[i].r, part_dpd[j].r);
      periodic_bc(&dr);
      r_ij = vmag(dr);

      // Check for cutoff distance
      if (r_ij < sys.r_c) {
        P.now += sys.a_ss*r_ij*(1-r_ij);
      }
    }
  }

  for (i = 0; i < sys.n_mon; i++) {
    // Contribution from monomer-monomer forces
    for (j = i+1; j < sys.n_mon; j++) {
      dr = vdist(part_mon[i].r, part_mon[j].r);
      periodic_bc(&dr);
      r_ij = vmag(dr);

      // Conservative force
      if (r_ij < sys.r_c) {
        P.now += sys.a_mm*r_ij*(1-r_ij);
      }

      // FENE spring force
      if (j == i+1) {
        fact = 1.0 - (r_ij - sys.r_eq) / sys.r_0;
        P.now += -sys.k_fene*r_ij*(r_ij-sys.r_eq) / fact;

        // There shouldn't be a case where r_ij > r_max because of bond breaks
        if (r_ij > sys.r_max) {
          printf("\nError: r_ij > r_max at monitor step %d\n", sys.monitor_step);
        }
      }
    }

    // Contribution from monomer-solvent forces
    for (j = 0; j < sys.n_dpd; j++) {
      dr = vdist(part_mon[i].r, part_dpd[j].r);
      periodic_bc(&dr);

      r_ij = vmag(dr);

      if (r_ij < sys.r_c) {
        P.now += sys.a_ms*r_ij*(1-r_ij);
      }
    }
  }

  P.now /= 3*sys.volume;
  P.now += sys.density*sys.temp;
}

void calc_re(void) {
  double rex, rey, rez;
  Vector dr;

  dr = vdist(part_mon[sys.n_mon-1].r, part_mon[0].r);
  rex = dr.x;
  rey = dr.y;
  rez = dr.z;

  RE2x.now = rex*rex;
  RE2y.now = rey*rey;
  RE2z.now = rez*rez;
  RE2.now = RE2x.now + RE2y.now + RE2z.now;
}

void calc_bond_length(void) {
  int i;
  Vector dr;

  for (i = 0; i < sys.n_mon-1; i++) {
     dr = vdist(part_mon[i].r, part_mon[i+1].r);
     BL.now += vmag(dr);
  }

  // Divide by number of bonds to get an average
  BL.now /= (sys.n_mon - 1);
}

void check_bond(int i) {
  double r_ij;
  Vector dr;

  sys.bond_break = 0;

  // If not the first monomer in the chain
  if (i != 0) {
    dr = vdist(part_mon[i].r, part_mon[i-1].r);
    periodic_bc(&dr);
    r_ij = vmag(dr);

    if (r_ij > sys.r_max) {
      sys.bond_break = 1;
    }
  }

  // If not the last monomer in the chain
  if (i != (sys.n_mon - 1)) {
    dr = vdist(part_mon[i].r, part_mon[i+1].r);
    periodic_bc(&dr);
    r_ij = vmag(dr);

    if (r_ij > sys.r_max) {
      sys.bond_break = 1;
    }
  }
}
