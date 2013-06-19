#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void sample(void) {
  int i;
  Vector dr;

  P.now = 0;
  BL.now = 0;
  Etot.now = sys.energy;

  // Easier to calculate pressure with monomer positions affected
  // by periodic boundary conditions
  calc_pressure();

  // Get an absolute position of the monomers without periodic boundaries
  for (i = 1; i < sys.n_mon; i++) {
    dr = vdist(part_mon[i].r, part_mon[i-1].r);
    periodic_bc_dr(&dr);

    part_mon[i].r.x = part_mon[i-1].r.x + dr.x;
    part_mon[i].r.y = part_mon[i-1].r.y + dr.y;
    part_mon[i].r.z = part_mon[i-1].r.z + dr.z;
  }

  calc_cm();
  calc_re();
  calc_rg();
  calc_bond_length();

  for (i = 0; i < sys.n_stats; i++) {
    sys.stats[i].sum += sys.stats[i].now;
    sys.stats[i].sumsq += sys.stats[i].now*sys.stats[i].now;
    sys.stats[i].num += 1;
  }

  // Revert the monomer positions
  for (i = 1; i < sys.n_mon; i++) {
    part_mon[i].r.x = part_mon[i].ro.x;
    part_mon[i].r.y = part_mon[i].ro.y;
    part_mon[i].r.z = part_mon[i].ro.z;
  }
}

void monitor(void) {
  sys.mon.energy[sys.monitor_step] = Etot.now;
  sys.mon.re2[sys.monitor_step] = RE2.now;
  sys.mon.rex[sys.monitor_step] = sqrt(RE2x.now);
  sys.mon.rey[sys.monitor_step] = sqrt(RE2y.now);
  sys.mon.rez[sys.monitor_step] = sqrt(RE2z.now);
  sys.mon.rg2[sys.monitor_step] = RG2.now;
  sys.mon.rgx[sys.monitor_step] = sqrt(RG2x.now);
  sys.mon.rgy[sys.monitor_step] = sqrt(RG2y.now);
  sys.mon.rgz[sys.monitor_step] = sqrt(RG2z.now);
  sys.mon.bond_length[sys.monitor_step] = BL.now;

  sys.monitor_step += 1;
}

void calc_cm(void) {
  int i;

  sys.mon.cmx[sys.monitor_step] = 0;
  sys.mon.cmy[sys.monitor_step] = 0;
  sys.mon.cmz[sys.monitor_step] = 0;

  for (i = 0; i < sys.n_mon; i++) {
    sys.mon.cmx[sys.monitor_step] += part_mon[i].r.x;
    sys.mon.cmy[sys.monitor_step] += part_mon[i].r.y;
    sys.mon.cmz[sys.monitor_step] += part_mon[i].r.z;
  }

  sys.mon.cmx[sys.monitor_step] /= sys.n_mon;
  sys.mon.cmy[sys.monitor_step] /= sys.n_mon;
  sys.mon.cmz[sys.monitor_step] /= sys.n_mon;
}

void calc_pressure(void) {
  // Calculate the pressure using the truncated virial expansion
  int i, j;
  double r_ij, fact;
  Vector dr;

  // Contribution from solvent-solvent forces
  for (i = 0; i < sys.n_solvent-1; i++) {
    for (j = i+1; j < sys.n_solvent; j++) {
      dr = vdist(part_dpd[i].r, part_dpd[j].r);
      periodic_bc_dr(&dr);
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
      periodic_bc_dr(&dr);
      r_ij = vmag(dr);

      // Conservative force
      if (r_ij < sys.r_c) {
        P.now += sys.a_mm*r_ij*(1-r_ij);
      }

      // FENE spring force
      if (j == i+1) {
        fact = 1.0 - (r_ij - sys.r_eq) / sys.r_0;
        P.now += -sys.k_fene*r_ij*(r_ij-sys.r_eq) / fact;
      }
    }

    // Contribution from monomer-solvent forces
    for (j = 0; j < sys.n_solvent; j++) {
      dr = vdist(part_mon[i].r, part_dpd[j].r);
      periodic_bc_dr(&dr);

      r_ij = vmag(dr);

      if (r_ij < sys.r_c) {
        P.now += sys.a_ms*r_ij*(1-r_ij);
      }
    }
  }

  // Contribution from wall and pore particles
  for (i = sys.n_solvent; i < sys.n_dpd; i++) {
    // Wall-solvent interaction
    for (j = 0; j < sys.n_solvent; j++) {
      dr = vdist(part_dpd[i].r, part_dpd[j].r);
      periodic_bc_dr(&dr);

      r_ij = vmag(dr);

      if (r_ij < sys.r_c) {
        P.now += sys.a_sw*r_ij*(1-r_ij);
      }
    }

    // Wall-monomer interaction
    for (j = 0; j < sys.n_mon; j++) {
      dr = vdist(part_dpd[i].r, part_mon[j].r);
      periodic_bc_dr(&dr);

      r_ij = vmag(dr);

      if (r_ij < sys.r_c) {
        P.now += sys.a_sw*r_ij*(1-r_ij);
      }
    }
  }

  P.now /= 3*sys.volume;
  P.now += sys.density_s*sys.temp;
}

void calc_re(void) {
  int i;
  Vector dr;

  dr = vdist(part_mon[sys.n_mon-1].r, part_mon[0].r);

  RE2x.now = dr.x * dr.x;
  RE2y.now = dr.y * dr.y;
  RE2z.now = dr.z * dr.z;
  RE2.now = RE2x.now + RE2y.now + RE2z.now;
}

void calc_rg(void) {
  int i;
  Vector dr;

  RG2x.now = 0;
  RG2y.now = 0;
  RG2z.now = 0;

  for (i = 0; i < sys.n_mon; i++) {
    dr.x = part_mon[i].r.x - sys.mon.cmx[sys.monitor_step];
    dr.y = part_mon[i].r.y - sys.mon.cmy[sys.monitor_step];
    dr.z = part_mon[i].r.z - sys.mon.cmz[sys.monitor_step];

    RG2x.now += dr.x * dr.x;
    RG2y.now += dr.y * dr.y;
    RG2z.now += dr.z * dr.z;
  }

  RG2x.now /= sys.n_mon;
  RG2y.now /= sys.n_mon;
  RG2z.now /= sys.n_mon;
  RG2.now = RG2x.now + RG2y.now + RG2z.now;
}

void calc_bond_length(void) {
  int i;
  Vector dr;

  for (i = 1; i < sys.n_mon; i++) {
     dr = vdist(part_mon[i].r, part_mon[i-1].r);
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
    periodic_bc_dr(&dr);
    r_ij = vmag(dr);

    if (r_ij > sys.r_max) {
      sys.bond_break = 1;
    }
  }

  // If not the last monomer in the chain
  if (i != (sys.n_mon - 1)) {
    dr = vdist(part_mon[i].r, part_mon[i+1].r);
    periodic_bc_dr(&dr);
    r_ij = vmag(dr);

    if (r_ij > sys.r_max) {
      sys.bond_break = 1;
    }
  }
}

void check_wall(Vector r) {
  sys.wall_overlap = 0;

  if (r.z >= sys.wall_min.z && r.z <= sys.wall_max.z) {
    if (r.x >= sys.pore_max.x || r.x <= sys.pore_min.x) {
      sys.wall_overlap = 1;
    } else if (r.y >= sys.pore_max.y || r.y <= sys.pore_min.y) {
      sys.wall_overlap = 1;
    }
  }
}

int check_pore(Vector r) {
  sys.pore_overlap = 0;

  if (r.z >= sys.pore_min.z && r.z <= sys.pore_max.z) {
    if ((r.x <= sys.pore_min.x && r.x >= sys.pore_min.x-(sys.n_layers*sys.r_wall))
      || (r.x >= sys.pore_max.x && r.x <= sys.pore_max.x+(sys.n_layers*sys.r_wall))) {
      if ((r.y <= sys.pore_min.y && r.y >= sys.pore_min.y-(sys.n_layers*sys.r_wall))
        || (r.y >= sys.pore_max.y && r.y <= sys.pore_max.y+(sys.n_layers*sys.r_wall))) {
        sys.pore_overlap = 1;
      }
    }
  }

  if (!sys.pore_overlap) {
    if ((r.x > sys.pore_min.x && r.x < sys.pore_max.x)
      && (r.y > sys.pore_min.y && r.y < sys.pore_min.y)
      && (r.z > sys.pore_min.z && r.z < sys.pore_max.z)) {
      return 1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}
