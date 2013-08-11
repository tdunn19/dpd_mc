#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void sample(void) {
  int i;
  Vector dr;

  // Get absolute position of the monomers without periodic boundary conditions
  for (i = 1; i < sys.n_mon; i++) {
    dr = vdist(part_mon[i].r, part_mon[i-1].r);
    periodic_bc_dr(&dr);

    part_mon[i].r.x = part_mon[i-1].r.x + dr.x;
    part_mon[i].r.y = part_mon[i-1].r.y + dr.y;
    part_mon[i].r.z = part_mon[i-1].r.z + dr.z;
  }

  calc_nseg();
  calc_q();
  calc_cm();
  calc_re();
  calc_rg();

  // Update measured quantities
  update_stats();
  update_monitor();

  // Revert the monomer positions
  for (i = 1; i < sys.n_mon; i++) {
    part_mon[i].r.x = part_mon[i].ro.x;
    part_mon[i].r.y = part_mon[i].ro.y;
    part_mon[i].r.z = part_mon[i].ro.z;
  }

}

void calc_nseg(void) {
  int i;

  sys.n_cis = 0;
  sys.n_trans = 0;

  for (i = 0; i < sys.n_mon; i++) {
    if (part_mon[i].r.z > sys.pore_max.z) {
      sys.n_trans++;
    } else if (part_mon[i].r.z < sys.pore_min.z) {
      sys.n_cis++;
    }
  }
}

void calc_q(void) {
  int i;

  sys.Q = 0;

  for (i = 0; i < sys.n_mon; i++) {
    if (part_mon[i].r.z > sys.pore_max.z) {
      sys.Q += 1.0;
    } else if (part_mon[i].r.z > sys.pore_min.z) {
      sys.Q += (part_mon[i].r.z-sys.pore_min.z) / sys.pore_length;
    }
  }

  sys.Q /= sys.n_mon;
}

void calc_cm(void) {
  int i;

  sys.cm.x = 0;
  sys.cm.y = 0;
  sys.cm.z = 0;

  for (i = 0; i < sys.n_mon; i++) {
    sys.cm.x += part_mon[i].r.x;
    sys.cm.y += part_mon[i].r.y;
    sys.cm.z += part_mon[i].r.z;
  }

  sys.cm.x /= sys.n_mon;
  sys.cm.y /= sys.n_mon;
  sys.cm.z /= sys.n_mon;

  sys.cm_cis.x = 0;
  sys.cm_cis.y = 0;
  sys.cm_cis.z = 0;

  if (sys.n_cis > 0) {
    for (i = sys.n_mon; i > sys.n_mon-sys.n_cis; i--) {
      sys.cm_cis.x += part_mon[i].r.x;
      sys.cm_cis.y += part_mon[i].r.y;
      sys.cm_cis.z += part_mon[i].r.z;
    }
    sys.cm_cis.x /= sys.n_cis;
    sys.cm_cis.y /= sys.n_cis;
    sys.cm_cis.z /= sys.n_cis;
  }

  sys.cm_trans.x = 0;
  sys.cm_trans.y = 0;
  sys.cm_trans.z = 0;

  if (sys.n_trans > 0) {
    for (i = 0; i < sys.n_trans; i++) {
      sys.cm_trans.x += part_mon[i].r.x;
      sys.cm_trans.y += part_mon[i].r.y;
      sys.cm_trans.z += part_mon[i].r.z;
    }
    sys.cm_trans.x /= sys.n_trans;
    sys.cm_trans.y /= sys.n_trans;
    sys.cm_trans.z /= sys.n_trans;
  }
}

void calc_re(void) {
  int i;
  Vector dr;

  dr = vdist(part_mon[sys.n_mon-1].r, part_mon[0].r);

  sys.re.x = dr.x;
  sys.re.y = dr.y;
  sys.re.z = dr.z;

  if (sys.n_cis > 0) {
    dr = vdist(part_mon[sys.n_mon-1].r, part_mon[sys.n_mon-sys.n_cis].r);

    sys.re_cis.x = dr.x;
    sys.re_cis.y = dr.y;
    sys.re_cis.z = dr.z;
  } else {
    sys.re_cis.x = 0;
    sys.re_cis.y = 0;
    sys.re_cis.z = 0;
  }

  if (sys.n_trans > 0) {
    dr = vdist(part_mon[0].r, part_mon[sys.n_trans+1].r);

    sys.re_trans.x = dr.x;
    sys.re_trans.y = dr.y;
    sys.re_trans.z = dr.z;
  } else {
    sys.re_trans.x = 0;
    sys.re_trans.y = 0;
    sys.re_trans.z = 0;
  }

}

void calc_rg(void) {
  int i;
  Vector dr;

  sys.rg2.x = 0;
  sys.rg2.y = 0;
  sys.rg2.z = 0;

  for (i = 0; i < sys.n_mon; i++) {
    dr.x = part_mon[i].r.x - sys.cm.x;
    dr.y = part_mon[i].r.y - sys.cm.y;
    dr.z = part_mon[i].r.z - sys.cm.z;

    sys.rg2.x += dr.x * dr.x;
    sys.rg2.y += dr.y * dr.y;
    sys.rg2.z += dr.z * dr.z;
  }

  sys.rg2.x /= sys.n_mon;
  sys.rg2.y /= sys.n_mon;
  sys.rg2.z /= sys.n_mon;

  sys.rg2_cis.x = 0;
  sys.rg2_cis.y = 0;
  sys.rg2_cis.z = 0;

  if (sys.n_cis > 0) {
    for (i = sys.n_mon; i > sys.n_mon-sys.n_cis; i--) {
      dr.x = part_mon[i].r.x - sys.cm_cis.x;
      dr.y = part_mon[i].r.y - sys.cm_cis.y;
      dr.z = part_mon[i].r.z - sys.cm_cis.z;

      sys.rg2_cis.x += dr.x * dr.x;
      sys.rg2_cis.y += dr.y * dr.y;
      sys.rg2_cis.z += dr.z * dr.z;
    }

    sys.rg2_cis.x /= sys.n_cis;
    sys.rg2_cis.y /= sys.n_cis;
    sys.rg2_cis.z /= sys.n_cis;
  }

  sys.rg2_trans.x = 0;
  sys.rg2_trans.y = 0;
  sys.rg2_trans.z = 0;

  if (sys.n_trans > 0) {
    for (i = 0; i < sys.n_trans; i++) {
      dr.x = part_mon[i].r.x - sys.cm_trans.x;
      dr.y = part_mon[i].r.y - sys.cm_trans.y;
      dr.z = part_mon[i].r.z - sys.cm_trans.z;

      sys.rg2_trans.x += dr.x * dr.x;
      sys.rg2_trans.y += dr.y * dr.y;
      sys.rg2_trans.z += dr.z * dr.z;
    }

    sys.rg2_trans.x /= sys.n_trans;
    sys.rg2_trans.y /= sys.n_trans;
    sys.rg2_trans.z /= sys.n_trans;
  }
}

void update_stats(void) {
  int i;

  Etot.now = sys.energy;
  RE2.now = sys.re.x*sys.re.x + sys.re.y*sys.re.y + sys.re.z*sys.re.z;
  RE2x.now = sys.re.x * sys.re.x;
  RE2y.now = sys.re.y * sys.re.y;
  RE2z.now = sys.re.z * sys.re.z;
  RG2.now = sys.rg2.x + sys.rg2.y + sys.rg2.z;
  RG2x.now = sys.rg2.x;
  RG2y.now = sys.rg2.y;
  RG2z.now = sys.rg2.z;

  for (i = 0; i < sys.n_stats; i++) {
    sys.stats[i].sum += sys.stats[i].now;
    sys.stats[i].sumsq += sys.stats[i].now*sys.stats[i].now;
    sys.stats[i].num += 1;
  }
}

void update_monitor(void) {
  sys.mon.energy[sys.monitor_step] = sys.energy;
  sys.mon.Q[sys.monitor_step] = sys.Q;
  sys.mon.cmx[sys.monitor_step] = sys.cm.x;
  sys.mon.cmy[sys.monitor_step] = sys.cm.y;
  sys.mon.cmz[sys.monitor_step] = sys.cm.z;
  sys.mon.re2[sys.monitor_step] = sys.re.x*sys.re.x + sys.re.y*sys.re.y + sys.re.z*sys.re.z;
  sys.mon.rex[sys.monitor_step] = sys.re.x;
  sys.mon.rey[sys.monitor_step] = sys.re.y;
  sys.mon.rez[sys.monitor_step] = sys.re.z;
  sys.mon.re2_cis[sys.monitor_step] = sys.re_cis.x*sys.re_cis.x + sys.re_cis.y*sys.re_cis.y + sys.re_cis.z*sys.re_cis.z;
  sys.mon.rex_cis[sys.monitor_step] = sys.re_cis.x;
  sys.mon.rey_cis[sys.monitor_step] = sys.re_cis.y;
  sys.mon.rez_cis[sys.monitor_step] = sys.re_cis.z;
  sys.mon.re2_trans[sys.monitor_step] = sys.re_trans.x*sys.re_trans.x + sys.re_trans.y*sys.re_trans.y + sys.re_trans.z*sys.re_trans.z;
  sys.mon.rex_trans[sys.monitor_step] = sys.re_trans.x;
  sys.mon.rey_trans[sys.monitor_step] = sys.re_trans.y;
  sys.mon.rez_trans[sys.monitor_step] = sys.re_trans.z;
  sys.mon.rg2[sys.monitor_step] = sys.rg2.x + sys.rg2.y + sys.rg2.z;
  sys.mon.rg2x[sys.monitor_step] = sys.rg2.x;
  sys.mon.rg2y[sys.monitor_step] = sys.rg2.y;
  sys.mon.rg2z[sys.monitor_step] = sys.rg2.z;
  sys.mon.rg2_cis[sys.monitor_step] = sys.rg2_cis.x + sys.rg2_cis.y + sys.rg2_cis.z;
  sys.mon.rg2x_cis[sys.monitor_step] = sys.rg2_cis.x;
  sys.mon.rg2y_cis[sys.monitor_step] = sys.rg2_cis.y;
  sys.mon.rg2z_cis[sys.monitor_step] = sys.rg2_cis.z;
  sys.mon.rg2_trans[sys.monitor_step] = sys.rg2_trans.x + sys.rg2_trans.y + sys.rg2_trans.z;
  sys.mon.rg2x_trans[sys.monitor_step] = sys.rg2_trans.x;
  sys.mon.rg2y_trans[sys.monitor_step] = sys.rg2_trans.y;
  sys.mon.rg2z_trans[sys.monitor_step] = sys.rg2_trans.z;

  sys.monitor_step += 1;
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

    if (r_ij > sys.r_max2) {
      sys.bond_break = 1;
    }
  }

  // If not the last monomer in the chain
  if (i != (sys.n_mon - 1)) {
    dr = vdist(part_mon[i].r, part_mon[i+1].r);
    periodic_bc_dr(&dr);
    r_ij = vmag(dr);

    if (r_ij > sys.r_max2) {
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

  // Pore wall overlap
  if (r.z >= sys.pore_min.z && r.z <= sys.pore_max.z) {
    if ((r.x <= sys.pore_min.x && r.x >= sys.pore_min.x-(sys.n_layers*sys.r_wall))
      || (r.x >= sys.pore_max.x && r.x <= sys.pore_max.x+(sys.n_layers*sys.r_wall))) {
      if (r.y <= sys.pore_max.y+(sys.n_layers*sys.r_wall)
        && r.y >= sys.pore_min.y-(sys.n_layers*sys.r_wall)) {
        sys.wall_overlap = 1;
      }
    }

    if ((r.y <= sys.pore_min.y && r.y >= sys.pore_min.y-(sys.n_layers*sys.r_wall))
      || (r.y >= sys.pore_max.y && r.y <= sys.pore_max.y+(sys.n_layers*sys.r_wall))) {
      if (r.x <= sys.pore_max.x+(sys.n_layers*sys.r_wall)
        && r.x >= sys.pore_min.x-(sys.n_layers*sys.r_wall)) {
        sys.wall_overlap = 1;
      }
    }
  }
}

int check_pore(Vector r) {
  if ((r.x > sys.pore_min.x && r.x < sys.pore_max.x)
    && (r.y > sys.pore_min.y && r.y < sys.pore_max.y)
    && (r.z > sys.pore_min.z && r.z < sys.pore_max.z)) {
    return 1;
  } else {
    return 0;
  }
}

void check_window() {
  sys.new_window = 0;

  calc_q();

  if (sys.Q < sys.Q_min + 1e-12) {
    sys.new_window = 1;
  } else if (sys.Q > sys.Q_max - 1e-12) {
    sys.new_window = 1;
  }
}

int check_side(Vector r) {
  if (r.z > sys.length.z/2) {
    return 1;
  } else {
    return 0;
  }
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
        P.now += sys.a_ms[part_dpd[j].side]*r_ij*(1-r_ij);
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
  P.now += sys.density_s;
  // P.now += sys.density_s*sys.temp;
}
