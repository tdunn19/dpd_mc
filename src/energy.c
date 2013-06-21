#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"


double calc_energy_dpd(int i) {
  int ix, iy, iz, jx, jy, jz, j, l, m, n;
  double E_ss = 0, E_ms = 0, E_sw = 0;
  Vector dr;

  // Determine cell number of particle i
  ix = (int) (part_dpd[i].r.x / sys.r_cell.x);
  iy = (int) (part_dpd[i].r.y / sys.r_cell.y);
  iz = (int) (part_dpd[i].r.z / sys.r_cell.z);

  // DPD-DPD interaction
  for (l = -1; l <= 1; l++) {
    for (m = -1; m <= 1; m++) {
      for (n = -1; n <= 1; n++) {
        // Determine nearest neighbor cell
        jx = mod(ix+l, sys.n_cell_1d.x);
        jy = mod(iy+m, sys.n_cell_1d.y);
        jz = mod(iz+n, sys.n_cell_1d.y);

        // First particle in the chain
        j = sys.hoc[jx][jy][jz];

        while (j != -1) {
          if (i != j) {
            dr = vdist(part_dpd[i].r, part_dpd[j].r);

            if (i < sys.n_solvent) {
              if (j < sys.n_solvent) {
                // Solvent-solvent interaction
                E_ss += energy_c(dr):
              } else {
                // Solvent-wall interaction
                E_sw += energy_c(dr);
              }
            } else {
              if (j < sys.n_solvent) {
                // Wall-solvent interaction
                E_sw += energy_c(dr);
              }
            }
          }
          // Next particle in the chain
          j = part_dpd[j].ll;
        }
      }
    }
  }
  E_ss *= sys.a_ss;
  E_sw *= sys.a_sw;

  // Solvent-monomer interaction
  for (j = 0; j < sys.n_mon; j++) {
    dr = vdist(part_dpd[i].r, part_mon[j].r);
    E_ms += energy_c(dr);
  }
  E_ms *= sys.a_ms;

  return E_ss + E_ms + E_sw;
}


double calc_energy_mon(int i) {
  int ix, iy, iz, jx, jy, jz, j, l, m, n;
  double E_ms = 0, E_mm = 0, E_fene = 0, E_mw = 0;
  Vector dr;

  // If not the first monomer in the chain
  if (i != 0) {
    E_fene += energy_fene(i, i-1);
  }
  // If not the last monomer in the chain
  if (i != (sys.n_mon-1)) {
    E_fene += energy_fene(i, i+1);
  }
  E_fene *= sys.k_fene / 2;

  // Check for bond break
  if (!sys.bond_break) {
    // Determine cell number of monomer i
    ix = (int) (part_mon[i].r.x / sys.r_cell.x);
    iy = (int) (part_mon[i].r.y / sys.r_cell.y);
    iz = (int) (part_mon[i].r.z / sys.r_cell.z);

    // Monomer-DPD interaction
    for (l = -1; l <= 1; l++) {
      for (m = -1; m <= 1; m++) {
        for (n = -1; n <= 1; n++) {
          // Determine nearest neighbor cell
          jx = mod(ix+l, sys.n_cell_1d.x);
          jy = mod(iy+m, sys.n_cell_1d.y);
          jz = mod(iz+n, sys.n_cell_1d.z);

          // First particle in the chain
          j = sys.hoc[jx][jy][jz];

          while (j != -1) {
            dr = vdist(part_mon[i].r, part_dpd[j].r);
            if (j < sys.n_solvent) {
              // Monomer-solvent interaction
              E_ms += energy_c(dr);
            } else {
              // Monomer-wall interaction
              E_mw += energy_c(dr);
            }
            // Next particle in the chain
            j = part_dpd[j].ll;
          }
        }
      }
    }
    E_ms *= sys.a_ms;
    E_mw *= sys.a_sw;

    // Monomer-monomer interaction
    for (j = 0; j < sys.n_mon; j++) {
      if (i != j) {
        dr = vdist(part_mon[i].r, part_mon[j].r);
        E_mm += energy_c(dr);
      }
    }
    E_mm *= sys.a_mm;
  }

  return E_ms + E_mm + E_fene;
}

void calc_energy_brute(void) {
  int i, j;
  double E_mm, E_ms, E_ss, E_fene, E_mw, E_sw;
  Vector dr;

  // Monomer energies
  for (i = 0; i < sys.n_mon; i++) {
    E_mm = 0;
    E_ms = 0;
    E_fene = 0;
    E_mw = 0;

    part_mon[i].Eo = part_mon[i].E;

    // If not the first monomer in the chain
    if (i != 0) {
        E_fene =+ energy_fene(i, i-1);
    }
    // If not the last monomer in the chain
    if (i != (sys.n_mon - 1)) {
        E_fene += energy_fene(i, i+1);
    }
    E_fene *= sys.k_fene / 2;

    for (j = 0; j < sys.n_mon; j++) {
      if (i != j) {
        dr = vdist(part_mon[i].r, part_mon[j].r);
        E_mm += energy_c(dr);
      }
    }
    E_mm *= sys.a_mm;

    for (j = 0; j < sys.n_solvent; j++) {
        dr = vdist(part_mon[i].r, part_dpd[j].r);
        E_ms += energy_c(dr);
    }
    E_ms *= sys.a_ms;

    for (j = sys.n_solvent; j < sys.n_dpd; j++) {
        dr = vdist(part_mon[i].r, part_dpd[j].r);
        E_mw += energy_c(dr);
    }
    E_mw *= sys.a_sw;

    part_mon[i].E = E_fene + E_mm + E_ms + E_mw;
  }

  // Solvent particle energies
  for (i = 0; i < sys.n_solvent; i++) {
    E_ss = 0;
    E_ms = 0;
    E_sw = 0;

    part_dpd[i].Eo = part_dpd[i].E;

    for (j = 0; j < sys.n_solvent; j++) {
      if (j != i) {
        dr = vdist(part_dpd[i].r, part_dpd[j].r);
        E_ss += energy_c(dr);
      }
    }
    E_ss *= sys.a_ss;

    for (j = sys.n_solvent; j < sys.n_dpd; j++) {
      dr = vdist(part_dpd[i].r, part_dpd[j].r);
      E_sw += energy_c(dr);
    }
    E_sw *= sys.a_sw;

    for (j = 0; j < sys.n_mon; j++) {
      dr = vdist(part_dpd[i].r, part_mon[j].r);
      E_ms += energy_c(dr);
    }
    E_ms *= sys.a_ms;

    part_dpd[i].E = E_ss + E_ms + E_sw;
  }
}

double energy_c(Vector dr) {
  double r_ij, E_ij;

  periodic_bc_dr(&dr);

  r_ij = vmag(dr);

  // Soft repulsive force
  if (r_ij < sys.r_c) {
    E_ij = (1 - r_ij/sys.r_c) * (1 - r_ij/sys.r_c) / 2;
  } else {
    E_ij = 0;
  }

  return E_ij;
}

double energy_fene(i, j) {
  double E, r_ij, r2;
  Vector dr;

  dr = vdist(part_mon[i].r, part_mon[j].r);
  periodic_bc_dr(&dr);
  r_ij = vmag(dr);

  if (r_ij <= sys.r_max) {
    r2 = sys.r_0 * sys.r_0;
    E = -1.0 * r2 * log(1 - (r_ij - sys.r_eq)*(r_ij - sys.r_eq)/r2);
    return E;
  }
  else {
    sys.bond_break = 1;
    return 0;
  }
}

double total_energy(void) {
  int i;
  double E;
  E = 0;

  for (i = 0; i < sys.n_solvent; i++) {
    E += part_dpd[i].E;
  }

  for (i = 0; i < sys.n_mon; i++) {
    E += part_mon[i].E;
  }

  return E;
}

