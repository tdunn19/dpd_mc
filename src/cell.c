/* cell.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dpd.h"

void new_list(void) {
  int i, ix, iy, iz, n;


  // Initialize head of chain for each cell.
  // -1 indiciates the end of the chain.
  for (ix = 0; ix < sys.n_cell; ix++) {
    for (iy = 0; iy < sys.n_cell; iy++) {
      for (iz = 0; iz < sys.n_cell; iz++) {
        sys.hoc[ix][iy][iz] = -1;
      }
    }
  }

  for (i = 0; i < sys.n_dpd; i++) {
    // Determine cell number of the particle
    ix = (int) part_dpd[i].r.x / sys.r_cell;
    iy = (int) part_dpd[i].r.y / sys.r_cell;
    iz = (int) part_dpd[i].r.z / sys.r_cell;

    // Link list the head of chain of cell i,j,k
    part_dpd[i].ll = sys.hoc[ix][iy][iz];

    // Make particle i the new head of chain
    sys.hoc[ix][iy][iz] = i;
  }
}

int check_cell(Vector r, Vector ro) {
  // Check for a new head of chain
  int ix, iy, iz, ixo, iyo, izo;

  ix = (int) r.x / sys.r_cell;
  iy = (int) r.y / sys.r_cell;
  iz = (int) r.z / sys.r_cell;

  ixo = (int) ro.x / sys.r_cell;
  iyo = (int) ro.y / sys.r_cell;
  izo = (int) ro.z / sys.r_cell;

  if (sys.hoc[ix][iy][iz] != sys.hoc[ixo][iyo][izo]) {
    // The cell has changed
    return 1;
  } else {
    return 0;
  }
}

double calc_energy_dpd(int i) {
  int ix, iy, iz, jx, jy, jz, j, l, m, n;
  double E_ss = 0, E_ms = 0;
  Vector dr;

  // Determine cell number of particle i
  ix = (int) part_dpd[i].r.x / sys.r_cell;
  iy = (int) part_dpd[i].r.y / sys.r_cell;
  iz = (int) part_dpd[i].r.z / sys.r_cell;

  // Solvent-solvent interaction
  for (l = -1; l <= 1; l++) {
    for (m = -1; m <= 1; m++) {
      for (n = -1; n <= 1; n++) {
        // Determine nearest neighbor cell
        jx = mod(ix+l, sys.n_cell);
        jy = mod(iy+m, sys.n_cell);
        jz = mod(iz+n, sys.n_cell);

        // First particle in the chain
        j = sys.hoc[jx][jy][jz];

        while (j != -1) {
          if (i != j) {
            dr = vdist(part_dpd[i].r, part_dpd[j].r);
            E_ss += energy_c(dr);
          }
          // Next particle in the chain
          j = part_dpd[j].ll;
        }
      }
    }
  }
  E_ss *= sys.a_ss;

  // Solvent-monomer interaction
  for (j = 0; j < sys.n_mon; j++) {
    dr = vdist(part_dpd[i].r, part_mon[j].r);
    E_ms += energy_c(dr);
  }
  E_ms *= sys.a_ms;

  return E_ss + E_ms;
}


double calc_energy_mon(int i) {
  int ix, iy, iz, jx, jy, jz, j, l, m, n;
  double E_ms = 0, E_mm = 0, E_fene = 0;;
  Vector dr;

  // If not the first monomer in the chain
  if (i != 0) {
    E_fene += energy_fene(i, i-1);
  }
  // If not the last monomer in the chain
  if (i != (sys.n_mon -1)) {
    E_fene += energy_fene(i, i+1);
  }
  E_fene *= sys.k_fene / 2;

  // Check for bond break
  if (!sys.bond_break) {
    // Determine cell number of monomer i
    ix = (int) part_mon[i].r.x / sys.r_cell;
    iy = (int) part_mon[i].r.y / sys.r_cell;
    iz = (int) part_mon[i].r.z / sys.r_cell;

    // Monomer-solvent interaction
    for (l = -1; l <= 1; l++) {
      for (m = -1; m <= 1; m++) {
        for (n = -1; n <= 1; n++) {
          // Determine nearest neighbor cell
          jx = mod(ix+l, sys.n_cell);
          jy = mod(iy+m, sys.n_cell);
          jz = mod(iz+n, sys.n_cell);

          // First particle in the chain
          j = sys.hoc[jx][jy][jz];

          while (j != -1) {
            dr = vdist(part_mon[i].r, part_dpd[j].r);
            E_ms += energy_c(dr);

            // Next particle in the chain
            j = part_dpd[j].ll;
          }
        }
      }
    }
    E_ms *= sys.a_ms;

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
