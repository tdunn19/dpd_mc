/* mc.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void monte_carlo(void) {
  int i, j, ix, iy, iz, ixo, iyo, izo;

  if (ran3() < sys.mc_ratio) {
    // A monomer was chosen
    sys.n_attempt_mon += 1;

    // Randomly move a monomer
    i = rand() % sys.n_mon;
    random_move_mon(i);

    if (!accept_move()) {
      part_mon[i].r.x = part_mon[i].ro.x;
      part_mon[i].r.y = part_mon[i].ro.y;
      part_mon[i].r.z = part_mon[i].ro.z;

      for (j = 0; j < sys.n_mon; j++) {
        part_mon[j].E = part_mon[j].Eo;
      }

      for (j = 0; j < sys.n_dpd; j++) {
        part_dpd[j].E = part_dpd[j].Eo;
      }
    } else {
      // Move was accepeted
      sys.n_accept_mon += 1;

      part_mon[i].ro.x = part_mon[i].r.x;
      part_mon[i].ro.y = part_mon[i].r.y;
      part_mon[i].ro.z = part_mon[i].r.z;

      for (j = 0; j < sys.n_mon; j++) {
        part_mon[j].Eo = part_mon[j].E;
      }

      for (j = 0; j < sys.n_dpd; j++) {
        part_dpd[j].Eo = part_dpd[j].E;
      }
    }
  } else {
    // A DPD particle was chosen
    sys.n_attempt_dpd += 1;

    // Randomly move a DPD particle
    i = rand() % sys.n_dpd;
    random_move_dpd(i);

    if (!accept_move()) {
      part_dpd[i].r.x = part_dpd[i].ro.x;
      part_dpd[i].r.y = part_dpd[i].ro.y;
      part_dpd[i].r.z = part_dpd[i].ro.z;

      for (j = 0; j < sys.n_dpd; j++) {
        part_dpd[j].E = part_dpd[j].Eo;
      }

      for (j = 0; j < sys.n_mon; j++) {
        part_mon[j].E = part_mon[j].Eo;
      }
    } else {
      // The move was accepted
      sys.n_accept_dpd += 1;

      part_dpd[i].ro.x = part_dpd[i].r.x;
      part_dpd[i].ro.y = part_dpd[i].r.y;
      part_dpd[i].ro.z = part_dpd[i].r.z;

      for (j = 0; j < sys.n_mon; j++) {
        part_mon[j].Eo = part_mon[j].E;
      }

      for (j = 0; j < sys.n_dpd; j++) {
        part_dpd[j].Eo = part_dpd[j].E;
      }

      // Generate a new list if the particle changed cells
      if (sys.calc_list && check_cell(part_dpd[i].r, part_dpd[i].ro)) {
        new_list();
      }
    }
  }
}

void random_move_dpd(int i) {
  int ix, ixo, iy, iyo, iz, izo, j, jx, jy, jz, dix, diy, diz, l, m, n;

  part_dpd[i].ro.x = part_dpd[i].r.x;
  part_dpd[i].ro.y = part_dpd[i].r.y;
  part_dpd[i].ro.z = part_dpd[i].r.z;

  // Displacement between -dr_max and +dr_max in each dimension
  part_dpd[i].r.x += sys.dr_max_dpd * (2*ran3() - 1);
  part_dpd[i].r.y += sys.dr_max_dpd * (2*ran3() - 1);
  part_dpd[i].r.z += sys.dr_max_dpd * (2*ran3() - 1);

  // Periodic boundary conditions
  if (part_dpd[i].r.x > sys.length) {
    part_dpd[i].r.x -= sys.length;
  } else if (part_dpd[i].r.x < 0) {
    part_dpd[i].r.x += sys.length;
  }

  if (part_dpd[i].r.y > sys.length) {
    part_dpd[i].r.y -= sys.length;
  } else if (part_dpd[i].r.y < 0) {
    part_dpd[i].r.y += sys.length;
  }

  if (part_dpd[i].r.z > sys.length) {
    part_dpd[i].r.z -= sys.length;
  } else if (part_dpd[i].r.z < 0) {
    part_dpd[i].r.z += sys.length;
  }

  if (sys.calc_list) {
    for (j = 0; j < sys.n_mon; j++) {
      part_mon[j].Eo = part_mon[j].E;
      part_mon[j].E = calc_energy_mon(j);
    }

    // Determine cell number of particle i
    ix = (int) part_dpd[i].r.x / sys.r_cell;
    iy = (int) part_dpd[i].r.y / sys.r_cell;
    iz = (int) part_dpd[i].r.z / sys.r_cell;

    for (l = -1; l <= 1; l++) {
      for (m = -1; m <= 1; m++) {
        for (n = -1; n <= 1; n++) {
          // Determine the cell
          jx = mod(ix+l, sys.n_cell);
          jy = mod(iy+m, sys.n_cell);
          jz = mod(iz+n, sys.n_cell);

          j = sys.hoc[jx][jy][jz];

          while (j != -1) {
            part_dpd[j].Eo = part_dpd[j].E;
            part_dpd[j].E = calc_energy_dpd(j);

            // Next particle in the linked list
            j = part_dpd[j].ll;
          }
        }
      }
    }

    if (check_cell(part_dpd[i].r, part_dpd[i].ro)) {
      // The particle entered a new cell
      // Must consider the neighbors of the old cell
      ixo = (int) part_dpd[i].ro.x / sys.r_cell;
      iyo = (int) part_dpd[i].ro.y / sys.r_cell;
      izo = (int) part_dpd[i].ro.z / sys.r_cell;

      dix = ix - ixo;
      diy = iy - iyo;
      diz = iz - izo;

      if (dix != 0) {
        jx = mod(ixo - dix, sys.n_cell);

        for (m = -1; m <= 1; m++) {
          for (n = -1; n <= 1; n++) {
            jy = mod(iyo+m, sys.n_cell);
            jz = mod(izo+n, sys.n_cell);

            j = sys.hoc_copy[jx][jy][jz];

            while (j != -1) {
              part_dpd[j].Eo = part_dpd[j].E;
              part_dpd[j].E = calc_energy_dpd(j);
              j = part_dpd[j].ll;
            }
            // To account for overcounting, set HOC to -1
            sys.hoc_copy[jx][jy][jz] = -1;
          }
        }
      }

      if (diy != 0) {
        jy = mod(iyo - diy, sys.n_cell);

        for (l = -1; l <= 1; l++) {
          for (n = -1; n <= 1; n++) {
            jx = mod(ixo+l, sys.n_cell);
            jz = mod(izo+n, sys.n_cell);

            j = sys.hoc_copy[jx][jy][jz];

            while (j != -1) {
              part_dpd[j].Eo = part_dpd[j].E;
              part_dpd[j].E = calc_energy_dpd(j);
              j = part_dpd[j].ll;
            }
            // To account for overcounting, set HOC to -1
            sys.hoc_copy[jx][jy][jz] = -1;
          }
        }
      }

      if (diz != 0) {
        jz = mod(izo - diz, sys.n_cell);

        for (l = -1; l <= 1; l++) {
          for (m = -1; m <= 1; m++) {
            jx = mod(ixo+l, sys.n_cell);
            jy = mod(iyo+m, sys.n_cell);

            j = sys.hoc_copy[jx][jy][jz];

            while (j != -1) {
              part_dpd[j].Eo = part_dpd[j].E;
              part_dpd[j].E = calc_energy_dpd(j);
              j = part_dpd[j].ll;
            }
          }
        }
      }
    }
  } else {
    calc_energy_brute();
  }
}

void random_move_mon(int i) {
  int ix, ixo, iy, iyo, iz, izo, j, jx, jy, jz, dix, diy, diz, l, m, n;

  part_mon[i].ro.x = part_mon[i].r.x;
  part_mon[i].ro.y = part_mon[i].r.y;
  part_mon[i].ro.z = part_mon[i].r.z;

  // Displacement between -dr_max and +dr_max in each dimension
  part_mon[i].r.x += sys.dr_max_mon * (2*ran3() - 1);
  part_mon[i].r.y += sys.dr_max_mon * (2*ran3() - 1);
  part_mon[i].r.z += sys.dr_max_mon * (2*ran3() - 1);

  // Periodic boundary conditions
  if (part_mon[i].r.x > sys.length) {
    part_mon[i].r.x -= sys.length;
  } else if (part_mon[i].r.x < 0) {
    part_mon[i].r.x += sys.length;
  }

  if (part_mon[i].r.y > sys.length) {
    part_mon[i].r.y -= sys.length;
  } else if (part_mon[i].r.y < 0) {
    part_mon[i].r.y += sys.length;
  }

  if (part_mon[i].r.z > sys.length) {
    part_mon[i].r.z -= sys.length;
  } else if (part_mon[i].r.z < 0) {
    part_mon[i].r.z += sys.length;
  }

  // Check for bond breaks
  check_bond(i);

  if (sys.calc_list && !sys.bond_break) {
    for (j = 0; j < sys.n_mon; j++) {
      part_mon[j].Eo = part_mon[j].E;
      part_mon[j].E = calc_energy_mon(j);
    }

    // Determine new energies for nearest neighbor cells
    ix = (int) part_mon[i].r.x / sys.r_cell;
    iy = (int) part_mon[i].r.y / sys.r_cell;
    iz = (int) part_mon[i].r.z / sys.r_cell;

    for (l = -1; l <= 1; l++) {
      for (m = -1; m <= 1; m++) {
        for (n = -1; n <= 1; n++) {
          // Determine the cell
          jx = mod(ix+l, sys.n_cell);
          jy = mod(iy+m, sys.n_cell);
          jz = mod(iz+n, sys.n_cell);

          j = sys.hoc[jx][jy][jz];

          while (j != -1) {
            part_dpd[j].Eo = part_dpd[j].E;
            part_dpd[j].E = calc_energy_dpd(j);

            // Next particle in the linked list
            j = part_dpd[j].ll;
          }
        }
      }
    }

    if (check_cell(part_mon[i].r, part_mon[i].ro)) {
      // The particle entered a new cell
      // Must consider the neighbors of the old cell
      ixo = (int) part_mon[i].ro.x / sys.r_cell;
      iyo = (int) part_mon[i].ro.y / sys.r_cell;
      izo = (int) part_mon[i].ro.z / sys.r_cell;

      dix = ix - ixo;
      diy = iy - iyo;
      diz = iz - izo;

      if (dix != 0) {
        jx = mod(ixo - dix, sys.n_cell);

        for (m = -1; m <= 1; m++) {
          for (n = -1; n <= 1; n++) {
            jy = mod(iyo+m, sys.n_cell);
            jz = mod(izo+n, sys.n_cell);

            j = sys.hoc_copy[jx][jy][jz];

            while (j != -1) {
              part_dpd[j].Eo = part_dpd[j].E;
              part_dpd[j].E = calc_energy_dpd(j);
              j = part_dpd[j].ll;
            }
            // To account for overcounting, set HOC to -1
            sys.hoc_copy[jx][jy][jz] = -1;
          }
        }
      }

      if (diy != 0) {
        jy = mod(iyo - diy, sys.n_cell);

        for (l = -1; l <= 1; l++) {
          for (n = -1; n <= 1; n++) {
            jx = mod(ixo+l, sys.n_cell);
            jz = mod(izo+n, sys.n_cell);

            j = sys.hoc_copy[jx][jy][jz];

            while (j != -1) {
              part_dpd[j].Eo = part_dpd[j].E;
              part_dpd[j].E = calc_energy_dpd(j);
              j = part_dpd[j].ll;
            }
            // To account for overcounting, set HOC to -1
            sys.hoc_copy[jx][jy][jz] = -1;
          }
        }
      }

      if (diz != 0) {
        jz = mod(izo - diz, sys.n_cell);

        for (l = -1; l <= 1; l++) {
          for (m = -1; m <= 1; m++) {
            jx = mod(ixo+l, sys.n_cell);
            jy = mod(iyo+m, sys.n_cell);

            j = sys.hoc_copy[jx][jy][jz];

            while (j != -1) {
              part_dpd[j].Eo = part_dpd[j].E;
              part_dpd[j].E = calc_energy_dpd(j);
              j = part_dpd[j].ll;
            }
          }
        }
      }
    }
  } else if (!sys.bond_break) {
    calc_energy_brute();
  }
}

int accept_move(void) {
  double En, Eo;

  if (sys.bond_break) {
    // Reject the movement
    return 0;
  } else {
    Eo = sys.energy;
    En = total_energy();

    if (Eo < En && ran3() > exp(-sys.temp*(En-Eo))) {
      // Reject the movement
      return 0;
    } else {
      // Accept the movement
      sys.energy = En;
      return 1;
    }
  }
}


