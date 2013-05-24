/* mc.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void monte_carlo(void) {
  int i, j, ix, iy, iz, ixo, iyo, izo;
  double Eo, En;

  // Store old energy
  Eo = sys.energy;
  sys.bond_break = 0;

  // Move a random particle
  i = rand() % (sys.n_dpd + sys.n_mon);
  if (i >= sys.n_dpd) {
    i = i % sys.n_dpd;
    random_move_mon(i);
  } else {
    random_move_dpd(i);
  }

  if (sys.bond_break) {
    // Reject the movement
    part_dpd[i].r.x = part_dpd[i].ro.x;
    part_dpd[i].r.y = part_dpd[i].ro.y;
    part_dpd[i].r.z = part_dpd[i].ro.z;

    for (j = 0; j < sys.n_dpd; j++) {
      part_dpd[j].E = part_dpd[j].Eo;
    }
  }
  else {
    // Calculate new energy
    En = total_energy();

    if (Eo < En && ran3() > exp(-sys.temp*(En-Eo))) {
      // Reject the movement
      part_dpd[i].r.x = part_dpd[i].ro.x;
      part_dpd[i].r.y = part_dpd[i].ro.y;
      part_dpd[i].r.z = part_dpd[i].ro.z;

      for (j = 0; j < sys.n_dpd; j++) {
        part_dpd[j].E = part_dpd[j].Eo;
      }
    } else {
      // Accept the movement
      sys.energy = En;

      if (sys.calc_list == 1) {
        // Check for new head of chain
        ix = (int) part_dpd[i].r.x / sys.r_cell;
        iy = (int) part_dpd[i].r.y / sys.r_cell;
        iz = (int) part_dpd[i].r.z / sys.r_cell;

        ixo = (int) part_dpd[i].ro.x / sys.r_cell;
        iyo = (int) part_dpd[i].ro.y / sys.r_cell;
        izo = (int) part_dpd[i].ro.z / sys.r_cell;

        if (sys.hoc[ix][iy][iz] != sys.hoc[ixo][iyo][izo]) {
          new_list();
        }
      }
    }
  }
}

void random_move_dpd(int i) {
  int ix, ixo, iy, iyo, iz, izo, j, jx, jy, jz, l, m, n;

  part_dpd[i].ro.x = part_dpd[i].r.x;
  part_dpd[i].ro.y = part_dpd[i].r.y;
  part_dpd[i].ro.z = part_dpd[i].r.z;

  // Displacement between -0.5 and +0.5 in each dimension
  part_dpd[i].r.x += ran3() - 0.5;
  part_dpd[i].r.y += ran3() - 0.5;
  part_dpd[i].r.z += ran3() - 0.5;

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

  // Calculate new energies of particles within two cells of particle i.
  if (sys.calc_list == 1) {
    // Determine cell number of particle i
    ix = (int) part_dpd[i].r.x / sys.r_cell;
    iy = (int) part_dpd[i].r.y / sys.r_cell;
    iz = (int) part_dpd[i].r.z / sys.r_cell;

    for (l = -2; l <= 2; l++) {
      for (m = -2; m <= 2; m++) {
        for (n = -2; n <= 2; n++) {
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

    for (j = 0; j < sys.n_mon; j++) {
      part_mon[j].Eo = part_mon[j].E;
      part_mon[j].E = calc_energy_mon(j);
    }

  } else {
    // Brute force calculation of energies
    for (j = 0; j < sys.n_dpd; j++) {
      part_dpd[j].Eo = part_dpd[j].E;
      part_dpd[j].E = calc_energy_dpd(j);
    }
  }
}

void random_move_mon(int i) {
  int ix, ixo, iy, iyo, iz, izo, j, jx, jy, jz, l, m, n;

  part_mon[i].ro.x = part_mon[i].r.x;
  part_mon[i].ro.y = part_mon[i].r.y;
  part_mon[i].ro.z = part_mon[i].r.z;

  // Displacement between -0.5 and +0.5 in each dimension
  part_mon[i].r.x += ran3() - 0.5;
  part_mon[i].r.y += ran3() - 0.5;
  part_mon[i].r.z += ran3() - 0.5;

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

  if (sys.calc_list == 1) {
    for (j = 0; j < sys.n_mon && !(sys.bond_break); j++) {
      part_mon[j].Eo = part_mon[j].E;
      part_mon[j].E = calc_energy_mon(j);
    }

    if (!sys.bond_break) {
      // Determine new energies for particles within two cells
      ix = (int) part_mon[i].r.x / sys.r_cell;
      iy = (int) part_mon[i].r.y / sys.r_cell;
      iz = (int) part_mon[i].r.z / sys.r_cell;

      for (l = -2; l <= 2; l++) {
        for (m = -2; m <= 2; m++) {
          for (n = -2; n <= 2; n++) {
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
    }
  } else {
    // Brute force calculation of energies
    for (j = 0; j < sys.n_mon; j++) {
      part_mon[j].Eo = part_mon[j].E;
      part_mon[j].E = calc_energy_mon(j);
    }

    for (j = 0; j < sys.n_dpd; j++) {
      part_mon[j].Eo = part_dpd[j].E;
      part_mon[j].E = calc_energy_dpd(j);
    }

  }
}


