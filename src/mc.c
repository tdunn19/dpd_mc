/* mc.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void monte_carlo(void) {
  int i, j;

  if (ran3() < sys.mc_ratio) {
    // A monomer was chosen
    sys.n_attempt_mon += 1;

    // Randomly move a monomer
    i = rand() % sys.n_mon;
    move_monomer(i);

    if (accept_move()) {
      // The move was accepted
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
    } else {
      // The move was rejected
      part_mon[i].r.x = part_mon[i].ro.x;
      part_mon[i].r.y = part_mon[i].ro.y;
      part_mon[i].r.z = part_mon[i].ro.z;

      for (j = 0; j < sys.n_mon; j++) {
        part_mon[j].E = part_mon[j].Eo;
      }

      for (j = 0; j < sys.n_solvent; j++) {
        part_dpd[j].E = part_dpd[j].Eo;
      }
    }
  } else {
    // A solvent particle was chosen
    sys.n_attempt_solvent += 1;

    // Randomly move a solvent particle
    i = rand() % sys.n_solvent;
    move_solvent(i);

    if (accept_move()) {
      // The move was accepted
      sys.n_accept_solvent += 1;

      part_dpd[i].ro.x = part_dpd[i].r.x;
      part_dpd[i].ro.y = part_dpd[i].r.y;
      part_dpd[i].ro.z = part_dpd[i].r.z;

      for (j = 0; j < sys.n_mon; j++) {
        part_mon[j].Eo = part_mon[j].E;
      }

      for (j = 0; j < sys.n_solvent; j++) {
        part_dpd[j].Eo = part_dpd[j].E;
      }
    } else {
      // The move was rejected
      part_dpd[i].r.x = part_dpd[i].ro.x;
      part_dpd[i].r.y = part_dpd[i].ro.y;
      part_dpd[i].r.z = part_dpd[i].ro.z;

      for (j = 0; j < sys.n_dpd; j++) {
        part_dpd[j].E = part_dpd[j].Eo;
      }

      for (j = 0; j < sys.n_mon; j++) {
        part_mon[j].E = part_mon[j].Eo;
      }
    }
  }

  // Generate a new list if the particle changed cells
  if (sys.calc_list && sys.new_cell) {
    new_list();
  }
}

void move_solvent(int i) {
  int ix, ixo, iy, iyo, iz, izo, j, jx, jy, jz, dix, diy, diz, l, m, n;
  double dE;
  Vector dr, dro;

  // Displacement between -dr_max and +dr_max in each dimension
  part_dpd[i].r.x += sys.dr_max_dpd * (2*ran3() - 1);
  part_dpd[i].r.y += sys.dr_max_dpd * (2*ran3() - 1);
  part_dpd[i].r.z += sys.dr_max_dpd * (2*ran3() - 1);

  // Periodic boundary conditions
  periodic_bc_r(&part_dpd[i].r);
  // Check for wall overlap
  check_wall(part_dpd[i].r);

  if (sys.calc_list && !sys.wall_overlap) {
    // Monomer-solvent interaction
    for (j = 0; j < sys.n_mon; j++) {
      dr = vdist(part_dpd[i].r, part_mon[j].r);
      dro = vdist(part_dpd[i].ro, part_mon[j].r);
      dE = sys.a_ms[part_dpd[i].side] * (energy_c(dr) - energy_c(dro));

      part_dpd[i].E += dE;
      part_mon[j].E += dE;
      sys.energy += dE;
    }

    // Determine cell number of particle i
    ix = (int) (part_dpd[i].r.x / sys.r_cell.x);
    iy = (int) (part_dpd[i].r.y / sys.r_cell.y);
    iz = (int) (part_dpd[i].r.z / sys.r_cell.z);

    for (l = -1; l <= 1; l++) {
      for (m = -1; m <= 1; m++) {
        for (n = -1; n <= 1; n++) {
          // Determine the neighbor cell
          jx = mod(ix+l, sys.n_cell_1d.x);
          jy = mod(iy+m, sys.n_cell_1d.y);
          jz = mod(iz+n, sys.n_cell_1d.z);

         j = sys.hoc[jx][jy][jz];

          while (j != -1) {
            if (i != j) {
              dr = vdist(part_dpd[i].r, part_dpd[j].r);
              dro = vdist(part_dpd[i].ro, part_dpd[j].r);
              dE = energy_c(dr) - energy_c(dro);

              if (j < sys.n_solvent) {
                // Solvent-solvent interaction
                dE *= sys.a_ss;
              } else {
                // Solvent-wall interaction
                dE *= sys.a_sw;
              }

              part_dpd[i].E += dE;
              part_dpd[j].E += dE;
              sys.energy += dE;
            }
            // Next particle in the linked list
            j = part_dpd[j].ll;
          }
        }
      }
    }

    check_cell(part_dpd[i].r, part_dpd[i].ro);

    if (sys.new_cell) {
      // The particle entered a new cell
      // Must consider the neighbors of the old cell
      ixo = (int) (part_dpd[i].ro.x / sys.r_cell.x);
      iyo = (int) (part_dpd[i].ro.y / sys.r_cell.y);
      izo = (int) (part_dpd[i].ro.z / sys.r_cell.z);

      dix = ix - ixo;
      diy = iy - iyo;
      diz = iz - izo;

      if (dix != 0) {
        jx = mod(ixo - dix, sys.n_cell_1d.x);

        for (m = -1; m <= 1; m++) {
          for (n = -1; n <= 1; n++) {
            jy = mod(iyo+m, sys.n_cell_1d.y);
            jz = mod(izo+n, sys.n_cell_1d.z);

            j = sys.hoc_copy[jx][jy][jz];

            while (j != -1) {
              dr = vdist(part_dpd[i].r, part_dpd[j].r);
              dro = vdist(part_dpd[i].ro, part_dpd[j].r);
              dE = energy_c(dr) - energy_c(dro);

              if (j < sys.n_solvent) {
                // Solvent-solvent interaction
                dE *= sys.a_ss;
              } else {
                // Solvent-wall interaction
                dE *= sys.a_sw;
              }

              part_dpd[i].E += dE;
              part_dpd[j].E += dE;
              sys.energy += dE;

              j = part_dpd[j].ll;
            }
            // To account for overcounting, set HOC to -1
            sys.hoc_copy[jx][jy][jz] = -1;
          }
        }
      }

      if (diy != 0) {
        jy = mod(iyo - diy, sys.n_cell_1d.y);

        for (l = -1; l <= 1; l++) {
          for (n = -1; n <= 1; n++) {
            jx = mod(ixo+l, sys.n_cell_1d.x);
            jz = mod(izo+n, sys.n_cell_1d.z);

            j = sys.hoc_copy[jx][jy][jz];

            while (j != -1) {
              dr = vdist(part_dpd[i].r, part_dpd[j].r);
              dro = vdist(part_dpd[i].ro, part_dpd[j].r);
              dE = energy_c(dr) - energy_c(dro);

              if (j < sys.n_solvent) {
                // Solvent-solvent interaction
                dE *= sys.a_ss;
              } else {
                // Solvent-wall interaction
                dE *= sys.a_sw;
              }

              part_dpd[i].E += dE;
              part_dpd[j].E += dE;
              sys.energy += dE;

              j = part_dpd[j].ll;
            }
            // To account for overcounting, set HOC to -1
            sys.hoc_copy[jx][jy][jz] = -1;
          }
        }
      }

      if (diz != 0) {
        jz = mod(izo - diz, sys.n_cell_1d.z);

        for (l = -1; l <= 1; l++) {
          for (m = -1; m <= 1; m++) {
            jx = mod(ixo+l, sys.n_cell_1d.x);
            jy = mod(iyo+m, sys.n_cell_1d.y);

            j = sys.hoc_copy[jx][jy][jz];

            while (j != -1) {
              dr = vdist(part_dpd[i].r, part_dpd[j].r);
              dro = vdist(part_dpd[i].ro, part_dpd[j].r);
              dE = energy_c(dr) - energy_c(dro);

             if (j < sys.n_solvent) {
                // Solvent-solvent interaction
                dE *= sys.a_ss;
              } else {
                // Solvent-wall interaction
                dE *= sys.a_sw;
              }

              part_dpd[i].E += dE;
              part_dpd[j].E += dE;
              sys.energy += dE;

              j = part_dpd[j].ll;
            }
          }
        }
      }
    }
  } else if (!sys.wall_overlap) {
    // Solvent-monomer interaction
    for (j = 0; j < sys.n_mon; j++) {
      dr = vdist(part_dpd[i].r, part_mon[j].r);
      dro = vdist(part_dpd[i].ro, part_mon[j].r);
      dE = sys.a_ms[part_dpd[j].side] * (energy_c(dr) - energy_c(dro));

      part_dpd[i].E += dE;
      part_mon[j].E += dE;
      sys.energy += dE;
    }

    // Solvent-solvent interaction
    for (j = 0; j < sys.n_solvent; j++) {
      if (i != j) {
        dr = vdist(part_dpd[i].r, part_dpd[j].r);
        dro = vdist(part_dpd[i].ro, part_dpd[j].r);
        dE = sys.a_ss * (energy_c(dr) - energy_c(dro));

        part_dpd[i].E += dE;
        part_dpd[j].E += dE;
        sys.energy += dE;
      }
    }

    // Solvent-wall interaction
    for (j = sys.n_solvent; j < sys.n_dpd; j++) {
      dr = vdist(part_dpd[i].r, part_dpd[j].r);
      dro = vdist(part_dpd[i].ro, part_dpd[j].r);
      dE = sys.a_sw * (energy_c(dr) - energy_c(dro));

      part_dpd[i].E += dE;
      part_dpd[j].E += dE;
      sys.energy += dE;
    }
  }
}

void move_monomer(int i) {
  int ix, ixo, iy, iyo, iz, izo, j, jx, jy, jz, dix, diy, diz, l, m, n;
  double dE;
  Vector dr, dro;

  // Displacement between -dr_max and +dr_max in each dimension
  part_mon[i].r.x += sys.dr_max_mon * (2*ran3() - 1);
  part_mon[i].r.y += sys.dr_max_mon * (2*ran3() - 1);
  part_mon[i].r.z += sys.dr_max_mon * (2*ran3() - 1);

  // Periodic boundary conditions
  periodic_bc_r(&part_mon[i].r);
  // Check for bond breaks
  check_bond(i);
  // Check for wall overlap
  check_wall(part_mon[i].r);
  // Check for new window
  check_window();

  if (sys.calc_list && !sys.bond_break && !sys.wall_overlap && !sys.new_window) {
    // Monomer-monomer interaction
    for (j = 0; j < sys.n_mon; j++) {
      if (i != j) {
        dr = vdist(part_mon[i].r, part_mon[j].r);
        dro = vdist(part_mon[i].ro, part_mon[j].r);
        dE = sys.a_mm * (energy_c(dr) - energy_c(dro));

        // If the monomers are bonded
        if (j == i-1) {
          dE += energy_fene(dr) - energy_fene(dro);
        } else if (j == i+1) {
          dE += energy_fene(dr) - energy_fene(dro);
        }

        part_mon[i].E += dE;
        part_mon[j].E += dE;
        sys.energy += dE;
      }
    }

    // Determine cell number of monomer i
    ix = (int) (part_mon[i].r.x / sys.r_cell.x);
    iy = (int) (part_mon[i].r.y / sys.r_cell.y);
    iz = (int) (part_mon[i].r.z / sys.r_cell.z);

    for (l = -1; l <= 1; l++) {
      for (m = -1; m <= 1; m++) {
        for (n = -1; n <= 1; n++) {
          // Determine the neighbor cell
          jx = mod(ix+l, sys.n_cell_1d.x);
          jy = mod(iy+m, sys.n_cell_1d.y);
          jz = mod(iz+n, sys.n_cell_1d.z);

          j = sys.hoc[jx][jy][jz];

          while (j != -1) {
            dr = vdist(part_mon[i].r, part_dpd[j].r);
            dro = vdist(part_mon[i].ro, part_dpd[j].r);
            dE = energy_c(dr) - energy_c(dro);

            if (j < sys.n_solvent) {
              // Monomer-solvent interaction
              dE *= sys.a_ms[part_dpd[j].side];
            } else {
              // Monomer-wall interaction
              dE *= sys.a_sw;
            }

            part_mon[i].E += dE;
            part_dpd[j].E += dE;
            sys.energy += dE;

            // Next particle in the linked list
            j = part_dpd[j].ll;
          }
        }
      }
    }

    check_cell(part_mon[i].r, part_mon[i].ro);
    if (sys.new_cell) {
      // The particle entered a new cell
      // Must consider the neighbors of the old cell
      ixo = (int) (part_mon[i].ro.x / sys.r_cell.x);
      iyo = (int) (part_mon[i].ro.y / sys.r_cell.y);
      izo = (int) (part_mon[i].ro.z / sys.r_cell.z);

      dix = ix - ixo;
      diy = iy - iyo;
      diz = iz - izo;

      if (dix != 0) {
        jx = mod(ixo - dix, sys.n_cell_1d.x);

        for (m = -1; m <= 1; m++) {
          for (n = -1; n <= 1; n++) {
            jy = mod(iyo+m, sys.n_cell_1d.y);
            jz = mod(izo+n, sys.n_cell_1d.z);

            j = sys.hoc_copy[jx][jy][jz];

            while (j != -1) {
              dr = vdist(part_mon[i].r, part_dpd[j].r);
              dro = vdist(part_mon[i].ro, part_dpd[j].r);
              dE = energy_c(dr) - energy_c(dro);

              if (j < sys.n_solvent) {
                // Monomer-solvent interaction
                dE *= sys.a_ms[part_dpd[j].side];
              } else {
                // Monomer-wall interaction
                dE *= sys.a_sw;
              }

              part_mon[i].E += dE;
              part_dpd[j].E += dE;
              sys.energy += dE;

              j = part_dpd[j].ll;
            }
            // To account for overcounting, set HOC to -1
            sys.hoc_copy[jx][jy][jz] = -1;
          }
        }
      }

      if (diy != 0) {
        jy = mod(iyo - diy, sys.n_cell_1d.y);

        for (l = -1; l <= 1; l++) {
          for (n = -1; n <= 1; n++) {
            jx = mod(ixo+l, sys.n_cell_1d.x);
            jz = mod(izo+n, sys.n_cell_1d.z);

            j = sys.hoc_copy[jx][jy][jz];

            while (j != -1) {
              dr = vdist(part_mon[i].r, part_dpd[j].r);
              dro = vdist(part_mon[i].ro, part_dpd[j].r);
              dE = energy_c(dr) - energy_c(dro);

              if (j < sys.n_solvent) {
                // Monomer-solvent interaction
                dE *= sys.a_ms[part_dpd[j].side];
              } else {
                // Monomer-wall interaction
                dE *= sys.a_sw;
              }

              part_mon[i].E += dE;
              part_dpd[j].E += dE;
              sys.energy += dE;

              j = part_dpd[j].ll;
            }
            // To account for overcounting, set HOC to -1
            sys.hoc_copy[jx][jy][jz] = -1;
          }
        }
      }

      if (diz != 0) {
        jz = mod(izo - diz, sys.n_cell_1d.z);

        for (l = -1; l <= 1; l++) {
          for (m = -1; m <= 1; m++) {
            jx = mod(ixo+l, sys.n_cell_1d.x);
            jy = mod(iyo+m, sys.n_cell_1d.y);

            j = sys.hoc_copy[jx][jy][jz];

            while (j != -1) {
              dr = vdist(part_mon[i].r, part_dpd[j].r);
              dro = vdist(part_mon[i].ro, part_dpd[j].r);
              dE = energy_c(dr) - energy_c(dro);

              if (j < sys.n_solvent) {
                // Monomer-solvent interaction
                dE *= sys.a_ms[part_dpd[j].side];
              } else {
                // Monomer-wall interaction
                dE *= sys.a_sw;
              }

              part_mon[i].E += dE;
              part_dpd[j].E += dE;
              sys.energy += dE;

              j = part_dpd[j].ll;
            }
          }
        }
      }
    }
  } else if (!sys.bond_break && !sys.wall_overlap && !sys.new_window) {
    // Monomer-monomer interaction
    for (j = 0; j < sys.n_mon; j++) {
      if (i != j) {
        dr = vdist(part_mon[i].r, part_mon[j].r);
        dro = vdist(part_mon[i].ro, part_mon[j].r);
        dE = sys.a_mm * (energy_c(dr) - energy_c(dro));

        // If the monomers are bonded
        if (j == i-1) {
          dE += energy_fene(dr) - energy_fene(dro);
        } else if (j == i+1) {
          dE += energy_fene(dr) - energy_fene(dro);
        }

        part_mon[i].E += dE;
        part_mon[j].E += dE;
        sys.energy += dE;
      }
    }

    // Monomer-solvent interaction
    for (j = 0; j < sys.n_solvent; j++) {
      dr = vdist(part_mon[i].r, part_dpd[j].r);
      dro = vdist(part_mon[i].ro, part_dpd[j].r);
      dE = sys.a_ms[part_dpd[j].side] * (energy_c(dr) - energy_c(dro));

      part_mon[i].E += dE;
      part_dpd[j].E += dE;
      sys.energy += dE;
    }

    // Monomer-wall interaction
    for (j = sys.n_solvent; j < sys.n_dpd; j++) {
      dr = vdist(part_mon[i].r, part_dpd[j].r);
      dro = vdist(part_mon[i].ro, part_dpd[j].r);
      dE = sys.a_sw * (energy_c(dr) - energy_c(dro));

      part_mon[i].E += dE;
      part_dpd[j].E += dE;
      sys.energy += dE;
    }
  }
}

int accept_move(void) {
  double dE;

  if (sys.bond_break || sys.wall_overlap || sys.new_window) {
    // Reject the movement and reset system flags
    sys.bond_break = 0;
    sys.wall_overlap = 0;
    sys.new_window = 0;

    sys.energy = sys.energy_old;

    return 0;
  } else {
    dE = sys.energy - sys.energy_old;

    if (dE > 0 && ran3() > exp(-dE)) {
      // Reject the movement
      sys.energy = sys.energy_old;
      return 0;
    } else {
      // Accept the movement
      sys.energy_old = sys.energy;
      return 1;
    }
  }
}
