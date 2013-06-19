/* mc.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "dpd.h"

void monte_carlo(void) {
  int i, j, ix, iy, iz, ixo, iyo, izo;

  if (ran3() < sys.mc_ratio) {
    // A monomer was chosen
<<<<<<< HEAD

=======
    sys.n_attempt_mon += 1;

    // Randomly move a monomer
    i = rand() % sys.n_mon;
    move_monomer(i);

    if (accept_move()) {
      // The move was accepeted
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

>>>>>>> polymer
      for (j = 0; j < sys.n_mon; j++) {
        part_mon[j].E = part_mon[j].Eo;
      }

      for (j = 0; j < sys.n_solvent; j++) {
        part_dpd[j].E = part_dpd[j].Eo;
      }
<<<<<<< HEAD

=======
>>>>>>> polymer
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

>>>>>>> polymer
      for (j = 0; j < sys.n_mon; j++) {
        part_mon[j].Eo = part_mon[j].E;
      }

      for (j = 0; j < sys.n_solvent; j++) {
        part_dpd[j].Eo = part_dpd[j].E;
      }
<<<<<<< HEAD
=======

>>>>>>> polymer
      // Generate a new list if the particle changed cells
      if (sys.calc_list && check_cell(part_dpd[i].r, part_dpd[i].ro)) {
        new_list();
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
}

void move_solvent(int i) {
  int ix, ixo, iy, iyo, iz, izo, j, jx, jy, jz, dix, diy, diz, l, m, n;
  double rij;
  // int ***hoc_copy = malloc(2*sizeof(sys.hoc));
  Vector dr;
  // int ***hoc_copy;
  // int ***hoc_copy = malloc(sizeof(pow((sys.n_cell+1)*sizeof(int),3)));

  part_dpd[i].ro.x = part_dpd[i].r.x;
  part_dpd[i].ro.y = part_dpd[i].r.y;
  part_dpd[i].ro.z = part_dpd[i].r.z;

  // Displacement between -dr_max and +dr_max in each dimension
  part_dpd[i].r.x += sys.dr_max_dpd * (2*ran3() - 1);
  part_dpd[i].r.y += sys.dr_max_dpd * (2*ran3() - 1);
  part_dpd[i].r.z += sys.dr_max_dpd * (2*ran3() - 1);

  // Periodic boundary conditions
  periodic_bc_r(&part_dpd[i].r);
  // Check for wall overlap
  check_wall(part_dpd[i].r);
  // Check for pore overlap
  check_pore(part_dpd[i].r);

  if (sys.calc_list && !sys.wall_overlap && !sys.pore_overlap) {
    for (j = 0; j < sys.n_mon; j++) {
      part_mon[j].Eo = part_mon[j].E;
      part_mon[j].E = calc_energy_mon(j);
      if (part_mon[j].E != part_mon[j].Eo) {
        dr = vdist(part_dpd[i].r,part_mon[j].r);
        periodic_bc(&dr);
        rij = vmag(dr);

        printf("part_mon[%d] changed, dr = (%lf,%lf,%lf) rij = %lf\n",
          j,dr.x,dr.y,dr.z,rij);
        part_mon[j].E = calc_energy_mon_debug(j);
      }
    }

    // Determine cell number of particle i
    ix = (int) (part_dpd[i].r.x / sys.r_cell.x);
    iy = (int) (part_dpd[i].r.y / sys.r_cell.y);
    iz = (int) (part_dpd[i].r.z / sys.r_cell.z);

    for (l = -1; l <= 1; l++) {
      for (m = -1; m <= 1; m++) {
        for (n = -1; n <= 1; n++) {
          // Determine the cell
          jx = mod(ix+l, sys.n_cell_1d.x);
          jy = mod(iy+m, sys.n_cell_1d.y);
          jz = mod(iz+n, sys.n_cell_1d.z);

          j = sys.hoc[jx][jy][jz];

          while (j != -1) {
            part_dpd[j].Eo = part_dpd[j].E;
            part_dpd[j].E = calc_energy_dpd(j);

            if (part_dpd[j].E != part_dpd[j].Eo) {
        dr = vdist(part_dpd[i].r,part_dpd[j].r);
        periodic_bc(&dr);
        rij = vmag(dr);
        if (rij > 1){
          dr = vdist(part_dpd[i].ro,part_dpd[j].r);
        periodic_bc(&dr);
        rij = vmag(dr);
        }
        printf("part_dpd[%d] changed, dr = (%lf,%lf,%lf) rij = %lf\n",
          j,dr.x,dr.y,dr.z,rij);
              part_dpd[j].E = calc_energy_dpd_debug(j,-2);
            }

            // Next particle in the linked list
            j = part_dpd[j].ll;
          }
        }
      }
    }

    if (check_cell(part_dpd[i].r, part_dpd[i].ro)) {
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

        printf("dix != 0:\n ");
        for (m = -1; m <= 1; m++) {
          for (n = -1; n <= 1; n++) {
            jy = mod(iyo+m, sys.n_cell_1d.y);
            jz = mod(izo+n, sys.n_cell_1d.z);

            j = sys.hoc_copy[jx][jy][jz];

            while (j >= 0) {
              part_dpd[j].Eo = part_dpd[j].E;
              part_dpd[j].E = calc_energy_dpd(j);
              if (part_dpd[j].E != part_dpd[j].Eo) {
        dr = vdist(part_dpd[i].r,part_dpd[j].r);
        periodic_bc(&dr);
        rij = vmag(dr);
        if (rij > 1){
          printf("\trij = %lf > 1, recalculating:\n", rij);
          dr = vdist(part_dpd[i].ro,part_dpd[j].r);
        periodic_bc(&dr);
        rij = vmag(dr);
        }
        printf("part_dpd[%d] changed, dr = (%lf,%lf,%lf) rij = %lf\n",
          j,dr.x,dr.y,dr.z,rij);
               part_dpd[j].E = calc_energy_dpd_debug(j,i);
             }
              j = part_dpd[j].ll;
            }
            // To account for overcounting, set HOC to -1
            sys.hoc_copy[jx][jy][jz] = -1;
          }
        }
      // printf("\n");
      }

      if (diy != 0) {
        jy = mod(iyo - diy, sys.n_cell_1d.y);

        for (l = -1; l <= 1; l++) {
          for (n = -1; n <= 1; n++) {
            jx = mod(ixo+l, sys.n_cell_1d.x);
            jz = mod(izo+n, sys.n_cell_1d.z);

            j = sys.hoc_copy[jx][jy][jz];

            while (j >= 0) {
              part_dpd[j].Eo = part_dpd[j].E;
              part_dpd[j].E = calc_energy_dpd(j);
              if (part_dpd[j].E != part_dpd[j].Eo) {
        dr = vdist(part_dpd[i].r,part_dpd[j].r);
        periodic_bc(&dr);
        rij = vmag(dr);
        if (rij > 1){
          dr = vdist(part_dpd[i].ro,part_dpd[j].r);
        periodic_bc(&dr);
        rij = vmag(dr);
        }
        printf("part_dpd[%d] changed, dr = (%lf,%lf,%lf) rij = %lf\n",
          j,dr.x,dr.y,dr.z,rij);
               part_dpd[j].E = calc_energy_dpd_debug(j,i);}
              j = part_dpd[j].ll;
            }
<<<<<<< HEAD
            // To account for overcounting, set HOC to -1 to indicate
            // this cell
=======
            // To account for overcounting, set HOC to -1
>>>>>>> polymer
            sys.hoc_copy[jx][jy][jz] = -1;
          }
        }
      // printf("\n");
      }

      if (diz != 0) {
        jz = mod(izo - diz, sys.n_cell_1d.z);

        for (l = -1; l <= 1; l++) {
          for (m = -1; m <= 1; m++) {
            jx = mod(ixo+l, sys.n_cell_1d.x);
            jy = mod(iyo+m, sys.n_cell_1d.y);

            j = sys.hoc_copy[jx][jy][jz];

            while (j >= 0) {
              part_dpd[j].Eo = part_dpd[j].E;
              part_dpd[j].E = calc_energy_dpd(j);
              if (part_dpd[j].E != part_dpd[j].Eo) {
        dr = vdist(part_dpd[i].r,part_dpd[j].r);
        periodic_bc(&dr);
        rij = vmag(dr);
        if (rij > 1){
          dr = vdist(part_dpd[i].ro,part_dpd[j].r);
        periodic_bc(&dr);
        rij = vmag(dr);
        }
        printf("part_dpd[%d] changed, dr = (%lf,%lf,%lf) rij = %lf\n",
          j,dr.x,dr.y,dr.z,rij);
               part_dpd[j].E = calc_energy_dpd_debug(j,i);}
              j = part_dpd[j].ll;
            }
          }
        }
      // printf("\n");
      }
    }
  } else if (!sys.wall_overlap && !sys.pore_overlap) {
    calc_energy_brute();
  }
}

void move_monomer(int i) {
  int ix, ixo, iy, iyo, iz, izo, j, jx, jy, jz, dix, diy, diz, l, m, n;

  part_mon[i].ro.x = part_mon[i].r.x;
  part_mon[i].ro.y = part_mon[i].r.y;
  part_mon[i].ro.z = part_mon[i].r.z;

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
  // Check for pore overlap
  check_pore(part_mon[i].r);

  if (sys.calc_list && !sys.bond_break && !sys.wall_overlap && !sys.pore_overlap) {
    for (j = 0; j < sys.n_mon; j++) {
      part_mon[j].Eo = part_mon[j].E;
      part_mon[j].E = calc_energy_mon(j);
 //     printf("mon[%d].E = %lf->%lf\n", j, part_mon[j].Eo, part_mon[j].E);
    }

    // Determine new energies for nearest neighbor cells
    ix = (int) (part_mon[i].r.x / sys.r_cell.x);
    iy = (int) (part_mon[i].r.y / sys.r_cell.y);
    iz = (int) (part_mon[i].r.z / sys.r_cell.z);

    for (l = -1; l <= 1; l++) {
      for (m = -1; m <= 1; m++) {
        for (n = -1; n <= 1; n++) {
          // Determine the cell
          jx = mod(ix+l, sys.n_cell_1d.x);
          jy = mod(iy+m, sys.n_cell_1d.y);
          jz = mod(iz+n, sys.n_cell_1d.z);

  //        printf("now calculating cell [%d][%d][%d]\n",jx,jy,jz);

          j = sys.hoc[jx][jy][jz];

          while (j != -1) {
            part_dpd[j].Eo = part_dpd[j].E;
            part_dpd[j].E = calc_energy_dpd(j);
            //printf("\tdpd[%d].E = %lf->%lf",j,part_dpd[j].Eo,part_dpd[j].E);
            // Next particle in the linked list
            j = part_dpd[j].ll;
          }
          //printf("\n");
        }
      }
    }

    if (check_cell(part_mon[i].r, part_mon[i].ro)) {
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
              part_dpd[j].Eo = part_dpd[j].E;
              part_dpd[j].E = calc_energy_dpd(j);
            //printf("\tdpd[%d].E = %lf->%lf",j,part_dpd[j].Eo,part_dpd[j].E);
              j = part_dpd[j].ll;
            }
            //printf("\n");
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
              part_dpd[j].Eo = part_dpd[j].E;
              part_dpd[j].E = calc_energy_dpd(j);
           // printf("\tdpd[%d].E = %lf->%lf",j,part_dpd[j].Eo,part_dpd[j].E);
              j = part_dpd[j].ll;
            }
            //printf("\n");
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
              part_dpd[j].Eo = part_dpd[j].E;
              part_dpd[j].E = calc_energy_dpd(j);
           // printf("\tdpd[%d].E = %lf->%lf",j,part_dpd[j].Eo,part_dpd[j].E);
              j = part_dpd[j].ll;
            }
            //printf("\n");
          }
        }
      }
    }
  } else if (!sys.bond_break && !sys.wall_overlap && !sys.pore_overlap) {
    calc_energy_brute();
  }
}

int accept_move(void) {
  double En, Eo;
  double randnum = ran3();

  if (sys.bond_break || sys.wall_overlap || sys.pore_overlap) {
    // Reject the movement
    sys.bond_break = 0;
    sys.wall_overlap = 0;
    sys.pore_overlap = 0;
    return 0;
  } else {
    printf("\nsys.Eo = %lf\n", sys.energy);
      Eo = sys.energy;
    En = total_energy();
    printf("\nsys.E = %lf\n", En);


    if (Eo < En && randnum > exp(-sys.temp*(En-Eo))) {
      // Reject the movement
      printf("REJECTED (%lf < %lf && %lf > exp(-%lf*(%lf-%lf) = %lf\n\n",
        Eo, En, randnum, sys.temp, En, Eo, exp(-sys.temp*(En-Eo)));
      return 0;
    } else {
      // Accept the movement
      printf("ACCEPTED (%lf > %lf) || %lf < exp(-%lf*(%lf-%lf) = %lf\n\n",
      Eo, En, randnum,sys.temp,En,Eo,exp(-sys.temp*(En-Eo)));
      sys.energy = En;
      return 1;
    }
  }
}


