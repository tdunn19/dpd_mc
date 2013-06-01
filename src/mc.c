/* mc.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "dpd.h"

void monte_carlo(void) {
  int i, j, ix, iy, iz, ixo, iyo, izo;

  // Choose a random particle
  i = rand() % (sys.n_dpd + sys.n_mon);

  if (i >= sys.n_dpd) {
    // A monomer was chosen
<<<<<<< HEAD

=======
    sys.n_attempt_mon += 1;
>>>>>>> polymer
    i = i % sys.n_dpd;
    //printf("monomer[%d] was chosen\n", i);
    random_move_mon(i);

    if (!accept_move()) {
    //  printf("MOVE REJECTED\n\n");
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
<<<<<<< HEAD
        //printf("MOVE ACCEPTED\n\n");
=======
      sys.n_accept_mon += 1;

>>>>>>> polymer
      for (j = 0; j < sys.n_mon; j++) {
        part_mon[j].Eo = part_mon[j].E;
      }

      for (j = 0; j < sys.n_dpd; j++) {
        part_dpd[j].Eo = part_dpd[j].E;
      }
<<<<<<< HEAD

=======
>>>>>>> polymer
    }
  } else {
    // A DPD particle was chosen
    sys.n_attempt_dpd += 1;
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
<<<<<<< HEAD
=======
      sys.n_accept_dpd += 1;

>>>>>>> polymer
      for (j = 0; j < sys.n_mon; j++) {
        part_mon[j].Eo = part_mon[j].E;
      }

      for (j = 0; j < sys.n_dpd; j++) {
        part_dpd[j].Eo = part_dpd[j].E;
      }
<<<<<<< HEAD
=======

>>>>>>> polymer
      // Generate a new list if the particle changed cells
      if (sys.calc_list && check_cell(part_dpd[i].r, part_dpd[i].ro)) {
        new_list();
      }
    }
  }
}

void random_move_dpd(int i) {
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
  printf("dpd[%d].r=(%lf,%lf,%lf)->(%lf,%lf,%lf)\n",i,
          part_dpd[i].ro.x,part_dpd[i].ro.y,part_dpd[i].ro.z,
          part_dpd[i].r.x,part_dpd[i].r.y,part_dpd[i].r.z);

  if (sys.calc_list) {
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
<<<<<<< HEAD
      // First, make a copy of the HOC array, which we can freely modify
      // memcpy(hoc_copy, sys.hoc, sizeof(hoc_copy));
      // memcpy(hoc_copy, sys.hoc, 2*sizeof(hoc_copy));
      // memcpy(hoc_copy, sys.hoc, sizeof(pow((sys.n_cell+1)*sizeof(int),3)));


=======
>>>>>>> polymer
      ixo = (int) part_dpd[i].ro.x / sys.r_cell;
      iyo = (int) part_dpd[i].ro.y / sys.r_cell;
      izo = (int) part_dpd[i].ro.z / sys.r_cell;

      dix = ix - ixo;
      diy = iy - iyo;
      diz = iz - izo;

      if (dix != 0) {
        jx = mod(ixo - dix, sys.n_cell);

        printf("dix != 0:\n ");
        for (m = -1; m <= 1; m++) {
          for (n = -1; n <= 1; n++) {
            jy = mod(iyo+m, sys.n_cell);
            jz = mod(izo+n, sys.n_cell);
        // printf("[%d][%d][%d], ",jx,jy,jz);

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
        jy = mod(iyo - diy, sys.n_cell);
        printf("diy != 0: \n");

        for (l = -1; l <= 1; l++) {
          for (n = -1; n <= 1; n++) {
            jx = mod(ixo+l, sys.n_cell);
            jz = mod(izo+n, sys.n_cell);
        // printf("[%d][%d][%d], ",jx,jy,jz);

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
        jz = mod(izo - diz, sys.n_cell);
        printf("diz != 0:\n ");

        for (l = -1; l <= 1; l++) {
          for (m = -1; m <= 1; m++) {
            jx = mod(ixo+l, sys.n_cell);
            jy = mod(iyo+m, sys.n_cell);
        // printf("[%d][%d][%d], ",jx,jy,jz);

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
  printf("mon[%d].r=(%lf,%lf,%lf)->(%lf,%lf,%lf)\n",i,
          part_mon[i].ro.x,part_mon[i].ro.y,part_mon[i].ro.z,
          part_mon[i].r.x,part_mon[i].r.y,part_mon[i].r.z);


  // Check for bond breaks
  check_bond(i);

  if (sys.calc_list && !sys.bond_break) {
//      printf("bond did not break\n");

      for (j = 0; j < sys.n_mon; j++) {
      part_mon[j].Eo = part_mon[j].E;
      part_mon[j].E = calc_energy_mon(j);
 //     printf("mon[%d].E = %lf->%lf\n", j, part_mon[j].Eo, part_mon[j].E);
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

      ixo = (int) part_mon[i].ro.x / sys.r_cell;
      iyo = (int) part_mon[i].ro.y / sys.r_cell;
      izo = (int) part_mon[i].ro.z / sys.r_cell;
      //printf("particle entered a new cell:\n");
      //printf("\t[%d][%d][%d] -> [%d][%d][%d]\n",ixo,iyo,izo,ix,iy,iz);

      dix = ix - ixo;
      diy = iy - iyo;
      diz = iz - izo;

      if (dix != 0) {
        jx = mod(ixo - dix, sys.n_cell);

        for (m = -1; m <= 1; m++) {
          for (n = -1; n <= 1; n++) {
            jy = mod(iyo+m, sys.n_cell);
            jz = mod(izo+n, sys.n_cell);
          //printf("now calculating cell [%d][%d][%d]\n",jx,jy,jz);

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
        jy = mod(iyo - diy, sys.n_cell);

        for (l = -1; l <= 1; l++) {
          for (n = -1; n <= 1; n++) {
            jx = mod(ixo+l, sys.n_cell);
            jz = mod(izo+n, sys.n_cell);
          //printf("now calculating cell [%d][%d][%d]\n",jx,jy,jz);

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
        jz = mod(izo - diz, sys.n_cell);

        for (l = -1; l <= 1; l++) {
          for (m = -1; m <= 1; m++) {
            jx = mod(ixo+l, sys.n_cell);
            jy = mod(iyo+m, sys.n_cell);
          //printf("now calculating cell [%d][%d][%d]\n",jx,jy,jz);

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
  } else if (!sys.bond_break) {
    calc_energy_brute();
  }
}

int accept_move(void) {
  double En, Eo;
  double randnum = ran3();

  if (sys.bond_break) {
    // Reject the movement
   printf("\nBOND BROKEN!\n\n");
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


