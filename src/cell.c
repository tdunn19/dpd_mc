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
        sys.hoc_copy[ix][iy][iz] = -1;
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
    sys.hoc_copy[ix][iy][iz] = i;
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

void periodic_bc(Vector *dr) {
  if ((*dr).x > sys.length/2) {
    (*dr).x -= sys.length;
  } else if ((*dr).x < -sys.length/2) {
    (*dr).x += sys.length;
  }

  if ((*dr).y > sys.length/2) {
    (*dr).y -= sys.length;
  } else if ((*dr).y < -sys.length/2) {
    (*dr).y += sys.length;
  }

  if ((*dr).z > sys.length/2) {
    (*dr).z -= sys.length;
  } else if ((*dr).z < -sys.length/2) {
    (*dr).z += sys.length;
  }
}
