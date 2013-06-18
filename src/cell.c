/* cell.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dpd.h"

void new_list(void) {
  int i, ix, iy, iz, n;

  // Initialize head of chain for each cell.
  // -1 indiciates the end of the chain.
  for (ix = 0; ix < sys.n_cell_1d.x; ix++) {
    for (iy = 0; iy < sys.n_cell_1d.y; iy++) {
      for (iz = 0; iz < sys.n_cell_1d.z; iz++) {
        sys.hoc[ix][iy][iz] = -1;
        sys.hoc_copy[ix][iy][iz] = -1;
      }
    }
  }

  for (i = 0; i < sys.n_dpd; i++) {
    // Determine cell number of the particle
    ix = (int) (part_dpd[i].r.x / sys.r_cell.x);
    iy = (int) (part_dpd[i].r.y / sys.r_cell.y);
    iz = (int) (part_dpd[i].r.z / sys.r_cell.z);

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

  ix = (int) (r.x / sys.r_cell.x);
  iy = (int) (r.y / sys.r_cell.y);
  iz = (int) (r.z / sys.r_cell.z);

  ixo = (int) (ro.x / sys.r_cell.x);
  iyo = (int) (ro.y / sys.r_cell.y);
  izo = (int) (ro.z / sys.r_cell.z);

  if (sys.hoc[ix][iy][iz] != sys.hoc[ixo][iyo][izo]) {
    // The cell has changed
    return 1;
  } else {
    return 0;
  }
}

void periodic_bc_r(Vector *r) {
  if ((*r).x > sys.length.x) {
    (*r).x -= sys.length.x;
  } else if ((*r).x < 0) {
    (*r).x += sys.length.x;
  }

  if ((*r).y > sys.length.y) {
    (*r).y -= sys.length.y;
  } else if ((*r).y < 0) {
    (*r).y += sys.length.y;
  }

  if ((*r).z > sys.length.z) {
    (*r).z -= sys.length.z;
  } else if ((*r).z < 0) {
    (*r).z += sys.length.z;
  }
}

void periodic_bc_dr(Vector *dr) {
  if ((*dr).x > sys.length.x/2) {
    (*dr).x -= sys.length.x;
  } else if ((*dr).x < -sys.length.x/2) {
    (*dr).x += sys.length.x;
  }

  if ((*dr).y > sys.length.y/2) {
    (*dr).y -= sys.length.y;
  } else if ((*dr).y < -sys.length.y/2) {
    (*dr).y += sys.length.y;
  }

  if ((*dr).z > sys.length.z/2) {
    (*dr).z -= sys.length.z;
  } else if ((*dr).z < -sys.length.z/2) {
    (*dr).z += sys.length.z;
  }
}
