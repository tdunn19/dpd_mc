#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"


double calc_energy_dpd(int i) {
  int ix, iy, iz, jx, jy, jz, j, l, m, n;
  double E_ss = 0, E_ms = 0;
  Vector dr;
  double Ec,rij;
  // Determine cell number of particle i
  ix = (int) part_dpd[i].r.x / sys.r_cell;
  iy = (int) part_dpd[i].r.y / sys.r_cell;
  iz = (int) part_dpd[i].r.z / sys.r_cell;

  // Solvent-solvent interaction
    printf("E_c(%d)=",i);
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
            Ec = energy_c(dr);
            periodic_bc(&dr);
            rij = vmag(dr);
            if (Ec > 0) {
              printf("%lf(%d,%lf)+",Ec,j,rij);
            }
            E_ss += Ec;
          }
          // Next particle in the chain
          j = part_dpd[j].ll;
        }
      }
    }
  }
  E_ss *= sys.a_ss;
  printf("\n\tE_ss = %lf\n", E_ss);

  // Solvent-monomer interaction
  for (j = 0; j < sys.n_mon; j++) {
    dr = vdist(part_dpd[i].r, part_mon[j].r);
    E_ms += energy_c(dr);
  }
  E_ms *= sys.a_ms;

  return E_ss + E_ms;
}

double calc_energy_dpd_debug(int i, int k) {
  int ix, iy, iz, jx, jy, jz, j, l, m, n;
  double E_ss = 0, E_ms = 0;
  Vector dr;
  double Ec,rij,Eco;
  // k is the particle that changed

  // Determine cell number of particle i
  ix = (int) part_dpd[i].r.x / sys.r_cell;
  iy = (int) part_dpd[i].r.y / sys.r_cell;
  iz = (int) part_dpd[i].r.z / sys.r_cell;

  printf("part_dpd[%d] = %lf->%lf (%lf,%lf,%lf) [%d][%d][%d]\n",
    i, part_dpd[i].Eo, part_dpd[i].E,part_dpd[i].r.x, part_dpd[i].r.y, part_dpd[i].r.z,
    ix,iy,iz);

            if (k > -1) {
              printf("\n\npart_dpd[%d].r=(%lf,%lf,%lf)\t",
                i,part_dpd[i].r.x,part_dpd[i].r.y,part_dpd[i].r.z);
              printf("part_dpd[%d].r=(%lf,%lf,%lf)->(%lf,%lf,%lf)\n",
                k,part_dpd[k].ro.x,part_dpd[k].ro.y,part_dpd[k].ro.z,
                part_dpd[k].r.x,part_dpd[k].r.y,part_dpd[k].r.z);
              dr = vdist(part_dpd[i].r,part_dpd[k].r);
              Ec = energy_c(dr);
              periodic_bc(&dr);
              rij = vmag(dr);
              printf("\trij = %lf (Ec = %lf)",rij,Ec);
              dr = vdist(part_dpd[i].r,part_dpd[k].ro);
              Eco = energy_c(dr);
              periodic_bc(&dr);
              rij = vmag(dr);
              printf("\trijo = %lf (Ec = %lf)",rij,Eco);
              printf("\ndpd[%d].E = %lf - %lf + %lf = %lf = %lf\n",
                i,part_dpd[i].Eo,Eco,Ec,part_dpd[i].Eo-Eco+Ec,part_dpd[i].E);
              }

  // Solvent-solvent interaction
  printf("E_c(%d)=",i);
  for (l = -1; l <= 1; l++) {
    for (m = -1; m <= 1; m++) {
      for (n = -1; n <= 1; n++) {
        // Determine nearest neighbor cell
        jx = mod(ix+l, sys.n_cell);
        jy = mod(iy+m, sys.n_cell);
        jz = mod(iz+n, sys.n_cell);
        // printf("\t[%d][%d][%d]: ", jx,jy,jz);
        // First particle in the chain
        j = sys.hoc[jx][jy][jz];

        while (j != -1) {
          if (i != j) {
            dr = vdist(part_dpd[i].r, part_dpd[j].r);
            // printf("%d")
            Ec = energy_c(dr);
            periodic_bc(&dr);
            rij = vmag(dr);
            if (Ec > 0) {
              printf("%lf(%d,%lf)+",Ec,j,rij);
            }
            E_ss += Ec;
          }
          // Next particle in the chain
          j = part_dpd[j].ll;
        }
      }
    }
  }
  E_ss *= sys.a_ss;
  printf("\n\tE_ss = %lf\t", E_ss);

  // Solvent-monomer interaction
  for (j = 0; j < sys.n_mon; j++) {
    dr = vdist(part_dpd[i].r, part_mon[j].r);
    E_ms += energy_c(dr);
  }
  E_ms *= sys.a_ms;
  printf("\tE_ms = %lf\n", E_ms);

  return E_ss + E_ms;
}

double calc_energy_mon_debug(int i) {
  int ix, iy, iz, jx, jy, jz, j, l, m, n;
  double E_ms = 0, E_mm = 0, E_fene = 0;
  Vector dr;

    ix = (int) part_mon[i].r.x / sys.r_cell;
    iy = (int) part_mon[i].r.y / sys.r_cell;
    iz = (int) part_mon[i].r.z / sys.r_cell;

  printf("part_mon[%d] = %lf->%lf (%lf,%lf,%lf) [%d][%d][%d]\n",
    i, part_mon[i].Eo, part_mon[i].E,part_mon[i].r.x, part_mon[i].r.y, part_mon[i].r.z,
    ix,iy,iz);

  // If not the first monomer in the chain
  if (i != 0) {
    E_fene += energy_fene(i, i-1);
  }
  // If not the last monomer in the chain
  if (i != (sys.n_mon - 1)) {
    E_fene += energy_fene(i, i+1);
  }
  E_fene *= sys.k_fene / 2;
  printf("\tE_fene = %lf\n", E_fene);

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
  printf("\tE_ms = %lf\n", E_ms);

    // Monomer-monomer interaction
    for (j = 0; j < sys.n_mon; j++) {
      if (i != j) {
        dr = vdist(part_mon[i].r, part_mon[j].r);
        E_mm += energy_c(dr);
      }
    }
    E_mm *= sys.a_mm;
  printf("\tE_mm = %lf\n", E_mm);
  }

  return E_ms + E_mm + E_fene;
}

double calc_energy_mon(int i) {
  int ix, iy, iz, jx, jy, jz, j, l, m, n;
  double E_ms = 0, E_mm = 0, E_fene = 0;
  Vector dr;

  // If not the first monomer in the chain
  if (i != 0) {
    E_fene += energy_fene(i, i-1);
  }
  // If not the last monomer in the chain
  if (i != (sys.n_mon - 1)) {
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

void calc_energy_brute(void) {
  int i, j;
  double E_mm, E_ms, E_ss, E_fene;
  Vector dr;

  // Monomer energies
  for (i = 0; i < sys.n_mon; i++) {
    E_mm = 0;
    E_ms = 0;
    E_fene = 0;

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

    for (j = 0; j < sys.n_dpd; j++) {
        dr = vdist(part_mon[i].r, part_dpd[j].r);
        E_ms += energy_c(dr);
    }
    E_ms *= sys.a_ms;

    part_mon[i].E = E_fene + E_mm + E_ms;
  }

  // DPD particle energies
  for (i = 0; i < sys.n_dpd; i++) {
    E_ss = 0;
    E_ms = 0;

    part_dpd[i].Eo = part_dpd[i].E;

    for (j = 0; j < sys.n_dpd; j++) {
      if (j != i) {
        dr = vdist(part_dpd[i].r, part_dpd[j].r);
        E_ss += energy_c(dr);
      }
    }
    E_ss *= sys.a_ss;

    for (j = 0; j < sys.n_mon; j++) {
      dr = vdist(part_dpd[i].r, part_mon[j].r);
      E_ms += energy_c(dr);
    }
    E_ms *= sys.a_ms;

    part_dpd[i].E = E_ss + E_ms;
  }
}

double energy_c(Vector dr) {
  double r_ij, E_ij;
  // Vector dro = dr;

  periodic_bc(&dr);
  // if (dro.x != dr.x || dro.y != dr.y || dro.z != dr.z) {
  //   printf("\n\ndr = %lf,%lf,%lf",dr.x,dr.y,dr.z);
  //   printf("\ndro = %lf,%lf,%lf\n\n",dro.x,dro.y,dro.z);
  // }

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
  periodic_bc(&dr);

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
  double dE_dpd = 0, dE_mon = 0;

  for (i = 0; i < sys.n_dpd; i++) {
    E += part_dpd[i].E;
    if (part_dpd[i].E != part_dpd[i].Eo) {
      printf("\tpart_dpd[%d].E = %lf->%lf (%lf,%lf,%lf)\n",i,part_dpd[i].Eo,part_dpd[i].E,
             part_dpd[i].r.x,part_dpd[i].r.y,part_dpd[i].r.z);
      dE_dpd += part_dpd[i].E - part_dpd[i].Eo;
    }
  }
  printf("\t\tdE_dpd = %lf\n", dE_dpd);

  for (i = 0; i < sys.n_mon; i++) {
    E += part_mon[i].E;
    if (part_mon[i].E != part_mon[i].Eo) {
      printf("\tpart_mon[%d].E = %lf->%lf (%lf,%lf,%lf)\n",i,part_mon[i].Eo,part_mon[i].E, part_mon[i].r.x,part_mon[i].r.y,part_mon[i].r.z);
      dE_mon += part_mon[i].E - part_mon[i].Eo;
    }
  }
  printf("\t\tdE_mon = %lf\n", dE_mon);

  return E;
}

