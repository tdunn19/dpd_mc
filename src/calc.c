#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void sample(void) {
    int i, j;
    double r_ij;
    Vector dr;

    P.now = 0;

    for (i = 0; i < sys.n_dpd; i++) {
        for (j = i+1; j < sys.n_dpd; j++) {
            dr = vdist(part_dpd[i].r, part_dpd[j].r);

            // Employ periodic boundary conditions
            if (dr.x > sys.length/2) {
                dr.x -= sys.length;
            } else if (dr.x < -sys.length/2) {
                dr.x += sys.length;
            }

            if (dr.y > sys.length/2) {
                dr.y -= sys.length;
            } else if (dr.y < -sys.length/2) {
                dr.y += sys.length;
            }

            if (dr.z > sys.length/2) {
                dr.z -= sys.length;
            } else if (dr.z < -sys.length/2) {
                dr.z += sys.length;
            }

            r_ij = vmag(dr);

            // Check for cutoff distance
            if (r_ij < sys.r_c) {
                P.now += sys.a_ij*r_ij*(1-r_ij);
            }
        }
    }

    P.now /= 3*sys.volume;
    P.now += sys.density*sys.temp;

    for (i = 0; i < sys.n_stats; i++) {
        sys.stats[i].sum += sys.stats[i].now;
        sys.stats[i].sumsq += sys.stats[i].now*sys.stats[i].now;
        sys.stats[i].num += 1;
    }
}

double calc_energy(int i) {
    int j;
    double E, r_ij;
    Vector dr;
    E = 0;

    for (j = 0; j < sys.n_dpd; j++) {
        if (j != i) {
            E += energy_ij(i, j);
        }
    }

    E *= sys.a_ij/2;

    return E;
}

double energy_ij(int i, int j) {
    double r_ij, E_ij;
    Vector dr;

    dr = vdist(part_dpd[i].r, part_dpd[j].r);

    // Periodic boundary conditions
    if (dr.x > sys.length/2) {
        dr.x -= sys.length/2;
    } else if (dr.x < -sys.length/2) {
        dr.x += sys.length/2;
    }

    if (dr.y > sys.length/2) {
        dr.y -= sys.length/2;
    } else if (dr.y < -sys.length/2) {
        dr.y += sys.length/2;
    }

    if (dr.z > sys.length/2) {
        dr.z -= sys.length/2;
    } else if (dr.z < -sys.length/2) {
        dr.z += sys.length/2;
    }

    r_ij = vmag(dr);

    // Soft repulsive force
    if (r_ij < sys.r_c) {
        E_ij = (1-r_ij/sys.r_c)*(1-r_ij/sys.r_c);
    } else {
        E_ij = 0;
    }

    return E_ij;
}

double total_energy(void) {
    int i;
    double E;
    E = 0;

    for (i = 0; i < sys.n_dpd; i++) {
        E += part_dpd[i].E;
    }

    return E;
}
