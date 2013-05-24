/* mc.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void monte_carlo(void) {
    int i, j, ix, iy, iz, ixo, iyo, izo;
    double Eo, En;

    // Choose a random particle
    i = rand() % sys.n_dpd;

    // Store old energy
    Eo = total_energy();
    // Displace particle i a random amount
    random_move(i);
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

void random_move(int i) {
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

        for (l =- 2; l <= 2; l++) {
            for (m =- 2; m <= 2; m++) {
                for (n =- 2; n <= 2; n++) {
                    // Determine the cell
                    jx = mod(ix+l, sys.n_cell);
                    jy = mod(iy+m, sys.n_cell);
                    jz = mod(iz+n, sys.n_cell);

                    j = sys.hoc[jx][jy][jz];

                    while (j != -1) {
                        part_dpd[j].E = calc_energy_list(j);
                        // Next particle in the linked list
                        j = part_dpd[j].ll;
                    }
                }
            }
        }
    } else {
        // Brute force calculation of energies
        for (j = 0; j < sys.n_dpd; j++) {
            part_dpd[j].E = calc_energy(j);
        }
    }
}


