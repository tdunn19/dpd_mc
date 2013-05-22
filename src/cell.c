/* cell.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dpd.h"

void new_list (void)
{
	int i,ix,iy,iz,n;


	// Initialize head of chain for each cell
	for (ix=0; ix<sys.n_cell; ix++)
	{
		for (iy=0; iy<sys.n_cell; iy++)
		{
			for (iz=0; iz<sys.n_cell; iz++)
			{
				// -1 indicates the end of the chain
				sys.hoc[ix][iy][iz] = -1;
			}
		}
	}

	for (i=0; i<sys.n_dpd; i++)
	{
		// Determine cell number of the particle
		ix = (int) part[i].r.x / sys.r_cell;
		iy = (int) part[i].r.y / sys.r_cell;
		iz = (int) part[i].r.z / sys.r_cell;

		// Link list the head of chain of cell i,j,k
		part[i].ll = sys.hoc[ix][iy][iz];

		// Make particle i the new head of chain
		sys.hoc[ix][iy][iz] = i;
	}
}

double calc_energy_list (int i)
{
	int ix,iy,iz,jx,jy,jz,j,l,m,n;
	double E;

	E = 0;

	// Determine cell number of particle i
	ix = (int) part[i].r.x / sys.r_cell;
	iy = (int) part[i].r.y / sys.r_cell;
	iz = (int) part[i].r.z / sys.r_cell;

	for (l=-1; l<=1; l++)
	{
		for (m=-1; m<=1; m++)
		{
			for (n=-1; n<=1; n++)
			{
				// Determine nearest neighbor cell
				jx = mod(ix+l,sys.n_cell);
				jy = mod(iy+m,sys.n_cell);
				jz = mod(iz+n,sys.n_cell);

				// First particle in the chain
				j = sys.hoc[jx][jy][jz];

				while (j!=-1)
				{
					if (i!=j)
					{
						E += energy_ij(i,j);
					}
					// Next particle in the chain
					j = part[j].ll;
				}
			}
		}
	}
}

