/* mc.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void monte_carlo (void)
{
	int i,j,ix,iy,iz,ixo,iyo,izo;
	double Eo, En;

	// Choose a random particle from 0 to n_dpd
	i = rand() % sys.n_dpd;

	// Store old energy
	Eo = total_energy();
	// Displace particle i a random amount
	random_move(i);
	// Calculate new energy
	En = total_energy();

	if (Eo < En && ran3() > exp(-sys.temp*(En-Eo)))
	{
		// Reject the movement
		part[i].r.x = part[i].ro.x;
		part[i].r.y = part[i].ro.y;
		part[i].r.z = part[i].ro.z;

		for (j=0; j<sys.n_dpd; j++)
		{
			part[j].E = part[j].Eo;
		}
	}
	else
	{
		// Accept the movement
		if (sys.calc_list == 1)
		{
			// Check for new head of chain
			ix = (int) part[i].r.x / sys.r_cell;
			iy = (int) part[i].r.y / sys.r_cell;
			iz = (int) part[i].r.z / sys.r_cell;

			ixo = (int) part[i].ro.x / sys.r_cell;
			iyo = (int) part[i].ro.y / sys.r_cell;
			izo = (int) part[i].ro.z / sys.r_cell;

			if (sys.hoc[ix][iy][iz] != sys.hoc[ixo][iyo][izo])
			{
				new_list();
			}
		}
	}
}

void random_move (int i)
{
	int ix,ixo,iy,iyo,iz,izo,j,jx,jy,jz,l,m,n;

	part[i].ro.x = part[i].r.x;
	part[i].ro.y = part[i].r.y;
	part[i].ro.z = part[i].r.z;

	// Displacement between -0.5 and +0.5 in each dimension
	part[i].r.x += ran3() - 0.5;
	part[i].r.y += ran3() - 0.5;
	part[i].r.z += ran3() - 0.5;

	// Periodic boundary conditions
	if (part[i].r.x > sys.length)
	{
		part[i].r.x -= sys.length;
	}
	else if (part[i].r.x < 0)
	{
		part[i].r.x += sys.length;
	}
	if (part[i].r.y > sys.length)
	{
		part[i].r.y -= sys.length;
	}
	else if (part[i].r.y < 0)
	{
		part[i].r.y += sys.length;
	}
	if (part[i].r.z > sys.length)
	{
		part[i].r.z -= sys.length;
	}
	else if (part[i].r.z < 0)
	{
		part[i].r.z += sys.length;
	}

	// Calculate new energies of particles within two cells of particle i.
	if (sys.calc_list == 1)
	{
		// Determine cell number of particle i
		ix = (int) part[i].r.x / sys.r_cell;
		iy = (int) part[i].r.y / sys.r_cell;
		iz = (int) part[i].r.z / sys.r_cell;

		for (l=-2; l<=2; l++)
		{
			for (m=-2; m<=2; m++)
			{
				for (n=-2; n<=2; n++)
				{
					// Determine the cell
					jx = mod(ix+l,sys.n_cell);
					jy = mod(iy+m,sys.n_cell);
					jz = mod(iz+n,sys.n_cell);

					j = sys.hoc[jx][jy][jz];

					while (j!=-1)
					{
						part[j].E = calc_energy_list(j);
						// Next particle in the linked list
						j = part[j].ll;
					}
				}
			}
		}
	}
	else
	{
		// Brute force calculation of energies
		for (j=0; j<sys.n_dpd; j++)
		{
			part[j].E = calc_energy(j);
		}
	}
}


