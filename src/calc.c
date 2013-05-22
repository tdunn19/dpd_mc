/* calc.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void sample (void)
{
	int i, j;
	double r_ij;
	Vector dr;

	P.now = 0;

	for (i=0; i<sys.n_dpd; i++)
	{
		for (j=i+1; j<sys.n_dpd; j++)
		{
			dr = vdist(part[i].r,part[j].r);

			// Employ periodic boundary conditions
			if (dr.x > sys.length/2)
			{
				dr.x -= sys.length;
			}
			else if (dr.x < -sys.length/2)
			{
				dr.x += sys.length;
			}
			if (dr.y > sys.length/2)
			{
				dr.y -= sys.length;
			}
			else if (dr.y < -sys.length/2)
			{
				dr.y += sys.length;
			}
			if (dr.z > sys.length/2)
			{
				dr.z -= sys.length;
			}
			else if (dr.z < -sys.length/2)
			{
				dr.z += sys.length;
			}

			r_ij = vmag(dr);

			// Check for cutoff distance
			if (r_ij < sys.r_c)
			{
				P.now += sys.a_ij*r_ij*(1-r_ij);
			}
		}
	}

	P.now /= 3*sys.volume;
	P.now += sys.density*sys.temp;

	for (i=0; i<sys.n_stats; i++)
	{
		sys.stats[i].sum += sys.stats[i].now;
		sys.stats[i].sumsq += sys.stats[i].now*sys.stats[i].now;
		sys.stats[i].num += 1;
	}
}

double calc_energy (int i)
{
	int j;
	double E, r_ij;
	Vector dr;
	E = 0;

	// Soft repulsive force
	for (j=0; j<sys.n_dpd; j++)
	{
		if (j!=i)
		{
			E += energy_ij(i,j);
		}
	}

	E *= sys.a_ij/2;

	return E;
}

double energy_ij (int i, int j)
{
	double r_ij, E_ij;
	Vector dr;

	dr = vdist(part[i].r,part[j].r);

	// Periodic boundary conditions
	if (dr.x > sys.length/2)
	{
		dr.x -= sys.length/2;
	}
	else if (dr.x < -sys.length/2)
	{
		dr.x += sys.length/2;
	}
	if (dr.y > sys.length/2)
	{
		dr.y -= sys.length/2;
	}
	else if (dr.y < -sys.length/2)
	{
		dr.y += sys.length/2;
	}
	if (dr.z > sys.length/2)
	{
		dr.z -= sys.length/2;
	}
	else if (dr.z < -sys.length/2)
	{
		dr.z += sys.length/2;
	}

	r_ij = vmag (dr);

	if (r_ij < sys.r_c)
	{
		E_ij = (1-r_ij/sys.r_c)*(1-r_ij/sys.r_c);
	}
	else
	{
		E_ij = 0;
	}

	return E_ij;
}

double total_energy (void)
{
	int i;
	double E;
	E = 0;

	for (i=0; i<sys.n_dpd; i++)
	{
		E += part[i].E;
	}

	return E;
}


// Maybe come back to this later if efficiency is an issue.
// Seems overly complicated right now. (May 21st, 2013)
/*
double total_energy_list (int i)
{
	// This will calculate new energies of the system after a MC move

	int ix,iy,iz,jx,jy,jz,j,l,m,n;

	// Determine cell number of particle i.
	ix = (int) part[i].r.x / sys.r_c;
	iy = (int) part[i].r.y / sys.r_c;
	iz = (int) part[i].r.z / sys.r_c;

	// Recalculate energies for the cells within two cell lengths to
	// account for the particle moving into a new cell.
	for (l=-2; l<=2; l++)
	{
		for (m=-2; m<=2; m++)
		{
			for (n=-2; n<=2; n++)
			{
				jx = mod(ix+l,sys.n_cell);
				jy = mod(iy+m,sys.n_cell);
				jz = mod(iz+n,sys.n_cell);

				// First particle in the chain:
				j = sys.hoc[jx][jy][jz];

				while (j!=-1)
				{
					if (i!=j)
					{

					}
				}

				j
			}
		}
	}
}
*/
