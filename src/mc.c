/* mc.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void monte_carlo (void)
{
	int i, j;
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
	// Else, accept the movement
}

void random_move (int i)
{
	int j;

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

	// Calculate new energies
	for (j=0; j<sys.n_dpd; j++)
	{
		part[j].Eo = part[j].E;
		if (sys.calc_list == 1)
		{
			part[j].E = calc_energy_list(j);
		}
		else
		{
			part[j].E = calc_energy(j);
		}
	}
}


