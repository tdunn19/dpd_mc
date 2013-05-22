/* init.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "dpd.h"

void initialize (void)
{
	input ();
	init_param();
	init_stats();
	write_log();
	setup_coords();
}

void init_param (void)
{
	int fact;
	int i, j;

	sys.n_dpd = sys.density * sys.volume;
	sys.length = pow(sys.volume,1.0/3.0);

	if (sys.calc_list == 1)
	{
		// Determine cell size and the number of cells
		fact = (int) (sys.length / sys.r_c + 1e-9);
		sys.r_cell = (double) sys.length / fact;
		sys.n_cell = (long) (sys.length / sys.r_cell + 1e-9);

		// Allocate memory for the head of chain array
		sys.hoc = (int ***) malloc((sys.n_cell+1)*sizeof(int **));
		for (i=0; i<sys.n_cell; i++)
		{
			sys.hoc[i] = (int **) malloc((sys.n_cell+1)*sizeof(int *));
			for (j=0; j<sys.n_cell; j++)
			{
				sys.hoc[i][j] = (int *) malloc((sys.n_cell+1)*sizeof(int));
			}
		}
	}
}

void init_stats (void)
{
	sys.n_stats = 1;
	sys.stats = (Stats *) calloc(sys.n_stats,sizeof(Stats));
	sys.stats[0].name = "Pressure            ";
}

void setup_coords (void)
{
	int i;

	part = (Particle *) calloc(sys.n_dpd,sizeof(Particle));

	for (i=0; i<sys.n_dpd; i++)
	{
		part[i].r.x = sys.length*ran3();
		part[i].r.y = sys.length*ran3();
		part[i].r.z = sys.length*ran3();
	}

	if (sys.calc_list == 1) new_list();

	for (i=0; i<sys.n_dpd; i++)
	{
		if (sys.calc_list == 1)
		{
			part[i].E = calc_energy_list(i);
		}
		else
		{
			part[i].E = calc_energy(i);
		}
	}
}
