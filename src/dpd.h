/* Header file for DPD Monte Carlo simulation program */


/* Inline substitutions */

#define P	sys.stats[0]

/* Structure defintions */

typedef struct vector_type
{
	double x, y, z;
} Vector;

typedef struct ivector_type
{
	int ix, iy, iz;
} Ivector;

typedef struct particle_type
{
	int ll;
	double E, Eo;
	Vector r, ro;
	// Ivector c, co;
} Particle;

typedef struct stats_type
{
	char *name;
	int num;
	double now, sum, sumsq, err;
} Stats;


typedef struct parameter_type
{
	int 	calc_list,
		***hoc;

	long 	freq_sample,
		iseed,
		nsteps,
	        n_cell,
		n_dpd,
		n_stats;

	double 	a_ij,
		density,
		length,
		r_c,
		r_cell,
		temp,
		volume;

	Stats *stats;
} System;


/* Global variables */

extern Particle *part;
extern System sys;


/* Global functions */

extern double 	calc_energy (int);
extern double 	calc_energy_list (int);
extern double	energy_ij (int, int);
extern void 	initialize (void);
extern void 	init_param (void);
extern void 	init_stats (void);
extern void 	input (void);
extern int	mod (int, int);
extern void 	monte_carlo (void);
extern void	new_list (void);
extern void 	output (void);
extern double 	ran3 (void);
extern void 	random_move (int);
extern void 	sample (void);
extern void 	setup_coords (void);
extern double	total_energy (void);
extern void 	write_log (void);
extern void	write_mon (void);
extern Vector 	vdist (Vector, Vector);
extern double 	vmag (Vector);

