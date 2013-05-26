/* Header file for DPD Monte Carlo simulation program */

// Inline substitutions

#define P   sys.stats[0]

/* Structure defintions */

typedef struct vector_type {
  double x, y, z;
} Vector;

typedef struct ivector_type {
  int ix, iy, iz;
} Ivector;

typedef struct particle_type {
  int ll;
  double E, Eo;
  Vector r, ro;
} Particle;

typedef struct stats_type {
  char *name;
  int num;
  double now, sum, sumsq, err;
} Stats;

typedef struct parameter_type {
  int
    bond_break,
    calc_list,
    freq_sample,
    iseed,
    ***hoc,
    nsteps,
    n_cell,
    n_dpd,
    n_mon,
    n_stats;

  double
    a_mm,
    a_ms,
    a_ss,
    density,
    energy,
    k_fene,
    length,
    r_c,
    r_cell,
    r_eq,
    r_max,
    r_0,
    temp,
    volume;

  Stats *stats;
} System;


/* Global variables */

extern Particle *part_dpd;
extern Particle *part_mon;
extern System sys;


/* Global functions */

extern int     accept_move(void);
extern double  calc_energy(int i);
extern double  calc_energy_dpd(int i);
extern double  calc_energy_mon(int i);
extern int     check_cell(Vector, Vector);
extern double  energy_c(Vector);
extern double  energy_fene(int i, int j);;
extern void    initialize(void);
extern void    init_param(void);
extern void    init_stats(void);
extern void    input(void);
extern int     mod(int, int);
extern void    monte_carlo(void);
extern void    new_list(void);
extern void    output(void);
extern double  ran3(void);
extern void    random_move_dpd(int i);
extern void    random_move_mon(int i);
extern void    sample(void);
extern void    setup_coords(void);
extern double  total_energy(void);
extern void    write_log(void);
extern void    write_mon(void);
extern Vector  vdist(Vector, Vector);
extern double  vmag(Vector);

