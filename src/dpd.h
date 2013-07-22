/* Header file for DPD Monte Carlo simulation program */

// Inline substitutions

#define P      sys.stats[0]
#define Etot   sys.stats[1]
#define RE2    sys.stats[2]
#define RE2x   sys.stats[3]
#define RE2y   sys.stats[4]
#define RE2z   sys.stats[5]
#define RG2    sys.stats[6]
#define RG2x   sys.stats[7]
#define RG2y   sys.stats[8]
#define RG2z   sys.stats[9]

/* Structure defintions */

typedef struct vector_type {
  double x, y, z;
} Vector;

typedef struct ivector_type {
  int x, y, z;
} Ivector;

typedef struct particle_type {
  int ll, side;
  double E, Eo;
  Vector r, ro;
} Particle;

typedef struct stats_type {
  char *name;
  int num;
  double now, sum, sumsq, err;
} Stats;

typedef struct monitor_type {
  double *energy, *Q,
    *re2, *rex, *rey, *rez,
    *re2_cis, *rex_cis, *rey_cis, *rez_cis,
    *re2_trans, *rex_trans, *rey_trans, *rez_trans,
    *rg2, *rg2x, *rg2y, *rg2z,
    *rg2_cis, *rg2x_cis, *rg2y_cis, *rg2z_cis,
    *rg2_trans, *rg2x_trans, *rg2y_trans, *rg2z_trans,
    *cmx, *cmy, *cmz;
} Monitor;

typedef struct parameter_type {
  int
    *bin_count,
    bond_break,
    calc_list,
    freq_sample,
    iQ_init,
    iseed,
    ***hoc,
    ***hoc_copy,
    monitor_step,
    n_accept_mon,
    n_accept_solvent,
    n_attempt_mon,
    n_attempt_solvent,
    n_bins,
    n_cis,
    n_dpd,
    n_layers,
    n_mon,
    n_pore,
    n_solvent,
    n_solvent_cis,
    n_solvent_trans,
    n_stats,
    n_steps,
    n_trans,
    n_wall,
    n_wins,
    new_cell,
    new_window,
    step,
    wall_overlap;

  double
    a_mm,
    a_ms,
    a_ss,
    a_sw,
    bin_width,
    bl_init,
    density_cis,
    density_trans,
    density_w,
    dr_max_dpd,
    dr_max_mon,
    energy,
    energy_old,
    k_fene,
    mc_ratio,
    pore_length,
    pore_radius,
    pore_volume,
    Q,
    Q_init,
    Q_max,
    Q_min,
    r_c,
    r_c2,
    r_eq,
    r_max,
    r_max2,
    r_pore,
    r_wall,
    r_0,
    volume,
    volume_cis,
    volume_trans,
    wall_volume,
    window_width;

  Vector
    cm,
    cm_cis,
    cm_trans,
    length,
    pore_max,
    pore_min,
    r_cell,
    re,
    re_cis,
    re_trans,
    rg2,
    rg2_cis,
    rg2_trans,
    wall_max,
    wall_min;

  Ivector n_cell_1d, n_pore_1d, n_wall_1d;
  Stats *stats;
  Monitor mon;
} System;


/* Global variables */

extern Particle *part_dpd;
extern Particle *part_mon;
extern System sys;
extern const double root_table[100001];

/* Global functions */

extern int     accept_move(void);
extern void    calc_energy_brute(void);
extern double  calc_energy_dpd(int i);
extern double  calc_energy_mon(int i);
extern void    calc_cm(void);
extern void    calc_nseg(void);
extern void    calc_pressure(void);
extern void    calc_q(void);
extern void    calc_re(void);
extern void    calc_rg(void);
extern void    check_bond(int i);
extern void    check_cell(Vector, Vector);
extern int     check_pore(Vector);
extern int     check_side(Vector);
extern void    check_wall(Vector);
extern double  energy_c(Vector);
extern double  energy_fene(Vector);
extern void    initialize(void);
extern void    init_energy(void);
extern void    init_monitor(void);
extern void    init_param(void);
extern void    init_polymer(void);
extern void    init_pore(void);
extern void    init_solvent(void);
extern void    init_stats(void);
extern void    init_wall(void);
extern void    input(void);
extern int     mod(int, int);
extern void    monte_carlo(void);
extern void    move_monomer(int i);
extern void    move_solvent(int i);
extern void    new_list(void);
extern void    periodic_bc_dr(Vector *);
extern void    periodic_bc_r(Vector *);
extern void    print_stats(void);
extern double  ran3(void);
extern double  root(double);
extern void    sample(void);
extern void    update_monitor(void);
extern void    update_stats(void);
extern void    write_log(void);
extern void    write_mon(void);
extern Vector  vdist(Vector, Vector);
extern double  vmag(Vector);

