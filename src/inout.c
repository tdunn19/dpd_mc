#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "dpd.h"

void input(void) {
  FILE *fp;
  if ((fp = fopen("dpd.inp", "r")) == NULL) {
    printf("Cannot open file: dpd.inp\n");
  } else {
    fscanf(fp, "%d%*s", &sys.n_mon);
    fscanf(fp, "%lf%*s", &sys.density_s);
    fscanf(fp, "%d%*s", &sys.calc_list);
    fscanf(fp, "%lf%*s", &sys.length.x);
    fscanf(fp, "%lf%*s", &sys.length.y);
    fscanf(fp, "%lf%*s", &sys.length.z);
    fscanf(fp, "%lf%*s", &sys.density_w);
    fscanf(fp, "%lf%*s", &sys.r_c);
    fscanf(fp, "%lf%*s", &sys.dr_max_dpd);
    fscanf(fp, "%lf%*s", &sys.dr_max_mon);
    fscanf(fp, "%lf%*s", &sys.a_mm);
    fscanf(fp, "%lf%*s", &sys.a_ms);
    fscanf(fp, "%lf%*s", &sys.a_ss);
    fscanf(fp, "%lf%*s", &sys.a_sw);
    fscanf(fp, "%d%*s", &sys.n_layers);
    fscanf(fp, "%d%*s", &sys.n_steps);
    fscanf(fp, "%lf%*s", &sys.mc_ratio);
    fscanf(fp, "%lf%*s", &sys.temp);
    fscanf(fp, "%d%*s", &sys.freq_sample);
    fscanf(fp, "%d%*s", &sys.freq_monitor);
    fscanf(fp, "%d%*s", &sys.iseed);
  }

  fclose(fp);
}

void write_log(void) {
  printf("\nDPD Monte Carlo Simulation Program");
  printf("\n\n");
  printf("n_mon       \t\t\t%10d\n", sys.n_mon);
  printf("density_s   \t\t\t%10.5lf\n", sys.density_s);
  printf("n_dpd       \t\t\t%10d\n", sys.n_dpd);
  printf("calc_list   \t\t\t%10d\n\n", sys.calc_list);
  printf("length_x    \t\t\t%10.5lf\n", sys.length.x);
  printf("length_y    \t\t\t%10.5lf\n", sys.length.y);
  printf("length_z    \t\t\t%10.5lf\n\n", sys.length.z);
  printf("volume      \t\t\t%10.5lf\n\n", sys.volume);
  printf("density_w   \t\t\t%10.5lf\n", sys.density_w);
  printf("n_layers    \t\t\t%10d\n", sys.n_layers);
  printf("n_wall      \t\t\t%10d\n\n", sys.n_wall);
  printf("r_c         \t\t\t%10.5lf\n", sys.r_c);
  printf("dr_max_dpd  \t\t\t%10.5lf\n", sys.dr_max_dpd);
  printf("dr_max_mon  \t\t\t%10.5lf\n\n", sys.dr_max_mon);
  printf("a_mm        \t\t\t%10.5lf\n", sys.a_mm);
  printf("a_ms        \t\t\t%10.5lf\n", sys.a_ms);
  printf("a_ss        \t\t\t%10.5lf\n", sys.a_ss);
  printf("a_sw        \t\t\t%10.5lf\n\n", sys.a_sw);
  printf("nsteps      \t\t\t%10d\n", sys.n_steps);
  printf("mc_ratio    \t\t\t%10.5lf\n", sys.mc_ratio);
  printf("temp        \t\t\t%10.5lf\n", sys.temp);
  printf("freq_sample \t\t\t%10d\n", sys.freq_sample);
  printf("freq_monitor\t\t\t%10d\n\n", sys.freq_monitor);
  printf("iseed       \t\t\t%10d\n", sys.iseed);
  printf("\n\n");
}

void final_stats(void) {
  int i;

  for (i = 0; i < sys.n_stats; i++) {
    sys.stats[i].sum /= sys.stats[i].num;
    sys.stats[i].sumsq /= sys.stats[i].num;
    sys.stats[i].err =
      sqrt(sys.stats[i].sumsq - sys.stats[i].sum*sys.stats[i].sum);
  }
}

void print_stats(void) {
  int i;
  double perc_dpd, perc_mon;

  printf("\nFinal averages\n\n");

  for (i = 0; i < sys.n_stats; i++) {
    printf("%s %12.8lf +/- %12.8lf\n",
      sys.stats[i].name, sys.stats[i].sum, sys.stats[i].err);
  }

  perc_dpd = (double) sys.n_accept_dpd / sys.n_attempt_dpd;
  perc_mon = (double) sys.n_accept_mon / sys.n_attempt_mon;

  printf("\n\nSelection and acceptance stats\n\n");
  printf("Particle    Moves    Accepted  Percent\n");
  printf("--------  ---------  --------  -------\n");
  printf("     DPD  %9d%10d%9.2lf\n",
    sys.n_attempt_dpd, sys.n_accept_dpd, perc_dpd);
  printf(" Monomer  %9d%10d%9.2lf\n",
    sys.n_attempt_mon, sys.n_accept_mon, perc_mon);
}

void write_mon(void) {
  int i;
  FILE *fp;

  if ((fp = fopen("energy.dat", "w")) == NULL) {
    printf("Cannot open file: energy.dat\n");
  } else {
    for (i = 0; i < sys.monitor_step; i++) {
      fprintf(fp, "%d  %lf\n", i, sys.mon.energy[i]);
    }
  }

  if ((fp = fopen("re.dat", "w")) == NULL) {
    printf("Cannot open file: re.dat\n");
  } else {
    for (i = 0; i < sys.monitor_step; i++) {
      fprintf(fp, "%d  %lf %lf %lf %lf\n",
       i, sys.mon.re2[i], sys.mon.rex[i], sys.mon.rey[i], sys.mon.rez[i]);
    }
  }

  if ((fp = fopen("rg.dat", "w")) == NULL) {
    printf("Cannot open file: rg.dat\n");
  } else {
    for (i = 0; i < sys.monitor_step; i++) {
      fprintf(fp, "%d  %lf %lf %lf %lf\n",
       i, sys.mon.rg2[i], sys.mon.rgx[i], sys.mon.rgy[i], sys.mon.rgz[i]);
    }
  }

  if ((fp = fopen("cm.dat", "w")) == NULL) {
    printf("Cannot open file: cm.dat\n");
  } else {
    for (i = 0; i < sys.monitor_step; i++) {
      fprintf(fp, "%d  %lf %lf %lf\n",
       i, sys.mon.cmx[i], sys.mon.cmy[i], sys.mon.cmz[i]);
    }
  }

  if ((fp = fopen("bond_length.dat", "w")) == NULL) {
    printf("Cannot open file: bond_length.dat\n");
  } else {
    for (i = 0; i < sys.monitor_step; i++) {
      fprintf(fp, "%d  %lf\n", i, sys.mon.bond_length[i]);
    }
  }
}
