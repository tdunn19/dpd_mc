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
    fscanf(fp, "%lf%*s", &sys.density);
    fscanf(fp, "%lf%*s", &sys.volume);
    fscanf(fp, "%lf%*s", &sys.r_c);
    fscanf(fp, "%lf%*s", &sys.a_mm);
    fscanf(fp, "%lf%*s", &sys.a_ms);
    fscanf(fp, "%lf%*s", &sys.a_ss);
    fscanf(fp, "%d%*s", &sys.nsteps);
    fscanf(fp, "%lf%*s", &sys.temp);
    fscanf(fp, "%d%*s", &sys.freq_sample);
    fscanf(fp, "%d%*s", &sys.freq_monitor);
    fscanf(fp, "%d%*s", &sys.calc_list);
    fscanf(fp, "%d%*s", &sys.iseed);
  }

  fclose(fp);
}

void write_log(void) {
  printf("\nDPD Monte Carlo Simulation Program");
  printf("\n\n");
  printf("n_mon       \t\t\t%10d\n", sys.n_mon);
  printf("density     \t\t\t%10.5lf\n", sys.density);
  printf("volume      \t\t\t%10.5lf\n", sys.volume);
  printf("n_dpd       \t\t\t%10d\n", sys.n_dpd);
  printf("r_c         \t\t\t%10.5lf\n", sys.r_c);
  printf("a_mm        \t\t\t%10.5lf\n", sys.a_mm);
  printf("a_ms        \t\t\t%10.5lf\n", sys.a_ms);
  printf("a_ss        \t\t\t%10.5lf\n", sys.a_ss);
  printf("nsteps      \t\t\t%10d\n", sys.nsteps);
  printf("temp        \t\t\t%10.5lf\n", sys.temp);
  printf("freq_sample \t\t\t%10d\n", sys.freq_sample);
  printf("freq_monitor\t\t\t%10d\n", sys.freq_monitor);
  printf("calc_list   \t\t\t%10d\n", sys.calc_list);
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

  printf("\nFinal averages\n\n");

  for (i = 0; i < sys.n_stats; i++) {
    printf("%s %12.8lf +/- %12.8lf\n",
      sys.stats[i].name, sys.stats[i].sum, sys.stats[i].err);
  }
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

  if ((fp = fopen("bond_length.dat", "w")) == NULL) {
    printf("Cannot open file: bond_length.dat\n");
  } else {
    for (i = 0; i < sys.monitor_step; i++) {
      fprintf(fp, "%d  %lf\n", i, sys.mon.bond_length[i]);
    }
  }
}
