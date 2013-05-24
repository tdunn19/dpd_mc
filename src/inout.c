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
        fscanf(fp, "%lf%*s", &sys.a_ij);
        fscanf(fp, "%d%*s", &sys.nsteps);
        fscanf(fp, "%lf%*s", &sys.temp);
        fscanf(fp, "%d%*s", &sys.freq_sample);
        fscanf(fp, "%d%*s", &sys.calc_list);
        fscanf(fp, "%d%*s", &sys.iseed);
    }

    fclose(fp);
}

void write_log(void) {
    printf("\nDPD Monte Carlo Simulation Program");
    printf("\n\n");
    printf("n_mon      \t\t\t%10d\n", sys.n_mon);
    printf("density    \t\t\t%10.5lf\n", sys.density);
    printf("volume     \t\t\t%10.5lf\n", sys.volume);
    printf("n_dpd      \t\t\t%10d\n", sys.n_dpd);
    printf("r_c        \t\t\t%10.5lf\n", sys.r_c);
    printf("a_ij       \t\t\t%10.5lf\n", sys.a_ij);
    printf("nsteps     \t\t\t%10d\n", sys.nsteps);
    printf("temp       \t\t\t%10.5lf\n", sys.temp);
    printf("freq_sample\t\t\t%10d\n", sys.freq_sample);
    printf("calc_list  \t\t\t%10d\n", sys.calc_list);
    printf("iseed      \t\t\t%10d\n", sys.iseed);
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

void output(void) {
    int i;

    printf("\nFinal averages\n\n");

    for (i = 0; i < sys.n_stats; i++) {
        printf("%s \t %lf +/- %lf\n",
            sys.stats[i].name, sys.stats[i].sum, sys.stats[i].err);
    }
}

void write_mon(void) {
    int i;
    FILE *fp;

    if ((fp = fopen("pressure.dat", "w")) == NULL) {
        printf("Cannot open file: pressure.dat\n");
    } else {
        fprintf(fp, "%lf +/- %lf\n", P.sum, P.err);
    }
}
