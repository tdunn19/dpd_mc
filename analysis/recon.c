#include <math.h>
#include <stdio.h> 
#include <string.h>
#include <stdlib.h>

int main()
{
  int i, j, k, nbin, nwin, nbintot, index, nbinhalf, symmetry;
  double delQ, M, x0, x1, x2, y0, y1, y2, c;

  double A[901];
  double B[901];
  double C[901];
  double Z[901];
  double F[10001];
  double PR[10001];

  int ibin[900][100];
  double PQ[900][100];

  scanf("%d",&nbin);
  scanf("%d",&nwin);
  scanf("%lf",&delQ);
  scanf("%d",&symmetry);

  nbinhalf = nbin/2;
  nbintot = 1.0/delQ + 0.5;

  for (i=1;i<=nwin;i++) {
    for (j=0;j<nbin;j++) {
      scanf("%d  %lf\n",&ibin[i][j],&PQ[i][j]);
    }
  }
  for (i=1;i<=nwin;i++) {
    M = 0.0;
    for (j=0;j<nbin;j++) M += PQ[i][j];
    for (j=0;j<nbin;j++) PQ[i][j] /= M;
    //printf("%d   %lf\n",i,M);
  }
  //exit(0);

  A[1] = 0;
  B[1] = -1.0;
  for (i=0;i<nbinhalf;i++)  B[1] += PQ[1][i]; 
  C[1] = 0.0;
  for (i=0;i<nbinhalf;i++)  C[1] += (PQ[1][i+nbinhalf] + PQ[2][i]);

  for (j=2;j<nwin;j++) {
    A[j] = 0.0;
    for (i=0;i<nbinhalf;i++)  A[j] += (PQ[j-1][i+nbinhalf] + PQ[j][i]);
    B[j] = -1.0;
    C[j] = 0.0;
    for (i=0;i<nbinhalf;i++)  C[j] += (PQ[j+1][i] + PQ[j][i+nbinhalf]);
  }

  Z[1] = 1.0;
  Z[2] = -1.0*B[1]/(B[1]+C[1]);

  for (i=2;i<nwin;i++) {
    Z[i+1] = Z[i]*Z[i-1]*A[i] + (Z[i]*Z[i-1]+Z[i]*Z[i]) * B[i];
    Z[i+1] /= Z[i-1]*A[i] + (Z[i-1]+Z[i]) * B[i] + (Z[i-1]+Z[i])*C[i];
    Z[i+1] *= -1.0;
  }


  for (j=0;j<nbintot;j++) PR[j] = 0.0;

  for (i=0;i<nbinhalf;i++) PR[i] = PQ[1][i];

  for (j=2;j<=nwin;j++) {
    for (i=0;i<nbinhalf;i++) {
      index = (j-1)*nbinhalf+i;
      PR[index] = PQ[j-1][i+nbinhalf] + PQ[j][i];
      PR[index] /= 1.0/Z[j-1] + 1.0/Z[j];
    }
  }
  for (i=nbinhalf;i<nbin;i++) {
    index = (nwin-1)*nbinhalf+i; 
    PR[index] = PQ[nwin][i];
    PR[index] /= 1.0/Z[nwin];
  }

  for (i=0;i<nbintot;i++) F[i] = -100.0;

  x0 = 0.5*delQ;
  y0 = -log(PR[0]);

  x1 = 1.5*delQ;
  y1 = -log(PR[1]);

  x2 = 2.5*delQ;
  y2 = -log(PR[2]);

  c = y0*(x1*x1*x2 - x1*x2*x2) - y1*(x0*x0*x2 - x0*x2*x2) + y2*(x0*x0*x1 - x0*x1*x1);
  c /= x1*x1*x2 - x1*x2*x2 - x0*x0*x2 + x0*x2*x2 + x0*x0*x1 - x0*x1*x1;

  if (symmetry) {
    for (i=0;i<nbintot/2;i++) {
      if (PR[i] > 0.0) {
        F[i] = -log(PR[i]) - c;
        F[nbintot-i-1] = F[i];
      }
    }
  } else {
    for (i=0;i<nbintot;i++) {
      if (PR[i] > 0.0) {
        F[i] = -log(PR[i]) - c;
      }
    }
  }

  printf("%lf  %lf\n",0.0,0.0);
  for (i=0;i<nbintot;i++) {
    if (F[i] > -99.0) printf("%lf  %lf\n",(i+0.5)*delQ,F[i]);
  }


}
