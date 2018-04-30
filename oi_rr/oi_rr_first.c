#include <math.h>
#include<stdio.h>


/* */ 
void oi_rr_first(int *no, 
                 double *innov,
                 double *SRinv,
                 double *vec) {
  int nov = no[0];
  int i,j,k;

  k=0;
  for (i=0;i<nov;i++) {
    vec[i]=0;
    for (j=0;j<nov;j++) {
      vec[i]=vec[i]+SRinv[k]*innov[j];
      k++;
    }
  }
  return;
}
