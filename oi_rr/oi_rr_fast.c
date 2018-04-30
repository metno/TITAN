#include <math.h>
#include<stdio.h>


/* bilinear interpolation from a coarser grid (cg) to a finer grid (fg) */
void oi_rr_fast(int *ng, 
                int *no, 
                double *xg,
                double *yg,
                double *zg,
                double *xo,
                double *yo,
                double *zo,
                double *Dh, 
                double *Dz,
                double *xb, 
                double *vec,
                double *xa) {
  int ngv = ng[0];
  int nov = no[0];
  double Dhv = Dh[0];
  double Dzv = Dz[0];
  int i,j;
  double g,Dh2,Dz2,hd2,vd2;

  Dh2=Dhv*Dhv;
  Dz2=Dzv*Dzv;
  for (i=0;i<ngv;i++) {
    xa[i]=xb[i];
  }
  for (j=0;j<nov;j++) { 
    for (i=0;i<ngv;i++) {
      hd2=((xg[i]-xo[j])*(xg[i]-xo[j])+(yg[i]-yo[j])*(yg[i]-yo[j])) / (1000.*1000.);
      vd2=(zg[i]-zo[j])*(zg[i]-zo[j]);
      g=exp(-0.5*(hd2/Dh2+vd2/Dz2));
      xa[i]=xa[i]+g*vec[j];
    }
  }
  for (i=0;i<ngv;i++) {
    if (xa[i]<0) xa[i]=0.;
  }

  return;
}
