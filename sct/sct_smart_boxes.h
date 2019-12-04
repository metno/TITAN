#include <assert.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void sct_smart_boxes(int *n, double *x, double *y, double *z, double *t,
                     int *nmax, int *nmin, int *nminprof, double *dzmin,
                     double *dhmin, double *dz, double *t2pos, double *t2neg,
                     double *eps2, int *flags, double *sct, double *rep,
                     int *boxids);

struct Box {
  int n;
  double *x;
  double *y;
  double *z;
  double *t;
  int *i;
};

struct BoxList {
  int n;
  struct Box *boxes;
};

void spatial_consistency_test(struct Box *currentBox, int *nminprof,
                              double *dzmin, double *dhmin, double *dz,
                              double *t2pos, double *t2neg, double *eps2,
                              int *flags, double *sct_out, double *rep_out);

int vertical_profile_optimizer(gsl_vector *input, struct Box *currentBox,
                               int nminprof, double dzmin, double *vp);

double basic_vertical_profile_optimizer_function(const gsl_vector *v,
                                                 void *data);
double vertical_profile_optimizer_function(const gsl_vector *v, void *data);

void basic_vertical_profile(int nz, double *z, double t0, double *t_out);

void vertical_profile(int nz, double *z, double t0, double gamma, double a,
                      double h0, double h1i, double *t_out);

struct BoxList control_box_division(int maxNumStationsInBox,
                                    int minNumStationsInBox,
                                    struct Box inputBox);

double compute_quantile(double quantile, double *array, int sizeArray);
double mean(const double *array, int sizeArray);
double max(const double *array, int sizeArray);
double min(const double *array, int sizeArray);
void print_vector(double *vector, int size);
void print_gsl_vector(gsl_vector *vector, int size);
void print_matrix(double **matrix, int rows, int columns);
void print_gsl_matrix(gsl_matrix *matrix, int rows, int columns);
void print_sub_gsl_matrix(gsl_matrix *matrix, int start, int stop);
gsl_matrix *inverse_matrix(const gsl_matrix *matrix);
