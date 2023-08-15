#ifndef __MATH_UTILS_INCLUDED__
#define __MATH_UTILS_INCLUDED__
double bessel_function(double x);

double scalar_product(double *u, double *v, int vec_len);

double **matrix_product(double **M, double **N, int m1, int m2, int n1, int n2);

double **pseudo_inv(double **M, int tm1, int m2);

double **invert(double **mat, int n);

double **transpose(double **M, int m1, int m2);

double primal_field(double hor_dist, double ver_dist);

double _Complex integral(double _Complex (*f)(double, double*), double *f_pars,
                        int par_num, double t0, double t1, double delta);

double **jacobian(double* (*f)(double*, double*), double *f_pars, double *point,
                  int in_dim, int out_dim, int par_dim, double delta);

double discrepancy(double *diff, double **R, int dim);
#endif
