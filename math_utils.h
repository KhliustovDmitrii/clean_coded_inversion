#ifndef __MATH_UTILS_INCLUDED__
#define __MATH_UTILS_INCLUDED__
long double bessel_function(long double x);

long double scalar_product(long double *u, long double *v, int vec_len);

long double **matrix_product(long double **M, long double **N, int m1, int m2, int n1, int n2);

long double **pseudo_inv(long double **M, int tm1, int m2);

long double **invert(long double **mat, int n);

long double **transpose(long double **M, int m1, int m2);

long double **cholesky(long double **M, int m);

long double primal_field(long double hor_dist, long double ver_dist);

long double _Complex integral(long double _Complex (*f)(long double, long double*), long double *f_pars,
                        int par_num, long double t0, long double t1, long double delta);

long double **jacobian(long double* (*f)(long double*, long double*), long double *f_pars, long double *point,
                  int in_dim, int out_dim, int par_dim, long double delta);

long double discrepancy(long double *diff, long double **R, int dim);
#endif
