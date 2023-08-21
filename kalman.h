#ifndef __KALMAN_H_INCLUDED__
#define __KALMAN_H_INCLUDED__

long double *kalman_sequential(long double *data, long double *observed,
                long double *x0,
                int data_dim, int obs_dim, int x_dim,
                long double* (*forward_fun)(long double*, long double*),
                long double **P0, long double **Q, long double **R,
                long double *lower_bounds, long double *upper_bounds,
                int num_iters, long double stop_val
                );


long double *kalman_step(long double *data, long double *observed, 
                long double *x0,
                int data_dim, int obs_dim, int x_dim, 
                long double* (*forward_fun)(long double*, long double*),
                long double **P0, long double **Q, long double **R,
                long double *lower_bounds, long double *upper_bounds,
                int num_iters, long double stop_val
                );

long double **kalman_gain(long double **P, long double **H, long double **R, 
                     int data_dim, int x_dim);

long double **update_cov(long double **P0, long double **K, long double **H, 
                    int data_dim, int x_dim);

#endif
