#ifndef __KALMAN_H_INCLUDED__
#define __KALMAN_H_INCLUDED__

double *kalman_step(double *data, double *observed, double *x0,
                int data_dim, int obs_dim, int x_dim, 
                double *forward_fun(double*, double*),
                double **P0, double **Q, double **R,
                double *lower_bounds, double *upper_bounds,
                int num_iters, double stop_val
                );

double **kalman_gain(double **P, double **H, double **R, 
                     int data_dim, int x_dim);

double **update_cov(double **P0, double **K, double **H, 
                    int data_dim, int x_dim);

#endif
