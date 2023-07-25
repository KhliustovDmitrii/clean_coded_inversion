#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 
 *This file contains code for Kalman filter
 */

double **kalman_step(double *data, double *observed, double *x0,
                int data_dim, int obs_dim, int x_dim, 
                double *forward_fun(double*, double*),
                double **P0, double **Q, double **R,
                double *lower_bounds, double *upper_bounds,
                int num_iters, double stop_val
                )
{
/*
 * This function computes Kalman inversion of data
 *Params:
 *        data: measured data to invert
 *        observed: known data for positioning etc.
 *        x0: initial estimate of hidden parameters
 *        data_dim: dimensionality of measurements
 *        obs_dim: dimensionality of observed 
 *        x_dim: number of hidden parameters to estimate
 *        forward_fun: function for solvin forward problem with signature
 *        f(x, pars)
 *        P0, Q, R: covariance of parameters, conseq steps and measurements
 *        *bounds: lower and upper bounds for parameters
 *        num_iters: max number of iterations
 *        stop_val: if discrepancy < stop_val, go to next point
 *Returns:
 *        result: vector, dim = 2*x_dim + 1,
 *        first elements contain estimated parameters,
 *        middle - discrepancy value, 
 *        last - uncertainty estimates for each parameter
 */
        int i, j;
        double discr;
        double *diff, *x, *result, *mod_val, *temp;
        double **H, **K, **P;
       
        x = (double *)malloc(x_dim*sizeof(double));
        for(i = 0; i < x_dim; i++)
                x[i] = x0[i];

        result = (double *)malloc((2*x_dim + 1)*sizeof(double));
        diff = (double *)malloc(x_dim*sizeof(double));
        temp = (double *)malloc(x_dim*sizeof(double));
        j = 0;
        s = stop_val+1;

        while(j<num_iters&&s>stop_val){
                H = jacobian(forward_fun, observed, x, x_dim, 
                             data_dim, obs_dim, 0.001);
                K = kalman_gain(P0, H, R, data_dim, x_dim);
                mod_val = forward_fun(x, observed, x_dim, obs_dim);

                for(i = 0; i < data_dim; i++)
                        diff[i] = data[i] - mod_val[i];
                for(i = 0; i < x_dim; i++)
                        temp[i] = scalar_product(K[i], diff);
                for(i = 0; i < x_dim; i++){
                        x[i] = min(x[i] + temp[i], upper_bound[i]);
                        x[i] = max(x[i], lower_bounds[i]);         
                }

                P = update_cov(P0, K, H, data_dim, x_dim);
                s = discrepancy(diff, R, data_dim);                
                j++;
        }
        for(i = 0; i < x_dim; i++){
                result[i] = x[i];
                result[x_dim + i + 1] = 1 - sqrt(P[i, i]/P0[i, i]);
        }
        result[x_dim] = s;
        return result;
}

double **kalman_gain(P, H, R, data_dim, x_dim)
{       
        int i, j;
        double **K, **HT;
        HT = transpose(H, data_dim, x_dim);
        K = matrix_product(P, HT, x_dim, x_dim, x_dim, data_dim);
        K = matrix_product(H, K, data_dim, x_dim, x_dim, data_dim);
        for(i = 0; i < data_dim; i++){
                for(j = 0; j < data_dim; j++)
                        K[i][j]+=R[i][j];
        }
        K = pseudo_inv(K, data_dim, data_dim);
        K = matrix_product(HT, K, x_dim, data_dim, data_dim, data_dim);
        K = matrix_product(P, K, x_dim, x_dim, x_dim, data_dim);
        return K;
}

double **update_cov(P0, K, H, data_dim, x_dim)
{
        int i, j;
        double **P;

        P = matrix_product(K, H, x_dim, data_dim, data_dim, x_dim);
        for(i = 0; i < x_dim; i++){
                for(j = 0; j < x_dim; j++){
                        P[i][j] = -P[i][j];
                        if(i==j)
                                P[i][j]+=1;
                }
        }
        P = matrix_product(P, P0);
        return P;
}
