#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kalman.h"
#include "math_utils.h"

#ifndef max
#define max(a, b)        (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a, b)        (((a) < (b)) ? (a) : (b))
#endif


/* 
 *This file contains code for Kalman filter
 */

long double *kalman_sequential(long double *data, long double *observed,
                               long double *x0,
                               int data_dim, int obs_dim, int x_dim,
                        long double* (*forward_fun)(long double*, long double*),
                        long double **P0, long double **Q, long double **R,
                        long double *lower_bounds, long double *upper_bounds,
                        int num_iters, long double stop_val)
{
/*
 * This function computes Kalman inversion of data
 * using sequential filter with square root variance decomposition
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
        int i, j, k, t;
        long double s, a, gamma;
        long double *diff, *x, *result, *mod_val, *phi, *K;
        long double **S, **temp, **H, **P;
        
        S = cholesky(P0, x_dim);
        x = (long double *)malloc(x_dim*sizeof(long double));
        for(i = 0; i < x_dim; i++)
                x[i] = x0[i];

        result = (long double *)malloc((2*x_dim + 1)*sizeof(long double));
        diff = (long double *)malloc(data_dim*sizeof(long double));
        temp = (long double **)malloc(x_dim*sizeof(long double *));
        phi = (long double *)malloc(x_dim*sizeof(long double));
        K = (long double *)malloc(x_dim*sizeof(long double));

        for(i = 0; i < x_dim; i++)
                temp[i] = (long double *)malloc(x_dim*sizeof(long double));
        j = 0;
        s = stop_val+1;
        
        while(j<num_iters&&s>stop_val){
                H = jacobian(forward_fun, observed, x, x_dim,
                             data_dim, obs_dim, 0.001);
                for(i = 0; i < data_dim; i++){
                        mod_val = forward_fun(x, observed);

                        for(k = 0; k < data_dim; k++)
                                diff[k] = data[k] - mod_val[k];

                        for(k = 0; k < x_dim; k++){
                                phi[k] = 0;
                                for(t = 0; t < x_dim; t++)
                                        phi[k]+=S[t][k]*H[i][t];
                        }

                        a = 0;
                        for(k = 0; k < x_dim; k++)
                                a+=phi[k]*phi[k];
                        a+=R[i][i];
                        a = 1/a;
                        
                        gamma = 1/(1 + sqrt(a*R[i][i])); 
                        for(k = 0; k < x_dim; k++){
                                for(t = 0; t < x_dim; t++){
                                        temp[k][t] = - a*gamma*phi[k]*phi[t];
                                        if(k==t){
                                                temp[k][t]+=1;
                                        } 
                                }
                        }
                        S = matrix_product(S, temp, x_dim, x_dim, x_dim, x_dim);
                        for(k = 0; k < x_dim; k++){
                                K[k] = 0;
                                for(t = 0; t < x_dim; t++)
                                        K[k]+=S[k][t]*phi[t];
                                K[k]*=a;
                        }
                     
                        for(k = 0; k < x_dim; k++)
                                x[k] += K[k]*diff[i];               
 
                }

                for(i = 0; i < x_dim; i++){
                        x[i] = min(x[i], upper_bounds[i]);
                        x[i] = max(x[i], lower_bounds[i]);
                }

                s = discrepancy(diff, R, data_dim);
                j++;
        }

        P = matrix_product(S, transpose(S, x_dim, x_dim),
                           x_dim, x_dim, x_dim, x_dim);

        for(i = 0; i < x_dim; i++){
                result[i] = x[i];
                result[x_dim + i + 1] = 1 - sqrt(P[i][i]/P0[i][i]);
        }
        result[x_dim] = s;
        free(diff);
        free(temp);
        free(H);
        free(K);
        free(S);
        free(x);
        free(P);
        free(mod_val);
        free(phi);
        return result;
}

long double *kalman_step(long double *data, long double *observed, long double *x0,
                int data_dim, int obs_dim, int x_dim, 
                long double* (*forward_fun)(long double*, long double*),
                long double **P0, long double **Q, long double **R,
                long double *lower_bounds, long double *upper_bounds,
                int num_iters, long double stop_val)
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
        long double discr, s;
        long double *diff, *x, *result, *mod_val, *temp;
        long double **H, **K, **P;
       
        x = (long double *)malloc(x_dim*sizeof(long double));
        for(i = 0; i < x_dim; i++)
                x[i] = x0[i];

        result = (long double *)malloc((2*x_dim + 1)*sizeof(long double));
        diff = (long double *)malloc(data_dim*sizeof(long double));
        temp = (long double *)malloc(x_dim*sizeof(long double));
        j = 0;
        s = stop_val+1;

        while(j<num_iters&&s>stop_val){
                H = jacobian(forward_fun, observed, x, x_dim, 
                             data_dim, obs_dim, 0.001);
                K = kalman_gain(P0, H, R, data_dim, x_dim);
                mod_val = forward_fun(x, observed);
                for(i = 0; i < data_dim; i++)
                        diff[i] = data[i] - mod_val[i];
                for(i = 0; i < x_dim; i++){
                        temp[i] = scalar_product(K[i], diff, data_dim);
                        if(num_iters >= 5 && j >= 3)
                                temp[i]/j;
                }
                for(i = 0; i < x_dim; i++){
                        x[i] = min(x[i] + temp[i], upper_bounds[i]);
                        x[i] = max(x[i], lower_bounds[i]);         
                }

                P = update_cov(P0, K, H, data_dim, x_dim);
                s = discrepancy(diff, R, data_dim);                
                j++;
        }
        for(i = 0; i < x_dim; i++){
                result[i] = x[i];
                result[x_dim + i + 1] = 1 - sqrt(P[i][i]/P0[i][i]);
        }
        result[x_dim] = s;
        free(diff);
        free(temp);
        free(H);
        free(K);
        free(P);
        free(x);
        free(mod_val);
        return result;
}

long double **kalman_gain(long double **P, long double **H, long double **R, 
                     int data_dim, int x_dim)
{       
        int i, j;
        long double **K, **HT;
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
        free(HT);
        return K;
}

long double **update_cov(long double **P0, long double **K, long double **H, 
                    int data_dim, int x_dim)
{
        int i, j;
        long double **P;

        P = matrix_product(K, H, x_dim, data_dim, data_dim, x_dim);
        for(i = 0; i < x_dim; i++){
                for(j = 0; j < x_dim; j++){
                        P[i][j] = -P[i][j];
                        if(i==j)
                                P[i][j]+=1;
                }
        }
        P = matrix_product(P, P0, x_dim, x_dim, x_dim, x_dim);
        return P;
}
