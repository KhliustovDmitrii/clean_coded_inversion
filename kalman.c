#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 
 *This file contains code for Kalman filter
 */

double **kalman(double **data, double **observed, double *x0,
                int num_points, int data_dim, int obs_dim, int x_dim, 
                double *(*forward_fun)(double*, double*),
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
 *        num_points: number of points in which to invert
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
 *        result: matrix, column_number = num_points, row_number = x_dim,
 *        in each column first rows contain estimated parameters, 
 *        last - uncertainty estimates for each
 */





}
